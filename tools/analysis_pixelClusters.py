# Functions to process pixelClusters dataframes

from tools.analysis_pixelHits import *
from opengate.logger import global_log
from .analysis_pixelHits import PIXEL_ID, TOA, ENERGY_keV, EVENTID
import pandas as pd

pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 400)
pd.set_option('display.max_rows', 1000)
pd.set_option('display.float_format',
              lambda x: f'{x:.9}')  # G4 steps are logged with f'{x:.3}'

# Cluster coordinates differ from pixel hit coordinates:
# PIX_X_ID, PIX_Y_ID are also in pixel indices but result from some processing function
# PIX_Z_ID result from another type of processing function, i.e. depth calculation
# Indices are fractional:
#  => 0 is center of lower left pixel (e.g. 0.5 is it's right edge)
#  => z-axis points towards the readout connected to the sensor
# PIX_X_mm, PIX_Y_mm, PIX_Z_mm are also global coordinates
#  => helps to directly compare with Allpix2 output text files
#  => fits Compton camera applications when detector position is used for reconstruction
pixelClusters_columns = [PIX_X_ID, PIX_Y_ID, PIXEL_ID, TOA, ENERGY_keV]  # TODO not used

def is_adjacent_1D(hit, cluster, n_pix):
    x1, y1 = get_pixID_2D(hit[PIXEL_ID], n_pix)
    return any(
        abs(x1 - x2) <= 1 and abs(y1 - y2) <= 1
        for x2, y2 in
        (get_pixID_2D(hit[PIXEL_ID], n_pix) for _, hit in cluster.iterrows())
    )

def is_adjacent_2D(hit, cluster, n_pix):
    x1, y1 = hit[PIX_X_ID], hit[PIX_Y_ID]
    return any(
        abs(x1 - x2) <= 1 and abs(y1 - y2) <= 1
        for _, h in cluster.iterrows()
        for x2, y2 in [(h[PIX_X_ID], h[PIX_Y_ID])]
    )

def get_is_adjacent(pixelHits_columns):
    if PIXEL_ID in pixelHits_columns:
        return is_adjacent_1D
    else:
        return is_adjacent_2D

def process_cluster_method1(cluster_df):
    cluster_total_energy = cluster_df[ENERGY_keV].sum()
    cluster_first_TOA = cluster_df[TOA].min()
    cluster_first_eventID = int(cluster_df[EVENTID].min())
    return pd.DataFrame({
        EVENTID: [cluster_first_eventID],
        ENERGY_keV: [cluster_total_energy],
        TOA: [cluster_first_TOA]
    })


def process_cluster_method2(cluster, n_pixels):
    """
    X and Y are in the sensor's local coordinates system, as in Allpix2
    => origin = center of the lower-left pixel
    """
    cluster_total_energy = cluster[ENERGY_keV].sum()
    cluster_first_TOA = cluster[TOA].min()
    cluster_first_eventID = int(cluster[EVENTID].min())

    if PIX_X_ID in cluster.columns and PIX_Y_ID in cluster.columns:
        pixX = cluster[PIX_X_ID].to_numpy()
        pixY = cluster[PIX_Y_ID].to_numpy()
    else:
        pixX, pixY = zip(*cluster[PIXEL_ID].apply(get_pixID_2D, args=(n_pixels,)))

    x = sum(pixX * cluster[ENERGY_keV]) / cluster_total_energy
    y = sum(pixY * cluster[ENERGY_keV]) / cluster_total_energy
    return pd.DataFrame({
        EVENTID: [cluster_first_eventID],
        PIX_X_ID: [x],
        PIX_Y_ID: [y],
        ENERGY_keV: [cluster_total_energy],
        TOA: [cluster_first_TOA]
    })


process_cluster_functions = {
    'm1': process_cluster_method1,
    'm2': process_cluster_method2
}


def new_clust(clust_list, cluster, hit, n_pixels, process_func, **kwargs):
    process_func = process_cluster_functions[process_func]
    clust_list.append(process_func(cluster, n_pixels, **kwargs))
    new_cluster_df = pd.DataFrame([hit])
    new_time_window_start = hit[TOA]
    return new_cluster_df, new_time_window_start


# TODO speed -> https://pandas.pydata.org/docs/user_guide/basics.html#iteration
def pixelHits2pixelClusters(pixelHits, npix, window_ns, f, **kwargs):
    stime = time.time()
    global_log.info(f"Offline [pixelClusters]: START")
    if not len(pixelHits):
        global_log.error(f"Empty pixel hits dataframe, probably no hit produced.")
        global_log.info(f"Offline [pixelClusters]: {get_stop_string(stime)}")
        return []
    else:
        global_log.debug(f"Input pixel hits dataframe with ({len(pixelHits)} entries)")

    pixelHits = pixelHits.sort_values(by=TOA)

    # 1st cluster & initialization
    clust = pd.DataFrame([pixelHits.iloc[0]])  # Initialize with first hit
    wst = pixelHits.iloc[0][TOA] # window start
    clusters = []

    is_adjacent = get_is_adjacent(pixelHits_columns)

    # Loop over hits
    for index, hit in pixelHits.iloc[1:].iterrows():
        if hit[TOA] - wst <= window_ns and is_adjacent(hit, clust, npix):
            clust = pd.concat([clust, hit.to_frame().T], ignore_index=True)
        else:
            clust, wst = new_clust(clusters, clust, hit, npix, f, **kwargs)

    # Last cluster
    new_clust(clusters, clust, hit, npix, f, **kwargs)

    df = pd.concat(clusters, ignore_index=True)
    global_log.debug(f"{len(clusters)} clusters")
    global_log_debug_df(df)
    global_log.info(f"Offline [pixelClusters]: {get_stop_string(stime)}")
    return df
