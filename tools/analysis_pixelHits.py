# Functions to process pixelHits dataframes

import os
import sys
import time

import pandas
import pandas as pd
import uproot
import matplotlib.pyplot as plt
from tools.utils import *
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator

from tools.utils import get_pixID
from opengate.logger import global_log

pandas.set_option('display.max_columns', 100)
pandas.set_option('display.width', 400)
pandas.set_option('display.max_rows', 1000)
pandas.set_option('display.float_format', lambda x: f'{x:.9}')  # G4 steps x:.3

PIXEL_ID = 'PixelID (int16)'
TOA = 'ToA (ns)'
ENERGY_keV = 'Energy (keV)'
PIX_X_ID = 'X'  # pixel X index (starts from 0, bottom left)
PIX_Y_ID = 'Y'  # pixel Y index (starts from 0, bottom left)
TOT = 'ToT'
pixelHits_columns = [PIX_X_ID, PIX_Y_ID, TOA, TOT, ENERGY_keV, PIXEL_ID]
EVENTID = 'EventID'
simulation_columns = [EVENTID]  # from Gate


def singles2pixelHits(file_path):
    if not os.path.isfile(file_path):
        sys.exit(f"{file_path} does not exist, probably no hit produced...")
    else:
        global_log.info(f"Offline [pixelHits]: START")
        global_log.debug(f"Input {file_path}")
    stime = time.time()
    singles = uproot.open(file_path)['Singles'].arrays(library='pd')
    global_log.debug(f"Number of singles: {len(singles)}")
    singles['HitUniqueVolumeID'] = singles['HitUniqueVolumeID'].astype(
        str).str.replace(r'0_', '', regex=True)
    singles.rename(columns={'HitUniqueVolumeID': PIXEL_ID}, inplace=True)
    singles[PIXEL_ID] = singles[PIXEL_ID].str.replace('pixel_param-', '', regex=False)
    singles[PIXEL_ID] = singles[PIXEL_ID].astype(int)
    singles.rename(columns={'TotalEnergyDeposit': ENERGY_keV}, inplace=True)
    singles[ENERGY_keV] = singles[ENERGY_keV] * 1e3  # Convert MeV to keV
    singles.rename(columns={'GlobalTime': TOA}, inplace=True)
    x, y = zip(*singles[PIXEL_ID].apply(get_pixID_2D, args=(256,)))
    singles[PIX_X_ID] = x
    singles[PIX_Y_ID] = y
    singles[TOT] = singles[ENERGY_keV] * 1e3  # TODO temporary
    global_log.debug(f"Number of pixel hits: {len(singles)}")
    global_log.info(f"Offline [pixelHits]: {get_stop_string(stime)}")
    return singles[simulation_columns + pixelHits_columns]


def pixelHits_fig_ax(pixelHits_df, n_pixels, fig, ax,
                     log_scale=[False, False, False]):
    df, np = pixelHits_df, n_pixels
    x, y = zip(*df[PIXEL_ID].apply(get_pixID_2D, args=(np,)))

    nc, ne, nt = [mcolors.LogNorm() if log else None for log in log_scale]

    hc = ax[0].hist2d(x, y, bins=[np] * 2, range=[[0, n_pixels]] * 2, norm=nc)
    cb = fig.colorbar(hc[3], ax=ax[0], label='Count')
    cb.locator = MaxNLocator(integer=True)
    cb.update_ticks()
    ax[0].set_title('Counts')

    he = ax[1].hist2d(x, y, bins=[np] * 2, weights=df[ENERGY_keV],
                      range=[[0, np]] * 2, norm=ne,
                      vmin=0.5 * df[ENERGY_keV].min() if not ne else None)
    fig.colorbar(he[3], ax=ax[1], label='Energy (keV)')
    ax[1].set_title('Energy')

    ht = ax[2].hist2d(x, y, bins=[np] * 2, weights=df[TOA],
                      range=[[0, np]] * 2, norm=nt,
                      vmin=0.9 * df[TOA].min() if not nt else None)
    fig.colorbar(ht[3], ax=ax[2], label='ToA (ns)')
    ax[2].set_title('Time')

    for a in ax:
        a.set_aspect('equal')
        a.set_xlabel('Pixel x')
        a.set_ylabel('Pixel y')
        a.xaxis.set_major_locator(MaxNLocator(integer=True))
        a.yaxis.set_major_locator(MaxNLocator(integer=True))

    return fig, ax


def plot_pixelHits_perEventID(pixelHits_df, n_pixels,
                              log_scale=[False, False, False]):
    unique_event_ids = pixelHits_df[EVENTID].unique()
    for event_id in unique_event_ids:
        fig, ax = plt.subplots(1, 3, figsize=(16, 4))
        df = pixelHits_df[pixelHits_df[EVENTID] == event_id]
        pixelHits_fig_ax(df, n_pixels, fig, ax, log_scale)
        plt.suptitle(f'Event ID: {event_id}')
        plt.tight_layout()
        plt.show()


def plot_pixelHits_comparison(pixelHits_df1, pixelHits_df2, n_pixels,
                              log_scale=[False, False, False]):
    fig, ax = plt.subplots(2, 3, figsize=(12, 6))
    pixelHits_fig_ax(pixelHits_df1, n_pixels, fig, ax[0], log_scale)
    pixelHits_fig_ax(pixelHits_df2, n_pixels, fig, ax[1], log_scale)
    plt.tight_layout()
    plt.show()


def plot_pixelHits_comparison_perEventID(pixelHits_df1, pixelHits_df2,
                                         n_pixels,
                                         log_scale=[False, False, False]):
    unique_event_ids = pixelHits_df1[EVENTID].unique()
    assert (unique_event_ids == pixelHits_df2[EVENTID].unique()).all()
    for event_id in unique_event_ids:
        df1 = pixelHits_df1[pixelHits_df1[EVENTID] == event_id]
        df2 = pixelHits_df2[pixelHits_df2[EVENTID] == event_id]
        fig, ax = plt.subplots(2, 3, figsize=(11, 6))
        pixelHits_fig_ax(df1, n_pixels, fig, ax[0], log_scale)
        pixelHits_fig_ax(df2, n_pixels, fig, ax[1], log_scale)
        plt.suptitle(f'Event ID: {event_id}')
        plt.tight_layout()
        plt.show()


def pixelHits2burdaman(pixelHits_df, out_path):
    # TODO set types correctly (else visu with TrackLab will not work)
    # => https://software.utef.cvut.cz/tracklab/manual/a01627.html
    # insert a column with 0s at the 3rd position
    pixelHits_df.insert(2, 'fTOA', 0)
    print(pixelHits_df)
    pixelHits_df.to_csv(out_path, header=False, index=False, sep='\t')

    # Dummy header
    # TODO replace values with NaNs
    custom_header = """# Start of measurement: 10/1/2017 17:34:41.8467094
# Start of measurement - unix time: 1506872081.846
# Chip ID: H3-W00036
# Readout IP address: 192.168.1.105
# Back-end location: Satigny, CH
# Detector mode: ToA & ToT
# Readout mode: Data-Driven Mode
# Bias voltage: 229.72V
# THL = 1570 (0.875V)
# Sensor temperature: 58.9°C
# Readout temperature: 42.9°C
# ------- Internal DAC values ---------------
# Ibias_Preamp_ON:\t128\t(1.208V)
# Ibias_Preamp_OFF:\t8\t(1.350V)
# VPreamp_NCAS:\t\t128\t(0.702V)
# Ikrum:\t\t15\t(1.128V)
# Vfbk:\t\t164\t(0.891V)
# Vthreshold_fine:\t505\t(0.877V)
# Vthreshold_coarse:\t7\t(0.875V)
# Ibias_DiscS1_ON:\t100\t(1.109V)
# Ibias_DiscS1_OFF:\t8\t(1.321V)
# Ibias_DiscS2_ON:\t128\t(0.396V)
# Ibias_DiscS2_OFF:\t8\t(0.256V)
# Ibias_PixelDAC:\t128\t(0.984V)
# Ibias_TPbufferIn:\t128\t(1.169V)
# Ibias_TPbufferOut:\t128\t(1.077V)
# VTP_coarse:\t\t128\t(0.693V)
# VTP_fine:\t\t256\t(0.724V)
# Ibias_CP_PLL:\t\t128\t(0.557V)
# PLL_Vcntrl:\t\t128\t(0.874V)
# BandGap output:\t--- \t(0.684V)
# BandGap_Temp:\t\t--- \t(0.733V)
# Ibias_dac:\t\t--- \t(1.241V)
# Ibias_dac_cas:\t\t--- \t(1.004V)
# DACs: \t128\t8\t128\t15\t164\t505\t7\t100\t8\t128\t8\t128\t128\t128\t128\t256\t128\t128
# DACs Scans: \t1.208V\t1.350V\t0.702V\t1.128V\t0.891V\t0.877V\t0.875V\t1.109V\t1.321V\t0.396V\t0.256V\t0.984V\t1.169V\t1.077V\t0.693V\t0.724V\t0.557V\t0.874V\t0.684V\t0.733V\t1.241V\t1.004V
# -----------------------------------------------------------------------------------------------------------------------------
"""

    # Read the CSV file and add the custom header
    with open(out_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    lines.insert(0, custom_header)

    # Write the modified content back to the file
    with open(out_path, 'w', encoding='utf-8') as file:
        file.writelines(lines)


def allpixTxt2pixelHit(text_file, n_pixels=256):
    global_log.info(f"Offline [pixelHits]: START")
    global_log.debug(f"Input {text_file}")
    # TODO adapt to different simulation chains

    stime = time.time()
    rows = []

    with open(text_file, "r") as file:
        event_id = None
        for line in file:
            line = line.strip()

            if line.startswith("==="):
                event_id = int(
                    line.split()[1]) - 1  # allpix adds 1 to event ID
                continue

            if line.startswith("---"):
                continue

            if line.startswith("PixelHit"):
                parts = line.split()
                x, y = int(parts[1].strip(',')), int(parts[2].strip(','))
                pixel_id = get_pixID(x, y, n_pixels=n_pixels)
                tot = float(parts[3].strip(','))
                toa = float(parts[4].strip(','))  # ToA from event start
                global_time = float(parts[5].strip(','))  # ToA from simu start
                position_x = float(parts[6].strip(','))
                position_y = float(parts[7].strip(','))
                position_z = float(parts[8].strip(','))

                rows.append({
                    EVENTID: event_id,
                    PIX_X_ID: x,
                    PIX_Y_ID: y,
                    PIXEL_ID: pixel_id,
                    TOT: tot,
                    ENERGY_keV: tot * 4.43 / 1000,
                    # TODO: adapt to qdc_resolution (on/off) in DefaultDigitizer
                    TOA: global_time
                })

    df = pd.DataFrame(rows, columns=simulation_columns + pixelHits_columns)
    global_log_debug_df(df)
    global_log.info(f"Offline [pixelHits]: {get_stop_string(stime)}")
    return df
