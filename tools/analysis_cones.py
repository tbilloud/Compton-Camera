# Utility functions to analyse output files
# Can be used in the main simulation script after sim.run() or offline (i.e. reading root files without simulation)
# E1 = energy deposited in the Compton scattering, as in CCMod paper

import os
import sys
import pandas
import uproot
from analysis_pixelHits import PIX_X_ID, PIX_Y_ID, EVENTID, ENERGY_keV, TOA
from tools.utils import *
from opengate.logger import global_log

pandas.set_option('display.max_columns', 100)
pandas.set_option('display.width', 400)
pandas.set_option('display.max_rows', 1000)
pandas.set_option('display.float_format', lambda x: f'{x:.9}')  # G4 step: .3

# TODO make order flexible (see below)
cones_columns = ['EventID', 'Apex_X', 'Apex_Y', 'Apex_Z', 'Direction_X',
                 'Direction_Y', 'Direction_Z', 'cosT', 'error']


# TODO: can be optimized using hits.keep_zero_edep = True in simulation settings
def gHits2cones_byEvtID(file_path, source_MeV):
    if not os.path.isfile(file_path):
        sys.exit(f"File {file_path} does not exist, probably no hit produced.")
    else:
        global_log.info(f"Offline [cones ghits]: START")
        global_log.debug(f"Input {file_path}")

    stime = time.time()
    hits = uproot.open(file_path)['Hits'].arrays(library='pd')
    grouped = hits.groupby('EventID')
    cones = []

    n_events_primary = 0
    n_events_full_edep = 0
    for eventid, grp in grouped:
        apex, direction, E1 = False, False, False
        # Sensor received primary gamma and it interacted TODO: is this correct with radioisotope source?
        if 1 in grp['TrackID'].values:
            n_events_primary += 1
            # All primary energy was deposited
            if round(grp['TotalEnergyDeposit'].sum(), 6) == source_MeV:
                # TODO round above is to avoid float precision issues
                n_events_full_edep += 1
                grp = grp.sort_values('GlobalTime')  # IMPORTANT !
                h = grp.iloc[0]
                # Gamma interacts via Compton, step has dE !=0 and is stored (recoil e- not tracked)
                if h['TrackID'] == 1 and grp['TrackID'].value_counts()[1] > 1:
                    # if value_counts()[1] == 1, TrackID 1 stopped at 1st step via photoelec (without prior Compton)
                    # TODO: what about rayleigh scattering and pair production?
                    apex = [h[f'PostPosition_{axis}'] for axis in 'XYZ']
                    direction = [-h[f'PostDirection_{axis}'] for axis in 'XYZ']
                    E1 = h['TotalEnergyDeposit']
                # Gamma interacts via Compton, step has dE = 0 and is not stored, but recoil e- tracked with TrackID=2
                # However I can't use direction of recoil e-... Need to go further
                elif h['TrackID'] == 2 and h['TrackCreatorProcess'] == 'compt':
                    apex = [h[f'PrePosition_{axis}'] for axis in 'XYZ']
                    E1 = h['KineticEnergy']

                    # Remove TrackID 2 and its descendants from group
                    def find_descendants(df, part_id):
                        descendants = set()
                        child = df[df['ParentID'] == part_id]['TrackID'].values
                        for child in child:
                            descendants.add(child)
                            descendants.update(find_descendants(df, child))
                        return descendants

                    desc_of_2 = find_descendants(grp, 2)
                    grp = grp[~grp['TrackID'].isin(desc_of_2.union({2}))]
                    h2 = grp.iloc[0]
                    # if post-Compton step of TrackID 1 has dE != 0, it is stored and is the next one in the
                    # time-sorted group, and it gives the direction
                    if h2['TrackID'] == 1:
                        direction = [-h2[f'PreDirection_{ax}'] for ax in 'XYZ']
                    # if not, there is a new track whose origin can be used to calculate the direction
                    else:
                        prepos = [h2[f'PrePosition_{axis}'] for axis in 'XYZ']
                        diff = np.array(apex) - np.array(prepos)
                        direction = (diff / np.linalg.norm(diff)).tolist()
        if apex:
            cosT = 1 - (0.511 * E1) / (source_MeV * (source_MeV - E1))
            cones.append([eventid] + apex + direction + [cosT] + [200])
            # TODO make order flexible

    df = pandas.DataFrame(cones, columns=cones_columns)
    global_log.debug(f"{n_events_primary} events with primary particle hitting sensor")
    global_log.debug(f"=> {n_events_full_edep} with full energy deposited in sensor")
    global_log.debug(f"  => {len(cones)} with at least one Compton interaction")
    global_log.info(f"Offline [cones ghits]: {len(cones)} cones")
    global_log_debug_df(df)
    global_log.info(f"Offline [cones ghits]: {get_stop_string(stime)}")
    return df


def pixelClusters2cones_byEvtID(pixelClusters, source_MeV, thickness_mm,
                                charge_speed_mm_ns, to_global=False):
    """
    Clusters have:
    - X/Y coordinates
    - ToA
    - ToT
    Cones need:
    - Apex (X,Y,Z)
    - Direction (X,Y,Z)
    - cosT
    - error

    to_global: use this to return cones in global coordinate system instead of local
    It is a list with [npix, sensor]:
     - npix is the number of pixels in one side of the sensor (256 for Timepix3)
     - sensor is an object with:
        sensor.size = list with x,y,z lengths in mm
        sensor.translation = list with x,y,z positions of the sensor's center in mm
        sensor.rotation = 3D rotation matrix
    """

    global_log.info(f"Offline [cones tpx]: START")
    global_log.debug(f"Input pixel clusters dataframe")
    stime = time.time()

    grouped = pixelClusters.groupby(EVENTID)
    grouped = [group for group in grouped if len(group[1]) == 2]

    cones = []

    for eventid, group in grouped:
        # 1) Distinguish compton vs photo-electric interactions
        group = group.sort_values(ENERGY_keV)
        Esum_MeV = 0.001 * (group.iloc[0][ENERGY_keV] + group.iloc[1][ENERGY_keV])
        if abs(Esum_MeV - source_MeV) < 0.1 and group.iloc[1][ENERGY_keV] > get_E1max(source_MeV):
            cl_photoel = group.iloc[1]
            cl_compton = group.iloc[0]
        else:
            continue

        # 2) Calculate depth difference
        dZ_mm = charge_speed_mm_ns * (cl_compton[TOA] - cl_photoel[TOA])
        dZ_frac = dZ_mm / thickness_mm

        # 3) Calculate absolute depth of Compton interaction (apex)
        z_compton = 0  # middle of sensor (in local fractional unit)
        # TODO or use cluster size/energy ?

        # 4) Complete 3D positions
        pos_compton = [cl_compton[PIX_X_ID], cl_compton[PIX_Y_ID], z_compton]
        pos_photoel = [cl_photoel[PIX_X_ID], cl_photoel[PIX_Y_ID], z_compton + dZ_frac]

        # 5) Construct cone
        E1_MeV = cl_compton[ENERGY_keV] / 1000
        cosT = 1 - (0.511 * E1_MeV) / (source_MeV * (source_MeV - E1_MeV))
        if to_global:
            npix, sensor = to_global
            apex = localFractional2globalCoordinates(pos_compton, sensor, npix)
            pos_photoel = localFractional2globalCoordinates(pos_photoel, sensor, npix)
            direction = np.array(apex) - np.array(pos_photoel)
        else:
            apex = pos_compton
            direction = np.array(apex) - np.array(pos_photoel)
        direction = (direction / np.linalg.norm(direction)).tolist()
        cones.append([eventid] + apex + direction + [cosT] + [200])
        # TODO make order flexible

    df = pandas.DataFrame(cones, columns=cones_columns)
    global_log.info(f"Offline [cones tpx]: {len(cones)} cones")
    global_log_debug_df(df)
    global_log.info(f"Offline [cones tpx]: {get_stop_string(stime)}")
    return df


def get_E1max(source_MeV):
    me = 0.511  # MeV
    E1max = source_MeV ** 2 / (source_MeV + me / 2)
    return E1max
