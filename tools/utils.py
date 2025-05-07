# Some utility functions
# WARNING: For print functions, make sure that dataframe columns are present in simulation settings (c.f. actor attribtues)

import time
import numpy as np
from opengate.logger import global_log
from opengate.utility import g4_units

um, mm, keV, MeV, deg, Bq, sec = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.s


# Prints hits like G4 steps are logged via sim.g4_verbose_level_tracking
# If pandas.set_option('display.float_format'...) is used in script calling the function, remove it
def print_hits_inG4format(hits_df):
    print(
        hits_df[['EventID', 'PostPosition_X', 'PostPosition_Y', 'PostPosition_Z',
                 'KineticEnergy', 'TotalEnergyDeposit',
                 'StepLength', 'TrackLength', 'HitUniqueVolumeID',
                 'ProcessDefinedStep', 'ParticleName', 'TrackID',
                 'ParentID', 'ParentParticleName',
                 'TrackCreatorProcess', 'TrackCreatorModelName'
                 ]])


# Prints only few relevant columns from hits tree
def print_hits_short(hits_df):
    print(hits_df[['EventID', 'TrackID', 'ParticleName', 'ParentID',
                   'ParentParticleName', 'KineticEnergy',
                   'TotalEnergyDeposit', 'ProcessDefinedStep',
                   'TrackCreatorProcess', 'GlobalTime',
                   'HitUniqueVolumeID']].to_string(index=False))


def print_hits_long(hits_df):
    print(hits_df[[
        'EventID', 'TrackID', 'ParticleName', 'ParentID', 'ParentParticleName',
        'KineticEnergy',
        'TotalEnergyDeposit', 'ProcessDefinedStep', 'TrackCreatorProcess',
        'PrePosition_X', 'PrePosition_Y', 'PrePosition_Z', 'PostPosition_X',
        'PostPosition_Y', 'PostPosition_Z',
        'PreDirection_X', 'PreDirection_Y', 'PreDirection_Z',
        'PostDirection_X', 'PostDirection_Y', 'PostDirection_Z'
    ]].to_string(index=False))


# Prints directional info from hits tree
# PreDirection and PostDirection seem to be the same very frequently but not always
def print_hits_direction(hits_df):
    print(hits_df[[
        'EventID', 'TrackID', 'ParticleName', 'ProcessDefinedStep',
        'PrePosition_X', 'PrePosition_Y', 'PrePosition_Z', 'PostPosition_X',
        'PostPosition_Y', 'PostPosition_Z',
        'PreDirection_X', 'PreDirection_Y', 'PreDirection_Z',
        'PostDirection_X', 'PostDirection_Y', 'PostDirection_Z'
    ]].to_string(index=False))


# Prints gamma interactions
def print_hits_gammas(hits_df):
    hits_df = hits_df[hits_df['ParticleName'] == 'gamma']
    print(hits_df[[
        'EventID', 'TrackID', 'ParentID', 'KineticEnergy',
        'TotalEnergyDeposit', 'ProcessDefinedStep', 'TrackCreatorProcess',
        'PostPosition_Z'
    ]].to_string(index=False))


# Prints time info
def print_hits_time(hits_df):
    print(hits_df[[
        'EventID', 'TrackID', 'ParticleName', 'GlobalTime', 'PreGlobalTime',
        'LocalTime', 'TimeFromBeginOfEvent',
        'TrackProperTime',
    ]].to_string(index=False))


# Prints processes info
def print_hits_processes(hits_df):
    print(hits_df[[
        'EventID', 'TrackID', 'ParticleName', 'ProcessDefinedStep',
        'TrackCreatorProcess', 'TrackCreatorModelName'
    ]].to_string(index=False))


def print_hits_inG4format_sortedByGlobalTime(hits_df):
    hits_df = hits_df.groupby('EventID').apply(
        lambda x: x.sort_values('GlobalTime'))
    print_hits_inG4format(hits_df)


def print_hits_long_sortedByGlobalTime(hits_df):
    hits_df = hits_df.groupby('EventID').apply(
        lambda x: x.sort_values('GlobalTime'))
    print_hits_long(hits_df)


def get_pixID(x, y, n_pixels=256):
    return x * n_pixels + y


def get_pixID_2D(pixel_id, n_pixels=256):
    x = pixel_id // n_pixels
    y = pixel_id % n_pixels
    return x, y


def sum_time_intervals(time_intervals):
    return sum([time_interval[1] - time_interval[0] for time_interval in
                time_intervals])


# Limit emission angle of source particles to the sensor area
def theta_phi(sensor, source):
    # TODO: only works if source and sensor have same x,y coordinates
    sensor_position = np.array(sensor.translation)
    source_position = np.array(source.position.translation)
    sensor_size = np.max(sensor.size[0:1])
    distance = np.linalg.norm(sensor_position - source_position) - sensor.size[
        2] / 2
    phi_deg = 180 - np.degrees(np.arctan(sensor_size / (2 * distance)))
    return [phi_deg * deg, 180 * deg], [0, 360 * deg]


def get_worldSize(sensor, source, margin=0.1):
    stype = source.position.type
    if stype not in ["point", "sphere", "box"]:
        raise ValueError(
            "Function get_worldSize() is only implemented for point/sphere/box sources")
    ssize = source.position.size if stype in ['point', 'box'] else [
                                                                       source.position.radius] * 3
    return [np.max([abs(st) + sz / 2, abs(sp) + ss / 2]) * (2 + margin) for
            st, sz, sp, ss in
            zip(sensor.translation, sensor.size, source.position.translation,
                ssize)]


def get_file_name(sim, doppler, fluo):
    sources_dict = sim.source_manager.sources
    if len(sources_dict) == 1:
        source = next(iter(sources_dict.values()))
    else:
        raise NotImplementedError(
            "Function get_file_name() is only implemented for one source")
    events = (f'{source.n}events' if source.n else
              f'{int(source.activity / Bq)}Bq_'
              f'{int(sum_time_intervals(sim.run_timing_intervals)) / sec}sec')
    energy = int(source.energy.mono / keV)
    return f'source{energy}keV_{events}_doppler{doppler}_fluo{fluo}.root'


def coordinateOrigin2arrayCenter(cp_array, vpitch, vsize):
    cp_array[:, 0] = cp_array[:, 0] + vpitch * vsize[0] / 2
    cp_array[:, 1] = cp_array[:, 1] + vpitch * vsize[1] / 2
    cp_array[:, 2] = cp_array[:, 2] + vpitch * vsize[2] / 2
    return cp_array


def coordinateOrigin2arrayCenter_df(df, vpitch, vsize):
    df['Apex_X'] = df['Apex_X'] + vpitch * vsize[0] / 2
    df['Apex_Y'] = df['Apex_Y'] + vpitch * vsize[1] / 2
    df['Apex_Z'] = df['Apex_Z'] + vpitch * vsize[2] / 2
    return df


def get_stop_string(stime):
    return f"STOP. Time: {time.time() - stime:.1f} seconds.\n" + '-' * 80


def global_log_debug_df(df):
    """
    Print preview (head) of the dataframe
    """
    if not df.empty:
        global_log.debug(f"Output preview:\n{df.head().to_string(index=False)}")


def localFractional2globalCoordinates(c, sensor, npix):
    pitch = sensor.size[0] / npix  # mm

    a = sensor.translation
    b = [-(npix / 2 - 0.5)] * 2 + [0]

    a, b, c = np.array(a), np.array(b), np.array(c)

    bc = b + c
    bc = np.dot(sensor.rotation, bc) * [pitch, pitch, sensor.size[2]]

    g = a + bc

    return g.tolist()


def global2localFractionalCoordinates(g, sensor, npix):
    pitch = sensor.size[0] / npix  # mm

    a = sensor.translation
    b = [-(npix / 2 - 0.5)] * 2 + [0]

    a, b, g = np.array(a), np.array(b), np.array(g)

    bc = -a + g
    bc = np.dot(sensor.rotation.T, bc) / [pitch, pitch, sensor.size[2]]

    c = bc - b

    return c.tolist()


def charge_speed_mm_ns(mobility_cm2_Vs, bias_V, thick_mm):
    """
    Assuming constant electric field (ohmic type sensors) and mobility, calculate speed of charges
    In mm per ns
    """

    efield = bias_V / (thick_mm / 10)  # [V/cm]
    elec_speed = mobility_cm2_Vs * efield  # [cm*cm/V/s] * [V/cm] => [cm/s]
    return elec_speed * 1e-8  # [mm/ns]
