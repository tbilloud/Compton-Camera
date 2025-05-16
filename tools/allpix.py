import subprocess
from scipy.spatial.transform import Rotation as R
import warnings
from tools.analysis_pixelHits import *


def run_allpix(sim,
               binary_path='allpix/allpix-squared/install-noG4/bin/',
               output_dir='allpix/', log_level='FATAL',
               config='default'):
    # ==========================
    # == INPUTS & INIT        ==
    # ==========================

    # Measure time spent with Allpix (run, and if used, weighting potential calculation)
    stime = time.time()

    # Fetch inputs
    hits_actor = sim.actor_manager.get_actor("Hits")
    hits_file = sim.output_dir + '/' + hits_actor.output_filename
    gateHits_df = uproot.open(hits_file)['Hits'].arrays(library='pd')
    sensor = sim.volume_manager.get_volume("sensor")
    source = sim.source_manager.get_source("source")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        angles = R.from_matrix(sensor.rotation).as_euler('xyz', degrees=True)

    # Prevent two competing visualizations
    if sim.visu is True:
        sys.exit("Allpix cannot be run with Gate visualization enabled")
    else:
        global_log.info(f"Offline [Allpix2]: START")
        global_log.debug(f"Input {hits_file}, {len(gateHits_df)} gHits")

    # Check that Gate geometry was adapted to Allpix
    try:
        pixel = sim.volume_manager.get_volume("pixel_param")
    except Exception:
        global_log.error(
            f"Pixels must be defined with RepeatParametrisedVolume() in Gate. Call the volume 'pixel'.")
        sys.exit()

    # Prepare the weighting potential file
    if config == 'precise':
        wp_fname = f"pitch{int(pixel.translation[0] * 1000)}um_thick{int(sensor.size[2] * 1000)}um"
    elif config == 'fast':
        if sensor.material != 'Silicon':
            global_log.error(
                f"The 'fast' configuration only works with Silicon sensors.")
            sys.exit()

    # ==========================
    # == PRODUCE CONFIG FILES ==
    # ==========================

    geometry_conf_content = f"""[0_0]
type = "detector_model"
position = {" ".join([f"{sensor.translation[i]}mm" for i in range(3)])}
orientation = {" ".join([f"{angles[i]}deg" for i in range(3)])}
    """

    detector_model_conf_content = f"""type = "hybrid"
geometry = "pixel"
number_of_pixels = {pixel.linear_repeat[0]} {pixel.linear_repeat[1]}
pixel_size = {pixel.translation[0]}mm {pixel.translation[1]}mm
sensor_thickness = {sensor.size[2]}mm
sensor_material = "{sensor.material}"
bump_sphere_radius = 9.0um
bump_cylinder_radius = 7.0um
bump_height = 20.0um
    """

    # TODO: allow setting important parameters in main gate script
    # TODO: is Jacoboni model valid for CdTe / GaAs ?
    # TODO: use CSADigitizer in 'precise' config
    # TODO: speed up precise config (e.g. charge groups etc)
    # TODO: deal with different units for pixelHits (energy/charge/bits for TOT, time/bits for TOA)
    configurations = {
        "fast": """
    [ElectricFieldReader]
    model = "linear"
    bias_voltage = -1000V # pixel side, - to collect electrons, + to collect holes
    [ProjectionPropagation] 
    # mobility model is Jacoboni 
    temperature = 293K
    integration_time = 1000s # default 25ns might stop charge propagation
    [PulseTransfer]
    # timestep = 1.6ns # 0.01ns by default, but Timepix3 clock is 1.6ns
    [DefaultDigitizer]
    threshold = 1e # 0e turns off ToA... 
    threshold_smearing = 0e
    electronics_noise = 0e
    """,
        "default": """ 
    [ElectricFieldReader]
    model = "constant"
    bias_voltage = -1000V # pixel side, - to collect electrons, + to collect holes
    [GenericPropagation]
    integration_time = 1us # default 25ns might stop propagation
    mobility_model = "constant"
    mobility_electron = 1000cm*cm/V/s
    mobility_hole = 100cm*cm/V/s
    [PulseTransfer]
    timestep = 1.6ns # 0.01ns by default, but Timepix3 clock is 1.6ns
    [DefaultDigitizer]
    threshold = 1e # 0e turns off ToA... 
    threshold_smearing = 0e
    electronics_noise = 0e
    """,
        "precise": f""" 
    [ElectricFieldReader]
    model="constant"
    bias_voltage=-1000V
    [WeightingPotentialReader]
    model = "mesh"
    file_name = "{wp_fname if config=='precise' else ''}_weightingpotential.apf"
    field_mapping = "PIXEL_FULL"
    [TransientPropagation]
    mobility_model = "constant"
    mobility_electron = 1000cm*cm/V/s
    mobility_hole = 500cm*cm/V/s # holes always propagated
    integration_time = 1us  # default 25ns might stop propagation
    timestep = 0.1ns
    distance = 0
    [PulseTransfer]
    [DefaultDigitizer]
    threshold = 1e # 0e turns off ToA... 
    threshold_smearing = 0e
    electronics_noise = 0e
    """
    }

    main_conf_content = f"""[Allpix]
log_level = {log_level}
log_format = "DEFAULT"
detectors_file = "geometry.conf"
number_of_events = {source.n if source.n else gateHits_df['EventID'].max() + 1}
model_paths = ["."]
output_directory = "."
random_seed = 1
[DepositionReader]
model = "root"
file_name = "../{hits_file}"
tree_name = "Hits"
detector_name_chars = 3
branch_names = ["EventID", "TotalEnergyDeposit", "GlobalTime", "Position_X", "Position_Y", "Position_Z", "HitUniqueVolumeID", "PDGCode", "TrackID", "ParentID"]
{configurations[config]}
[TextWriter]
include = "PixelHit"
    """

    with open(output_dir + 'geometry.conf', 'w') as geometry_conf_file:
        geometry_conf_file.write(geometry_conf_content)

    with open(output_dir + 'detector_model.conf',
              'w') as detector_model_conf_file:
        detector_model_conf_file.write(detector_model_conf_content)

    with open(output_dir + 'main.conf', 'w') as main_conf_file:
        main_conf_file.write(main_conf_content)

    # ===========================
    # === RUN BINARIES        ===
    # ===========================

    if config == 'precise':
        if os.path.isfile(f'allpix/{wp_fname}_weightingpotential.apf'):
            global_log.info(f"Offline [Allpix2]: Using {wp_fname + '_weightingpotential.apf'}")
        else:
            global_log.warning(f"Weighting potential file not found. Generating it...")
            subprocess.run([binary_path + 'generate_potential', '--model',
                            output_dir + 'detector_model.conf', '--output',
                            f'allpix/{wp_fname}', '-v', log_level],
                           check=True)

    subprocess.run([binary_path + 'allpix', '-c', output_dir + 'main.conf'], check=True)

    event_time_offset_flag = False
    if source.n:
        global_log.warning(
            f"Use source.activity (not source.n) for realistic TOA. Adding 1us/event.")
        event_time_offset_flag = True

    global_log.info(f"Offline [Allpix2]: {get_stop_string(stime)}")

    return event_time_offset_flag


# TODO: I've seen negative ToT values in data.txt
def gHits2allpix2pixelHits(sim, npix,
                           binary_path='allpix/allpix-squared/install-noG4/bin/',
                           config='default',
                           log_level='FATAL'):
    time_offset = run_allpix(sim, binary_path, output_dir='allpix/',
                             log_level=log_level, config=config)
    pixelHits = allpixTxt2pixelHit('allpix/data.txt', n_pixels=npix)
    if time_offset: pixelHits[TOA] += pixelHits.groupby(EVENTID).ngroup() * 1e3
    return pixelHits
    # Lines starting with PixelHit in data.txt have:
    # PixelHit X_ID, Y_ID, TOT, TOA, global_time, X_global, Y_global, Z_global
