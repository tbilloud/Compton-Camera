from pathlib import Path
import opengate_core
from opengate.managers import Simulation
from opengate.utility import g4_units
from tools.analysis_cones import gHits2cones_byEvtID, pixelClusters2cones_byEvtID
from tools.analysis_pixelClusters import *
from tools.point_source_validation import *
from tools.allpix import *
from tools.utils_opengate import setup_pixels, set_fluorescence

um, mm, keV, MeV, deg, Bq, sec = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.s

if __name__ == "__main__":
    check_gate_version()
    sim, sim.output_dir = Simulation(), "output"
    um, mm, keV, MeV, deg, Bq, ms, sec = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.ms, g4_units.s
    sim.volume_manager.add_material_database('tools/GateMaterials.db')
    sim.random_engine, sim.random_seed = "MersenneTwister", 1
    sim.visu = False
    sim.verbose_level = 'DEBUG' # DEBUG for data preview, INFO for algo timing only

    # ===========================
    # ==   GEOMETRY            ==
    # ===========================
    npix, pitch, thickness = 256, 55 * um, 1 * mm
    sim.world.material = "Vacuum"
    sensor = sim.add_volume("Box", "sensor")
    sensor.material = "cadmium_telluride" # or 'Silicon'
    sensor.size = [npix * pitch, npix * pitch, thickness]
    sensor.translation = [0 * mm, 0 * mm, 5 * mm]
    # sensor.rotation = R.from_euler('xyz', [0,45,0], degrees=True).as_matrix()
    setup_pixels(sim, npix, sensor, pitch, thickness)

    ## ===========================
    ## ==  PHYSICS              ==
    ## ===========================
    sim.physics_manager.physics_list_name = 'G4EmLivermorePhysics' # for Doppler effect
    set_fluorescence(sim)

    ## =============================
    ## == ACTORS                  ==
    ## =============================
    hits = sim.add_actor('DigitizerHitsCollectionActor', 'Hits')
    hits.attached_to = sensor.name
    hits.authorize_repeated_volumes = True
    hits.attributes = opengate_core.GateDigiAttributeManager.GetInstance().GetAvailableDigiAttributeNames()
    hits.output_filename = 'hits.root'
    # hits.keep_zero_edep = True # TODO compatible with gHits2cones_byEventID ?

    ## ============================
    ## == SOURCE                 ==
    ## ============================
    source = sim.add_source("GenericSource", "source")
    source.n = 108
    # source.activity, sim.run_timing_intervals = 100_000 * Bq, [[0, 2 * ms]]
    source.particle = "gamma"
    source.energy.mono = 140 * keV
    source.position.translation = [0 * mm, 0 * mm, -5 * mm]
    source.direction.type, source.direction.momentum = "momentum", [0, 0, 1]
    # source.direction.theta, source.direction.phi = theta_phi(sensor, source)
    sim.world.size = get_worldSize(sensor, source, margin=10)

    ## ============================
    ## ==  RUN                   ==
    ## ============================
    sim.run()

    ## ============================
    ## ==  OFFLINE ANALYSIS      ==
    ## ============================

    hits_path = Path(sim.output_dir) / hits.output_filename
    hits_df = uproot.open(hits_path)['Hits'].arrays(library='pd')

    # ################# PIXEL HITS ########################
    pixelHits = gHits2allpix2pixelHits(sim, npix, config='default', log_level='FATAL')

    # ################# PIXEL CLUSTERS ####################
    pixelClusters = pixelHits2pixelClusters(pixelHits, npix=npix, window_ns=100, f='m2')

    # #################### CONES ##########################
    # =======> GROUND TRUTH <=======
    ctruth = gHits2cones_byEvtID(hits_path, source.energy.mono)
    ctruth.to_csv(Path(sim.output_dir) / 'cones_truth.csv', index=False)
    # # =========> TIMEPIX <==========
    spd = charge_speed_mm_ns(mobility_cm2_Vs=1000, bias_V=1000, thick_mm=sensor.size[2])
    ctpx = pixelClusters2cones_byEvtID(pixelClusters,
                                            source_MeV=source.energy.mono,
                                            thickness_mm=thickness,
                                            charge_speed_mm_ns=spd,
                                            to_global=[npix,sensor]  # for global coord
                                            )
    ctpx.to_csv(Path(sim.output_dir) / 'cones_timepix.csv', index=False)

    # ########## VALIDATION WITH POINT SOURCE #############
    sp, vp, vs = source.position.translation, 0.1, (256, 256, 256)
    sth = valid_psource(ctruth, src_pos=sp, vpitch=vp, vsize=vs, plot_seq=0, plot_stk=1)
    stpx = valid_psource(ctpx, src_pos=sp, vpitch=vp, vsize=vs, plot_seq=0, plot_stk=1)