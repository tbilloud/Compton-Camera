from pathlib import Path
import opengate_core
from opengate.managers import Simulation
from opengate.utility import g4_units
from tools.analysis_cones import gHits2cones_byEvtID, pixelClusters2cones_byEvtID
from tools.analysis_pixelHits import *
from tools.analysis_pixelClusters import pixelHits2pixelClusters
from tools.point_source_validation import *
from tools.reco_backprojection import *
from tools.utils_opengate import setup_pixels, theta_phi

if __name__ == "__main__":
    sim, sim.output_dir = Simulation(), "output"
    um, mm, keV, MeV, deg, Bq, sec = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.s
    sim.volume_manager.add_material_database('tools/GateMaterials.db')
    sim.random_engine, sim.random_seed = "MersenneTwister", 1
    sim.visu = False
    sim.verbose_level = 'DEBUG'

    # ===========================
    # ==   GEOMETRY            ==
    # ===========================
    npix, pitch, thickness = 256, 55 * um, 1 * mm
    sim.world.material = "Vacuum"
    sensor = sim.add_volume("Box", "sensor")
    sensor.material = "cadmium_telluride"
    sensor.size = [npix * pitch, npix * pitch, thickness]
    sensor.translation = [0 * um, 0 * um, 5 * mm]
    setup_pixels(sim, npix, sensor, pitch, thickness)

    ## ===========================
    ## ==  PHYSICS              ==
    ## ===========================
    # sim.physics_manager.physics_list_name = 'G4EmLivermorePhysics'
    sim.physics_manager.global_production_cuts.gamma = 1 * um
    sim.physics_manager.global_production_cuts.electron = 1 * um
    sim.physics_manager.em_parameters.update(
        {'fluo': True, 'pixe': True, 'deexcitation_ignore_cut': False,
         'auger': True, 'auger_cascade': True})

    ## =============================
    ## == ACTORS                  ==
    ## =============================
    hits = sim.add_actor('DigitizerHitsCollectionActor', 'Hits')
    hits.attached_to = sensor.name
    hits.authorize_repeated_volumes = True
    hits.attributes = opengate_core.GateDigiAttributeManager.GetInstance().GetAvailableDigiAttributeNames()
    hits.output_filename = 'hits.root'
    singles = sim.add_actor("DigitizerAdderActor", "Singles")
    singles.authorize_repeated_volumes = True
    singles.input_digi_collection = "Hits"
    singles.policy = "EnergyWeightedCentroidPosition"
    singles.output_filename = 'singles.root'

    ## ============================
    ## == SOURCE                 ==
    ## ============================
    source = sim.add_source("GenericSource", "source")
    source.activity, sim.run_timing_intervals = 10 * Bq, [[0, 10 * sec]]
    source.particle = "gamma"
    source.energy.mono = 140 * keV
    source.position.translation = [0 * mm, 0 * mm, -5 * mm]
    source.direction.theta, source.direction.phi = theta_phi(sensor, source)
    sim.world.size = get_worldSize(sensor, source, margin=0.1)

    ##=====================================================
    ##   RUN
    ##=====================================================
    sim.run()

    ##=====================================================
    ##   ANALYSIS AND RECONSTRUCTION
    ##=====================================================

    # INPUTS
    hits_path = Path(sim.output_dir) / hits.output_filename
    singles_path = Path(sim.output_dir) / singles.output_filename

    # ################# PIXEL HITS ########################
    pixelHits = singles2pixelHits(singles_path)
    # plot_pixelHits_perEventID(pixelHits,n_pixels=npix,log_scale=[False, False, True])

    # ################# PIXEL CLUSTERS ####################
    pixelClusters = pixelHits2pixelClusters(pixelHits, npix=npix, window_ns=100, f='m2')

    # #################### CONES ##########################
    # =======> GROUND TRUTH <=======
    cn1 = gHits2cones_byEvtID(hits_path, source.energy.mono)
    # # =========> TIMEPIX <==========
    sp = charge_speed_mm_ns(mobility_cm2_Vs=1000, bias_V=1000, thick_mm=sensor.size[2])
    cn2 = pixelClusters2cones_byEvtID(pixelClusters,
                                            source_MeV=source.energy.mono,
                                            thickness_mm=thickness,
                                            charge_speed_mm_ns=sp,
                                            to_global=[npix,sensor]  # for global coord
                                            )

    # ########## VALIDATION WITH POINT SOURCE #############
    sp = source.position.translation
    vs = (256, 256, 256)
    stk1 = valid_psource(cn1, src_pos=sp, vpitch=0.1, vsize=vs, plot_seq=0, plot_stk=1)
    stk2 = valid_psource(cn2, src_pos=sp, vpitch=0.1, vsize=vs, plot_seq=0, plot_stk=1)
    # plot_stack_napari(sk, vsize = vs) # TODO check with np and cp

    # ################## RECONSTRUCTION ####################
    d = {'size': sensor.size, 'position': sensor.translation}
    vol = reco_bp(cn1, vpitch=0.1, vsize=vs, det=d)
    # TODO: add 3d plot with matplotlib and optionally napari
