import sys
import opengate_core
from opengate.managers import Simulation
from opengate.geometry.volumes import *
from tools.analysis_pixelClusters import *
from tools.point_source_validation import *
from tools.reco_backprojection import *
from tools.allpix import *

if __name__ == "__main__":
    sim, sim.output_dir = Simulation(), "output"
    um, mm, keV, MeV, deg, Bq, ms, sec = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.ms, g4_units.s
    sim.volume_manager.add_material_database('tools/GateMaterials.db')
    sim.random_engine, sim.random_seed = "MersenneTwister", 1
    sim.visu = False
    sim.verbose_level = 'INFO' # DEBUG for data preview, INFO for algo timing only

    # ===========================
    # ==   GEOMETRY            ==
    # ===========================
    npix, pitch, thickness = 256, 55 * um, 1 * mm
    sim.world.material = "Vacuum"
    sensor = sim.add_volume("Box", "sensor")
    sensor.material = "cadmium_telluride"
    sensor.size = [npix * pitch, npix * pitch, thickness]
    sensor.translation = [0 * mm, 0 * mm, 5 * mm]
    # sensor.rotation = R.from_euler('xyz', [0,45,0], degrees=True).as_matrix()
    if not sim.visu: # because 256 x 256 pixels are too heavy for visualization
        pixel = sim.add_volume("Box", "pixel")
        pixel.mother, pixel.size = sensor.name, [pitch, pitch, thickness]
        pixel.material = sensor.material
        par = RepeatParametrisedVolume(repeated_volume=pixel)
        par.linear_repeat, par.translation = [npix, npix, 1], [pitch, pitch, 0]
        sim.volume_manager.add_volume(par)
    else:
        global_log.warning("Simulation was ran without pixels, analysis will not work.")

    ## ===========================
    ## ==  PHYSICS              ==
    ## ===========================
    doppler = False
    fluo = False
    if doppler: sim.physics_manager.physics_list_name = 'G4EmLivermorePhysics'
    if fluo:
        sim.physics_manager.global_production_cuts.gamma = 1 * um
        sim.physics_manager.global_production_cuts.electron = 100 * um
    sim.physics_manager.em_parameters.update(
        {'fluo': fluo, 'pixe': fluo, 'deexcitation_ignore_cut': False,
         'auger': fluo, 'auger_cascade': fluo})
    # TODO: deexcitation_ignore_cut impacts number of hits, and depends on cuts

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
    source.n = 1000
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
    pixelHits = gHits2allpix2pixelHits(sim, npix, config='fast', log_level='FATAL')
    # TODO: check warning about (mobility?) model not adapted to CdTe

    # ################# PIXEL CLUSTERS ####################
    pixelClusters = pixelHits2pixelClusters(pixelHits, npix=npix, window_ns=100, f='m2')

    # #################### CONES ##########################
    # =======> GROUND TRUTH <=======
    ctruth = gHits2cones_byEvtID(hits_path, source.energy.mono)
    # # =========> TIMEPIX <==========
    sp = charge_speed_mm_ns(mobility_cm2_Vs=1000, bias_V=1000, thick_mm=sensor.size[2])
    ctpx = pixelClusters2cones_byEvtID(pixelClusters,
                                            source_MeV=source.energy.mono,
                                            thickness_mm=thickness,
                                            charge_speed_mm_ns=sp,
                                            to_global=[npix,sensor]  # for global coord
                                            )

    # ########## VALIDATION WITH POINT SOURCE #############
    sp, vp, vs = source.position.translation, 0.1, (256, 256, 256)
    sth = valid_psource(ctruth, src_pos=sp, vpitch=vp, vsize=vs, plot_seq=0, plot_stk=1)
    stpx = valid_psource(ctpx, src_pos=sp, vpitch=vp, vsize=vs, plot_seq=0, plot_stk=1)

    # ################## RECONSTRUCTION ####################
    # TODO