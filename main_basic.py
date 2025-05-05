import opengate_core
from opengate.managers import Simulation
from opengate.geometry.volumes import *
from tools.analysis_basics import *
from tools.point_source_validation import *
from tools.reco_backprojection import *
from tools.allpix import *

if __name__ == "__main__":
    sim, sim.output_dir = Simulation(), "output"
    um, mm, keV, MeV, deg, Bq, sec = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.s
    sim.volume_manager.add_material_database('tools/GateMaterials.db')
    sim.random_engine, sim.random_seed = "MersenneTwister", 1
    sim.visu = False

    # ===========================
    # ==   GEOMETRY            ==
    # ===========================
    npix, pitch, thickness = 256, 55 * um, 1 * mm
    sim.world.material = "Vacuum"
    sensor = sim.add_volume("Box", "sensor")
    sensor.material = "cadmium_telluride"
    sensor.size = [npix * pitch, npix * pitch, thickness]
    sensor.translation = [0 * um, 0 * um, 5 * mm]
    pixel = sim.add_volume("Box", "pixel")
    pixel.mother, pixel.size = sensor.name, [pitch, pitch, thickness]
    pixel.material = sensor.material
    par = RepeatParametrisedVolume(repeated_volume=pixel)
    par.linear_repeat, par.translation = [npix, npix, 1], [pitch, pitch, 0]
    sim.volume_manager.add_volume(par)

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

    # BASICS
    hits_path = Path(sim.output_dir) / hits.output_filename
    singles_path = Path(sim.output_dir) / singles.output_filename
    analyse_hits(hits_path)
    analyse_singles(singles_path)

    # PIXEL HITS (optional)
    # pixelHits = singles2pixelHits(singles_path)
    # pixelHits = gHits2allpix2pixelHits(sim, npix)
    # plot_pixelHits_perEventID(pixelHits,n_pixels=npix,log_scale=[False, False, True])

    # CONES
    cones = gHits2cones_byEvtID(hits_path, source.energy.mono)

    # POINT SOURCE VALIDATION
    sp = source.position.translation
    validate_psource(cones, source_pos=sp, vpitch=0.1, vsize=(256, 256, 256),
                     plot_seq=False, plot_stack=True, plot_napari=True)

    # RECONSTRUCTION
    d = {'size': sensor.size, 'position': sensor.translation}
    reco_bp(cones, vpitch=0.1, vsize=(256, 256, 256), napari=True, det=d)
