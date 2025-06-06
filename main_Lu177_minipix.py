from pathlib import Path
import opengate_core
from opengate.managers import Simulation
from opengate.utility import g4_units
from tools.analysis_cones import gHits2cones_byEvtID, pixelClusters2cones_byEvtID
from tools.analysis_pixelClusters import *
from tools.point_source_validation import *
from tools.allpix import *
from tools.utils_opengate import setup_pixels, theta_phi, get_isotope_data, \
    set_fluorescence
from tools.utils_plot import plot_hitsNclusters

um, mm, keV, MeV, deg, Bq, sec = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.s

if __name__ == "__main__":
    check_gate_version()
    sim, sim.output_dir = Simulation(), "output"
    um, mm, keV, MeV, deg, Bq, ms, sec, min, hour, day, year = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.ms, g4_units.s, g4_units.min, g4_units.hour, g4_units.day, g4_units.year
    kBq, MBq = 1_000 * Bq, 1_000_000 * Bq
    sim.volume_manager.add_material_database('tools/GateMaterials.db')
    sim.random_engine, sim.random_seed = "MersenneTwister", 1
    sim.visu = False
    sim.verbose_level = 'DEBUG' # DEBUG for data preview, INFO for algo timing only

    # ===========================
    # ==   GEOMETRY            ==
    # ===========================
    npix, pitch, thickness = 256, 55 * um, 1 * mm
    sim.world.material = "Air"
    sensor = sim.add_volume("Box", "sensor")
    sensor.material = "cadmium_telluride" # or 'Silicon'
    sensor.size = [npix * pitch, npix * pitch, thickness]
    sensor.translation = [0 * mm, 0 * mm, 10 * mm]
    setup_pixels(sim, npix, sensor, pitch, thickness)

    ## ===========================
    ## ==  PHYSICS              ==
    ## ===========================
    sim.physics_manager.physics_list_name = 'G4EmLivermorePhysics' # For Doppler effect
    set_fluorescence(sim) # For fluorescence and Auger electrons
    sim.physics_manager.enable_decay = True # For radioactive sources

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
    source.activity, sim.run_timing_intervals = 250 * kBq, [[0, 10 * ms]] # TODO: file size wrong above 100ms
    source.particle, source.half_life = 'ion 71 177', 6.65 * day # Lu177
    global_log.debug(get_isotope_data(source))
    source.direction.theta, source.direction.phi = theta_phi(sensor, source)
    sim.world.size = get_worldSize(sensor, source, margin=10)

    ## ============================
    ## ==  RUN                   ==
    ## ============================
    sim.run()

    ## ============================
    ## ==  OFFLINE ANALYSIS      ==
    ## ============================

    hits_path = Path(sim.output_dir) / hits.output_filename
    singles_path = Path(sim.output_dir) / singles.output_filename

    # ################# PIXEL HITS ########################
    #pixelHits = gHits2allpix2pixelHits(sim, npix, config='precise', log_level='FATAL')
    #pixelHits.to_csv(f'output/pixelHits_{int(source.activity/kBq)}kBq_{int(sim.run_timing_intervals[0][1]/ms)}ms.csv', index=False)
    pixelHits = singles2pixelHits(singles_path)

    # ################# PIXEL CLUSTERS ####################
    pixelClusters = pixelHits2pixelClusters(pixelHits, npix=npix, window_ns=100, f='m2')

    plot_hitsNclusters(pixelHits, pixelClusters, max_keV=300)
