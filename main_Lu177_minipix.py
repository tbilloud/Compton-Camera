from pathlib import Path
import opengate_core
from opengate.managers import Simulation
from opengate.utility import g4_units
from tools.analysis_cones import gHits2cones_byEvtID, pixelClusters2cones_byEvtID
from tools.analysis_pixelClusters import *
from tools.point_source_validation import *
from tools.allpix import *
from tools.utils_opengate import setup_pixels, theta_phi

um, mm, keV, MeV, deg, Bq, sec = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.s

if __name__ == "__main__":
    check_gate_version()
    sim, sim.output_dir = Simulation(), "output"
    um, mm, keV, MeV, deg, Bq, ms, sec, min, hour, day, year = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.ms, g4_units.s, g4_units.min, g4_units.hour, g4_units.day, g4_units.year
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
    doppler = True
    fluo = True
    if doppler: sim.physics_manager.physics_list_name = 'G4EmLivermorePhysics'
    if fluo:
        sim.physics_manager.global_production_cuts.gamma = 1 * um
        sim.physics_manager.global_production_cuts.electron = 100 * um
    sim.physics_manager.em_parameters.update(
        {'fluo': fluo, 'pixe': fluo, 'deexcitation_ignore_cut': False,
         'auger': fluo, 'auger_cascade': fluo})
    # TODO: deexcitation_ignore_cut impacts number of hits, and depends on cuts
    sim.physics_manager.enable_decay = True # for radioactive sources

    ## =============================
    ## == ACTORS                  ==
    ## =============================
    hits = sim.add_actor('DigitizerHitsCollectionActor', 'Hits')
    hits.attached_to = sensor.name
    hits.authorize_repeated_volumes = True
    hits.attributes = opengate_core.GateDigiAttributeManager.GetInstance().GetAvailableDigiAttributeNames()
    hits.output_filename = 'hits.root'

    ## ============================
    ## == SOURCE                 ==
    ## ============================
    source = sim.add_source("GenericSource", "source")
    source.activity, sim.run_timing_intervals = 1_000_000 * Bq, [[0, 100 * ms]]
    source.particle, source.half_life = 'ion 71 177', 6.65 * day # Lu177
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
    hits_df = uproot.open(hits_path)['Hits'].arrays(library='pd')

    # ################# PIXEL HITS ########################
    pixelHits = gHits2allpix2pixelHits(sim, npix, config='default', log_level='FATAL')

    # ################# PIXEL CLUSTERS ####################
    pixelClusters = pixelHits2pixelClusters(pixelHits, npix=npix, window_ns=100, f='m2')
    plt.hist(pixelClusters['Energy (keV)'], bins=240, range=(10, 250))
    plt.xlabel('Energy (keV)')
    plt.ylabel('Count')
    plt.show()

