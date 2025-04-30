# Utility functions to analyse output files
# Can be used in the main simulation script after sim.run() or offline (i.e. reading root files without simulation)
import os
import sys
import pandas
import uproot
import SimpleITK as sitk
import matplotlib.pyplot as plt
from pandas import Series
from tools.utils import *
from opengate.utility import g4_units

pandas.set_option('display.max_columns', 100)
pandas.set_option('display.width', 400)
pandas.set_option('display.max_rows', 1000)
pandas.set_option('display.float_format', lambda x: f'{x:.3}')  # G4 steps are logged with f'{x:.3}'


# List of all possible attributes for Hits:
# ['Direction', 'EventDirection', 'EventID', 'EventKineticEnergy', 'EventPosition', 'GlobalTime', 'HitUniqueVolumeID',
# 'KineticEnergy', 'LocalTime', 'PDGCode', 'ParentID', 'ParentParticleName', 'ParticleName', 'ParticleType', 'Position',
# 'PostDirection', 'PostKineticEnergy', 'PostPosition', 'PostPositionLocal', 'PostStepUniqueVolumeID', 'PostStepVolumeCopyNo',
# 'PreDirection', 'PreDirectionLocal', 'PreGlobalTime', 'PreKineticEnergy', 'PrePosition', 'PrePositionLocal', 'PreStepUniqueVolumeID',
# 'PreStepVolumeCopyNo', 'ProcessDefinedStep', 'RunID', 'StepLength', 'ThreadID', 'TimeFromBeginOfEvent', 'TotalEnergyDeposit',
# 'TrackCreatorModelIndex', 'TrackCreatorModelName', 'TrackCreatorProcess', 'TrackID', 'TrackLength', 'TrackProperTime',
# 'TrackVertexKineticEnergy', 'TrackVertexMomentumDirection', 'TrackVertexPosition', 'TrackVolumeCopyNo', 'TrackVolumeInstanceID',
# 'TrackVolumeName', 'UnscatteredPrimaryFlag', 'Weight']
# Obtained with print(opengate_core.GateDigiAttributeManager.GetInstance().GetAvailableDigiAttributeNames())
def analyse_hits(file_path):
    tree_hits = uproot.open(file_path)['Hits']
    print(int(tree_hits.num_entries), 'entries in tree Hits')
    # print(tree.keys())
    hits = tree_hits.arrays(library='pd', entry_stop=None)  # None to read all entries
    # print('Number of events', hits['EventID'].nunique())
    # print_hits_short(hits)
    # print_hits_short_sortedByGlobalTime(hits)
    # print(hits.to_string(index=False))
    # print(hits.groupby('EventID').first().to_string(index=False))
    # print(hits[hits['ParticleName'] == 'gamma'].to_string(index=False))
    # print(hits[hits['TrackCreatorProcess'] == 'compt'].to_string(index=False))
    # print(hits['TotalEnergyDeposit'].sum())
    # print(hits[hits['TrackID']==2]['TotalEnergyDeposit'].sum())
    # print(hits.groupby('EventID').first())
    # print_hits_inG4format(hits)
    # print_hits_inG4format_sortedByGlobalTime(hits)
    return hits


# List of available attributes for Singles:
#  EventID  TrackID  ParentID  ParentParticleName  ParticleName  KineticEnergy  TotalEnergyDeposit  TrackCreatorProcess
#  ProcessDefinedStep     Position_X     Position_Y    Position_Z  PreStepUniqueVolumeID  PostPosition_X  PostPosition_Y
#  PostPosition_Z GlobalTime
def analyse_singles(file_path):
    tree_singles = uproot.open(file_path)['Singles']
    print(int(tree_singles.num_entries), 'entries in tree Singles')
    singles = tree_singles.arrays(library='pd', entry_stop=None)  # None to read all entries
    # print(singles[['EventID','TotalEnergyDeposit','KineticEnergy','HitUniqueVolumeID']].to_string(index=False))
    # print(singles[singles['TrackCreatorProcess'] == 'compt'].to_string(index=False))
    # print(Series(singles['PreStepUniqueVolumeID'].to_numpy()).value_counts(normalize=True) * 100,'\n')  # !! entry_stop = None  !!
    # singles = singles.loc[:, ~singles.columns.str.contains('Position', case=False)]
    return singles


def plot_DigitizerProjectionActor(sim):
    Bq, sec = g4_units.Bq, g4_units.s
    proj = sim.actor_manager.get_actor("Projection")
    file_path = sim.output_dir + '/' + proj.output_filename  # Replace with the actual path to your .mhd file
    image = sitk.ReadImage(file_path)
    im = sitk.GetArrayFromImage(image)[0, :, :].astype(int)
    plt.imshow(im, cmap='gray',vmax=2)
    source = sim.source_manager.get_source("source")
    events = f'{source.n} events' if source.n else f'{int(source.activity / Bq)}Bq {sum_time_intervals(sim.run_timing_intervals) / sec} sec'
    plt.title(f'{source.particle} {source.energy.mono} MeV \n {events}')
    cbar = plt.colorbar(label='number of pixel hits summed over all events')
    cbar.set_ticks(np.arange(np.min(im), np.max(im) + 1))
    plt.show()


def plot_hits_TotalEnergyDeposit(file_path, bins=100):
    if not os.path.isfile(file_path):
        sys.exit(f"File {file_path} does not exist, probably no hit produced...")
    hits = uproot.open(file_path)['Hits'].arrays(library='pd')  # None to read all entries
    plt.hist(hits['TotalEnergyDeposit'], bins=bins)
    plt.xlabel('TotalEnergyDeposit [MeV]')
    plt.ylabel('Counts')
    plt.show()


def plot_hits_TotalEnergyDeposit_sumPerEvent(file_path, bins=100):
    if not os.path.isfile(file_path):
        sys.exit(f"File {file_path} does not exist, probably no hit produced...")
    hits = uproot.open(file_path)['Hits'].arrays(library='pd')  # None to read all entries
    plt.hist(hits.groupby('EventID')['TotalEnergyDeposit'].sum(), bins=bins)
    plt.xlabel('Sum of TotalEnergyDeposit per EventID [MeV]')
    plt.ylabel('Counts')
    plt.show()
