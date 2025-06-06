# Reconstruct images from cone data
# Data can be measured or from Gate/Allpix simulation
# Can be used with an offline venv (and should be on MacOS): see README.md

import sys
import pandas as pd
import matplotlib.pyplot as plt
import logging

try:
    from opengate.logger import global_log
    global_log.setLevel(logging.DEBUG)
except ImportError:
    global_log = logging.getLogger("dummy")
    global_log.addHandler(logging.NullHandler())

# ===========================
# ==   INPUT PIXEL HITS    ==
# ===========================

# SIMULATION
pixelHits = pd.read_csv('output/pixelHits_250kBq_100ms.csv')

# MEASUREMENT
from tools.analysis_pixelHits import pixet2pixelHit
file_t3pa = '/home/billoud/DATA/data-driven/Lu177_300s.t3pa'
pixelHits_meas = pixet2pixelHit(file_t3pa, calib='./minipix/', max_rows=10000)

# ===========================
# ==   PIXEL CLUSTERS      ==
# ===========================
from tools.analysis_pixelClusters import pixelHits2pixelClusters

# SIMULATION
pixelClusters = pixelHits2pixelClusters(pixelHits, npix=256, window_ns=100, f='meas_calib')
plt.hist(pixelClusters['Energy (keV)'], bins=300, range=(0, 300))
plt.xlabel('Energy (keV)'),plt.ylabel('Count'),plt.title('Single Clusters SIMULATED')
plt.show()

# MEASUREMENT
# print(pixelHits_meas)
pixelClusters_meas = pixelHits2pixelClusters(pixelHits_meas, npix=256, window_ns=100, f='meas_calib')
plt.hist(pixelClusters_meas['Energy (keV)'], bins=300, range=(0, 300))
plt.xlabel('Energy (keV)'),plt.ylabel('Count'),plt.title('Single Clusters MEASURED')
plt.show()

sys.exit()

# ===========================
# ==   INPUT CONES         ==
# ===========================

# Simulation or experiment parameters
src_pos = [0 , 0, -5]
npix, pitch, thickness = 256, 0.055, 1
sensorsize = [npix * pitch, npix * pitch, thickness]
sensortranslation = [0 , 0 , 5 ]

# Image parameters
vp, vs = 0.1, (256, 256, 256)

# Read cones from CSV file
cones = pd.read_csv('output/cones_comb.csv')
print(cones)

# ===========================
# = POINT SOURCE VALIDATION =
# ===========================
# # => Only in case a point source was used
# from tools.point_source_validation import valid_psource
# valid_psource(cones, src_pos=src_pos, vpitch=vp, vsize=vs, plot_seq=0, plot_stk=1)

# ===========================
# == 3D RECONSTRUCTION     ==
# ===========================
# Using back-projection
from tools.reco_backprojection import reco_bp
d = {'size': sensorsize, 'position': sensortranslation}
vol = reco_bp(cones, vpitch=0.1, vsize=vs, det=d)

# ===========================
# == 3D VISUALIZATION      ==
# ===========================
# Using napari
# On MacOS, the 'offline' venv is necessary, since opengate is not compatible with napari
from tools.napari_plot import plot_reconstruction_napari
import napari
d = {'size': sensorsize, 'position': sensortranslation}
plot_reconstruction_napari(vol, vs, vp, detector=d)