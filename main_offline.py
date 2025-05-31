# Reconstruct images from cone data
# Data can be measured or from Gate/Allpix simulation
# Can be used with an offline venv (and should be on MacOS): see README.md

# ===========================
# ==   INPUT               ==
# ===========================

# Simulation or experiment parameters
src_pos = [0 , 0, -5]
npix, pitch, thickness = 256, 0.055, 1
sensorsize = [npix * pitch, npix * pitch, thickness]
sensortranslation = [0 , 0 , 5 ]

# Image parameters
vp, vs = 0.1, (256, 256, 256)

# Read cones from CSV file
import pandas as pd
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