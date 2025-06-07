import sys
import pandas as pd
import matplotlib.pyplot as plt
import logging
from tools.analysis_cones import pixelClusters2cones_byEvtID
from tools.utils_plot import plot_hitsNclusters
from tools.analysis_pixelClusters import pixelHits2pixelClusters
from tools.point_source_validation import valid_psource

try:
    from opengate.logger import global_log
    global_log.setLevel(logging.DEBUG)
except ImportError:
    global_log = logging.getLogger("dummy")
    global_log.addHandler(logging.NullHandler())


pixelHits = pd.read_csv('output/pixelHits_250kBq_100ms.csv')
pixelClusters = pixelHits2pixelClusters(pixelHits, npix=256, window_ns=100, f='meas_calib')
plot_hitsNclusters(pixelHits, pixelClusters, max_keV=300)

# Simulation or experiment parameters
src_pos = [0 , 0, -5]
npix, pitch, thickness = 256, 0.055, 1
sensorsize = [npix * pitch, npix * pitch, thickness]
sensortranslation = [0 , 0 , 5 ]

cones = pixelClusters2cones_byEvtID(pixelClusters,
                                   source_MeV=0.092,
                                   thickness_mm=thickness,
                                   charge_speed_mm_ns=spd,
                                   to_global=[npix, sensor]  # for global coord
                                   )
vp, vs = 0.1, (256, 256, 256)
valid_psource(cones, src_pos=src_pos, vpitch=vp, vsize=vs, plot_seq=0, plot_stk=1)