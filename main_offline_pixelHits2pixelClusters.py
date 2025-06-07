# Reconstruct images from cone data
# Data can be measured or from Gate/Allpix simulation
# Can be used with an offline venv (and should be on MacOS): see README.md

import sys
import pandas as pd
import matplotlib.pyplot as plt
import logging
from tools.utils_plot import plot_hitsNclusters

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

# MEASUREMENT
pixelClusters_meas = pixelHits2pixelClusters(pixelHits_meas, npix=256, window_ns=100, f='meas_calib')

plot_hitsNclusters(pixelHits_meas, pixelClusters_meas, max_keV=300)
