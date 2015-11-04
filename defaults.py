"""
defaults.py

Default simulator constants. Not included are channel-dependent spectrograph
parameters, which need to be given as C constants. Those are listed in
"armparams.h"
"""
import numpy as np

DEFAULT_DWT = 2.0 # [Angstrom] @MZ 2 not 10
DEFAULT_NS = 30 # @CM 60 was too high
DEFAULT_NWF = 10    # default nw factor
DEFAULT_DW = 0.001    # default dw factor -> for 5500 A, R>~820000
DEFAULT_NR = 1e8 # default number of rays for MonteCarlo sampling
NORM_CCD_VAL = 40000.0
SIM_LIST = ['B', 'C', 'S', 'L', 'T', 'F', 'W', 'P', 'U', 'X', '0']
DATA_DTYPE = np.float64
RESOLUTION = 82000.0
MEM_LIM = 5e7
DMAX_SNR = 200.0
TELLURIC_LINES = "CH4 CO2 H2O N2O O2 O3"

# noise defaults
sig_d = 0.05    # [e-/s] dark current
sig_rn = 5.0   # [e-] readout noise
gain = 1.00     # [ADU/e-] HARPS value is 0.75
time = 600.0     # [s] integration time
bias = 300      # [DN] for now, arbitrary value
qe = 0.70       # % quantum efficiency
