import os
from physconst import *

DEFAULT_PATH = os.path.join(os.getcwd(),'output')
ALLOW_FORMAT = ['pdf','png','jpeg','tiff']
DEFAULT_FORMAT = 'pdf'
DEFAULT_FIGURE_SIZE = (10,10)
DEFAULT_LW = 2
DEFAULT_LEGEND_POS = (0.5, 0., 0.5, 0.5)
DEFAULT_LEGEND_LOC = 'best'
SIGMA_COND = 1e-3*e0
SIGMA_VAL = 3e-4*e0
