import os
from physconst import *
## io settings
DEFAULT_PATH = os.path.join(os.getcwd(),'output')
ALLOW_FORMAT = ['pdf','png','jpeg','tiff']
DEFAULT_FORMAT = 'pdf'
## plot settings
DEFAULT_FIGURE_SIZE = (10,10)
DEFAULT_LW = 2
DEFAULT_LEGEND_POS = (0.5, 0., 0.5, 0.5)
DEFAULT_LEGEND_LOC = 'best'
DEFAULT_EBPLOT_XLABEL = '$B$ [T]'
DEFAULT_EBPLOT_YLABEL = '$E [eV]$'
DEFAULT_NBPLOT_XLABEL = '$B$ [T]'
DEFAULT_NBPLOT_YLABEL = '$n$ [1/m$^2$]'
DEFAULT_DOSBPLOT_XLABEL = '$B$ [T]'
DEFAULT_DOSBPLOT_YLABEL = 'DOS [a.u.]'
# model settings
SIGMA_COND = 1e-3*e0
SIGMA_VAL = 3e-4*e0
