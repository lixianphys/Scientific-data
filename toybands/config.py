import os
import matplotlib.pyplot as plt
from physconst import *
## IO settings
DEFAULT_PATH = os.path.join(os.getcwd(),'output')
ALLOW_FORMAT = ['pdf','png','jpeg','tiff']
DEFAULT_FORMAT = 'pdf'
## plot settings
DEFAULT_FIGURE_SIZE = (10,10)
DEFAULT_LW = 2
DEFAULT_LEGEND_POS = (0.5, 0., 0.5, 0.5)
DEFAULT_LEGEND_LOC = 'best'
DEFAULT_XLABEL = '$B$ [T]'
DEFAULT_EBPLOT_YLABEL = '$E [eV]$'
DEFAULT_NBPLOT_YLABEL = '$n$ [1/m$^2$]'
DEFAULT_DOSMAP_YLABEL = '$n$ [1/m$^2$]'
DEFAULT_DOSBPLOT_YLABEL = 'DOS [a.u.]'

plt.rc('lines', lw=1, color='k')  # thicker black lines
plt.rc('grid', c='0.5', ls='-', lw=0.5)  # solid gray grid lines
plt.rc('savefig', dpi=600)  # higher res outputs
plt.rc("font", size=20, family='arial', weight='light')
plt.rc("axes", labelsize=20, titlesize=20, linewidth=1)
plt.rc("xtick", direction='in', labelsize=20)
plt.rc('xtick.major', size=10, pad=7)
plt.rc('xtick.minor', size=5, pad=7, visible=True)
plt.rc("ytick", direction='in', labelsize=20)
plt.rc('ytick.major', size=10, pad=7)
plt.rc('ytick.minor', size=5, pad=7, visible=True)
plt.rc("legend", fontsize=20)
plt.rcParams['figure.constrained_layout.use'] = True

# model settings
SIGMA_COND = 1e-3*e0
SIGMA_VAL = 3e-4*e0
D_PARAM = 5.2e-21*e0
# D_PARAM = 0
