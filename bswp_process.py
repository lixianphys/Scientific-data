# library
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from data_import import import_col
from bswp_plot import bswp_plot

# data import
dir = r'\\system01\QT\Lixian\data\QC 0492 0h\QC 0492 0H-exc500mv\bsweep\gate1p2to0V (02)'
os.chdir(dir)
fname = sorted(glob.glob('*gate1p2to0V (02).dat')) # sort the filenames by their value
print(fname)
cols = (0,4,5,-1)
skiprows = 40
Dic_data = import_col(dir,fname,cols,skiprows,'bs',-1,1)
# define physical constants
e = 1.6021766208E-19
h = 6.62607015E-34
# data processing
# choose dataset for Rxx - temp dependence plot
Volt_pick = [0.3,0.6,0.84,0,1.2]
add_path = r'\\system01\QT\Lixian\data\QC 0492 0h\QC 0492 0H-exc500mv\gatesweep\m6to4V (03).dat'
fig_title = 'bsweep_at_{}V_exc500mV'.format(Volt_pick)
# data plot
fig = bswp_plot(Dic_data,fname,Volt_pick,add_path,(15,18),40,(0,-4),fig_title,-1)
fig[0].savefig(fig_title+'.pdf')

