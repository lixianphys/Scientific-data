# data import
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from data_import import import_col
from gswp_plot import gswp_plot

e = 1.6021766208E-19
h = 6.62607015E-34

# data import
dir = r'\\system01\QT\Lixian\data\QC 0492 0h\QC 0492 0H-difTemp\gatesweep\B-10T-Extn1mV_I_38_Rxx_BE_Rxy_BC_Gate4tom6V (02)'
os.chdir(dir)
fname = sorted(glob.glob('TempControl*(01).dat')) # sort the filenames by their value
print(fname)
cols = (0,-4,-2,1)
skiprows = 45
Dic_data = import_col(dir,fname,cols,skiprows,'gs',1,1)
# choose dataset for Rxx - temp dependence plot
Volt_pick = [-6,-3,-2,-1,0,2,4]
fig_title = 'gatesweep_at_diffTemp_10T'
# data plot
fig = gswp_plot(Dic_data,fname,Volt_pick,(15,18),fig_title)
fig[0].savefig(fig_title+'.pdf')