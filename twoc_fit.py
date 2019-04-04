import numpy as np
import glob
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import os
from pandas_scidata import pd_import,pd_output


# physical constant
e0 = 1.6021766208E-19
h0 = 6.62607015E-34


# set your directory
search_dir = r'C:\LAB\local data\qc460_gr_we01\raw_data\QC0461GR WE I\B_sweep\Exc100mV_Bsweep_Top_gate_m0p2_0p1_m0p7V_Back_gate_gnd (05)'


files = list(filter(os.path.isfile,glob.glob(os.path.join(search_dir,'*.dat'))))
files.sort(key=lambda x:os.path.getmtime(x))

# set other parameters for input
Ref = 4958.2
spr = 22
ucols = (0,1,2,3)
nms = ['bf','uxx','uxy','curr']
step = np.linspace(-0.2,-0.7,6)
# import data
pd_data = pd_import(files,spr,ucols,nms,step,Ref)

#for t in step:
#        plt.plot(pd_data.loc[pd_data['tgate'] == t,'bf'][10:],pd_data.loc[pd_data['tgate'] == t,'rxy'][10:])


# Write structured data
#WorkingFolder = r'C:\LAB\local data\qc460_gr_we01\tg_bsweep'
#filename = 'pd_tg_bsweep03.txt'
#pd_output(pd_data,WorkingFolder,filename,step)
t = step[3]
bf = pd_data.loc[(pd_data['tgate'] == t) & (0.2>pd_data['bf']) & (pd_data['bf']>0),'bf']
hl = pd_data.loc[(pd_data['tgate'] == t) & (0.2>pd_data['bf']) & (pd_data['bf']>0),'rxy']
def func(x,a,b):
    return a+b*x

fitParams,fitCovariances = curve_fit(func,bf,hl)
plt.plot(bf,hl,"b^",bf,func(bf,*fitParams),"r-")
print(fitParams)
plt.show()

# model from PHYSICAL REVIEW B 81, 054515 (2010)
#fitfunc = lambda p, x: -x*((p[3]-p[1]*p[0]**2)+p[0]**2*p[2]**2*x**2*(p[3]-p[1]))/((p[0]*p[1]+p[3])**2+p[0]**2*p[2]**2*x**2*(p[3]-p[1])**2)
#errfunc = lambda p, x, y: abs(fitfunc(p,x)-y)
#p0 = [1,1,1,1]
#p1,success = optimize.leastsq(errfunc,p0[:],args = (bf,hl))
#print(p1)
#plt.plot(bf,hl,"b^",bf,fitfunc(p0,bf),"r-")
#plt.show()










