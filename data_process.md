## data processing after import
```python
# data import
import data_import
import numpy as np
import matplotlib.pyplot as plt
from data_import import fname
from data_import import Dic_data

e = 1.6021766208E-19
h = 6.62607015E-34
# data process

# choose dataset for Rxx - temp dependence plot
Voltpick = [-6,-1.2,0.5,4]

# plot
newf = plt.figure(figsize=(15,18)) # Create a new figure
ax_rxx = newf.add_subplot(221)
ax_rxy = newf.add_subplot(222)
ax_temp = newf.add_subplot(223)
ax_tlog = newf.add_subplot(224)

# plot Rxx-temp dependence
for volt in Voltpick:
       index = 0
       Temp = [None]*len(fname)
       R = [None]*len(fname)
       for f in fname:
                index_pick = min(range(len(Dic_data[f].gVolt)), key=lambda k: abs(Dic_data[f].gVolt[k]-volt))
                #print('len(Dic_data[f].gVolt)={}'.format(len(Dic_data[f].gVolt)))
                #print('index_pick={}'.format(index_pick))
                Temp[index] = np.mean(Dic_data[f].tmp)
                R[index] = Dic_data[f].rxx[index_pick]
                index+=1
       # plot every single curve in outer loop
       ax_temp.scatter(Temp,R,label='Vg = %02.2fV'% volt)
       # plot in in Temp**power semilog
       ax_tlog.scatter([x**-0.5 for x in Temp],R,label='Vg = %02.2fV'% volt)
       ax_tlog.set_yscale('log')

ax_temp.set_xlabel('$Temp (K)$')
ax_temp.set_ylabel('Rxx ($\Omega$)')
ax_temp.legend()
ax_tlog.set_xlabel('$Temp^{-1/2} (K)$')
ax_tlog.set_ylabel('Rxx ($\Omega$)')
ax_tlog.legend()

# plot Rxx
fname = sorted(fname,key = lambda x: np.mean(Dic_data[x].tmp)) # sort fname by its temp value
for f in fname:
       ax_rxx.plot(Dic_data[f].gVolt,Dic_data[f].rxx,label='Temp = %02.1fK'% np.mean(Dic_data[f].tmp))
ax_rxx.set_xlabel('Gate Voltage (V)')
ax_rxx.set_ylabel('Rxx ($\Omega$)')
ax_rxx.legend()
for volt in Voltpick:
        ax_rxx.axvline(x=volt, linestyle=':', color='c')

# plot Rxy
for f in fname:
       ax_rxy.plot(Dic_data[f].gVolt,Dic_data[f].rxy,label='Temp = %02.1fK'% np.mean(Dic_data[f].tmp))
ax_rxy.set_xlabel('Gate Voltage (V)')
ax_rxy.set_ylabel('Rxy ($\Omega$)')
ax_rxy.legend()
for volt in Voltpick:
        ax_rxy.axvline(x = volt, linestyle = ':', color = 'c')
for i in range(1,5):
        ax_rxy.axhline(y = h/i/e**2,linestyle = ':',color = 'm')
plt.show()
newf.savefig('plot_metallike.pdf')
```
