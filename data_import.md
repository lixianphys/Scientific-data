## import data from multiple files with names of a specific pattern
```
import os
import glob
import numpy as np
# path = input()
path = '//system01\QT\Lixian\data\QC 0492 0h\QC 0492 0H-difTemp\gatesweep\B-10T-Extn1mV_I_38_Rxx_BE_Rxy_BC_Gate4tom6V (02)'
os.chdir(path)

fname = sorted(glob.glob('TempControl*.dat')) # sort the filenames by their value
print('dictionary keys:')
Dic_data = {}

class sData:  # define a class type as a container
            def __init__(self, Name, Gatevoltage, Temp, Rxx,Rxy):
                self.name = Name
                self.gVolt = Gatevoltage
                self.tmp = Temp
                self.rxx = Rxx
                self.rxy = Rxy

for f in fname:
            read_data = np.loadtxt(f,skiprows=32,usecols=(0,3,7,-2))
            read_data = np.array(read_data)
            Dic_data[f] = sData(f,read_data[:,0],read_data[:,1],read_data[:,2],read_data[:,3])
            print(f)

```
