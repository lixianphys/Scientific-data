import os
import glob
import numpy as np
# path = input()
path = 'YOUR PATH HERE'
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

