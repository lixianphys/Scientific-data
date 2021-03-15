import numpy as np
from SciData.toybands.classes import Band, System
from SciData.toybands.functions import *
from SciData.physconst import *
import pdb

bdA = Band(density=1e15, is_cond=True, is_dirac=True, gfactor=28, M=0.1, vf=1e6)
bdB = Band(density=-3e15, is_cond=False, is_dirac=True, gfactor=28, M=0.1, vf=1e6)
bdC = Band(density=2e15, is_cond=True, is_dirac=False, gfactor=10, meff=0.1, spin=1)
bdD = Band(density=-3e15, is_cond=False, is_dirac=False, gfactor=10, meff=0.2, spin=-1)

e_list, B, Nmax, angle_in_deg, sigma = [0.1, 0.2, 0.3], 10, 2, 0, 0.1 * e0
E, B, sigma, angle_in_deg, e_lls  = 0.5,10,0.1,0,[0.5,0.7,0.8,0.9,0.5]
# print(e_density_of_state(E, B, sigma, angle_in_deg, e_lls))
# print(e_idos_gen(e_list, B, sigma, angle_in_deg, e_lls))

# print(bdA.__dict__)
# print(bdC.__dir__())


# for index,band in enumerate([bdA,bdC,bdB,bdD]):
#     print(index)
#     print([idos/1e15 for idos in band.cal_dos_b(e_list, B, Nmax, angle_in_deg, sigma)])




# onebandsys = System(bdA)
# twobandsys = System(bdA, bdB)
# threebandsys = System(bdA, bdB, bdC)
fourbandsys = System(bdA, bdB, bdC, bdD)

# print(onebandsys.dos_gen([0.1, 0.2, 0.3], 10, 2, 0, 0.1))
# print(twobandsys.dos_gen([0.1, 0.2, 0.3], 10, 2, 0, 0.1))
# print(threebandsys.bands)
a = fourbandsys.get_info()
print(fourbandsys.tot_density())
print(fourbandsys.mu(e_list, B, Nmax, angle_in_deg, sigma))
# fourbandsys.del_band()
print(a.iloc[3].dropna().to_dict())
