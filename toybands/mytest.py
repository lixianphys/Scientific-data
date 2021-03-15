import numpy as np
from classes import Band, System

e0 = 1.6021766208e-19  # The elementary charge
h0 = 6.62607015e-34  # The Planck's constant
hbar = h0 / np.pi / 2
muB = -9.284764e-24
me = 9.109e-31  # The static mass of electron

bdA = Band(density=1e15, is_cond=True, is_dirac=True, gfactor=28, M=0.1, vf=1e6)
bdB = Band(density=3e15, is_cond=False, is_dirac=True, gfactor=28, M=0.1, vf=1e6)
bdC = Band(density=2e15, is_cond=True, is_dirac=False, gfactor=10, meff=0.1, spin=1)
bdD = Band(density=3e15, is_cond=False, is_dirac=False, gfactor=10, meff=0.2, spin=-1)
onebandsys = System(bdA)
twobandsys = System(bdA, bdB)
threebandsys = System(bdA, bdB, bdC)
fourbandsys = System(bdA, bdB, bdC, bdD)
# print(fourbandsys.mu([0.1, 0.2, 0.3], 10, 2, 0, 0.1 * e0))

print(onebandsys.dos_gen([0.1, 0.2, 0.3], 10, 2, 0, 0.1 * e0))
print(twobandsys.dos_gen([0.1, 0.2, 0.3], 10, 2, 0, 0.1 * e0))
print(threebandsys.bands)
print(fourbandsys.bands)
