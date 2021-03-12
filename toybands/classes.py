'''
All the classes
'''
import sys
import os

import numpy as np
import pandas as pd



from toybands.functions import lldirac_gen, llconv_gen

class Band():
    def __init___(self,Ebb,is_cond,is_dirac,gfactor,*kwargs):
        self.Ebb = Ebb
        self.is_cond = is_cond
        self.is_dirac = is_dirac
        self.gfactor = gfactor
        for key, value in kwargs.items():
            if key == 'M' and is_dirac:
                self.M = value
            elif key == 'vf' and is_dirac:
                self.vf = value
            elif key == 'meff' and not is_dirac:
                self.meff = value
            elif key == 'spin' and not is_dirac:
                self.spin = value
            else:
                raise ValueError('your input is invalid, initialization of Band failed')
    def cal_energy(self,b_list,Nmax,angle_in_deg):
        if not isinstance(b_list,list):
            raise ValueError(f'b_list = {b_list} is not list')
        if not isinstance(Nmax,int):
            raise ValueError(f'Nmax={Nmax} is not integer')
        if self.is_dirac:
            return pd.DataFrame([{f'{N}':map(lambda x: lldirac_gen(x, x*np.cos(angle_in_deg), N, self.is_cond, self.gfactor, M, vf), blist)}])
        else:
            return pd.DataFrame([{f'{N}':map(lambda x: llconv_gen(x, x*np.cos(angle_in_deg), N, spin, gfactor, meff), blist)}])
    




