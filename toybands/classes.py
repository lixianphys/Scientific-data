"""
All the classes
"""
import sys
import os
import warnings

import numpy as np
import pandas as pd
from functools import reduce

from physconst import *
from toybands.functions import *
from toybands.config import *


class Band:
    def __init__(self, density, is_cond, is_dirac, gfactor, M, vf, meff, spin):
        self.density = abs(density)
        self.is_cond = is_cond
        self.is_dirac = is_dirac
        self.gfactor = gfactor
        self.M = M if M is not None else None
        self.vf = vf
        self.meff = meff if meff is not None else None
        self.spin = spin

        if is_dirac and is_cond:
            self.Ebb = -hbar * self.vf * (4 * np.pi * self.density) ** 0.5
        elif is_dirac and not is_cond:
            self.Ebb = hbar * self.vf * (4 * np.pi * self.density) ** 0.5
        elif not is_dirac and is_cond:
            self.Ebb = -(hbar ** 2) * self.density/ np.pi / self.meff/me / 2
        elif not is_dirac and not is_cond:
            self.Ebb = (hbar ** 2) * self.density / np.pi / self.meff/me / 2

    def cal_energy(self, b_list, Nmax, angle_in_deg):
        if not isinstance(b_list, list):
            raise TypeError(f"b_list = {b_list} is not list")
        if not isinstance(Nmax, int):
            raise TypeError(f"Nmax={Nmax} is not integer")
        if self.is_dirac:
            ll_dict = {}
            for N in range(Nmax):
                ll_dict.update(
                    {
                        f"#{N}": [
                            self.Ebb
                            + lldirac_gen(
                                x,
                                x * np.cos(angle_in_deg * np.pi / 180),
                                N,
                                self.is_cond,
                                self.gfactor,
                                self.M,
                                self.vf,
                            )
                            for x in b_list
                        ]
                    }
                )
            return pd.DataFrame.from_dict(ll_dict)
        else:
            ll_dict = {}
            for N in range(Nmax):
                ll_dict.update(
                    {
                        f"#{N}": [
                            self.Ebb
                            + llconv_gen(
                                x,
                                x * np.cos(angle_in_deg * np.pi / 180),
                                N,
                                self.is_cond,
                                self.spin,
                                self.gfactor,
                                self.meff,
                            )
                            for x in b_list
                        ]
                    }
                )
            return pd.DataFrame.from_dict(ll_dict)

    def cal_dos_b(self, e_list, B, Nmax, angle_in_deg, sigma):
        if not isinstance(e_list, list):
            raise ValueError(f"e_list = {e_list} is not a list")
        else:
            if self.is_dirac and self.is_cond:
                e_lls = [
                    self.Ebb
                    + lldirac_gen(
                        B,
                        B * np.cos(angle_in_deg * np.pi / 180),
                        N,
                        self.is_cond,
                        self.gfactor,
                        self.M,
                        self.vf,
                    )
                    for N in range(Nmax)
                ]
                return e_idos_gen(e_list, B, sigma, angle_in_deg, e_lls)
            elif self.is_dirac and not self.is_cond:
                h_lls = [
                    self.Ebb
                    + lldirac_gen(
                        B,
                        B * np.cos(angle_in_deg * np.pi / 180),
                        N,
                        self.is_cond,
                        self.gfactor,
                        self.M,
                        self.vf,
                    )
                    for N in range(Nmax)
                ]
                return h_idos_gen(e_list, B, sigma, angle_in_deg, h_lls)
            elif not self.is_dirac and self.is_cond:
                e_lls = [
                    self.Ebb
                    + llconv_gen(
                        B,
                        B * np.cos(angle_in_deg * np.pi / 180),
                        N,
                        self.is_cond,
                        self.spin,
                        self.gfactor,
                        self.meff,
                    )
                    for N in range(Nmax)
                ]
                return e_idos_gen(e_list, B, sigma, angle_in_deg, e_lls)
            elif not self.is_dirac and not self.is_cond:
                h_lls = [
                    self.Ebb
                    + llconv_gen(
                        B,
                        B * np.cos(angle_in_deg * np.pi / 180),
                        N,
                        self.is_cond,
                        self.spin,
                        self.gfactor,
                        self.meff,
                    )
                    for N in range(Nmax)
                ]
                return h_idos_gen(e_list, B, sigma, angle_in_deg, h_lls)

    def print_band(self):
        band_dict = self.__dict__
        print("---------------------------")
        print(self)
        [print(f"{key}={band_dict.get(key)}") for key in band_dict.keys()]
        print("---------------------------")

    def get_info(self):
        return self.__dict__


class System:
    def __init__(self, *args):
        print("You are building a multi-band toy model with toybands")
        self.bands = []
        if len(args) > 0:
            for arg in args:
                if isinstance(arg, Band):
                    self.bands.append(arg)
                else:
                    raise ValueError(f"{arg} is not an instance of Band")
        else:
            warnings.warn(f"Initialization of an empty System")

    def get_info(self):
        band_info = []
        for band in self.bands:
            band_dict = band.__dict__
            band_info.append(band_dict)
        return pd.DataFrame(band_info)

    def set_all_density(self,den_list):
        if not len(self.bands) == len(den_list):
            sys.stderr.write(f'The number of bands {len(self.bands)} does not match up with the input densities {len(den_list)}')
        for band,den in zip(self.bands,den_list):
            band.density = abs(den)
            if band.is_dirac and band.is_cond:
                band.Ebb = -hbar * band.vf * (4 * np.pi * band.density) ** 0.5
            elif band.is_dirac and not band.is_cond:
                band.Ebb = hbar * band.vf * (4 * np.pi * band.density) ** 0.5
            elif not band.is_dirac and band.is_cond:
                band.Ebb = -(hbar ** 2) * band.density/ np.pi / band.meff/me / 2
            elif not band.is_dirac and not band.is_cond:
                band.Ebb = (hbar ** 2) * band.density / np.pi / band.meff/me / 2

    def add_band(self, band):
        if not isinstance(band, Band):
            raise TypeError(f"{band} is not Band object")
        elif band in self.bands:
            raise ValueError(f"{band} is already a part of System")
        else:
            self.bands.append(band)

    def del_band(self, band):
        if not isinstance(band, Band):
            raise TypeError(f"{band} is not Band object")
        elif band not in self.bands:
            raise ValueError(f"{band} is not a part of System")
        else:
            self.bands.remove(band)

    def tot_density(self):
        if self.bands:
            bands_den = []
            for band in self.bands:
                if band.is_cond:
                    bands_den.append(band.density)
                else:
                    bands_den.append(-band.density)
            return sum(bands_den)
        else:
            return 0

    def dos_gen(self, e_list, B, Nmax, angle_in_deg, sigma):
        if not isinstance(e_list, list):
            raise TypeError(f"{e_list} is not a list")
        elif not self.bands:
            raise ValueError(f"No band added into the system")
        else:
            return reduce(
                (lambda x, y: add_list(x, y)),
                [
                    band.cal_dos_b(e_list, B, Nmax, angle_in_deg, sigma)
                    for band in self.bands
                ],
            )

    def mu(self, e_list, B, Nmax, angle_in_deg, sigma):
        return np.interp(
            x=self.tot_density(),
            xp=self.dos_gen(e_list, B, Nmax, angle_in_deg, sigma),
            fp=e_list,
        )
    def databdl_write_csv(self,filename,bfrange,y_databdl,indicator,plotrange=None):
        if len(filename.split('.'))>1:
            filename = filename.split('.')[-2]
        path = os.path.join(DEFAULT_PATH,filename+'.csv')
        if not any(isinstance(el, list) for el in y_databdl):
            raise TypeError(f'y_databdl is not nested list')
        if indicator == 'enplot':
            columns = ['B','E','N','Band']
        elif indicator == 'denplot':
            columns = ['B','den','N','Band']
        elif indicator == 'simu':
            columns = ['B','den','N','Band','system']
        else:
            sys.stderr.write(f'Invalid indicator {indicator}, use enplot, denplot or simu instead')
            exit()
        df = pd.DataFrame(columns=columns)
        if indicator in ['enplot','denplot']:
            for n_band, y_data in enumerate(y_databdl):
                for n_ll, y in enumerate(y_data):
                    df_toappend = pd.DataFrame(np.transpose(np.array([bfrange,y,[n_ll]*len(y),[n_band]*len(y)])),columns=columns)
                    df = df.append(df_toappend,ignore_index=True)
            df.to_csv(path,mode='w',index=False)
        elif indicator == 'simu':
            if self.bands:
                _den = list(map(lambda x:x.density,self.bands))
            else:
                sys.stderr.write('There is no band in this system')
                exit()
            for n_band, y_data in enumerate(y_databdl):
                for n_ll, y in enumerate(y_data):
                    bf_p,y_p = bfrange,y
                    str_den = ' '.join(["{:e}".format(den) for den in _den])
                    if plotrange is not None:
                        bf_p = extract_list(bfrange,[yy>plotrange[0] and yy<plotrange[1] for yy in y])
                        y_p = extract_list(y,[yy>plotrange[0] and yy<plotrange[1] for yy in y])
                    df_toappend = pd.DataFrame(np.transpose(np.array([bf_p,y_p,[n_ll]*len(y_p),[n_band]*len(y_p),[str_den]*len(y_p)])),columns=columns)
                    df = df.append(df_toappend,ignore_index=True)
            if os.path.isfile(path):
                df.to_csv(path,mode='a',index=False,header=False)
            else:
                df.to_csv(path,mode='a',index=False)
