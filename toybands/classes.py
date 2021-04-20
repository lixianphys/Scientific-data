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
from toybands.config import DEFAULT_AUTONAME,DEFAULT_PATH


class Band:
    def __init__(self, density, is_cond, is_dirac, gfactor,vf, dparam, meff, spin):      
        self.density = abs(density)
        self.is_cond = is_cond
        self.is_dirac = is_dirac
        self.gfactor = gfactor
        self.vf = vf
        self.dparam = dparam
        self.meff = meff
        self.spin = spin
        self.active = True
        self.Ebb = den2en(abs(density),is_dirac,is_cond,vf,dparam,meff)
   
    def get(self, attr):
        if attr in ['density','is_cond','is_dirac','M','vf','meff','spin','Ebb']:
            return self.__dict__.get(attr)
        else:
            sys.stderr.write(f'{attr} is not an attribute for Band. Use density,is_cond,is_dirac,M,vf,meff,spin instead\n')
    def set_den(self,value):
        if not isinstance(value,(int,float)):
            sys.stderr.write(f'value need to be a number\n')
        self.density = value
        self.Ebb = den2en(self.density,self.is_dirac,self.is_cond,self.vf,self.dparam,self.meff)
        return None

    def disable(self):
        self.active = False 
    def enable(self):
        self.active = True
    def status(self):
        return self.active

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
                                self.vf,
                                self.dparam
                            )
                            for x in b_list
                        ]
                    }
                )
            return pd.DataFrame.from_dict(ll_dict)
        else:
            ll_dict = {}
            # spin-degenerate band with spin-up and spin-down branch
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
                    })
            return pd.DataFrame.from_dict(ll_dict)

    def cal_idos_b(self, e_list, B, Nmax, angle_in_deg, sigma, dirac_parity=True):
        if not isinstance(e_list, list):
            raise ValueError(f"e_list = {e_list} is not a list")
        if self.is_dirac and self.is_cond:
            e_lls = [
                self.Ebb
                + lldirac_gen(
                    B,
                    B * np.cos(angle_in_deg * np.pi / 180),
                    N,
                    self.is_cond,
                    self.gfactor,
                    self.vf,
                    self.dparam
                )
                for N in range(Nmax)
            ]
            return e_idos_gen(e_list, B, sigma, angle_in_deg, e_lls,compensate_on=dirac_parity)
        elif self.is_dirac and not self.is_cond:
            h_lls = [
                self.Ebb
                + lldirac_gen(
                    B,
                    B * np.cos(angle_in_deg * np.pi / 180),
                    N,
                    self.is_cond,
                    self.gfactor,
                    self.vf,
                    self.dparam
                )
                for N in range(Nmax)
            ]
            return h_idos_gen(e_list, B, sigma, angle_in_deg, h_lls,compensate_on=dirac_parity)
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
            return e_idos_gen(e_list, B, sigma, angle_in_deg, e_lls,compensate_on=False)
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
            return h_idos_gen(e_list, B, sigma, angle_in_deg, h_lls,compensate_on=False)

    def cal_dos(self, E, B, Nmax, angle_in_deg, sigma,dirac_parity=True):
        if self.is_dirac and self.is_cond:
            e_lls = [
                self.Ebb
                + lldirac_gen(
                    B,
                    B * np.cos(angle_in_deg * np.pi / 180),
                    N,
                    self.is_cond,
                    self.gfactor,
                    self.vf,
                    self.dparam
                )
                for N in range(Nmax)
            ]
            return e_density_of_state(E, B, sigma, angle_in_deg, e_lls,compensate_on=dirac_parity)
        elif self.is_dirac and not self.is_cond:
            h_lls = [
                self.Ebb
                + lldirac_gen(
                    B,
                    B * np.cos(angle_in_deg * np.pi / 180),
                    N,
                    self.is_cond,
                    self.gfactor,
                    self.vf,
                    self.dparam
                )
                for N in range(Nmax)
            ]
            return h_density_of_state(E, B, sigma, angle_in_deg, h_lls,compensate_on=dirac_parity)
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
            return e_density_of_state(E, B, sigma, angle_in_deg, e_lls,compensate_on=False)
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
            return h_density_of_state(E, B, sigma, angle_in_deg, h_lls,compensate_on=False)

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

    def get_info(self):
        band_info = []
        for band in self.bands:
            band_dict = band.__dict__
            band_info.append(band_dict)
        return pd.DataFrame(band_info)

    def set_all_density(self,den_list):
        if not len(self.bands) == len(den_list):
            num_active_band = len(self.bands)
            sys.stderr.write(f'The number of active bands {num_active_band} does not match up with the input densities {len(den_list)}\n')
        for band,den in zip(self.bands,den_list):
            band.set_den(abs(den))
    
    def get_band(self,band):
        if len(self.bands) == 0:
            sys.stdout.write('No bands in this system')
            return None
        if band == 'a':
            return [band for band in self.bands if band.status()]
        elif isinstance(band,int):
            if band < len(self.bands) and band >= 0:
                return self.bands[band]
            else:
                sys.stderr.write(f'#{band} band does not exist. You only have {len(self.bands)} bands\n')
        else:
            sys.stderr.write(f'Invalid {band} use integer as index or a to choose all\n') 

    def add_band(self, band):
        if not isinstance(band, Band):
            raise TypeError(f"{band} is not Band object\n")
        elif band in self.bands:
            raise ValueError(f"{band} is already a part of System\n")
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
        if self.get_band('a'):
            tot_den = 0
            for band in self.get_band('a'):
                if band.get('is_cond'):
                    tot_den+=band.get('density')
                else:
                    tot_den-=band.get('density')
            return tot_den
        else:
            return 0

    def dos_gen(self, e_list, B, Nmax, angle_in_deg, sigma_list,dirac_parity=True):
        if not isinstance(e_list, list):
            raise TypeError(f"{e_list} is not a list")
        elif not self.get_band('a'):
            raise ValueError(f"No active band found in the system")
        else:
            return reduce(
                (lambda x, y: add_list(x, y)),
                [
                    band.cal_idos_b(e_list, B, Nmax, angle_in_deg, sigma_list[idx],dirac_parity)
                    for idx,band in enumerate(self.get_band('a'))
                ],
            )

    def mu(self, e_list, B, Nmax, angle_in_deg, sigma_list,dirac_parity=True):
        return np.interp(
            x=self.tot_density(),
            xp=self.dos_gen(e_list, B, Nmax, angle_in_deg, sigma_list,dirac_parity),
            fp=e_list,
        )
    def databdl_write_csv(self,filename,bfrange,y_databdl,indicator,plotrange=None):
        if filename is None:
            filename = DEFAULT_AUTONAME
        elif len(filename.split('.'))>1:
            filename = filename.split('.')[-2]
        path = os.path.join(mkdir(filename),filename+'.csv')

        if indicator == 'enplot':
            columns = ['B','E','N','Band']
        elif indicator == 'denplot':
            columns = ['B','den','N','Band']
        elif indicator == 'simu':
            columns = ['B','den','N','Band','System([band density])']
        elif indicator == 'dos':
            columns = ['B','dos_at_mu','Band']
        elif indicator == 'dosm':
            columns = ['B','dos_at_mu','System([band density])']
        else:
            raise ValueError(f'Invalid indicator {indicator}, use enplot, denplot, simu, dos, dosm instead')
        
        df = pd.DataFrame(columns=columns)
        
        if indicator in ['enplot','denplot']:
            for n_band, y_data in enumerate(y_databdl):
                for n_ll, y in enumerate(y_data):
                    df_toappend = pd.DataFrame(np.transpose(np.array([bfrange,y,[n_ll]*len(y),[n_band]*len(y)])),columns=columns)
                    df = df.append(df_toappend,ignore_index=True)
            df.to_csv(path,mode='w',index=False)
        elif indicator == 'simu':
            if self.get_band('a'):
                _den = list(map(lambda x:x.density,self.get_band('a')))
            else:
                sys.stderr.write('There is no active band in this system\n')
                exit()
            str_den = ' '.join(["{:e}".format(den) for den in _den])
            for n_band, y_data in enumerate(y_databdl):
                for n_ll, y in enumerate(y_data):
                    bf_p,y_p = bfrange,y
                    if plotrange is not None:
                        bf_p = extract_list(bfrange,[yy>plotrange[0] and yy<plotrange[1] for yy in y])
                        y_p = extract_list(y,[yy>plotrange[0] and yy<plotrange[1] for yy in y])
                    df_toappend = pd.DataFrame(np.transpose(np.array([bf_p,y_p,[n_ll]*len(y_p),[n_band]*len(y_p),[str_den]*len(y_p)])),columns=columns)
                    df = df.append(df_toappend,ignore_index=True)
            if os.path.isfile(path):
                df.to_csv(path,mode='a',index=False,header=False)
            else:
                df.to_csv(path,mode='a',index=False)

        elif indicator == 'dos':
            for n_band, y in enumerate(y_databdl):
                df_toappend = pd.DataFrame(np.transpose(np.array([bfrange,y,[n_band]*len(y)])),columns=columns)
                df = df.append(df_toappend,ignore_index=True)
            df.to_csv(path, mode='w', index=False)

        elif indicator == 'dosm':
            if self.get_band('a'):
                _den = list(map(lambda x:x.density,self.get_band('a')))
            else:
                sys.stderr.write('There is no active band in this system\n')
                exit()
            str_den = ' '.join(["{:e}".format(den) for den in _den])

            df_toappend = pd.DataFrame(np.transpose(np.array([bfrange,y_databdl,[str_den]*len(y_databdl)])),columns=columns)
            df = df.append(df_toappend,ignore_index=True)
            if os.path.isfile(path):
                df.to_csv(path,mode='a', index=False, header=False)
            else:
                df.to_csv(path,mode='a',index=False)
    
    def ydata_gen(self, nmax, angle_in_deg, bfrange, enrange, indicator,scond, sval, dirac_parity=True):
        if indicator == 'enplot':
            return [
            [
                band.cal_energy(bfrange, nmax, angle_in_deg)[
                    f"#{N}"
                ].tolist()
                for N in range(nmax)
            ]
            for band in self.get_band('a')]
        elif indicator == 'denplot':
            IDOS = [
            self.dos_gen(
                enrange, B, nmax, angle_in_deg, [scond if band.get('is_cond') else sval for band in self.get_band('a')],dirac_parity)
            for B in bfrange]
            return [
            [
                [
                    np.interp(x=x, xp=enrange, fp=IDOS[index])
                    for index, x in enumerate(
                        band.cal_energy(bfrange, nmax, angle_in_deg)[
                            f"#{N}"
                        ].tolist()
                    )
                ]
                for N in range(nmax)
            ]
            for band in self.get_band('a')]
        elif indicator == 'dos':
            mus = [self.mu(enrange, B, nmax, angle_in_deg, [scond if band.get('is_cond') else sval for band in self.get_band('a')],dirac_parity) for B in bfrange]
            return [[band.cal_dos(mu_at_B, B, nmax, angle_in_deg, scond if band.get('is_cond') else sval,dirac_parity) for B, mu_at_B in zip(bfrange,mus)] for band in self.get_band('a')]
        else:
            raise ValueError(f'Invalid indicator {indicator}, use enplot, denplot, dos instead\n')