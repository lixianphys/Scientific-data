"""
All the classes
"""
import sys
import os
import warnings

import numpy as np
import pandas as pd
from functools import reduce

from functions import *

e0 = 1.6021766208e-19  # The elementary charge
h0 = 6.62607015e-34  # The Planck's constant
hbar = h0 / np.pi / 2
muB = -9.284764e-24
me = 9.109e-31  # The static mass of electron


class Band:
    def __init__(self, density, is_cond, is_dirac, gfactor, **kwargs):
        self.density = density
        self.is_cond = is_cond
        self.is_dirac = is_dirac
        self.gfactor = gfactor

        for key, value in kwargs.items():
            if key == "M" and is_dirac:
                self.M = value
            elif key == "vf" and is_dirac:
                self.vf = value
            elif key == "meff" and not is_dirac:
                self.meff = value
            elif key == "spin" and not is_dirac:
                self.spin = value
            elif (key == "M" or key == "vf") and not is_dirac:
                raise ValueError(
                    "Conventional bands do not have M and vf as parameters, use meff and spin"
                )
            elif (key == "meff" or key == "spin") and is_dirac:
                raise ValueError(
                    "Dirac bands do not have meff and spin as parameters, use vf and M"
                )
            elif key not in ["M", "vf", "meff", "spin"]:
                raise ValueError(f"Invalid argument {key}")

        if is_dirac:
            self.Ebb = -hbar * self.vf * (4 * np.pi * self.density) ** 0.5
        else:
            self.Ebb = -(hbar ** 2) * self.density * np.pi / (self.meff * me) / 2

    def cal_energy(self, b_list, Nmax, angle_in_deg):
        if not isinstance(b_list, list):
            raise TypeError(f"b_list = {b_list} is not list")
        if not isinstance(Nmax, int):
            raise TypeError(f"Nmax={Nmax} is not integer")
        if self.is_dirac:
            return pd.DataFrame(
                [
                    {
                        f"{N}": map(
                            lambda x: self.Ebb
                            + lldirac_gen(
                                x,
                                x * np.cos(angle_in_deg),
                                N,
                                self.is_cond,
                                self.gfactor,
                                self.M,
                                self.vf,
                            ),
                            b_list,
                        )
                    }
                    for N in range(Nmax)
                ]
            )
        else:
            return pd.DataFrame(
                [
                    {
                        f"{N}": map(
                            lambda x: self.Ebb
                            + llconv_gen(
                                x,
                                x * np.cos(angle_in_deg),
                                N,
                                self.spin,
                                self.gfactor,
                                self.meff,
                            ),
                            b_list,
                        )
                    }
                    for N in range(Nmax)
                ]
            )

    def cal_dos_b(self, e_list, B, Nmax, angle_in_deg, sigma):
        if not isinstance(e_list, list):
            raise ValueError(f"e_list = {e_list} is not a list")
        else:
            if self.is_dirac and self.is_cond:
                e_lls = [
                    self.Ebb
                    + lldirac_gen(
                        B,
                        B * np.cos(angle_in_deg),
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
                        B * np.cos(angle_in_deg),
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
                        B * np.cos(angle_in_deg),
                        N,
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
                        B * np.cos(angle_in_deg),
                        N,
                        self.spin,
                        self.gfactor,
                        self.meff,
                    )
                    for N in range(Nmax)
                ]
                return h_idos_gen(e_list, B, sigma, angle_in_deg, h_lls)


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

    def get_band(self):
        return self.bands

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
            return reduce((lambda x, y: x + y), [band.den for band in self.bands])
        else:
            return 0

    def dos_gen(self, e_list, B, Nmax, angle_in_deg, sigma):
        if not isinstance(e_list, list):
            raise ValueError(f"{e_list} is not a list")
        elif not self.bands:
            raise ValueError(f"No band added")
        else:
            # print(
            #     [
            #         band
            #         for band in self.bands
            #     ]
            # )
            return reduce(
                (lambda x, y: add_list(x, y)),
                [
                    band.cal_dos_b(e_list, B, Nmax, angle_in_deg, sigma)
                    for band in self.bands
                ],
            )

    def mu(self, e_list, B, Nmax, angle_in_deg, sigma):

        return np.interp(
            x=self.tot_density,
            xp=self.dos_gen(e_list, B, Nmax, angle_in_deg, sigma),
            fp=e_list,
        )
