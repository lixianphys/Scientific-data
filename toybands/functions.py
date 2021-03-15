"""
All the functions
"""
import os
import sys
import pdb

import numpy as np
import pandas as pd

from functools import reduce
from scipy.integrate import quad

from SciData.physconst import *




def lldirac_gen(B, B_perp, N, is_cond, gfactor, M, vf):
    """Calculate the energy of Landau level in Dirac dispersion with a Zeeman term
    Arguments:
    B: Total magnetic field
    B_perp: The perpendicular component of B
    N: Landau index
    is_cond: True for conduction and False for valence band
    gfactor: g-factor
    vf: Fermi-velocity
    Return:
    Energy
    """
    if not isinstance(is_cond, bool):
        raise ValueError(f"your input is_cond ={is_cond} is not boolen")
    if N < 0 or not isinstance(N, int):
        raise ValueError(f"your input N = {N} should be an integer no less than zero")
    if gfactor < 0:
        raise ValueError(f"your input gfactor = {gfactor}<0")
    alpha = 1 if is_cond else -1
    return (
        alpha
        * (2 * e0 * hbar * vf ** 2 * B_perp * N + (gfactor * muB * B) ** 2 + M ** 2)
        ** 0.5
    )


def llconv_gen(B, B_perp, N, spin, gfactor, meff):
    """Calculate the energy of Landau level in conventional dispersion with a Zeeman term
    Arguments:
    B: Total magnetic field
    B_perp: The perpendicular component of B
    s: spin indicator
    meff: effective mass in units of the rest mass of electron me
    gfactor: g-factor
    Return:
    Energy
    """
    if not spin in [1, -1]:
        raise ValueError(f"your input spin ={spin} is neither 1 or -1")
    return (N + 0.5) * hbar * e0 * B_perp / me / meff + spin * gfactor * muB * B / 2


def _e_integral(fun, ymin, y_list, args):
    if not isinstance(y_list,list):
        raise ValueError("ylist is not a list")
    if not isinstance(args, tuple):
        raise ValueError("args must be a tuple")
    output = []
    pre_integral = 0
    yint = abs(y_list[0]-y_list[1])
    for index, y in enumerate(y_list):
        if y > ymin + yint and index == 0:
            result = quad(fun, ymin, y, args=args)[0]
        elif y > ymin + yint and index > 0:
            result = quad(fun, y - yint, y, args=args)[0]
        elif y <= ymin:
            result = 0
        else:
            result = quad(fun, ymin, y, args=args)[0]
        result = result + pre_integral
        pre_integral = result
        output.append(result)
    return output


def _h_integral(fun, ymax, y_list, args):
    if not isinstance(y_list,list):
        raise ValueError("ylist is not a list")
    if not isinstance(args, tuple):
        raise ValueError("args must be a tuple")
    output = []
    pre_integral = 0
    yint = abs(y_list[0]-y_list[1])
    for index, y in enumerate(sorted(y_list,reverse=True)):
        if y < ymax - yint and index == 0:
            result = quad(fun, y, ymax, args=args)[0]
        elif y < ymax - yint and index > 0:
            result = quad(fun, y, y + yint, args=args)[0]
        elif y >= ymax:
            result = 0
        else:
            result = quad(fun, y, ymax, args=args)[0]
        result = result + pre_integral
        pre_integral = result
        output.append(result)
    output.reverse()
    return output


def e_density_of_state(E, B, sigma, angle_in_deg, e_lls):
    """Calculate the density of state at a set of certain chemical potential and magnetic field for top/bottom surface states
    Arguments:
    E: position of chemical potential
    B: total magnetic field
    sigma: broadening of Landau level by assuming a Gaussian-shape distribution around the central energy
    angle: the angle of magnetic field with the normal of sample plane
    llenergy_top_surface: energy of Landau levels from top surface state
    llenergy_bottom_surface: energy of Landau levels from bottom surface state
    Return:
    Density of state from all the bands at (E,B)
    """
    # degeneracy of Landau levels at a certain field
    lldegeneracy = B * np.cos(angle_in_deg * np.pi / 180) * e0 / h0

    # both top and bottom surfaces right at the chemical potential
    return reduce(
        (lambda x, y: x + y),
        [
            lldegeneracy
            * np.exp(-0.5 * (E - e_ll) ** 2 / sigma ** 2)
            / sigma
            / (2 * np.pi) ** 0.5
            for e_ll in e_lls
        ],
    )


def h_density_of_state(E, B, sigma, angle_in_deg, h_lls):
    """Calculate the density of state at a set of certain chemical potential and magnetic field for top/bottom surface states
    Arguments:
    E: position of chemical potential
    B: total magnetic field
    sigma: broadening of Landau level by assuming a Gaussian-shape distribution around the central energy
    angle: the angle of magnetic field with the normal of sample plane
    llenergy_top_surface: energy of Landau levels from top surface state
    llenergy_bottom_surface: energy of Landau levels from bottom surface state
    Return:
    Density of state from all the bands at (E,B)
    """
    # degeneracy of Landau levels at a certain field
    lldegeneracy = B * np.cos(angle_in_deg * np.pi / 180) * e0 / h0

    # both top and bottom surfaces right at the chemical potential
    return reduce(
        (lambda x, y: x + y),
        [
            -lldegeneracy
            * np.exp(-0.5 * (E - h_ll) ** 2 / sigma ** 2)
            / sigma
            / (2 * np.pi) ** 0.5
            for h_ll in h_lls
        ],
    )


def e_idos_gen(e_list, B, sigma, angle_in_deg, e_lls):
    lowest_ll_eng = min(e_lls)
    output = _e_integral(
            e_density_of_state,
            lowest_ll_eng - 3 * sigma,
            e_list,
            (B, sigma, angle_in_deg, e_lls),
        )
    return output


def h_idos_gen(e_list, B, sigma, angle_in_deg, h_lls):
    highest_ll_eng = max(h_lls)
    output=_h_integral(
            h_density_of_state,
            highest_ll_eng + 3 * sigma,
            e_list,
            (B, sigma, angle_in_deg, h_lls),
        )
    return output


def add_list(a, b):

    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        output = []
        for aa, bb in zip(a, b):
            output.append(aa + bb)
        return output
    elif not isinstance(a, list) or not isinstance(b, list):
        raise TypeError("Input must be a list")
    elif len(a) != len(b):
        raise ValueError("Two lists must be of the same length")


