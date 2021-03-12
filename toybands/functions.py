import os
import sys

import numpy as np
import pandas as pd



def lldirac_gen(B, B_perp, N, is_cond, gfactor, M, vf):

    """ Calculate the energy of Landau level in Dirac dispersion with a Zeeman term
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
    if not isinstance(alpha,bool):
        raise ValueError(f'your input alpha ={alpha} is not boolen')
    if N<0 or not isinstance(N,int):
        raise ValueError(f'your input N = {N} should be an integer no less than zero')
    if gfactor<0:
        raise ValueError(f'your input gfactor = {gfactor}<0')
    return alpha*(2 * e0 * hbar * vf ** 2 * B_perp * N + (gfactor * muB * B) ** 2+M**2) ** 0.5

def llconv_gen(B, B_perp, N, spin, gfactor, meff):
    """ Calculate the energy of Landau level in conventional dispersion with a Zeeman term
    Arguments:
    B: Total magnetic field
    B_perp: The perpendicular component of B
    s: spin indicator
    meff: effective mass in units of the rest mass of electron me
    gfactor: g-factor
    Return:
    Energy
    """
    if not spin in [1,-1]:
        raise ValueError(f'your input spin ={spin} is neither 1 or -1')
    return (N + 0.5) * hbar * e0 * B_perp / me / meff + spin*gfactor * muB * B / 2