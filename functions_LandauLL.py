# Standard library imports
import os

# Third party imports
import matplotlib.pyplot as plt
import numpy as np

# Local application import
from physconst import *



def LL_energy(B_tot, B_perp, N, gfactor, vf):  # N = 0,1,2,3...
    return (2 * e0 * hbar * vf ** 2 * B_perp * N + (gfactor * muB * B_tot) ** 2) ** 0.5


def LL_fillfactor(fermi_energy, LLenergy: list):
    index = 0
    while fermi_energy > LLenergy[index]:
        if index == len(LLenergy) - 1:
            break
        else:
            index = index + 1
    return index - 0.5


def density2energy(density, vf):
    return hbar * vf * (4 * np.pi * density) ** 0.5


def llenergy_dirac(B, B_perp, N, vf=1e6, gfactor=28):
    if N > 0:
        return (2 * e0 * hbar * vf ** 2 * B_perp * N + (gfactor * muB * B) ** 2) ** 0.5
    elif N == 0:
        return -gfactor * muB * B
    else:
        return -(2 * e0 * hbar * vf ** 2 * B_perp * (-N) + (gfactor * muB * B) ** 2) ** 0.5


def llenergy_conv(B, B_perp, N, s, meff, gfactor=6):
    if s == 1:
        return (N + 0.5) * hbar * e0 * B_perp / me / meff + gfactor * muB * B / 2
    else:
        return (N + 0.5) * hbar * e0 * B_perp / me / meff - gfactor * muB * B / 2


def llenergy_generator(Ets, Ebs, Evp, B, angle, meff, Nmax=30, vf=1e6, gfactor=28):
    llenergy_top_surface = [Ets + llenergy_dirac(B, B * np.cos(angle * np.pi / 180), N, vf, gfactor) for N in
                            range(Nmax - 1)]
    llenergy_bottom_surface = [Ebs + llenergy_dirac(B, B * np.cos(angle * np.pi / 180), N, vf, gfactor) for N in
                               range(Nmax - 1)]
    llenergy_vps_up = [Evp + llenergy_conv(B, B * np.cos(angle * np.pi / 180), N, 1, meff) for N in range(Nmax - 1)]
    llenergy_vps_down = [Evp + llenergy_conv(B, B * np.cos(angle * np.pi / 180), N, -1, meff) for N in range(Nmax - 1)]

    return llenergy_top_surface, llenergy_bottom_surface, llenergy_vps_up, llenergy_vps_down


# def Gaussian_broadening(energy,sigma,energy_center,lldegeneracy):
#     return lldegeneracy*np.exp(-0.5*(energy-energy_center)**2/sigma**2)/sigma/(2*np.pi)**0.5

def electron_density_of_state(E, B, sigma, angle, llenergy_top_surface, llenergy_bottom_surface):
    electron_density_of_state = []

    lldegeneracy = B * np.cos(angle * np.pi / 180) * e0 / h0

    for ll in llenergy_top_surface + llenergy_bottom_surface:
        electron_density_of_state.append(
            lldegeneracy * np.exp(-0.5 * (E - ll) ** 2 / sigma ** 2) / sigma / (2 * np.pi) ** 0.5)

    if all([llenergy_top_surface, llenergy_bottom_surface]):
        compensate = 0.5 * lldegeneracy * np.exp(-0.5 * (E - min(llenergy_top_surface)) ** 2 / sigma ** 2) / sigma / (
                    2 * np.pi) ** 0.5 + 0.5 * lldegeneracy * np.exp(
            -0.5 * (E - min(llenergy_bottom_surface)) ** 2 / sigma ** 2) / sigma / (
                                 2 * np.pi) ** 0.5  # DOS from 0LL should be half of other LLs.
    elif llenergy_top_surface:
        compensate = 0.5 * lldegeneracy * np.exp(-0.5 * (E - min(llenergy_top_surface)) ** 2 / sigma ** 2) / sigma / (
                    2 * np.pi) ** 0.5
    elif llenergy_bottom_surface:
        compensate = 0.5 * lldegeneracy * np.exp(
            -0.5 * (E - min(llenergy_bottom_surface)) ** 2 / sigma ** 2) / sigma / (2 * np.pi) ** 0.5
    else:
        print('Laudau levels can not be both empty!')
    return sum(electron_density_of_state) - compensate


def hole_density_of_state(E, B, sigma, angle, llenergy_vps_up, llenergy_vps_down):
    hole_density_of_state = []

    lldegeneracy = B * np.cos(angle * np.pi / 180) * e0 / h0

    for ll in llenergy_vps_up + llenergy_vps_down:
        hole_density_of_state.append(
            -lldegeneracy * np.exp(-0.5 * (E - ll) ** 2 / (3 * sigma) ** 2) / (3 * sigma) / (2 * np.pi) ** 0.5)

    return sum(hole_density_of_state)


def Integral_electron_DOS(E, B, sigma, angle, llenergy_top_surface, llenergy_bottom_surface):
    lowest_energy = min(llenergy_top_surface + llenergy_bottom_surface)
    # it is crucial to integrate from lowest_energy-3*sigma to take into account the broadening effect.
    if E > lowest_energy:
        result, _ = quad(electron_density_of_state, lowest_energy - sigma, E,
                         args=(B, sigma, angle, llenergy_top_surface, llenergy_bottom_surface))
        return result
    else:
        return 0


def fastIntegral_electron_DOS(Energy, B, sigma, angle, llenergy_top_surface, llenergy_bottom_surface):
    lowest_energy = min(llenergy_top_surface + llenergy_bottom_surface)
    output = []
    last_result = 0
    energy_interval = Energy[1] - Energy[0]
    for index, energy in enumerate(Energy):
        if all([energy > lowest_energy - 3 * sigma + energy_interval, index == 0]):
            result, _ = quad(electron_density_of_state, lowest_energy - 3 * sigma, energy,
                             args=(B, sigma, angle, llenergy_top_surface, llenergy_bottom_surface))
        elif all([energy > lowest_energy - 3 * sigma + energy_interval, index > 0]):
            result, _ = quad(electron_density_of_state, energy - energy_interval, energy,
                             args=(B, sigma, angle, llenergy_top_surface, llenergy_bottom_surface))
        elif energy <= lowest_energy - 3 * sigma:
            result = 0
        else:
            result, _ = quad(electron_density_of_state, lowest_energy - 3 * sigma, energy,
                             args=(B, sigma, angle, llenergy_top_surface, llenergy_bottom_surface))

        result = last_result + result
        output.append(result)
        last_result = result
    return output


def Integral_hole_DOS(E, B, sigma, angle, llenergy_vps_up, llenergy_vps_down):
    highest_energy = max(llenergy_vps_up + llenergy_vps_down)
    # it is crucial to integrate to highest_energy+3*sigma to take into account the broadening effect.
    if E < highest_energy:
        result, _ = quad(hole_density_of_state, highest_energy + 3 * sigma, E,
                         args=(B, sigma, angle, llenergy_vps_up, llenergy_vps_down))
        return result
    else:
        return 0


def fastIntegral_hole_DOS(Energy, B, sigma, angle, llenergy_vps_up, llenergy_vps_down):
    highest_energy = max(llenergy_vps_up + llenergy_vps_down)

    output = []
    last_result = 0
    energy_interval = Energy[1] - Energy[0]
    Ennergy_r = Energy.tolist()
    Ennergy_r.reverse()
    for index, energy in enumerate(Ennergy_r):
        if all([energy < highest_energy + 3 * sigma - energy_interval, index == 0]):
            result, _ = quad(hole_density_of_state, energy, highest_energy + 3 * sigma,
                             args=(B, sigma, angle, llenergy_vps_up, llenergy_vps_down))
        elif all([energy < highest_energy + 3 * sigma - energy_interval, index > 0]):
            result, _ = quad(hole_density_of_state, energy, energy + energy_interval,
                             args=(B, sigma, angle, llenergy_vps_up, llenergy_vps_down))
        elif energy >= highest_energy + 3 * sigma:
            result = 0
        else:
            result, _ = quad(hole_density_of_state, energy, highest_energy + 3 * sigma,
                             args=(B, sigma, angle, llenergy_vps_up, llenergy_vps_down))

        result = last_result + result
        output.append(result)
        last_result = result
    output.reverse()
    return output


def find_energy_bydensity(target_density, B, IDOS_B, energy):
    return np.interp(x=target_density, xp=IDOS_B, fp=energy)




# define a container for passing on frequently used global arguments
class TBLLsimu():
    def __init__(self, vf, gfactor, sigma, meff):
        self.vf = vf
        self.gfactor = gfactor
        self.sigma = sigma
        self.meff = meff

    def __repr__(self):
        return f'TBLLsimu(vf = {self.vf},gfactor = {self.gfactor},sigma = {self.sigma},meff = {self.meff}))'

    def get_ll_en(self, angle, Brange, Nmax, den_top, den_bot, den_vps=None, threeband=False):

        vf = self.vf
        sigma = self.sigma
        gfactor = self.gfactor
        meff = self.meff

        Ets = -hbar * vf * (4 * np.pi * den_top) ** 0.5
        Ebs = -hbar * vf * (4 * np.pi * den_bot) ** 0.5

        LLenergy_top_surface = []
        LLenergy_bottom_surface = []

        if threeband:
            Evp = -hbar ** 2 * den_vps * np.pi / (meff * me) / 2
            LLenergy_vps_up = []
            LLenergy_vps_down = []

            for B in Brange:
                llenergy_top_surface, llenergy_bottom_surface, llenergy_vps_up, llenergy_vps_down = llenergy_generator(
                    Ets, Ebs, Evp, B, angle, meff, Nmax, vf, gfactor)
                LLenergy_top_surface.append(llenergy_top_surface)
                LLenergy_bottom_surface.append(llenergy_bottom_surface)
                LLenergy_vps_up.append(llenergy_vps_up)
                LLenergy_vps_down.append(llenergy_vps_down)
            return LLenergy_top_surface, LLenergy_bottom_surface, LLenergy_vps_up, LLenergy_vps_down

        else:

            for B in Brange:
                llenergy_top_surface, llenergy_bottom_surface, _, _ = llenergy_generator(Ets, Ebs, 0, B, angle, meff,
                                                                                         Nmax, vf, gfactor)
                LLenergy_top_surface.append(llenergy_top_surface)
                LLenergy_bottom_surface.append(llenergy_bottom_surface)

            return LLenergy_top_surface, LLenergy_bottom_surface

    def get_ll_den(self, angle, Brange, Erange, Nmax, den_top, den_bot, den_vps=None, threeband=False):

        vf = self.vf
        sigma = self.sigma
        gfactor = self.gfactor
        meff = self.meff

        N_top_surface = []
        N_bottom_surface = []
        if threeband:

            LLenergy_top_surface, LLenergy_bottom_surface, LLenergy_vps_up, LLenergy_vps_down = self.get_ll_en(angle,
                                                                                                               Brange,
                                                                                                               Nmax,
                                                                                                               den_top,
                                                                                                               den_bot,
                                                                                                               den_vps,
                                                                                                               threeband)
            IDOS = self.IDOS_generator(angle, Brange, Erange, LLenergy_top_surface, LLenergy_bottom_surface,
                                       LLenergy_vps_up, LLenergy_vps_down)
            N_vps_up = []
            N_vps_down = []

            for LLenergy, LLn in zip(
                    [LLenergy_top_surface, LLenergy_bottom_surface, LLenergy_vps_up, LLenergy_vps_down],
                    [N_top_surface, N_bottom_surface, N_vps_up, N_vps_down]):
                for index, ll in enumerate(LLenergy):
                    IDOS_B = IDOS[index]
                    LLn.append([np.interp(x=x, xp=Erange, fp=IDOS_B) for x in ll])

            return N_top_surface, N_bottom_surface, N_vps_up, N_vps_down

        else:

            LLenergy_top_surface, LLenergy_bottom_surface = self.get_ll_en(angle, Brange, Nmax, den_top, den_bot,
                                                                           den_vps, threeband)
            IDOS = self.IDOS_generator(angle, Brange, Erange, LLenergy_top_surface, LLenergy_bottom_surface)
            for LLenergy, LLn in zip([LLenergy_top_surface, LLenergy_bottom_surface],
                                     [N_top_surface, N_bottom_surface]):
                for index, ll in enumerate(LLenergy):
                    IDOS_B = IDOS[index]
                    LLn.append([np.interp(x=x, xp=Erange, fp=IDOS_B) for x in ll])

            return N_top_surface, N_bottom_surface

    def plot_ll_en(self, angle, Brange, Erange, Nmax, den_top, den_bot, den_vps=None, threeband=False):

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)

        if threeband:

            LL_ts, LL_bs, LL_vpsup, LL_vpsdown = self.get_ll_en(angle, Brange, Nmax, den_top, den_bot, den_vps,
                                                                threeband)
            IDOS = self.IDOS_generator(angle, Brange, Erange, LL_ts, LL_bs, LL_vpsup, LL_vpsdown)
            ax.plot(Brange,
                    [find_energy_bydensity(den_top + den_bot - den_vps, B, IDOS_B, Erange) * 1e3 / e0 for B, IDOS_B in
                     zip(Brange, IDOS)], linewidth=1, color='k')

            for ll_ts in np.transpose(LL_ts):
                ax.plot(Brange, ll_ts * 1e3 / e0, 'r-')
            for ll_bs in np.transpose(LL_bs):
                ax.plot(Brange, ll_bs * 1e3 / e0, 'b-')
            for ll_vpsup in np.transpose(LL_vpsup):
                ax.plot(Brange, ll_vpsup * 1e3 / e0, 'k-')
            for ll_vpsdown in np.transpose(LL_vpsdown):
                ax.plot(Brange, ll_vpsdown * 1e3 / e0, 'y-')

        else:
            LL_ts, LL_bs = self.get_ll_en(angle, Brange, Nmax, den_top, den_bot, den_vps, threeband)
            IDOS = self.IDOS_generator(angle, Brange, Erange, LL_ts, LL_bs)
            #                     [ax.plot(Brange,[find_energy_bydensity(den,B,IDOS_B,Erange)*1e3/e0 for B,IDOS_B in zip(Brange,IDOS)],linewidth=1,color='g',linestyle='--') for den in np.linspace(1e15,1e16,10)]
            ax.plot(Brange, [find_energy_bydensity(den_top + den_bot, B, IDOS_B, Erange) * 1e3 / e0 for B, IDOS_B in
                             zip(Brange, IDOS)], linewidth=1, color='k')
            for ll_ts in np.transpose(LL_ts):
                ax.plot(Brange, ll_ts * 1e3 / e0, 'r-')
            for ll_bs in np.transpose(LL_bs):
                ax.plot(Brange, ll_bs * 1e3 / e0, 'b-')

        ax.set_ylim([-10, 30])
        ax.set_xlabel('B (T)')
        ax.set_ylabel('Energy (meV)')
        return fig, ax

    def plot_ll_den(self, angle, Brange, Erange, Nmax, den_top, den_bot, den_vps=None, threeband=False):

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        if threeband:
            N_ts, N_bs, N_vpsup, N_vpsdown = self.get_ll_den(angle, Brange, Erange, Nmax, den_top, den_bot, den_vps,
                                                             threeband)
            ax.axhline(y=(den_top + den_bot - den_vps) / 1e15, linewidth=1)
            for n_ts in np.transpose(N_ts):
                ax.plot(Brange, n_ts / 1e15, 'r-')
            for n_bs in np.transpose(N_bs):
                ax.plot(Brange, n_bs / 1e15, 'b-')
            for n_vpsup in np.transpose(N_vpsup):
                ax.plot(Brange, n_vpsup / 1e15, 'k-')
            for n_vpsdown in np.transpose(N_vpsdown):
                ax.plot(Brange, n_vpsdown / 1e15, 'y-')

        else:
            N_ts, N_bs = self.get_ll_den(angle, Brange, Erange, Nmax, den_top, den_bot, den_vps, threeband)
            ax.axhline(y=(den_top + den_bot) / 1e15, linewidth=1)
            for n_ts in np.transpose(N_ts):
                ax.plot(Brange, n_ts / 1e15, 'r-')
            for n_bs in np.transpose(N_bs):
                ax.plot(Brange, n_bs / 1e15, 'b-')
        ax.set_ylim(-4, 4)
        ax.set_xlabel('B (T)')
        ax.set_ylabel('Density ($10^{11}cm^{-2}$)')
        return fig, ax

    def IDOS_generator(self, angle, Brange, Erange, LLenergy_top_surface, LLenergy_bottom_surface, LLenergy_vps_up=None,
                       LLenergy_vps_down=None):
        sigma = self.sigma
        IDOS = []

        if all([LLenergy_vps_up, LLenergy_vps_down]):
            for B, llenergy_top_surface, llenergy_bottom_surface, llenergy_vps_up, llenergy_vps_down in zip(Brange,
                                                                                                            LLenergy_top_surface,
                                                                                                            LLenergy_bottom_surface,
                                                                                                            LLenergy_vps_up,
                                                                                                            LLenergy_vps_down):
                IDOS_B = [x + y for x, y in zip(
                    fastIntegral_electron_DOS(Erange, B, sigma, angle, llenergy_top_surface, llenergy_bottom_surface),
                    fastIntegral_hole_DOS(Erange, B, sigma, angle, llenergy_vps_up, llenergy_vps_down))]
                IDOS.append(IDOS_B)
        else:
            for B, llenergy_top_surface, llenergy_bottom_surface in zip(Brange, LLenergy_top_surface,
                                                                        LLenergy_bottom_surface):
                IDOS_B = fastIntegral_electron_DOS(Erange, B, sigma, angle, llenergy_top_surface,
                                                   llenergy_bottom_surface)
                IDOS.append(IDOS_B)

        return IDOS

    def plot_DOS(self, angle, Bfield, Erange, Nmax, den_top, den_bot, den_vps=None, threeband=False):

        vf = self.vf
        sigma = self.sigma
        gfactor = self.gfactor
        meff = self.meff

        Ets = -hbar * vf * (4 * np.pi * den_top) ** 0.5
        Ebs = -hbar * vf * (4 * np.pi * den_bot) ** 0.5
        lldegeneracy = Bfield * np.cos(angle * np.pi / 180) * e0 / h0

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)

        if threeband:
            Evp = -hbar ** 2 * den_vps * np.pi / (meff * me) / 2
            llenergy_top_surface, llenergy_bottom_surface, llenergy_vps_up, llenergy_vps_down = llenergy_generator(Ets,
                                                                                                                   Ebs,
                                                                                                                   Evp,
                                                                                                                   Bfield,
                                                                                                                   angle,
                                                                                                                   meff,
                                                                                                                   Nmax,
                                                                                                                   vf,
                                                                                                                   gfactor)

            ax.fill_between(1e3 * Erange / e0,
                            [electron_density_of_state(E, Bfield, sigma, angle, llenergy_top_surface, []) for E in
                             Erange], 0, color='r', alpha=0.3)
            ax.fill_between(1e3 * Erange / e0,
                            [electron_density_of_state(E, Bfield, sigma, angle, [], llenergy_bottom_surface) for E in
                             Erange], 0, color='b', alpha=0.3)
            ax.fill_between(1e3 * Erange / e0,
                            [hole_density_of_state(E, Bfield, sigma, angle, llenergy_vps_up, []) for E in Erange], 0,
                            color='k', alpha=0.3)
            ax.fill_between(1e3 * Erange / e0,
                            [hole_density_of_state(E, Bfield, sigma, angle, [], llenergy_vps_down) for E in Erange], 0,
                            color='y', alpha=0.3)


        else:

            llenergy_top_surface, llenergy_bottom_surface, _, _ = llenergy_generator(Ets, Ebs, 0, Bfield, angle, meff,
                                                                                     Nmax, vf, gfactor)
            ax.fill_between(1e3 * Erange / e0,
                            [electron_density_of_state(E, Bfield, sigma, angle, llenergy_top_surface, []) for E in
                             Erange], 0, color='r', alpha=0.3)
            ax.fill_between(1e3 * Erange / e0,
                            [electron_density_of_state(E, Bfield, sigma, angle, [], llenergy_bottom_surface) for E in
                             Erange], 0, color='b', alpha=0.3)

        ax.set_xlabel('Energy (meV)')
        return fig, ax

    def plot_muDOS(self, angle, Brange, Erange, Nmax, den_top, den_bot, den_vps=None, threeband=False):

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        sigma = self.sigma

        if threeband:

            LL_ts, LL_bs, LL_vpsup, LL_vpsdown = self.get_ll_en(angle, Brange, Nmax, den_top, den_bot, den_vps,
                                                                threeband)
            IDOS = self.IDOS_generator(angle, Brange, Erange, LL_ts, LL_bs, LL_vpsup, LL_vpsdown)
            mu = [find_energy_bydensity(den_top + den_bot - den_vps, B, IDOS_B, Erange) for B, IDOS_B in
                  zip(Brange, IDOS)]
            ax.plot(Brange, [electron_density_of_state(E, B, sigma, angle, ll_ts, []) for E, B, ll_ts in
                             zip(mu, Brange, LL_ts)], linewidth=1, color='r')
            ax.plot(Brange, [electron_density_of_state(E, B, sigma, angle, [], ll_bs) for E, B, ll_bs in
                             zip(mu, Brange, LL_bs)], linewidth=1, color='b')
            ax.plot(Brange,
                    [hole_density_of_state(E, B, sigma, angle, ll_up, []) for E, B, ll_up in zip(mu, Brange, LL_vpsup)],
                    linewidth=1, color='k')
            ax.plot(Brange, [hole_density_of_state(E, B, sigma, angle, [], ll_down) for E, B, ll_down in
                             zip(mu, Brange, LL_vpsdown)], linewidth=1, color='y')

        else:
            LL_ts, LL_bs = self.get_ll_en(angle, Brange, Nmax, den_top, den_bot, den_vps, threeband)
            IDOS = self.IDOS_generator(angle, Brange, Erange, LL_ts, LL_bs)
            mu = [find_energy_bydensity(den_top + den_bot, B, IDOS_B, Erange) for B, IDOS_B in zip(Brange, IDOS)]
            ax.plot(Brange, [electron_density_of_state(E, B, sigma, angle, ll_ts, []) for E, B, ll_ts in
                             zip(mu, Brange, LL_ts)], linewidth=1, color='r')
            ax.plot(Brange, [electron_density_of_state(E, B, sigma, angle, [], ll_bs) for E, B, ll_bs in
                             zip(mu, Brange, LL_bs)], linewidth=1, color='b')
            ax.plot(Brange, [electron_density_of_state(E, B, sigma, angle, ll_ts, ll_bs) for E, B, ll_ts, ll_bs in
                             zip(mu, Brange, LL_ts, LL_bs)], linewidth=2, color='k', linestyle='--')

        return fig, ax