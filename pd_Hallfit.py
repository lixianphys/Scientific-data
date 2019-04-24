import numpy as np
import glob
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import os


# define custom function
def H1st_ft(Bf,Rxx,Rxy,AspRatio=3,threshold = 25):
    def func_one(x, a, b):
        return a + b * x
    e0 = 1.6021766208E-19
    h0 = 6.62607015E-34
    try:
        fitParams, fitCovariances = curve_fit(func_one, Bf, Rxy)
    except:
        print('The fitting program failed')
        mobility = 0
        density = 0
    else:
        dev = np.sqrt(np.diag(fitCovariances))
        if sum(dev) <= threshold:
            density = 1 / fitParams[1] / e0 / 1e4
            rxx0 = Rxx.tolist()[list(map(abs,Bf.tolist())).index(min(map(abs,Bf.tolist())))]
            mobility = AspRatio / density / e0 / rxx0
            plt.plot(Bf, Rxy, "b-", Bf, func_one(Bf, *fitParams), "r-")
            plt.show()
        else:
            print('The fitting results is not acceptable, fitCov is {}'.format(fitCovariances))
            plt.plot(Bf, Rxy, "b-", Bf, func_one(Bf, *fitParams), "r-")
            mobility = 0
            density = 0
    return density,mobility

def H2nd_ft(Bf,Rxx,Rxy,AspRatio=3):
    def func_two(x, n1, m1, n2, m2):
        return e0 * x * (n1 * m1 ** 2 / (1 + m1 ** 2 * x ** 2) + n2 * m2 ** 2 / (1 + m2 ** 2 * x ** 2)) #model from PHYSICAL REVIEW B 95, 115126 (2017)
    e0 = 1.6021766208E-19
    h0 = 6.62607015E-34
    sxy = Rxy / ((Rxx / AspRatio) ** 2 + Rxy ** 2)
    try:
            popt, pcov = curve_fit(func_two, Bf, sxy, bounds=((2e13, 1, -1e16, 0.1), (2e15, 100, -1e14, 10)))
            plt.plot(Bf, sxy, "b-", Bf, func_two(Bf, *popt), "r-")
            plt.show()
            print('The fitCov is {}'.format(pcov))
    except:
            print('The fitting program failed')
    return  popt


