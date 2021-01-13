# Standard library imports
import os

# Third party imports
import matplotlib.pyplot as plt
import pandas as pd

# Local application import
# from Scientific_data import *
import SciData
from SciData import *
from functions import *
from physconst import *


# plot default
plt.rc('lines', lw=1, color='k')  # thicker black lines
plt.rc('grid', c='0.5', ls='-', lw=0.5)  # solid gray grid lines
plt.rc('savefig', dpi=600)  # higher res outputs
plt.rc("font", size=20, family='arial', weight='light')
plt.rc("axes", labelsize=20, titlesize=20, linewidth=1)
plt.rc("xtick", direction='in', labelsize=20)
plt.rc('xtick.major', size=10, pad=7)
plt.rc('xtick.minor', size=5, pad=7, visible=True)
plt.rc("ytick", direction='in', labelsize=20)
plt.rc('ytick.major', size=10, pad=7)
plt.rc('ytick.minor', size=5, pad=7, visible=True)
plt.rc("legend", fontsize=20)
plt.rc('lines', linewidth=2)

# interface
print('welcome')
print('You are using the init.py located in ', __file__)
print('You are using the SciData package located in ', SciData.__file__)