import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as mplcolors

# define the norm in plotting
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



def make_n_colors(n,cmap,vstart,vend):
    if not isinstance(n,int):
        raise TypeError('n must be an integer')
    if not isinstance(cmap,str):
        raise TypeError('cmap should be the name of an existing cmap instance')
    if vstart>1 or vstart<0:
        raise ValueError('vstart should be in [0.0,1.0]')
    if vend>1 or vend<0:
        raise ValueError('vend should be in [0.0,1.0]')
    try:
        cmap = mplcm.get_cmap(cmap)
    except:
        raise ValueError('failed to fetch colormap {cmap}')
    colors_list = np.linspace(vstart,vend,n).tolist()
    colors = [cmap(x) for x in colors_list]
    return colors

def make_1d_E_B_plots(bfrange,y_databdl,colors,figsize=(10,10),filename=None,linewidth=2):
    if not isinstance(bfrange,list) or not isinstance(y_databdl,list) or not isinstance(colors,list):
        raise TypeError(f'either x_databdl or y_databdl or colors is not list')
    if not len(y_databdl)==len(colors):
        raise ValueError(f'y_databdl, colors are not of the same length')
    if not any(isinstance(el, list) for el in y_databdl):
        raise TypeError(f'y_databdl is not nested list')
    if not isinstance(figsize,tuple):
        raise TypeError(f'figsize should be a tuple like (10,10)')
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    for y_data, color in zip(y_databdl,colors):
        for y in y_data:
            ax.plot(bfrange,y,linewidth=linewidth,color=color)
    ax.set_xlabel('$B$ [T]')
    ax.set_ylabel('$E$ [eV]')
    plt.savefig('test.pdf')


