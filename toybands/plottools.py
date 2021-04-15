import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as mplcolors
import pdb
from matplotlib.backends.backend_pdf import PdfPages

from physconst import *
from utils import flattenList, div
from toybands.functions import extract_list,add_list, extents, mkdir
from toybands.config import *


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


def make_1d_E_B_plots(bfrange,y_databdl,colors, mu_pos = None,enrange=None,figsize=DEFAULT_FIGURE_SIZE,linewidth=DEFAULT_LW,ax=None,plotrange = None,legend=True):
    if not isinstance(bfrange,list) or not isinstance(y_databdl,list) or not isinstance(colors,list):
        raise TypeError(f'either x_databdl or y_databdl or colors is not list')
    if not len(y_databdl)==len(colors):
        raise ValueError(f'y_databdl, colors are not of the same length')
    if not any(isinstance(el, list) for el in y_databdl):
        raise TypeError(f'y_databdl is not nested list')
    if not isinstance(figsize,tuple):
        raise TypeError(f'figsize should be a tuple like (10,10)')
    if ax is None:
        ax = make_canvas(figsize=figsize)
    for n_band, (y_data, color) in enumerate(zip(y_databdl,colors)):
        for y in y_data:
            line, = ax.plot(bfrange,div(y,e0),linewidth=linewidth,color=color)
        line.set_label(f'Band{n_band}')
    if legend:
        ax.legend(loc=DEFAULT_LEGEND_LOC,bbox_to_anchor=DEFAULT_LEGEND_POS)
    if mu_pos and len(mu_pos)==len(bfrange):
        ax.plot(bfrange,div(mu_pos,e0),linewidth=linewidth,color='k')
    if enrange is not None:
        ax.set_ylim(min(enrange)/e0,max(enrange)/e0)

    ax.set_xlabel(DEFAULT_XLABEL)
    ax.set_ylabel(DEFAULT_EBPLOT_YLABEL)
    return ax




def make_1d_den_B_plots(bfrange,y_databdl, colors, tot_den = None,enrange=None,figsize=DEFAULT_FIGURE_SIZE,linewidth=DEFAULT_LW, ax=None,plotrange = None,legend=True):
    if not isinstance(bfrange,list) or not isinstance(y_databdl,list) or not isinstance(colors,list):
        raise TypeError(f'either x_databdl or y_databdl or colors is not list')
    if not len(y_databdl)==len(colors):
        raise ValueError(f'y_databdl, colors are not of the same length')
    if not any(isinstance(el, list) for el in y_databdl):
        raise TypeError(f'y_databdl is not nested list')
    if not isinstance(figsize,tuple):
        raise TypeError(f'figsize should be a tuple like (10,10)')
    if ax is None:
        ax = make_canvas(figsize=figsize)
    if plotrange is None:
        for n_band, (y_data, color) in enumerate(zip(y_databdl,colors)):
            for n_ll, y in enumerate(y_data):
                line, = ax.plot(bfrange,y,linewidth=linewidth,color=color)
            line.set_label(f'Band{n_band}')
    else:
        for y_data, color in zip(y_databdl, colors):
            for y in y_data:
                ax.scatter(extract_list(bfrange,[yy>plotrange[0] and yy<plotrange[1] for yy in y]),extract_list(y,[yy>plotrange[0] and yy<plotrange[1] for yy in y]), color=color)
    if legend:
        ax.legend(loc=DEFAULT_LEGEND_LOC,bbox_to_anchor=DEFAULT_LEGEND_POS)
    if tot_den is not None:
        ax.axhline(y = tot_den, linewidth=linewidth, color='k')
    if enrange is not None:
        ax.set_ylim(min(enrange)/e0,max(enrange)/e0)
    
    ax.set_xlabel(DEFAULT_XLABEL)
    ax.set_ylabel(DEFAULT_NBPLOT_YLABEL)
    return ax

def make_1d_dos_B_plot(bfrange,y_databdl,colors,figsize=DEFAULT_FIGURE_SIZE,linewidth=DEFAULT_LW, ax=None,plotrange = None,legend=True):
    if ax is None:
        ax = make_canvas(figsize=figsize)
    y_tot = [0]*len(bfrange)
    for n_band, (y_data, color) in enumerate(zip(y_databdl,colors)):
        line, = ax.plot(bfrange,y_data,color=color)
        line.set_label(f'Band{n_band}')
        y_tot = add_list(y_tot,y_data)
    if legend:
        ax.legend(loc=DEFAULT_LEGEND_LOC,bbox_to_anchor=DEFAULT_LEGEND_POS)
    ax.plot(bfrange,y_tot,linestyle='--',color='k')
    ax.set_xlabel(DEFAULT_XLABEL)
    ax.set_ylabel(DEFAULT_DOSBPLOT_YLABEL)
    return ax

def make_2d_dos_map(bfrange,y_tot,y_databdl,cmap,figsize=DEFAULT_FIGURE_SIZE,ax=None,legend=False):
    if ax is None:
        ax = make_canvas(figsize=figsize)
    y_array = np.array(y_databdl)
    # make extents define the bottom as the last row
    if y_tot[0]>y_tot[-1]:
        y_tot.reverse()
    ax.imshow(y_array, aspect='auto', interpolation='none', extent=(extents(bfrange) + extents(y_tot)), origin = 'upper')
    ax.set_xlabel(DEFAULT_XLABEL)
    ax.set_ylabel(DEFAULT_DOSMAP_YLABEL)
    ax.set_ylim(min(y_tot),max(y_tot))
    
    return ax

def draw_band(system,figsize=DEFAULT_FIGURE_SIZE,legend=True):
    colors = make_n_colors(len(system.get_band('a')), DEFAULT_CMAP, DEFAULT_CMAP_VMIN, DEFAULT_CMAP_VMAX)
    k = np.linspace(-0.5,0.5,101)*1e9
    ax = make_canvas(figsize=figsize)
    for n_band,(band,color) in enumerate(zip(system.get_band('a'),colors)):
        if band.is_dirac and band.is_cond:
            enk = [(band.Ebb-kk*hbar*band.vf+band.dparam*kk**2)/e0 for kk in k if kk<0]+[(band.Ebb+kk*hbar*band.vf+band.dparam*kk**2)/e0 for kk in k if kk>=0]
        elif band.is_dirac and not band.is_cond:
            enk = [(band.Ebb+kk*hbar*band.vf-band.dparam*kk**2)/e0 for kk in k if kk<0]+[(band.Ebb-kk*hbar*band.vf-band.dparam*kk**2)/e0 for kk in k if kk>=0]
        elif not band.is_dirac and band.is_cond:
            enk = [(band.Ebb+(hbar*kk)**2/2/band.meff/me)/e0 for kk in k]
        elif not band.is_dirac and not band.is_cond:
            enk = [(band.Ebb-(hbar*kk)**2/2/band.meff/me)/e0 for kk in k]
        line, = ax.plot(k/1e9,enk,color=color,linewidth = DEFAULT_LW)
        line.set_label(f'Band{n_band}')
    if legend:
        ax.legend(loc=DEFAULT_LEGEND_LOC,bbox_to_anchor=DEFAULT_LEGEND_POS)
    ax.axhline(y=0, linestyle = '--', linewidth = DEFAULT_LW)
    ax.set_xlabel('k [nm$^{-1}$]')
    ax.set_ylabel('E (eV)')

def plot_from_csv(path,ax=None,cmap=None,legend=True):
    if not isinstance(path,str):
        sys.stderr.write(f'the path {path} is not a string\n')
    if not os.path.isfile(path):
        sys.stderr.write(f'the file {path} does not exist\n')
    if not path.endswith('.csv'):
        sys.stderr.write(f'the file {path} is not a csv file\n')
    try:
        df = pd.read_csv(path)
    except:
        sys.stderr.write(f'Failed to read the csv file\n')
    if cmap is None:
        cmap = DEFAULT_CMAP
    colors = make_n_colors(len(df.Band.unique()),cmap,DEFAULT_CMAP_VMIN,DEFAULT_CMAP_VMAX)
    if ax is None:
        ax = make_canvas()
    if 'N' in df.columns:
        x,y,N,iBand = [],[],df.iloc[0].N,df.iloc[0].Band
    else:
        x,y,iBand = [],[],df.iloc[0].Band
    
    # define plot parameter and type
    if 'System([band density])' in df.columns:
        ind, plottype = 'den','scatter'
    elif 'den' in df.columns:
        ind,plottype = 'den','plot'
    elif 'E' in df.columns:
        ind,plottype = 'E','plot'
    elif 'dos_at_mu' in df.columns:
        ind,plottype = 'dos_at_mu','plot'
    else:
        sys.stderr.write('This file is not readable by toybands\n')
    
    # plot
    for i in range(len(df)):
        if 'N' in df.columns:
            if df.iloc[i].N != N or df.iloc[i].Band != iBand:
                N = df.iloc[i].N
                if plottype == 'plot':
                    if ind == 'E':
                        line, = ax.plot(x, div(y,e0), color=colors[int(iBand)],linewidth=DEFAULT_LW)
                    else:
                        line, = ax.plot(x, y, color=colors[int(iBand)],linewidth=DEFAULT_LW)
                else:
                    ax.scatter(x, y, color=colors[int(iBand)])
                if df.iloc[i].Band != iBand:
                    if plottype == 'plot':
                        line.set_label(f'Band{int(iBand)}')
                    iBand = df.iloc[i].Band
                x,y = [],[]
            x.append(df.iloc[i].B)
            y.append(df.iloc[i][ind])

        else:
            if df.iloc[i].Band != iBand:
                line, = ax.plot(x, y, color=colors[int(iBand)],linewidth=DEFAULT_LW)
                line.set_label(f'Band{int(iBand)}')
                iBand = df.iloc[i].Band
                x,y = [],[]
                x.append(df.iloc[i].B)
                y.append(df.iloc[i][ind])
            else:
                x.append(df.iloc[i].B)
                y.append(df.iloc[i][ind])
    # plot the stored last curve
    if 'System([band density])' in df.columns:
        ax.scatter(x, y, color=colors[int(iBand)]) 
    else:        
        line, = ax.plot(x, y, color=colors[int(iBand)],linewidth=DEFAULT_LW)
        line.set_label(f'Band{int(iBand)}')
    
    # label and legend
    ax.set_xlabel(DEFAULT_XLABEL)
    if ind == 'den':
        ax.set_ylabel(DEFAULT_NBPLOT_YLABEL)
    elif ind == 'E':
        ax.set_ylabel(DEFAULT_EBPLOT_YLABEL)
    elif ind == 'dos_at_mu':
        ax.set_ylabel(DEFAULT_DOSBPLOT_YLABEL)
    if legend and not 'System([band density])' in df.columns:
        ax.legend(loc=DEFAULT_LEGEND_LOC,bbox_to_anchor=DEFAULT_LEGEND_POS)
    
    super_save(filename=path.split('/')[-1].split('.')[-2]+'-replot')


def make_canvas(figsize=DEFAULT_FIGURE_SIZE):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    return ax

def super_save(filename=None,path=None):
    if filename is None:
        filename = DEFAULT_AUTONAME
    if path is None:
        path = mkdir(filename)
    if not os.path.isdir(path):
        os.mkdir(path)
        sys.stdout.write(f'path created under {path}')
    if len(filename.split('.'))>1:
        fmt = filename.split('.')[-1]
        if fmt in ALLOW_FORMAT:
            plt.savefig(os.path.join(path,filename))
        else:
            plt.savefig(os.path.join(path,filename.split('.')[0]+'.'+DEFAULT_FORMAT))
            sys.stderr.write(f'The format {fmt} is not supported, save as {DEFAULT_FORMAT} file instead\n')
    else:
        plt.savefig(os.path.join(path, filename+'.'+DEFAULT_FORMAT))
        sys.stdout.write(f'By default, saved as a {DEFAULT_FORMAT} file\n')

def pdf_save(filename=None,path=None,end=True):
    if filename is None:
        filename = '[auto]default.pdf'
    if path is None:
        path = mkdir(filename)
    if len(filename.split('.'))>1:
        filename = filename.split('.')[-2]
    filename = filename+'.pdf'
    if not os.path.isdir(path):
        os.mkdir(path)
        sys.stdout.write(f'path created under {path}')
    with PdfPages(os.path.join(path,filename)) as pdf:
        pdf.savefig()
        if end:
            plt.close()

def make_slices(den_list,numofsteps):
    if not isinstance(den_list, list):
        sys.stderr.write('Input is not list\n')
    all_density  = []
    for start,end in zip(den_list[::2],den_list[1::2]):
        all_density.append(np.linspace(start,end,numofsteps).tolist())
    all_density = np.array(all_density)
    output = np.transpose(all_density)
    return output
