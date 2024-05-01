#Copyright 2021 Lixian WANG. All Rights Reserved.
# Standard library
import os

# Third party
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# Local application
from .physconst import *
from .functions import *

__all__ = ['Databs','Datags','Datamap','DataX']

class Datajungle:
    """ Parent Class for Generic data type
    Auguments:
    directory: list of filenames
    step: step values
    ucols: extracted columns from original source file
    spr: skipped rows in the header of source file
    ref: reference resistance in series
    AspRatio: aspect ratio of Hall bar, set to 3 by default

    Return:

    CHILDREN CLASS:
    Databs, Datags, Datafc
    """
    def __init__(self,directory,step,ucols,nms,spr):
        self.dir = directory # directories for all the files in current folder.
        self.step = step # values of parameter B
        self.ucols = ucols # choose columns to import
        self.nms = nms # to name chosen columns
        self.spr = spr  # to skip rows in the header of .dat file



class Databs(Datajungle):
    """Inherent from Class Datajungle
    METHODS:
    getdata: return a panda.DataFrame type data
    hallfit: linear hall fit and return density and mobility in a DataFrame format
    plotdata: plot magnetic field sweep type data in a specific way
    plotfc: plot fan chart"""

    def __init__(self, directory, step, ucols, spr, ref, nms=['bf', 'curr', 'uxx', 'uxy'],AspRatio=3):
        super().__init__(directory,step,ucols,nms,spr)
        self.ref = ref # reference resistance in series
        self.AspRatio = AspRatio  # Aspect ratio of Hall bar structure

    def __str__(self):
        return ', '.join(['Databs', ', '.join(['{key} = {value}'.format(key = key, value = self.__dict__[key]) for key in ['AspRatio','ucols','nms','spr','ref']])])

    def __repr__(self):
        pass

    def getdata(self):
        databundle = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None, encoding='unicode_escape')
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['gate'] = self.step[i]
            # databundle = databundle.append(data) # Deprecated in Pandas 1.4.0 and above
            databundle = pd.concat([databundle, data],ignore_index=False)
        return databundle

    def hallfit(self,fitrange):
        Dens = []
        Mob = []
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None, encoding='unicode_escape')
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            bf_fit = data['bf'][(data['bf']<fitrange[1])&(data['bf']>fitrange[0])]
            rxx_fit = data['rxx'][(data['bf']<fitrange[1])&(data['bf']>fitrange[0])]
            rxy_fit = data['rxy'][(data['bf']<fitrange[1])&(data['bf']>fitrange[0])]
            dens,mob = H1st_ft(bf_fit,rxx_fit,rxy_fit,AspRatio=AspRatio)
            Dens.append(dens)
            Mob.append(mob)
        FitRes = pd.DataFrame({'gate':self.step,'dens':Dens,'mob':Mob})
        return FitRes

    def plotdata(self,label_value = '$V_g$={:02.2f}V' ):
        font = {'family' : 'normal','weight' : 'normal','size'   : 15}
        matplotlib.rc('font', **font)
        jet= plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0,1,len(self.dir))))
        plt.figure(figsize=(16,16))
        ax_rxy = plt.subplot(2,1,1)
        ax_rxx = plt.subplot(2,2,3)
        ax_sxy = plt.subplot(2,2,4)
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            line_color = next(colors)
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None, encoding='unicode_escape')
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            ax_rxx.plot(data.bf,data['rxx'],color = line_color,label=label_value.format(self.step[i]))
            ax_rxy.plot(data.bf,data['rxy'],color = line_color,label=label_value.format(self.step[i]))
            ax_sxy.plot(data.bf,data['sxy'],color = line_color,label=label_value.format(self.step[i]))
        ax_rxx.set_xlabel(r'$B_{field}(T)$',fontsize = 18)
        ax_rxx.set_ylabel(r'$R_{xx}(\Omega)$',fontsize = 18)
        ax_rxy.set_xlabel(r'$B_{field}(T)$',fontsize = 18)
        ax_rxy.set_ylabel(r'$R_{xy}(\Omega)$',fontsize = 18)
        ax_sxy.set_xlabel(r'$B_{field}(T)$',fontsize = 18)
        ax_sxy.set_ylabel(r'$\sigma_{xy}(e^2/h)$',fontsize = 18)
        ax_sxy.set_ylim([-10,10])
        ax_rxy.legend(loc = 'best')
        for mark in range(-5,6):
                   ax_sxy.axhline(y=mark,linestyle=':',color='c')
        for mark in [-5,-4,-3,-2,-1,1,2,3,4,5]:
                   ax_rxy.axhline(y=h0/e0**2/mark,linestyle=':',color='c')
        return ax_rxx,ax_rxy,ax_sxy


class Datags(Datajungle):
    '''Inherent from Class Datajungle
    METHODS:
    getdata: return a panda.DataFrame type data
    plotdata: plot gate sweep type data in a specific way
    plotfc: plot fan chart'''

    def __init__(self, directory, step, ucols, spr, ref, nms=['gate', 'curr', 'uxx', 'uxy'], AspRatio=3):
        super().__init__(directory,step,ucols,nms,spr)
        self.ref = ref # reference resistance in series
        self.AspRatio = AspRatio  # Aspect ratio of Hall bar structure

    def __str__(self):
        return ', '.join(['Datags', ', '.join(['{key} = {value}'.format(key = key, value = self.__dict__[key]) for key in ['AspRatio','ucols','nms','spr','ref']])])

    def __repr__(self):
        pass

    def getdata(self):
        databundle = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None, encoding='unicode_escape')
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['bf'] = self.step[i]
            # databundle = databundle.append(data) # Deprecated in Pandas 1.4.0 and above
            databundle = pd.concat([databundle,data],ignore_index=False)
        return databundle
    def plotdata(self,label_value='$B$={:02.2f}T'):
        font = {'family' : 'normal','weight' : 'normal','size' : 15}
        matplotlib.rc('font', **font)
        jet= plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0,1,len(self.dir))))
        plt.figure(figsize=(16,16))
        ax_rxy = plt.subplot(2,1,1)
        ax_rxx = plt.subplot(2,2,3)
        ax_sxy = plt.subplot(2,2,4)
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            line_color = next(colors)
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None, encoding='unicode_escape')
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            ax_rxx.plot(data.gate,data['rxx'],color = line_color,label=label_value.format(self.step[i]))
            ax_rxy.plot(data.gate,data['rxy'],color = line_color,label=label_value.format(self.step[i]))
            ax_sxy.plot(data.gate,data['sxy'],color = line_color,label=label_value.format(self.step[i]))
        ax_rxx.set_xlabel(r'$V_g(V)$',fontsize = 18)
        ax_rxx.set_ylabel(r'$R_{xx}(\Omega)$',fontsize = 18)
        ax_rxy.set_xlabel(r'$V_g(V)$',fontsize = 18)
        ax_rxy.set_ylabel(r'$R_{xy}(\Omega)$',fontsize = 18)
        ax_sxy.set_xlabel(r'$V_g(V)$',fontsize = 18)
        ax_sxy.set_ylabel(r'$\sigma_{xy}(e^2/h)$',fontsize = 18)
        ax_sxy.set_ylim([-10,10])
        ax_rxy.legend(loc = 'best')
        for mark in range(-5,6):
                   ax_sxy.axhline(y=mark,linestyle=':',color='c')
        for mark in [-5,-4,-3,-2,-1,1,2,3,4,5]:
                   ax_rxy.axhline(y=h0/e0**2/mark,linestyle=':',color='c')
        return ax_rxx,ax_rxy,ax_sxy

class Datamap(Datajungle):
    '''Inherent from Class Datajungle
    METHODS:
    getdata: return x, y and a 2D array with z-value
    plotmap: plot 2D mapping'''

    def __init__(self,directory,step,ucols,spr,ref,nms=['v1','curr','uxx','uxy'],AspRatio=3):
        super().__init__(directory,step,ucols,nms,spr)
        self.ref = ref # reference resistance in series
        self.AspRatio = AspRatio  # Aspect ratio of Hall bar structure

    def __str__(self):
        return ', '.join(['Datagmap', ', '.join(['{key} = {value}'.format(key = key, value = self.__dict__[key]) for key in ['AspRatio','ucols','nms','spr','ref']])])

    def __repr__(self):
        pass

    def getdata(self):
        databundle = pd.DataFrame()
        diffsxy2D_v1 = pd.DataFrame()
        diffsxy2D_v2 = pd.DataFrame()
        rxx2D = pd.DataFrame()
        rxy2D = pd.DataFrame()
        sxy2D = pd.DataFrame()
        sxx2D = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None, encoding='unicode_escape')
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['v2'] = self.step[i]
            rxx2D = pd.concat([rxx2D,data['rxx']],axis=1)
            # rxy2D = rxy2D.append(data['rxy'])
            rxy2D = pd.concat([rxy2D,data['rxy']],axis=1)
            # sxy2D = sxy2D.append(data['sxy'])
            sxy2D = pd.concat([sxy2D,data['sxy']],axis=1)
            # sxx2D = sxx2D.append(data['sxx'])
            sxx2D = pd.concat([sxx2D,data['sxx']],axis=1)
            # databundle = databundle.append(data)
            databundle = pd.concat([databundle,data],ignore_index=False)

        diffsxy2D_v1 = sxy2D.diff(axis=0)/abs(databundle['v1'].unique()[0]-databundle['v1'].unique()[1])
        diffsxy2D_v2 = sxy2D.diff(axis=1)/abs(self.step[0]-self.step[1])

        datafc = {'v1':databundle['v1'].unique(),'v2':self.step,'dv1':diffsxy2D_v1,'dv2':diffsxy2D_v2,'rxx2d':rxx2D,'rxy2d':rxy2D,'sxy2d':sxy2D,'sxx2d':sxx2D}
        return datafc, databundle

    def plotmap(self,vmin1,vmax1,vmin2,vmax2,cmap='terrain'): # plot gate-mapping
        fc,_ = self.getdata()
        v1 = fc['v1']
        v2 = fc['v2']
        diffsxy2D_v1 = fc['dv1']
        diffsxy2D_v2 = fc['dv2']

        fig = plt.figure(figsize=(8,12),constrained_layout=True)
        gs = fig.add_gridspec(2, 1)
        ax1,ax2 = [fig.add_subplot(x) for x in gs]
        dv1plot = ax1.imshow(diffsxy2D_v1,aspect='auto', interpolation='none',extent=extents(v1.tolist()) + extents(v2.tolist()), origin='lower',cmap=cmap,vmin=vmin1, vmax=vmax1)
        cbaxes1 = fig.add_axes([.8, .9, 0.15, 0.01])
        cb1 = fig.colorbar(dv1plot, cax=cbaxes1, orientation='horizontal', pad=10)
        cb1.ax.set_xlabel('$d\sigma_{xy}/dV_{1}$')
        dv2plot = ax2.imshow(diffsxy2D_v2,aspect='auto', interpolation='none',extent=extents(v1.tolist()) + extents(v2.tolist()), origin='lower',cmap=cmap,vmin=vmin2, vmax=vmax2)
        cbaxes2 = fig.add_axes([.8, .4, 0.15, 0.01])
        cb2 = fig.colorbar(dv2plot, cax=cbaxes2, orientation='horizontal', pad=10)
        cb2.ax.set_xlabel('$d\sigma_{xy}/dV_{2}$')
        ax1.set_ylabel('$V_{2}$')
        ax1.set_xlabel('$V_{1}$')
        ax2.set_ylabel('$V_{2}$')
        ax2.set_xlabel('$V_{1}$')
        return [ax1,ax2]

class DataX(Datajungle):
    '''Inherent from Class Datajungle
    METHODS:
    getdata: return databundle'''

    def __init__(self,directory,step,ucols,spr,nms):
        super().__init__(directory,step,ucols,nms,spr)

    def __str__(self):
        return ', '.join(['DataX', ', '.join(['{key} = {value}'.format(key = key, value = self.__dict__[key]) for key in ['ucols','nms','spr','step']])])

    def __repr__(self):
        pass

    def getdata(self):
        databundle = pd.DataFrame()
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None, encoding='unicode_escape')
            data['x'] = self.step[i]
            # databundle = databundle.append(data) # Deprecated in Pandas 1.4.0 and above
            databundle = pd.concat([databundle,data],ignore_index=False)
        return databundle

def denCal_single(data_formatted, AspRatio, bf_range_to_fit, gate_range_to_fit, residual_field_in_T, call=False):
    # PURPOSE: calculate the carier density/mobility by the low-field Hall and transverse resistance
    # INPUT: databs | type class Databs() or Datafc()
    # OUTPUT: return dens/mob | type list
    if isinstance(data_formatted, (Databs, Datags)):
        data = data_formatted.getdata()
    elif isinstance(data_formatted, Datamap):
        _, data = data_formatted.getdata()
    else:
        raise TypeError('Error: Wrong input data type')

    gates = gate_range_to_fit

    bf_range = bf_range_to_fit  # range to fit

    dens = []
    mob = []

    for gate in gates:
        data_p = df_range(df_range(data, 'bf', bf_range), 'gate', [gate - 0.005, gate + 0.005])  # data
        para1, para2 = H1st_ft(data_p.bf - residual_field_in_T, data_p.rxx, data_p.rxy, AspRatio=AspRatio,
                               threshold=1000)
        dens.append(para1 / 1e11)
        mob.append(para2)
    return dens, mob


def denCal_double(data_formatted, AspRatio, bf_range_to_fit, gate_range_to_fit, residual_field_in_T, alternating=False, call=False):
    # PURPOSE: calculate the carier density/mobility by a two-carrier model
    # INPUT: databs | type class Databs, Datags and Datamap
    # OUTPUT: return dens/mob | type list
    if isinstance(data_formatted, (Databs, Datags)):
        data = data_formatted.getdata()
    elif isinstance(data_formatted, Datamap):
        _, data = data_formatted.getdata()
    else:
        raise TypeError('Error: Wrong input data type')

    gates = gate_range_to_fit
    bf_range = bf_range_to_fit  # range to fit
    ndens, nmob, pdens, pmob = [], [], [], []

    if call == True:
        fig = plt.figure(figsize=(5, 5), constrained_layout=True)
        ax1 = fig.add_subplot(111)
    for n, gate in enumerate(gates):
        if alternating == True:
            residual_field_in_T = -residual_field_in_T
        data_p = df_range(df_range(data, 'bf', bf_range), 'gate', [gate - 0.005, gate + 0.005])  # data
        para, fun = twocarrierfit(data_p.bf - residual_field_in_T, data_p.rxy)
        if call == True:
            ax1.plot(data_p.bf - residual_field_in_T, data_p.rxy, c='r')
            ax1.plot(data_p.bf - residual_field_in_T, fun, c='k')
        ndens.append(para[0] / 1e15)
        nmob.append(para[1] * 1e4)
        pdens.append(para[2] / 1e15)
        pmob.append(para[3] * 1e4)
    return ndens, nmob, pdens, pmob




def main():
    ''' main function '''
    print('running SciData in main program')

if __name__ == '__main__':
    main()
