#Copyright 2021 Lixian WANG. All Rights Reserved.
# Standard library imports
import os

# Third party imports
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# Local application import
from physconst import *
from functions import *



class Datajungle:
    ''' Parent Class for Generic data type
    Auguments:
    directory: list of filenames
    step: step values
    ucols: extracted columns from original source file
    spr: skipped rows in the header of source file 
    ref: reference resistance in series
    AspRatio: aspect ratio of Hall bar, set to 3 by default

    Return:

    CHILDREN CLASS:
    Databs, Datags'''
    def __init__(self,directory,step,ucols,nms,spr,ref,AspRatio=3):
        self.dir = directory
        self.step = step
        self.AspRatio = AspRatio
        self.ucols = ucols
        self.nms = nms
        self.spr = spr
        self.ref = ref


class Databs(Datajungle):
    '''Inherent from Class Datajungle
    METHODS:
    getdata: return a panda.DataFrame type data
    hallfit: linear hall fit and return density and mobility in a DataFrame format
    plotdata: plot magnetic field sweep type data in a specific way'''
    def __str__(self):
        return ', '.join(['Databs', ', '.join(['{key} = {value}'.format(key = key, value = self.__dict__[key]) for key in ['AspRatio','ucols','nms','spr','ref']])])

    def __repr__(self):
        pass


    def getdata(self):
        databundle = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['gate'] = self.step[i]
            databundle = databundle.append(data)
        return databundle

    def hallfit(self,fitrange):
        Dens = []
        Mob = []
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
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

    def plotdata(self,label_value = '$U_g$={:02.2f}V' ):
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
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
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
    plotdata: plot gate sweep type data in a specific way'''
    def __str__(self):
        return ', '.join(['Datags', ', '.join(['{key} = {value}'.format(key = key, value = self.__dict__[key]) for key in ['AspRatio','ucols','nms','spr','ref']])])

    def __repr__(self):
        pass


    def getdata(self):
        databundle = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['bf'] = self.step[i]
            databundle = databundle.append(data)
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
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            ax_rxx.plot(data.gate,data['rxx'],color = line_color,label=label_value.format(self.step[i]))
            ax_rxy.plot(data.gate,data['rxy'],color = line_color,label=label_value.format(self.step[i]))
            ax_sxy.plot(data.gate,data['sxy'],color = line_color,label=label_value.format(self.step[i]))
        ax_rxx.set_xlabel(r'$U_g(V)$',fontsize = 18)
        ax_rxx.set_ylabel(r'$R_{xx}(\Omega)$',fontsize = 18)
        ax_rxy.set_xlabel(r'$U_g(V)$',fontsize = 18)
        ax_rxy.set_ylabel(r'$R_{xy}(\Omega)$',fontsize = 18)
        ax_sxy.set_xlabel(r'$U_g(V)$',fontsize = 18)
        ax_sxy.set_ylabel(r'$\sigma_{xy}(e^2/h)$',fontsize = 18)
        ax_sxy.set_ylim([-10,10])
        ax_rxy.legend(loc = 'best')
        for mark in range(-5,6):
                   ax_sxy.axhline(y=mark,linestyle=':',color='c')
        for mark in [-5,-4,-3,-2,-1,1,2,3,4,5]:
                   ax_rxy.axhline(y=h0/e0**2/mark,linestyle=':',color='c')
        return ax_rxx,ax_rxy,ax_sxy

class Datafc(Datajungle):
    '''Inherent from Class Datajungle
    METHODS:
    plotfc: plot fanchart
    getdata: return x, y and a 2D array with z-value
    plotbs: extract magnetic field sweeps
    plotgs: plot single sxy vs gate sweeps'''
    def __str__(self):
        return ', '.join(['Datafc', ', '.join(['{key} = {value}'.format(key = key, value = self.__dict__[key]) for key in ['AspRatio','ucols','nms','spr','ref']])])

    def __repr__(self):
        pass



    def plotfc(self,vm,vmx,cmap='inferno'):
        diffsxy2D = pd.DataFrame()
        databundle = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['diffsxy'] = data['sxy'].diff()/(data['gate'][0]-data['gate'][1])
            data['bf'] = self.step[i]
            diffsxy2D = diffsxy2D.append(data['diffsxy'].dropna())
            databundle = databundle.append(data)
        x = databundle['gate'].unique()

        y = self.step
        fig = plt.figure(figsize=(15,8))
        ax1 = fig.add_subplot(111)
        plt.imshow(diffsxy2D,aspect='auto', interpolation='none',extent=extents(x.tolist()) + extents(y), origin='lower',cmap=cmap,vmin=vm, vmax=vmx)
        ax1.set_ylabel('B(T)')
        ax1.set_xlabel('$U_{tg}(V)$')
        return ax1



    def getdata(self):
        databundle = pd.DataFrame()
        diffsxy2D = pd.DataFrame()
        rxx2D = pd.DataFrame()
        rxy2D = pd.DataFrame()
        sxy2D = pd.DataFrame()
        sxx2D = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['diffsxy'] = data['sxy'].diff()/(data['gate'][0]-data['gate'][1])
            data['bf'] = self.step[i]
            diffsxy2D = diffsxy2D.append(data['diffsxy'].dropna())
            rxx2D = rxx2D.append(data['rxx'])
            rxy2D = rxy2D.append(data['rxy'])
            sxy2D = sxy2D.append(data['sxy'])
            sxx2D = sxx2D.append(data['sxx'])
            databundle = databundle.append(data)
        datafc = {'x':databundle['gate'].unique(),'y':self.step,'z1':diffsxy2D,'z2':rxx2D,'z3':rxy2D,'z4':sxy2D,'z5':sxx2D}
        return datafc, databundle

    def plotbs(self,gate_list):
        databundle = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        font = {'family' : 'normal','weight' : 'normal','size'   : 15}
        matplotlib.rc('font', **font)
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['bf'] = self.step[i]
            databundle = databundle.append(data)
        fig = plt.figure(figsize=(10,10))
        ax1 = plt.subplot(2,1,1)
        ax2 = plt.subplot(2,1,2)
        jet= plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0,1,len(gate_list))))
        for i in range(len(gate_list)):
            line_color = next(colors)
            x = databundle['bf'].unique()
            y1 = databundle[databundle['gate']==gate_list[i]].rxx
            y2 = databundle[databundle['gate']==gate_list[i]].rxy
            ax1.plot(x,y1,color=line_color,label = '$U_g$ = {:02.2f}V'.format(gate_list[i]))
            ax2.plot(x,y2,color=line_color,label = '$U_g$ = {:02.2f}V'.format(gate_list[i]))
        ax1.set_xlabel(r'B(T)')
        ax1.set_ylabel(r'$R_{xx}(\Omega)$')
        ax1.legend()
        ax2.set_xlabel(r'B(T)')
        ax2.set_ylabel(r'$R_{xy}(\Omega)$')
        ax2.legend()
        return ax1,ax2

    def plotsxy(self,b_list):
        databundle = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        font = {'family' : 'normal','weight' : 'normal','size'   : 15}
        matplotlib.rc('font', **font)
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/AspRatio/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['bf'] = self.step[i]
            databundle = databundle.append(data)
        fig = plt.figure(figsize=(12,5))
        ax1 = plt.subplot(1,1,1)
        jet= plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0,1,len(b_list))))
        for i in range(len(b_list)):
            line_color = next(colors)
            x = databundle['gate'].unique()
            y1 = databundle[databundle['bf']==b_list[i]].sxy/e0**2*h0
            ax1.plot(x,y1,color=line_color,label = 'B= {:02.2f}T'.format(b_list[i]))
        ax1.set_xlabel(r'$U_g(V)$')
        ax1.set_ylabel(r'$\sigma_{xy}(e^2/h)$')
        ax1.legend()
        return ax1



def main():
    ''' main function '''
    print('running SciData in main program')

if __name__ == '__main__':
    main()
