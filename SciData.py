import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pd_Hallfit import H1st_ft
e0 = 1.6021766208E-19 # The elementary charge
h0 = 6.62607015E-34 # The Planck's constant

class Datajungle:
    ''' Parent Class for Generic data type
    ATTRIBUTES:
    directory: list of file names
    step: step values
    ucols and spr : usecols and skiprows for panda.DataFrame type
    Ref: reference resistance in series
    AspRatio: aspect ratio of Hall bar, set to 3 by default
    CHILDREN CLASS:
    Databs, Datags'''
    def __init__(self,directory,step,ucols,spr,Ref,AspRatio=3):
        self.dir = directory
        self.step = step
        self.AspRatio = AspRatio
        self.ucols = ucols
        self.spr = spr
        self.ref = Ref

        
class Databs(Datajungle):
    '''Inherent from Class Datajungle
    METHODS:
    getdata: return a panda.DataFrame type data
    hallfit: linear hall fit and return density and mobility in a DataFrame format
    plotdata: plot magnetic field sweep type data in a specific way'''
    def getdata(self,nms = ['bf','curr','uxx','uxy']):
        databundle = pd.DataFrame()
        Ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['gate'] = self.step[i]
            databundle = databundle.append(data)
        return databundle

    def hallfit(self,fitrange,nms = ['bf','curr','uxx','uxy']):
        Dens = []
        Mob = []
        Ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            bf_fit = data['bf'][(data['bf']<fitrange[1])&(data['bf']>fitrange[0])]
            rxx_fit = data['rxx'][(data['bf']<fitrange[1])&(data['bf']>fitrange[0])]
            rxy_fit = data['rxy'][(data['bf']<fitrange[1])&(data['bf']>fitrange[0])]
            dens,mob = H1st_ft(bf_fit,rxx_fit,rxy_fit)
            Dens.append(dens)
            Mob.append(mob)
        FitRes = pd.DataFrame({'gate':self.step,'dens':Dens,'mob':Mob})
        return FitRes
    
    def plotdata(self,nms = ['bf','curr','uxx','uxy']):
        font = {'family' : 'normal','weight' : 'normal','size'   : 15}
        matplotlib.rc('font', **font)
        jet= plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0,1,len(self.dir))))
        plt.figure(figsize=(16,16))
        ax_rxy = plt.subplot(2,1,1)
        ax_rxx = plt.subplot(2,2,3)
        ax_sxy = plt.subplot(2,2,4)
        Ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            line_color = next(colors)
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            ax_rxx.plot(data.bf,data['rxx'],color = line_color,label='$U_g$={:02.2f}V'.format(self.step[i]))
            ax_rxy.plot(data.bf,data['rxy'],color = line_color,label='$U_g$={:02.2f}V'.format(self.step[i]))
            ax_sxy.plot(data.bf,data['sxy'],color = line_color,label='$U_g$={:02.2f}V'.format(self.step[i]))
        ax_rxx.set_xlabel('$B_{field}(T)$',fontsize = 18)
        ax_rxx.set_ylabel('$R_{xx}(\Omega)$',fontsize = 18)
        ax_rxy.set_xlabel('$B_{field}(T)$',fontsize = 18)
        ax_rxy.set_ylabel('$R_{xy}(\Omega)$',fontsize = 18)
        ax_sxy.set_xlabel('$B_{field}(T)$',fontsize = 18)
        ax_sxy.set_ylabel('$\sigma_{xy}(e^2/h)$',fontsize = 18)
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
    def getdata(self,nms = ['gate','curr','uxx','uxy']):
        databundle = pd.DataFrame()s
        Ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['bf'] = self.step[i]
            databundle = databundle.append(data)
        return databundle
    def plotdata(self,nms = ['gate','curr','uxx','uxy']):
        font = {'family' : 'normal','weight' : 'normal','size' : 15}
        matplotlib.rc('font', **font)
        jet= plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0,1,len(self.dir))))
        plt.figure(figsize=(16,16))
        ax_rxy = plt.subplot(2,1,1)
        ax_rxx = plt.subplot(2,2,3)
        ax_sxy = plt.subplot(2,2,4)
        Ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            line_color = next(colors)
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            ax_rxx.plot(data.gate,data['rxx'],color = line_color,label='$B$={:02.2f}T'.format(self.step[i]))
            ax_rxy.plot(data.gate,data['rxy'],color = line_color,label='$B$={:02.2f}T'.format(self.step[i]))
            ax_sxy.plot(data.gate,data['sxy'],color = line_color,label='$B$={:02.2f}T'.format(self.step[i]))
        ax_rxx.set_xlabel('$U_g(V)$',fontsize = 18)
        ax_rxx.set_ylabel('$R_{xx}(\Omega)$',fontsize = 18)
        ax_rxy.set_xlabel('$U_g(V)$',fontsize = 18)
        ax_rxy.set_ylabel('$R_{xy}(\Omega)$',fontsize = 18)
        ax_sxy.set_xlabel('$U_g(V)$',fontsize = 18)
        ax_sxy.set_ylabel('$\sigma_{xy}(e^2/h)$',fontsize = 18)
        ax_sxy.set_ylim([-10,10])
        ax_rxy.legend(loc = 'best')
        for mark in range(-5,6):
                   ax_sxy.axhline(y=mark,linestyle=':',color='c')
        for mark in [-5,-4,-3,-2,-1,1,2,3,4,5]:
                   ax_rxy.axhline(y=h0/e0**2/mark,linestyle=':',color='c')
        return ax_rxx,ax_rxy,ax_sxy

