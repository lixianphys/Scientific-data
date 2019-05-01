#!/usr/bin/env python
# coding: utf-8

# In[7]:


class Datajungle:
    def __init__(self,directory,step,ucols,spr,AspRatio=3):
        self.dir = directory
        self.step = step
        self.AspRatio = AspRatio
        self.ucols = ucols
        self.spr = spr
    
class Databs(Datajungle):
    def getdata(self,nms = ['bf','curr','uxx','uxy']):
        databundle = panda.DataFrame()
        for i in range(len(directory)):
            data = pd.read_csv(directory[i], sep="\t",skiprows=spr, usecols=ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['gate'] = step[i]
            databundle = databundle.append(data)
        return databundle

    def hallfit(self,fitrange):
        Dens = []
        Mob = []
        for i in range(len(directory)):
            data = pd.read_csv(directory[i], sep="\t",skiprows=spr, usecols=ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            bf_fit = Ddata['bf'][(Ddata['bf']<fitrange[1])&(Ddata['bf']>fitrange[0])]
            rxx_fit = Ddata['rxx'][(Ddata['bf']<fitrange[1])&(Ddata['bf']>fitrange[0])]
            rxy_fit = Ddata['rxy'][(Ddata['bf']<fitrange[1])&(Ddata['bf']>fitrange[0])]
            dens,mob = H1st_ft(bf_fit,rxx_fit,rxy_fit)
            Dens = Dens.append(dens)
            Mob = Mob.append(mob)
        FitRes = panda.DataFrame({'gate':step,'dens':Dens,'mob':Mob})
        return FitRes
    
    def plotdata(self):
        font = {'family' : 'normal','weight' : 'normal','size'   : 15}
        matplotlib.rc('font', **font)
        jet= plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0,1,len(directory))))
        plt.figure(figsize=(16,16))
        ax_rxy = plt.subplot(2,1,1)
        ax_rxx = plt.subplot(2,2,3)
        ax_sxy = plt.subplot(2,2,4)
        for i in range(len(directory)):
            line_color = next(colors)
            data = pd.read_csv(directory[i], sep="\t",skiprows=spr, usecols=ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            ax_rxx.plot(Ddata.bf,data['rxx'],color = line_color)
            ax_rxy.plot(Ddata.bf,data['rxy'],color = line_color)
            ax_sxy.plot(Ddata.bf,data['sxy'],color = line_color)
        ax_rxx.set_xlabel('$B_{field}(T)$',fontsize = 18)
        ax_rxx.set_ylabel('$R_{xx}(\Omega)$',fontsize = 18)
        ax_rxy.set_xlabel('$B_{field}(T)$',fontsize = 18)
        ax_rxy.set_ylabel('$R_{xy}(\Omega)$',fontsize = 18)
        ax_sxy.set_xlabel('$B_{field}(T)$',fontsize = 18)
        ax_sxy.set_ylabel('$\sigma_{xy}(e^2/h)$',fontsize = 18)
        ax_sxy.set_ylim([-10,10])
        for mark in np.linspace(3,10,8):
                   ax_sxy.axhline(y=mark,linestyle=':',color='c')
                   ax_rxy.axhline(y=h0/e0**2/mark,linestyle=':',color='c')
        plt.show()
        return ax_rxx,ax_rxy,ax_sxy
        

    
class Datags(Datajungle):
    def getdata(self,ucols,spr,AspRatio,nms = ['gate','curr','uxx','uxy']):
        databundle = panda.DataFrame()
        for i in range(len(directory)):
            data = pd.read_csv(directory[i], sep="\t",skiprows=spr, usecols=ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['bf'] = step[i]
            databundle = databundle.append(data)
        return databundle
    def plotdata(self):
        font = {'family' : 'normal','weight' : 'normal','size' : 15}
        matplotlib.rc('font', **font)
        jet= plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0,1,len(directory))))
        plt.figure(figsize=(16,16))
        ax_rxy = plt.subplot(2,1,1)
        ax_rxx = plt.subplot(2,2,3)
        ax_sxy = plt.subplot(2,2,4)
        for i in range(len(directory)):
            line_color = next(colors)
            data = pd.read_csv(directory[i], sep="\t",skiprows=spr, usecols=ucols, names=nms, header=None)
            data['rxx'] = data.uxx/data.curr*Ref
            data['rxy'] = data.uxy/data.curr*Ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            ax_rxx.plot(data.gate,data['rxx'],color = line_color)
            ax_rxy.plot(data.gate,data['rxy'],color = line_color)
            ax_sxy.plot(data.gate,data['sxy'],color = line_color)
        ax_rxx.set_xlabel('$U_g(V)$',fontsize = 18)
        ax_rxx.set_ylabel('$R_{xx}(\Omega)$',fontsize = 18)
        ax_rxy.set_xlabel('$U_g(V)$',fontsize = 18)
        ax_rxy.set_ylabel('$R_{xy}(\Omega)$',fontsize = 18)
        ax_sxy.set_xlabel('$U_g(V)$',fontsize = 18)
        ax_sxy.set_ylabel('$\sigma_{xy}(e^2/h)$',fontsize = 18)
        ax_sxy.set_ylim([-10,10])
        for mark in np.linspace(3,10,8):
                   ax_sxy.axhline(y=mark,linestyle=':',color='c')
                   ax_rxy.axhline(y=h0/e0**2/mark,linestyle=':',color='c')
        plt.show()
        return ax_rxx,ax_rxy,ax_sxy


# In[ ]:




