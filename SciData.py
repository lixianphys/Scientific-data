import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pd_Hallfit import H1st_ft

e0 = 1.6021766208E-19 # The elementary charge
h0 = 6.62607015E-34 # The Planck's constant

def extents(f):
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

def plot_esh(handle,list_esh):
    '''
    handle: 'matplotlib.axes._subplots.AxesSubplot'
    list_esh: list[int]
    no return
    '''
    for esh in list_esh:
            handle.axhline(y=esh,linestyle=':',color='c')

def dir2fnm(directory):
    import glob
    import os
    filename = list(filter(os.path.isfile,glob.glob(os.path.join(directory,'*.dat'))))
    filename.sort(key=lambda x:os.path.getmtime(x))
    return filename

def interp_user(x,y,n_interp):
#     PARAM:
#     x: list[float]
#     y: list[float]
#     n_interp : int
#     RETURN:
#     x_vals: equally spaced sequence (x_vals_i)
#     yinterp: interpolated y_i(x_vals_i) sequence
    import numpy as np
    xmin = min(x)
    xmax = max(x)
    x.sort() # make sure your x is sorted in ascending order
    x_vals = np.linspace(xmin,xmax,n_interp)
    yinterp = np.interp(x_vals,x,y)
    return x_vals,yinterp

# FFT analysis
def FFT_bs(x,y):
#    PARAM:
#    x: time series
#    y: time-dependent variable
#    RETURN:
#    frq: frequency
#    Y: FFT outputed frequency series
    Y = np.fft.fft(y)/len(y)
    Y = Y[np.arange(len(y)/2,dtype = int)]
    k = np.arange(len(y)/2)
    FS = len(y)/(max(x)-min(x))
    frq = k*FS/len(y)
    return frq,Y

def diffz_df(dataframe,axes,z_vec):
    ''' Parameters:
    ===================
        dataframe (pd.DataFrame): operated object
        axes (list): [key for another axis in plots,key for the diff axis]
        z_vec(str) key for z directional z

        Returns:
    ====================
        a_rest(list): vector in rest axis
        a_diff(list): vector in diff axis
        z_array: yielded z array in shape of (len(a_rest),len(a_diff)-1)
        z_df: same content with z_array but in DataFrame format (useful in recursive calls)
    '''
    after_diff = dataframe.sort_values(by=axes).diff()
    rest_dim = axes[0]
    diff_dim = axes[1]

    a_rest = sorted(dataframe[rest_dim].unique())
    a_diff = sorted(dataframe[diff_dim].unique())

    z_values = after_diff[z_vec].tolist()
    z_array = np.zeros([len(a_rest),len(a_diff[1:])])
    z_list = []
    for x_i,x in enumerate(a_rest):
        for y_i,y in enumerate(a_diff[1:]):
                z_array[x_i,y_i] = z_values[len(a_diff)*x_i+y_i+1]
                z_dict = {rest_dim:x,diff_dim:y,z_vec:z_values[len(a_diff)*x_i+y_i+1]}
                z_list.append(z_dict)
    z_df = pd.DataFrame(z_list)
    print('The output array is in shape {}\nwith x_ass of length of {} and y_diff of length of {}'.format(z_array.shape,len(x_ass),len(y_diff)))
    return a_rest,a_diff,z_array,z_df

def fc_interp(x_vec,y_vec,z_df,mult_factor=3):
    ''' Parameters:
    ===================
        x_vec (list): vector in x axis
        y_vec (list): vector in y axis
        z_df(DataFrame): DataFrame data to be interpolated
        mult_factor: determine the intensity of interpolation (default=3)
        Returns:
    ====================
        grid_z: yielded z array in shape of (len(x_vec)*mult_factor,len(y_vec)*mult_factor)
    '''
    from scipy.interpolate import griddata
    values = z_df.values
    points = np.zeros((len(values),2))
    for x_i,x in enumerate(x_vec):
        for y_i,y in enumerate(y_vec):
            points[x_i*len(y_vec)+y_i,0] = x
            points[x_i*len(y_vec)+y_i,1] = y
    # grids for interpolation
    grid_x,grid_y = np.mgrid[x_vec[0]:x_vec[-1]:complex(0,len(x_vec)*mult_factor), y_vec[0]:y_vec[-1]:complex(0,len(y_vec)*mult_factor)]
    # interpolation
    grid_z = griddata(points, values, (grid_x, grid_y), method='nearest')
    return grid_z

class Datajungle:
    ''' Parent Class for Generic data type
    ATTRIBUTES:
    directory: list of file names
    step: step values
    ucols and spr : usecols and skiprows for panda.DataFrame type
    ref: reference resistance in series
    AspRatio: aspect ratio of Hall bar, set to 3 by default
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


    # @abstractmethod
    # def getdata(cls):

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
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
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
            dens,mob = H1st_ft(bf_fit,rxx_fit,rxy_fit)
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
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            ax_rxx.plot(data.bf,data['rxx'],color = line_color,label=label_value.format(self.step[i]))
            ax_rxy.plot(data.bf,data['rxy'],color = line_color,label=label_value.format(self.step[i]))
            ax_sxy.plot(data.bf,data['sxy'],color = line_color,label=label_value.format(self.step[i]))
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
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
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
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)/e0**2*h0
            ax_rxx.plot(data.gate,data['rxx'],color = line_color,label=label_value.format(self.step[i]))
            ax_rxy.plot(data.gate,data['rxy'],color = line_color,label=label_value.format(self.step[i]))
            ax_sxy.plot(data.gate,data['sxy'],color = line_color,label=label_value.format(self.step[i]))
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



    def plotfc(self,vm,vmx):
        diffsxy2D = pd.DataFrame()
        databundle = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['diffsxy'] = data['sxy'].diff()/(data['gate'][0]-data['gate'][1])
            data['bf'] = self.step[i]
            diffsxy2D = diffsxy2D.append(data['diffsxy'].dropna())
            databundle = databundle.append(data)
        x = databundle['gate'].unique()

        y = self.step
        fig = plt.figure(figsize=(15,8))
        ax1 = fig.add_subplot(111)
        plt.imshow(diffsxy2D,aspect='auto', interpolation='bilinear',extent=extents(x.tolist()) + extents(y), origin='lower',cmap='jet',vmin=vm, vmax=vmx)
        ax1.set_ylabel('B(T)')
        ax1.set_xlabel('$U_{tg}(V)$')
        return ax1

    def getdata(self):
        databundle = pd.DataFrame()
        diffsxy2D = pd.DataFrame()
        rxx2D = pd.DataFrame()
        ref = self.ref
        AspRatio = self.AspRatio
        for i in range(len(self.dir)):
            data = pd.read_csv(self.dir[i], sep="\t",skiprows=self.spr, usecols=self.ucols, names=self.nms, header=None)
            data['rxx'] = data.uxx/data.curr*ref
            data['rxy'] = data.uxy/data.curr*ref
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['sxy'] = data['rxy']/((data['rxx']/AspRatio)**2+data['rxy']**2)
            data['diffsxy'] = data['sxy'].diff()/(data['gate'][0]-data['gate'][1])
            data['bf'] = self.step[i]
            diffsxy2D = diffsxy2D.append(data['diffsxy'].dropna())
            rxx2D = rxx2D.append(data['rxx'])
            databundle = databundle.append(data)
        datafc = {'x':databundle['gate'].unique(),'y':self.step,'z1':diffsxy2D,'z2':rxx2D}
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
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)
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
        ax1.set_xlabel('B(T)')
        ax1.set_ylabel('$R_{xx}(\Omega)$')
        ax1.legend()
        ax2.set_xlabel('B(T)')
        ax2.set_ylabel('$R_{xy}(\Omega)$')
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
            data['sxx'] = data['rxx']/((data['rxx']/AspRatio)**2+data['rxy']**2)
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
        plot_esh(ax1,range(-10,11))
        ax1.set_xlabel('$U_g(V)$')
        ax1.set_ylabel('$\sigma_{xy}(e^2/h)$')
        ax1.legend()
        return ax1



def main():
    ''' main function '''
    print('running SciData in main program')

if __name__ == '__main__':
    main()
