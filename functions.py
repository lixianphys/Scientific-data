# Copyright 2021 Lixian WANG. All Rights Reserved.
# Standard library imports
import os
import re

# Third party imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
# Local application import
from physconst import *
# General use

def getnumber(fnm):
    fnm_strip = fnm.strip('.dat').split('/')[-1] if '/' in fnm else fnm.strip('.dat').split('\\')[-1]
    # num_str = fnm_strip.replace('T', '').replace(',', '.').replace('DCV', '').replace('B-Field', '').replace('V', '').replace('m','-').replace('p','.').replace('(01)','')
    num_str = re.sub(r"DCV|set|B-Field|Heater|Ch4|power|%|V|T|\(\d\d\)", "", fnm_strip.split('_')[-1])
    num_str = re.sub(r"(,)|(p)", ".", num_str)
    num_str = re.sub(r"(m)", "-", num_str)
    num = float(num_str)
    return num


def dir2fnm(directory, sort_by_fnm=False):
    # TODO: add a custmizable feature in string pattern recognization
    '''
    Convert a directory to a list of filenames contained inside
    :param directory: directory
    :return: a list of filenames sorted by the creation time in ascend order
    '''
    import glob
    import os
    filename = list(filter(os.path.isfile, glob.glob(os.path.join(directory, '*.dat'))))

    if sort_by_fnm:
        filename.sort(key=lambda x: getnumber(x))
    else:
        filename.sort(key=lambda x: os.path.getmtime(x))
    return filename


def read_file(directory, sort_by_fnm=False):
    '''
    Extract the numbers in the filenames of a batch of files and output them in a list
    '''
    filenames = dir2fnm(directory, sort_by_fnm=sort_by_fnm)
    num_list = [0] * len(filenames)
    for i, fnm in enumerate(filenames):
        num_list[i] = float(getnumber(fnm))
    return num_list


def pos_neg(num):
    if num > 0:
        return 1
    else:
        return -1

def is_close(num_list: list, match_num: float, precision=1e-6) -> bool:
    '''
    Target the nearest number in list of numbers (num_list)
    '''
    return [abs(num - match_num) < precision for num in num_list]


def df_range(df, column, col_range):
    return df[(df[column] > col_range[0]) & (df[column] < col_range[1])]


def range_pick(yourlist, lb, ub):
    lb_set = yourlist > lb
    ub_set = yourlist < ub
    chosen_set = [all([x, y]) for x, y in zip(lb_set, ub_set)]
    return chosen_set


# Calculation

def H1st_ft(Bf,Rxx,Rxy,AspRatio=3,threshold = 25, fitpara_output=False):

    '''
    Linear fit model for Hall analysis

    :param Bf:
    :param Rxx:
    :param Rxy:
    :param AspRatio:
    :param threshold:
    :return:
    '''
    def func_one(x, a, b):
        return a + b * x
    e0 = 1.6021766208E-19
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
        else:
            print('The fitting results is not acceptable, fitCov is {}'.format(fitCovariances))
            plt.plot(Bf, Rxy, "b-", Bf, func_one(Bf, *fitParams), "r-")
            mobility = 0
            density = 0
    if fitpara_output:
        return density,mobility,fitParams
    else:
        return density,mobility


def H2nd_ft(Bf, Rxx, Rxy, AspRatio=3):
    '''
    Two carrier (electron-hole) model for Hall analysis

    :param Bf:
    :param Rxx:
    :param Rxy:
    :param AspRatio:
    :return:
    '''
    e0 = 1.6021766208E-19

    def func_two(x, n1, m1, n2, m2):
        return e0 * x * (n1 * m1 ** 2 / (1 + m1 ** 2 * x ** 2) + n2 * m2 ** 2 / (
                    1 + m2 ** 2 * x ** 2))  # model from PHYSICAL REVIEW B 95, 115126 (2017)

    sxy = Rxy / ((Rxx / AspRatio) ** 2 + Rxy ** 2)
    try:
        popt, pcov = curve_fit(func_two, Bf, sxy, bounds=((1e14, 5, -1e16, 0), (5e15, 15, -1e14, 3)))
        # plt.plot(Bf, sxy, "b-", Bf, func_two(Bf, *popt), "r-")
        return popt, func_two(Bf, *popt)
    except:
        print('The fitting program failed')

def twocarrierfit(Bf, Rxy):
    '''
    Two carrier (electron-hole) model for Hall analysis
    :param Bf:
    :param Rxy:
    :return:
    '''
    e0 = 1.6021766208E-19

    def func(x, n1, m1, n2, m2):
        return  -((n2*m2**2-n1*m1**2)+m2**2*m1**2*x**2*(n2-n1))*x/e0/((n2*m2+n1*m1)**2+m2**2*m1**2*x**2*(n2-n1)**2)  # Reference: Li, Cai-Zhen, et al. ACS nano 10.6 (2016): 6020-6028.

    try:
        popt, pcov = curve_fit(func, Bf, Rxy, bounds=((1e14, 5, 1e14, 0), (5e15, 15, 1e16, 3)))
        # plt.plot(Bf, sxy, "b-", Bf, func_two(Bf, *popt), "r-")
        return popt, func(Bf, *popt)
    except:
        print('The two carrier fit failed')


def cutout_bkgd(x, y):
    '''
    Remove the smooth background of a function y = f(x) by subtracting a polynomial function matching the shape of f(x)
    properly

    :param x: variable
    :param y: function

    :return: the y values after removing the background
    '''
    from scipy.optimize import curve_fit

    def func(x, a, b, c, d, e, f, g):
        return a + b * x + c * x ** 2 + d * x ** 3 + e * x ** 4 + f * x ** 5 + g * x ** 6

    try:
        fit_params, _ = curve_fit(func, x, y)
        y_bkgd = fit_params[0] + fit_params[1] * x + fit_params[2] * x ** 2 + fit_params[3] * x ** 3 + fit_params[
            4] * x ** 4 + fit_params[5] * x ** 5 + fit_params[6] * x ** 6
        y_signal = y - y_bkgd
    except:
        y_signal = None
        print('In called cutout_bkgd function, the processing of polynomial fit failed, the return is None')
    return y_signal


def interp_user(x, y, n_interp):
    '''
    :param x: list[float]
    :param y: list[float]
    :param n_interp: int

    :return:
    :x_vals: equally spaced sequence (x_vals_i)
    :yinterp: interpolated y_i(x_vals_i) sequence
    '''
    xmin = min(x)
    xmax = max(x)
    ziplist = zip(x, y)
    sortedziplist = sorted(ziplist)
    x_ascd = [element for element, _ in sortedziplist]  # make sure your x is sorted in ascending order
    y_ascd = [element for _, element in sortedziplist]
    x_vals = np.linspace(xmin, xmax, n_interp)
    yinterp = np.interp(x_vals, x_ascd, y_ascd)
    return x_vals, yinterp


def FFT_bs(x, y):
    '''
    Get the FFT result of a function y = f(x)

    :param x: time series
    :param y: time-dependent variable

    :return:
    :frq: frequency
    :Y: FFT amplitude
    '''

    Y = np.fft.fft(y) / len(y)
    Y = Y[np.arange(len(y) / 2, dtype=int)]
    k = np.arange(len(y) / 2)
    FS = len(y) / (max(x) - min(x))
    frq = k * FS / len(y)
    return frq, Y


def diffz_df(dataframe, axes, z_vec, check_output=False):
    '''
    Perform a one-dimensional diff operation on a 2d data.

    :param dataframe: original 2d data
    :param axes: [key for another axis in plots,key for the diff axis]
    :param z_vec: key for z axis
    :param check_output: print the

    :return:
    :a_rest(list): vector in rest axis
    :a_diff(list): vector in diff axis
    :z_array: yielded z array in shape of (len(a_rest),len(a_diff)-1)
    :z_df: same content with z_array but in DataFrame format (useful in recursive calls)
    '''

    after_diff = dataframe.sort_values(by=axes).diff()
    rest_dim = axes[0]
    diff_dim = axes[1]

    ax_rest = sorted(dataframe[rest_dim].unique())  # get the x-vector
    ax_diff = sorted(dataframe[diff_dim].unique())  # get the y-vector

    z_values = after_diff[z_vec].tolist()
    z_array = np.zeros([len(ax_rest), len(ax_diff[1:])])
    z_list = []
    for x_i, x in enumerate(ax_rest):
        for y_i, y in enumerate(ax_diff[1:]):
            z_array[x_i, y_i] = z_values[len(ax_diff) * x_i + y_i + 1]
            z_dict = {rest_dim: x, diff_dim: y, z_vec: z_values[len(ax_diff) * x_i + y_i + 1]}
            z_list.append(z_dict)
    z_df = pd.DataFrame(z_list)
    if check_output:
        print('The output array is in shape {}\nwith ax_rest of length of {} and ax_diff of length of {}'.format(
            z_array.shape, len(ax_rest), len(ax_diff)))
    else:
        pass
    return ax_rest, ax_diff, z_array, z_df


def fc_interp(x_vec, y_vec, z_df, diff=True, mult_factor=3):
    '''
    Interpolate 2d data z_df

    :param x_vec: vector in x axis
    :param y_vec: vector in y axis
    :param z_df: DataFrame data to be interpolated
    :param diff: True=interpolate 1st differential data (size-1), False = interpolate normal size data
    :param mult_factor: determine how dense the interpolation could be performed. [size of output] = mult_factor*[size of input] (default=3)

    :return:
    :grid_z: yielded z array in shape of (len(x_vec)*mult_factor,len(y_vec)*mult_factor)

    '''

    from scipy.interpolate import griddata
    values = z_df.values
    points = np.zeros((len(values), 2))
    if diff:
        for x_i, x in enumerate(x_vec):
            for y_i, y in enumerate(y_vec[:-2], 1):
                points[x_i * (len(y_vec) - 1) + y_i, 0] = x
                points[x_i * (len(y_vec) - 1) + y_i, 1] = y
    else:
        for x_i, x in enumerate(x_vec):
            for y_i, y in enumerate(y_vec):
                points[x_i * len(y_vec) + y_i, 0] = x
                points[x_i * len(y_vec) + y_i, 1] = y
    # grids for interpolation
    grid_x, grid_y = np.mgrid[x_vec[0]:x_vec[-1]:complex(0, len(x_vec) * mult_factor),
                     y_vec[0]:y_vec[-1]:complex(0, len(y_vec) * mult_factor)]
    # interpolation
    grid_z = griddata(points, values, (grid_x, grid_y), method='nearest')
    return grid_z


# Plotting

def quickplot(path, num_plot, PhyQty, ref, skiprows, nms, ucols, AspRatio=3):
    '''
    Quick plot for multiple files containing the same type of data
    Arguments:
    path: The directory of multiple files
    PhyQty: The physical quantities to be plotted('bf','gate','rxx','rxy','sxx','sxy' and etc)
    ref: reference resistor
    skiprows: skipped rows from the header
    nms: Names for all used columns
    ucols: Used columns
    AspRatio: The aspect ratio of the Hall bar. Default is 3

    Return:
    the handle of axes to facilitate further adjustment if necessary
    '''

    fig = plt.figure(figsize=(10, 5 * len(PhyQty)))
    fnm = dir2fnm(path)
    jet = plt.get_cmap('jet')
    colors = iter(jet(np.linspace(0, 1, len(fnm))))
    if len(PhyQty) > 1:
        gs = fig.add_gridspec(len(PhyQty), 1)
        plots_ax = [fig.add_subplot(x) for x in gs]
    else:
        plots_ax = fig.add_subplot(111)

    for file in fnm:
        color = next(colors)
        data = pd.read_csv(file, sep="\t", skiprows=skiprows, usecols=ucols, names=nms, header=None)
        data['rxx'] = data.uxx / data.curr * ref
        data['rxy'] = data.uxy / data.curr * ref
        data['sxx'] = data['rxx'] / AspRatio / ((data['rxx'] / AspRatio) ** 2 + data['rxy'] ** 2) / e0 ** 2 * h0
        data['sxy'] = data['rxy'] / ((data['rxx'] / AspRatio) ** 2 + data['rxy'] ** 2) / e0 ** 2 * h0
        if len(PhyQty) > 1:
            for index, phyqty in enumerate(PhyQty):
                plot_ax = plots_ax[index]
                plot_ax.plot(data.x, data[phyqty], color=color)
        else:
            plots_ax.plot(data.x, data[PhyQty[0]], color=color)
    return plots_ax


def extents(f):
    '''
    Calculate the extent parameter for mapping-like plots:
    extent=[horizontal_min,horizontal_max,vertical_min,vertical_max]
    :param f: list type x-vector or y-vector
    :return:
    '''
    delta = f[1] - f[0]
    return [f[0] - delta / 2, f[-1] + delta / 2]


## interactive plotting

def plot_fc_analysis(datafc, label, vmin=-0.005, vmax=0, equal_spaced=True, bgortg=True, axis_diff='gate', zoom_in=[]):
    from ipywidgets import interactive, FloatSlider, Dropdown
    fc, data = datafc.getdata()
    diffsxy2D = fc['z1']
    x = fc['x']
    y = fc['y']
    bf = data['bf'].apply(lambda x: round(x, 4)).unique()
    bf_step = abs(round(np.mean(np.diff(bf)), 4))
    gates = data['gate'].apply(lambda x: round(x, 4)).unique()
    gate_step = abs(np.mean(np.diff(gates)))
    x_bybf, y_bybf, diffsxy2d_bybf, diffsxy_bybf = diffz_df(data, ['gate', 'bf'], 'sxy')

    def plot_animation(uplim, gate):

        fig = plt.figure(figsize=(15, 12))
        ax1 = plt.subplot2grid((5, 5), (2, 0), colspan=4, rowspan=3)

        if equal_spaced:
            if axis_diff == 'gate':
                ax1.imshow(diffsxy2D, aspect='auto', interpolation='none',
                           extent=extents(x.tolist()) + extents(y.tolist()), origin='lower', cmap='inferno', vmin=vmin,
                           vmax=vmax)
            else:
                ax1.imshow(diffsxy2d_bybf.T, aspect='auto', interpolation='none',
                           extent=extents(x_bybf) + extents(y_bybf), origin='lower', cmap='inferno', vmin=vmin,
                           vmax=vmax)
        elif axis_diff == 'gate':
            interp_diffsxy = fc_interp(y, x, data.diffsxy.dropna())
            ax1.imshow(interp_diffsxy, aspect='auto', interpolation='none', extent=extents(x) + extents(y),
                       origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
        else:
            interp_diffsxy = fc_interp(x_bybf, y_bybf, diffsxy_bybf.sxy)
            ax1.imshow(interp_diffsxy.T, aspect='auto', interpolation='none', extent=extents(x_bybf) + extents(y_bybf),
                       origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)

        if bgortg:
            ax1.set_xlabel('$U_{tg}(V)$')
        else:
            ax1.set_xlabel('$U_{bg}(V)$')

        ax1.set_ylabel('B(T)')
        if not zoom_in:
            ax1.set_ylim([min(bf), max(bf)])
            ax1.set_xlim([min(gates), max(gates)])
        else:
            ax1.set_xlim([zoom_in[0], zoom_in[1]])
            ax1.set_ylim([zoom_in[2], zoom_in[3]])
        ax1.axhline(y=uplim, linestyle='--', color='b', linewidth=2)
        ax1.axvline(x=gate, linestyle='--', color='g', linewidth=2)
        ax2 = plt.subplot2grid((5, 5), (0, 0), colspan=4, rowspan=2)
        ax3 = ax2.twinx()

        if equal_spaced:
            if axis_diff == 'gate':
                data_p = df_range(data, 'bf', [uplim - bf_step / 2, uplim + bf_step / 2])
            else:
                data_p = df_range(data, 'bf', [uplim - bf_step / 2, uplim + bf_step / 2])
                data_pbybf = df_range(diffsxy_bybf, 'bf', [uplim - bf_step / 2, uplim + bf_step / 2])

        elif axis_diff == 'gate':
            data_p = df_range(data, 'bf', [uplim - 0.001, uplim + 0.001])
        else:
            data_p = df_range(data, 'bf', [uplim - 0.001, uplim + 0.001])
            data_pbybf = df_range(diffsxy_bybf, 'bf', [uplim - 0.001, uplim + 0.001])

        ax2.plot(data_p.gate, data_p.sxy / e0 ** 2 * h0, 'k-', linewidth=3, label=r'$\sigma_{xy}(e^2/h)$')

        if axis_diff == 'gate':
            ax3.plot(data_p.gate[1:], -data_p.diffsxy[1:], 'r-', linewidth=2, label=r'$-d\sigma_{xy}/dU_{tg}$')
            ax3.set_ylabel(r'$-d\sigma_{xy}/dU_{tg}$')
        else:
            ax3.plot(x_bybf[1:], -data_pbybf.sxy[1:], 'r-', linewidth=2, label=r'$-d\sigma_{xy}/dB$')
            ax3.set_ylabel(r'$-d\sigma_{xy}/dB$')

        ax2.set_xlim([min(gates), max(gates)])
        ax2.set_ylim([-5, 15])
        ax2.legend(bbox_to_anchor=(0, 0.8), loc='center left')
        ax3.legend(bbox_to_anchor=(0.15, 0.8), loc='center left')
        ax2.set_ylabel(r'$\sigma_{xy}/(e^2/h)$')

        ax2.axhline(y=0, linestyle=':', color='c', linewidth=2)
        [ax2.axhline(y=yi, linestyle=':', color='y', linewidth=2) for yi in range(1, 15)]
        [ax2.axhline(y=-yi, linestyle=':', color='g', linewidth=2) for yi in range(1, 15)]
        data_gp = df_range(data, 'gate', [gate - gate_step / 2, gate + gate_step / 2])
        ax4 = plt.subplot2grid((5, 5), (2, 4), colspan=1, rowspan=3)
        ax4.plot(data_gp.rxy, data_gp.bf, 'k-', linewidth=3, label='$r_{xy}$')
        ax4.set_xlim([min(data_gp.rxy) - 500, max(data_gp.rxy) + 500])
        ax4.set_ylim([min(bf), max(bf)])
        ax4.set_xlabel(r'$R_{xy}(\Omega)$')
        [ax4.axvline(x=h0 / e0 ** 2 / xi, linestyle=':', color='y', linewidth=2) for xi in range(1, 10)]
        [ax4.axvline(x=-h0 / e0 ** 2 / xi, linestyle=':', color='g', linewidth=2) for xi in range(1, 10)]
        props = dict(boxstyle='round', fc='b', alpha=0.5)

        if bgortg:
            textstr = ''.join((r'$U_{bg} = $', label))
        else:
            textstr = ''.join((r'$U_{tg} = $', label))

        ax2.text(0.85, 0.9, textstr, transform=ax2.transAxes, fontsize=14, verticalalignment='top', bbox=props,
                 color='w')

    if equal_spaced:
        return interactive(plot_animation,
                           uplim=FloatSlider(min=min(bf), max=max(bf), step=bf_step, continuous_update=False),
                           gate=FloatSlider(min=min(gates), max=max(gates), step=gate_step, continuous_update=False)
                           )
    else:
        return interactive(plot_animation, uplim=Dropdown(
            options=y,
            value=y[1],
            description='field(T):',
            disabled=False,
        ),
                           gate=FloatSlider(min=min(gates), max=max(gates), step=gate_step, continuous_update=False))


def plot_fftmap(datafc, vmin=0, vmax=25, tgorbg=True, bf_range=[0.25, 1]):
    from ipywidgets import interactive, FloatSlider
    _, data = datafc.getdata()
    gates = data['gate'].apply(lambda x: round(x, 3)).unique()
    gate_step = abs(np.mean(np.diff(gates)))
    ## extract pieces of data

    data_p = df_range(data, 'bf', bf_range)
    fft2d = np.zeros([len(gates), (len(data_p) // len(gates) + 1) // 2])
    ## obtain fft2d values in 2d array format
    for index, gate in enumerate(gates):
        data_pp = df_range(data_p, 'gate', [gate - gate_step / 2, gate + gate_step / 2])
        x_vals, yinterp = interp_user(1 / data_pp.bf.values, cutout_bkgd(1. / data_pp.bf.values, data_pp.rxx.values),
                                      len(data_pp.bf))  # interpolation if applicable
        frq, Y = FFT_bs(x_vals, yinterp)
        fft2d[index, :] = (abs(Y) / np.mean(
            abs(Y))) ** 2  # normalized amplitude of fft and the power square is for color coding.

    ## transform frequency into 2D electron/hole density
    n2d = e0 * frq / h0 / 1e15
    x = n2d.tolist()
    y = [round(x, 3) for x in gates.tolist()]

    def plot_animation(volt_slice):

        fig = plt.figure(figsize=(14, 10))
        ax1 = plt.subplot2grid((5, 5), (0, 0), colspan=3, rowspan=3)
        pos = ax1.imshow(fft2d, aspect='auto', interpolation='none', extent=(extents(x) + extents(y)), origin='cool',
                         cmap='inferno', vmin=vmin, vmax=vmax)
        ax1.set_ylim([min(gates), max(gates)])
        if tgorbg:
            ax1.set_ylabel(r'$U_{tg}$ (V)')
        else:
            ax1.set_ylabel(r'$U_{bg}$ (V)')
        ax1.axhline(y=volt_slice, color='w', linestyle=':', linewidth=2)
        ax2 = plt.subplot2grid((5, 5), (3, 0), colspan=3, rowspan=2)
        ax2.plot(n2d, fft2d[y.index(volt_slice), :], '-x')
        ax2.set_ylim([vmin, vmax])
        ax2.set_xlabel('$n_{2d}$ in $10^{11} cm^{-2}$')
        ax2.set_ylabel('FFT (a.u.)')
        ax3 = plt.subplot2grid((5, 5), (0, 3), colspan=2, rowspan=5)
        data_pp = df_range(data_p, 'gate', [volt_slice - gate_step / 2, volt_slice + gate_step / 2])
        ax3.plot(1 / data_pp.bf, cutout_bkgd(1 / data_pp.bf.values, data_pp.rxx.values), 'r-x', linewidth=1,
                 label='raw data - background')
        ax4 = plt.twinx(ax3)
        ax4.plot(1 / data_pp.bf, data_pp.rxx.values, 'b-x', linewidth=1, label='raw data')
        ax3.set_xlabel(r'$B^{-1} (T^{-1})$')
        #         ax3.yaxis.tick_right()
        #         ax3.yaxis.set_label_position('right')
        ax3.set_ylabel(r'$R_{xx} (\Omega)$')
        #         ax3.set_xlim([0.25,0.6])
        ax3.legend(loc='lower right')
        ax4.legend(loc='upper right')
        fig.tight_layout()

    return interactive(plot_animation,
                       volt_slice=FloatSlider(min=min(gates), max=max(gates), step=gate_step, continuous_update=False))
