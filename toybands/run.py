import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pdb
import xml.etree.ElementTree as ET
from tqdm import tqdm
from datetime import datetime

from physconst import e0
from utils import flattenList, div
from toybands.functions import *
from toybands.classes import *
from toybands.plottools import (make_n_colors, make_1d_E_B_plots, make_1d_den_B_plots, make_1d_dos_B_plot, make_2d_dos_map, super_save, make_slices,make_canvas,draw_band)
from toybands.config import *

def multi_floats(value):
    values = value.split()
    values = map(float, values)
    return list(values)


def run():
    my_parser = argparse.ArgumentParser(
        prog="run", description="Run your calculation here by throwing arguments into the <toybands> model"
    )
    plot_group = my_parser.add_mutually_exclusive_group(required=True)
    den_group = my_parser.add_mutually_exclusive_group(required=False)
    plot_group.add_argument(
        "-enplot",
        action="store_true",
        help="plot the energy versus bfield (yes/no)",
    )

    plot_group.add_argument(
        "-denplot",
        action="store_true",
        help="plot the density versus bfield (yes/no)",
    )

    plot_group.add_argument(
        "-simu",
        action="store_true",
        help="dynamically generate relationship between the density and the bfield at steps of input density (yes/no)",
    )

    plot_group.add_argument("-dos",action="store_true",help="plot dos versus bfield (yes/no)")
    plot_group.add_argument("-dosm",action="store_true",help="map dos at different (B,n) (yes/no)")

    den_group.add_argument(
        "--allden",
        action="store",
        type = multi_floats,
        help="densities for each band: start1 end1 start2 end2 .... ",
    )

    den_group.add_argument(
        "-loadden",
        action="store",
        type = str,
        help="file to load the densities",
    )

    my_parser.add_argument(
        "-nos",
        action="store",
        type= int,
        help="number of steps in the simulation ",
    )

    my_parser.add_argument(
        "-dir",
        action="store",
        help="relative output directory",
    )

    my_parser.add_argument(
        "-fnm",
        action="store",
        help="filename",
    )

    my_parser.add_argument(
        "--enrange",
        action="store",
        type=float,
        nargs=3,
        help="energy range: start end numofpoints",
    )

    my_parser.add_argument(
        "--bfrange",
        action="store",
        type=float,
        nargs=3,
        help="magnetic field range: start end numofpoints",
    )

    my_parser.add_argument(
        "-nmax",
        action="store",
        type=int,
        default=20,
        help="number of Landau levels involved (default=20)",
    )

    my_parser.add_argument(
        "-angle",
        action="store",
        type=float,
        default=0,
        help="angle in degree made with the sample plane norm by the external field (default=0)",
    )

    my_parser.add_argument(
        "-scond",
        action="store",
        type=float,
        help="LL broadening parameter sigma at B=1T for all conduction bands (meV,default value see config.py)",
    )

    my_parser.add_argument(
        "-sval",
        action="store",
        type=float,
        help="LL broadening parameter sigma at B=1T for all valence bands (meV,default value see config.py)",
    )

    my_parser.add_argument(
        "-parity",
        action="store_true",
        help="to decide if reducing the DOS of 0LL to half of others",
    )

    args = my_parser.parse_args()
    print('Your inputs:\n')
    print(vars(args))
    return args


def enplot(args,newsystem,bfrange,enrange,colors = None):
    if args.nmax is not None and args.angle is not None:
        y_databdl = newsystem.ydata_gen(args.nmax,args.angle, bfrange,enrange, 'enplot',SIGMA_COND,SIGMA_VAL,args.parity)
        if colors is None:
            colors = make_n_colors(len(newsystem.get_band('a')), DEFAULT_CMAP, DEFAULT_CMAP_VMIN, DEFAULT_CMAP_VMAX)
        mu_pos = [
            newsystem.mu(
                enrange,
                B,
                args.nmax,
                args.angle,
                [SIGMA_COND if band.get('is_cond') else SIGMA_VAL for band in newsystem.get_band('a')],
                args.parity
            )
            for B in bfrange
        ]
        ax = make_1d_E_B_plots(bfrange, y_databdl, colors, mu_pos)
        newsystem.databdl_write_csv(args.fnm,bfrange,y_databdl,'enplot')
        super_save(args.fnm, args.dir)
        system_stamp_csv(args)
    else:
        sys.stderr.write("The arguments -nmax and -angle are needed\n")

def denplot(args,newsystem,bfrange,enrange,colors = None):
    if args.nmax is not None and args.angle is not None:

        # bundle data from each Landau level originating from each band
        y_databdl = newsystem.ydata_gen(args.nmax,args.angle,bfrange, enrange, 'denplot',SIGMA_COND,SIGMA_VAL,args.parity)
        if colors is None:
            colors = make_n_colors(len(newsystem.get_band('a')), DEFAULT_CMAP, DEFAULT_CMAP_VMIN, DEFAULT_CMAP_VMAX)
        tot_den = newsystem.tot_density()
        make_1d_den_B_plots(bfrange, y_databdl, colors, tot_den)
        newsystem.databdl_write_csv(args.fnm,bfrange,y_databdl,'denplot')
        system_stamp_csv(args)
        super_save(args.fnm, args.dir)
    else:
        sys.stderr.write('The arguments -nmax and -angle are needed\n')

def simu(args,newsystem,bfrange,enrange,colors = None):
    if args.nmax is not None and args.angle is not None:
        if args.allden is not None and args.nos is not None:
            den_slices = make_slices(args.allden,args.nos)
        elif args.loadden is not None:
            den_slices = read_den_from_csv(args.loadden)
        else:
            sys.stderr.write('The arguments -nmax and -angle and -allden and -nos are needed\n')
        tot_den_list = [sum(den_slice) for den_slice in den_slices]
        tot_den_int = np.mean([abs(tot_den_list[i]-tot_den_list[i+1]) for i in range(len(tot_den_list)-1)])
        ax = make_canvas()
        if colors is None:
            colors = make_n_colors(len(newsystem.get_band('a')),DEFAULT_CMAP,DEFAULT_CMAP_VMIN,DEFAULT_CMAP_VMAX)
        for den_slice in tqdm(den_slices):
            # colors_p = color will only generate a new pointer to the same list
            colors_p = [color for color in colors]
            newsystem.set_all_density(den_slice)
            tot_den = newsystem.tot_density()
            IDOS = [newsystem.dos_gen(enrange, B, args.nmax, args.angle, [SIGMA_COND if band.get('is_cond') else SIGMA_VAL for band in newsystem.get_band('a')],args.parity) for B in bfrange]
            # disable the band with zero density
            idx_disabled = []
            # if any band has a density <=0 then disable it
            for idx,den in enumerate(den_slice):
                if den <= 0:
                    newsystem.get_band(idx).disable()
                    idx_disabled.append(idx)
                    colors_p[idx] = None
            colors_p = [color for color in colors_p if color is not None]
            y_databdl = [[[np.interp(x=x, xp=enrange, fp=IDOS[index]) for index, x in enumerate(band.cal_energy(bfrange,args.nmax,args.angle)[f'#{N}'].tolist())] for N in range(args.nmax)] for band in newsystem.get_band('a')]
            plotrange = [tot_den-0.5*tot_den_int,tot_den+0.5*tot_den_int]
            ax = make_1d_den_B_plots(bfrange,y_databdl,colors_p,ax=ax,plotrange=plotrange,legend=False)
            newsystem.databdl_write_csv(args.fnm,bfrange,y_databdl,'simu',plotrange=plotrange)
            # enable the disabled band again for next loop
            if idx_disabled:
                for idx in idx_disabled:
                    newsystem.get_band(idx).enable()
        system_stamp_csv(args)    
        super_save(args.fnm,args.dir)

    else:
        sys.stderr.write('The arguments -nmax and -angle are needed\n')


def dos_at_mu(args,newsystem,bfrange,enrange,colors=None):
    if args.nmax is not None and args.angle is not None:
        if colors is None:
            colors = make_n_colors(len(newsystem.get_band('a')),DEFAULT_CMAP,DEFAULT_CMAP_VMIN,DEFAULT_CMAP_VMAX)
        y_databdl = newsystem.ydata_gen(args.nmax,args.angle, bfrange,enrange, 'dos',SIGMA_COND,SIGMA_VAL,args.parity)
        make_1d_dos_B_plot(bfrange,y_databdl,colors)
        newsystem.databdl_write_csv(args.fnm,bfrange,y_databdl,'dos')
        system_stamp_csv(args)
        super_save(args.fnm,args.dir)
    else:
        sys.stderr.write('The arguments -nmax and -angle are needed\n')

def dos_map(args,newsystem,bfrange,enrange,cmap=None):
    if args.nmax is not None and args.angle is not None:
        if cmap is None:
            cmap = DEFAULT_CMAP
        if args.allden is not None and args.nos is not None:
            den_slices = make_slices(args.allden,args.nos)
        elif args.loadden is not None:
            den_slices = read_den_from_csv(args.loadden)
        else:
            sys.stderr.write('The arguments -nmax and -angle and -allden and -nos are needed\n')
        ax = make_canvas()
        # data storage
        y_databdl = []
        y_tot = []
        for den_slice in tqdm(den_slices):
            newsystem.set_all_density(den_slice)
            tot_den = newsystem.tot_density()
            y_tot.append(tot_den)
            # disable the band with zero density or below
            idx_disabled = []
            # if any band has a density <=0 then disable it
            for idx,den in enumerate(den_slice):
                if den <= 0:
                    newsystem.get_band(idx).disable()
                    idx_disabled.append(idx)
            # calculate the chemical potential for the new system        
            mus = [newsystem.mu(enrange, B, args.nmax, args.angle, [SIGMA_COND if band.get('is_cond') else SIGMA_VAL for band in newsystem.get_band('a')],args.parity) for B in bfrange]
            # calculate the dos for each band at each chemical potential along the B field axis
            to_append =  [sum([(1 if band.get('is_cond') else -1)*band.cal_dos(mu_at_B, B, args.nmax, args.angle, SIGMA_COND if band.get('is_cond') else SIGMA_VAL,args.parity) for band in newsystem.get_band('a')]) for B,mu_at_B in zip(bfrange,mus)]
            y_databdl.append(to_append)
            newsystem.databdl_write_csv(args.fnm,bfrange,to_append,'dosm')
            # enable the disabled band again for next loop
            if idx_disabled:
                for idx in idx_disabled:
                    newsystem.get_band(idx).enable()
        make_2d_dos_map(bfrange,y_tot,y_databdl,cmap,ax=ax)
        system_stamp_csv(args)
        super_save(args.fnm,args.dir)
    else:
        sys.stderr.write('The arguments -nmax and -angle are needed\n')

def draw(args,newsystem):
    draw_band(newsystem)
    super_save(args.fnm+'_bandsketch',os.path.join(DEFAULT_PATH,args.fnm))


def loadsys(df):
    load_system = System()
    for i in range(len(df)):
        dt = df.iloc[i].to_dict()
        load_band = Band(
            density=dt["density"],
            is_cond=dt["is_cond"],
            is_dirac=dt["is_dirac"],
            gfactor=dt["gfactor"],
            meff=dt["meff"],
            spin=dt["spin"],
            dparam = dt["dparam"],
            vf=dt["vf"],
        )
        load_system.add_band(load_band)
    return load_system

def output_xml(args):
    now = datetime.now()
    data = ET.Element('data')

    time = ET.SubElement(data,'time')
    time_ymd = ET.SubElement(time,'ymd')
    time_ymd.set('name','yyyy/mm/dd')
    time_ymd.text = str(now.strftime("%Y/%m/%d"))
    time_hms = ET.SubElement(time,'hms')
    time_hms.set('name','hh/mm/ss')
    time_hms.text = str(now.strftime("%H:%M:%S"))

    configs = ET.SubElement(data,'configs')
    sig_cond = ET.SubElement(configs,'sigma_cond')
    sig_val = ET.SubElement(configs,'sigma_val')
    sig_cond.set('name','SIGMA_COND')
    sig_val.set('name','SIGMA_VAL')
    sig_cond.text = str(SIGMA_COND)
    sig_val.text = str(SIGMA_VAL)
    
    items = ET.SubElement(data,'args')
    for key,arg in vars(args).items():
        item = ET.SubElement(items,'arg')
        item.set('name',str(key))
        item.text = str(arg)
    if args.fnm:
        myfile = open(os.path.join(mkdir(args.fnm),args.fnm.split('.')[0]+'_args.xml'),'w')
    else:
        myfile = open(os.path.join(mkdir(DEFAULT_AUTONAME),DEFAULT_AUTONAME+'_args.xml'),'w')
    myfile.write(et_pretty(data))

if __name__ == "__main__":
    args = run()
    output_xml(args)
    if args.scond is not None:
        SIGMA_COND = args.scond*e0*1e-3
    if args.sval is not None:
        SIGMA_VAL = args.sval*e0*1e-3
    if os.path.isfile("system.json"):
        df = pd.read_json("system.json")
        newsystem = loadsys(df)
        enrange = list(
            np.linspace(
                args.enrange[0] * e0, args.enrange[1] * e0, int(args.enrange[2])
            )
        )
        bfrange = list(
            np.linspace(args.bfrange[0], args.bfrange[1], int(args.bfrange[2]))
        )
        if args.enplot:
            enplot(args,newsystem,bfrange,enrange)
        if args.denplot:
            denplot(args,newsystem,bfrange,enrange)
        if args.simu:
            simu(args,newsystem,bfrange,enrange)
        if args.dos:
            dos_at_mu(args,newsystem,bfrange,enrange)
        if args.dosm:
            dos_map(args,newsystem,bfrange,enrange)
        draw(args,newsystem)
    else:
        sys.stderr.write("no system (system.json) exist\n")