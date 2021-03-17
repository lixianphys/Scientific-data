import argparse
import os
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from toybands.functions import *
from toybands.classes import *
from toybands.plottools import make_n_colors, make_1d_E_B_plots


def run():
    my_parser = argparse.ArgumentParser(
        prog="run", description="A band model to play with"
    )

    my_parser.add_argument(
        "-enplot",
        action="store_true",
        help="plot the energy versus bfield",
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
        action="store_const",
        const=20,
        help="number of Landau levels involved",
    )

    my_parser.add_argument(
        "-angle",
        action="store_const",
        const=0,
        help="angle in degree made with the sample plane norm by the external field",
    )


    args = my_parser.parse_args()
    print(vars(args))
    return args


if __name__ == "__main__":
    run()
    if os.path.isfile("system.json"):
        df = pd.read_json("system.json")
        newsystem = System()
        for i in range(len(df)):
            dt = df.iloc[i].to_dict()
            newband = Band(
                density=dt["density"],
                is_cond=dt["is_cond"],
                is_dirac=dt["is_dirac"],
                gfactor=dt["gfactor"],
                meff=dt["meff"],
                spin=dt["spin"],
                M=dt["M"],
                vf=dt["vf"],
            )
            newsystem.add_band(newband)
        args = run()
        enrange = list(
            np.linspace(
                int(args.enrange[0]*e0), int(args.enrange[1]*e0), int(args.enrange[2]*e0)
            )
        )
        bfrange = list(
            np.linspace(
                int(args.bfrange[0]), int(args.bfrange[1]), int(args.bfrange[2])
            )
        )
        if args.enplot:
            if args.nmax is not None and args.angle is not None:
                y_databdl = [[[x/e0 for x in band.cal_energy(bfrange,args.nmax,args.angle)[f'#{N}'].tolist()] for N in range(args.nmax)]  for band in newsystem.bands]
                colors = make_n_colors(len(y_databdl),'jet',0.1,0.9)
                make_1d_E_B_plots(bfrange,y_databdl,colors)
            else:
                sys.stderr.write('The argument -nmax and -angle is needed')
                exit(1)
    else:
        sys.stderr.write("no system (system.json) exist")