import argparse
import os
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from toybands.functions import *
from toybands.classes import *


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
        enrange = np.linspace(
            int(args.enrange[0]), int(args.enrange[1]), int(args.enrange[2])
        )
        bfrange = np.linspace(
            int(args.bfrange[0]), int(args.bfrange[1]), int(args.bfrange[2])
        )
        fig = plt.figure(figsize=(18, 18))
        ax = fig.add_subplot(111)
        for band in newsystem.bands:
            df = band.cal_energy(list(bfrange), 20, 0)
            print(df)
        # plt.savefig("test.pdf")
    else:
        print("no system exist")