import argparse
import os
import sys

from toybands.functions import *
from toybands.classes import *


def addband():
    my_parser = argparse.ArgumentParser(
        prog="addband", description="Add band to a system"
    )

    my_parser.add_argument(
        "density",
        metavar="density",
        type=float,
        help="density to create a band (float,1/m2)",
    )
    my_parser.add_argument(
        "-is_cond",
        action="store_true",
        help="conduction/valence band (bool)",
    )
    my_parser.add_argument(
        "-is_dirac",
        action="store_true",
        help="Dirac/Conventional band (bool)",
    )
    my_parser.add_argument(
        "-gfactor",
        metavar="gfactor",
        type=float,
        action="store",
        help="gfactor to create a band (float)",
    )
    my_parser.add_argument(
        "-M",
        metavar="M",
        type=float,
        action="store",
        help="M to create a band (float,eV)",
    )
    my_parser.add_argument(
        "-vf",
        metavar="vf",
        type=float,
        action="store",
        help="M to create a band (float,m/s)",
    )
    my_parser.add_argument(
        "-meff",
        metavar="meff",
        type=float,
        action="store",
        help="meff to create a band (float,me)",
    )
    my_parser.add_argument(
        "-spin",
        metavar="spin",
        type=int,
        action="store",
        help="spin to create a band (float)",
        choices=[-1, 1],
    )

    args = my_parser.parse_args()
    newband = Band(
        density=args.density,
        is_cond=args.is_cond,
        is_dirac=args.is_dirac,
        gfactor=args.gfactor,
        M=args.M,
        vf=args.vf,
        meff=args.meff,
        spin=args.spin,
    )
    print("This band added")
    newband.print_band()
    return newband


if __name__ == "__main__":
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
        inputband = addband()
        newsystem.add_band(inputband)
        df = newsystem.get_info()
        df.to_json(os.path.join(os.getcwd(), "system.json"))
        print(df)
    else:
        newsystem = System()
        inputband = addband()
        newsystem.add_band(inputband)
        df = newsystem.get_info()
        df.to_json(os.path.join(os.getcwd(), "system.json"))
        print(df)
