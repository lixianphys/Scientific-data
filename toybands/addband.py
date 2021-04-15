import argparse
import os
import sys

from toybands.functions import *
from toybands.classes import *


def addband():
    my_parser = argparse.ArgumentParser(
        prog="addband", description="Add band to a system"
    )
    group = my_parser.add_mutually_exclusive_group(required=True)

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
        type=float,
        action="store",
        help="gfactor to create a band (float)",
    )
    group.add_argument(
        "-dp",
        "-diracpara",
        type=float,
        nargs=2,
        help="(vf,dparam) parameters to create a Dirac-like band (float:m/s,float:meV nm^2)",
    )
    group.add_argument(
        "-cp",
        "-convpara",
        type= float,
        nargs = 2,
        help="(meff,spin) parameters to create a conventional band (float:me,int:-1,0,+1)",
    )


    args = my_parser.parse_args()
    if args.dp is None and args.cp is not None:
        newband = Band(
            density=args.density,
            is_cond=args.is_cond,
            is_dirac=args.is_dirac,
            gfactor=args.gfactor,
            vf=None,
            dparam=None,
            meff=args.cp[0] if isinstance(args.cp[0],(int, float, complex)) else None,
            spin=args.cp[1] if isinstance(args.cp[1],(int, float, complex)) else None,
        )
    elif args.dp is not None and args.cp is None:
        newband = Band(
            density=args.density,
            is_cond=args.is_cond,
            is_dirac=args.is_dirac,
            gfactor=args.gfactor,
            vf=args.dp[0] if isinstance(args.dp[0],(int, float, complex)) else None,
            dparam=args.dp[1]*e0*1e-21 if isinstance(args.dp[1],(int, float, complex)) else None,
            meff= None,
            spin= None,
        )
        
    print("This band added")
    newband.print_band()
    return newband


if __name__ == "__main__":
    if os.path.isfile("system.json"):
        df = pd.read_json("system.json")
        # rebuild the system from JSON file
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
                dparam = dt["dparam"],
                vf=dt["vf"],
            )
            newsystem.add_band(newband)
        inputband = addband()
        newsystem.add_band(inputband)
        df = newsystem.get_info()
        df.to_json(os.path.join(os.getcwd(), "system.json"))
        pretty_print(df)
    else:
        newsystem = System()
        inputband = addband()
        newsystem.add_band(inputband)
        df = newsystem.get_info()
        df.to_json(os.path.join(os.getcwd(), "system.json"))
        pretty_print(df)
