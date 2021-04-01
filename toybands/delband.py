import argparse
import os
import sys

import pandas as pd
from toybands.functions import pretty_print


def get_args():
    my_parser = argparse.ArgumentParser(
        prog="delband", description="Delete a band or all bands"
    )
    my_parser.add_argument(
    "-i",
    action="store",
    type=str,
    help="index number of band to delete"
    )
    args = my_parser.parse_args()
    return args

def no_ui(df,index):
    if index == 'e':
        sys.stderr.write('Exit...\n')
    elif index == 'all':
        df.drop(df.index).to_json("system.json")
        print('No bands in your system, please add more using addband.py')
    else:
        df.drop([int(index)]).to_json("system.json")

def with_ui(df):
    if not df.empty:
        print("Your current system:")
        pretty_print(df)
        index = input("which band to delete? input the index number: Or press 'e' to exit, 'all' to clear up\n")
        if index == 'e':
            sys.stderr.write('Exit...\n')
        elif index == 'all':
            df.drop(df.index).to_json("system.json")
            print('No bands in your system, please add more using addband.py')
        else:
            print("Deleted No.", index, "band")
            df.drop([int(index)]).to_json("system.json")
            if not df.drop([int(index)]).empty:
                print("Now, your system:")
                pretty_print(df.drop([int(index)]))
            else:
                print('No bands in your system, please add more using addband.py\n')
    else:
        print('No bands in your system, please add more using addband.py\n')

if __name__ == "__main__":
    if os.path.isfile("system.json"):
        df = pd.read_json("system.json")
        args = get_args()
        if args.i is not None:
            no_ui(df,args.i)
        else:
            with_ui(df)
    else:
        sys.stderr.write("No system exist\n")
        exit()
# once delete a band not at the bottom (end of the df), then the end of the df is NAN in density/Ebb
