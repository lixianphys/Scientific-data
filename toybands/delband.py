import argparse
import os
import sys

import pandas as pd
from toybands.functions import pretty_print


if __name__ == "__main__":
    if os.path.isfile("system.json"):
        df = pd.read_json("system.json")
        if not df.empty:
            print("Your current system:")
            pretty_print(df)
            index = input("which band to delete? input the index number: Or press 'e' to exit\n")
            if index == 'e':
                sys.stderr.write('Exit...')
                exit()
            print("Deleted No.", index, "band")
            df.drop([int(index)]).to_json("system.json")
            if not df.drop([int(index)]).empty:
                print("Your current system:")
                pretty_print(df.drop([int(index)]))
            else:
                print('No bands in your system, please add more using addband.py')
        else:
            print('No bands in your system, please add more using addband.py')
    else:
        sys.stderr.write("No system exist")
        exit()
