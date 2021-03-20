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
            exit()
        else:
            print('No bands in your system, please add more using addband.py')
    else:
        sys.stderr.write("No system exist")
        exit()