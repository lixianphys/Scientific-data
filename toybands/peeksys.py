import argparse
import os
import sys

import pandas as pd
from toybands.functions import pretty_print

def peek():
    my_parser = argparse.ArgumentParser(
        prog="peek", description="Peek into system"
    )
    df = pd.read_json("system.json")
    if not df.empty:
        print("Your current system:")
        pretty_print(df)
        exit()
    else:
        print('No bands in your system, please add more using addband.py')

if __name__ == "__main__":
    if os.path.isfile("system.json"):
        peek()
    else:
        sys.stderr.write("No system exist")
        exit()