import argparse
import os
import sys

import pandas as pd


if __name__ == "__main__":
    if os.path.isfile("system.json"):
        df = pd.read_json("system.json")
        print(df)
        order = input("which band to delete, input the number:")
        print("going to delete no.", order, "band")
        df.drop([int(order)])
        df.drop([int(order)]).to_json(os.path.join(os.getcwd(), "system.json"))
        print("Your current system:")
        print(df.drop([int(order)]))
    else:
        print("no system exist")
