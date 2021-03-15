import argparse
import os
import sys

from SciData.toybands.functions import *
from SciData.toybands.classes import *

my_parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
my_parser.version = '1.0'   
def addband():
    my_parser = argparse.ArgumentParser(prog='addband',
                                        description='A band model to play with')

    my_parser.add_argument('density',
                        metavar='density',
                        type=float,
                        help='density to create a band (float)')
    my_parser.add_argument('is_cond',
                        metavar='is_cond',
                        type=bool,
                        help='conduction/valence band (bool)')
    my_parser.add_argument('is_dirac',
                        metavar='is_dirac',
                        type=bool,
                        help='Dirac/Conventional band (bool)')
    my_parser.add_argument('gfactor',
                        metavar='gfactor',
                        type=float,
                        help='gfactor to create a band (float)')

    args = my_parser.parse_args()
    newband = Band(density=args.density, is_cond=args.is_cond, is_dirac=args.is_dirac, gfactor=args.gfactor, M=0.1, vf=1e6)
    print('This band added')
    newband.print_band()
    return newband
 
if __name__ == "__main__":
    if os.path.isfile('system.json'):
        df = pd.read_json('system.json')
        newsystem = System()
        for i in range(len(df)):
            dt = df.iloc[i].dropna().to_dict()
            if 'meff' in dt.keys():
                newband = Band(density=dt.density,is_cond=dt.is_cond, is_dirac=dt.is_dirac, gfactor=dt.gfactor,meff = dt.meff, spin = dt.spin)
            elif 'M' in dt.keys():
                newband = Band(density=dt.density,is_cond=dt.is_cond, is_dirac=dt.is_dirac, gfactor=dt.gfactor,M = dt.M, vf = dt.vf)
            newsystem.add_band(newband)
        inputband = addband()
        newsystem.add_band(inputband)
        print(newsystem.get_info())
        df.append(pd.DataFrame.from_dict(inputband.get_info()))
        df.to_json('system.json')
    else:
        newsystem = System()    
        inputband = addband()
        newsystem.add_band(inputband)
        df = newsystem.get_info()
        print(df)
        df.to_json('system.json')
        