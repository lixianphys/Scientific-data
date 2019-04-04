import pandas as pd
import os


# physical constant
e0 = 1.6021766208E-19
h0 = 6.62607015E-34


# define pandas import
def pd_import(files,skiprows,usecols,names,step,Ref):
            super_data = pd.DataFrame()
            for f in files:
                    data = pd.read_csv(f, sep="\t",skiprows=skiprows, usecols=usecols, names=names, header=None)
                    data['rxx'] = data.uxx / data.curr * Ref
                    data['rxy'] = data.uxy / data.curr * Ref
                    data['sxx'] = data.rxx / (data.rxx ** 2 + data.rxy ** 2)
                    data['sxy'] = data.rxy / (data.rxx ** 2 + data.rxy ** 2)
                    data['tgate'] = step[files.index(f)]
                    data['source'] = f
                    super_data = super_data.append(data)
            return super_data

# define pandas output
def pd_output(pd_data,WorkingFolder, filename, step):
    if not os.path.exists(WorkingFolder):
        os.makedirs(WorkingFolder, 0o700)

    with open(os.path.join(WorkingFolder, filename), 'w') as outfile:
        pd_data.to_string(outfile)

    desc = filename.strip('.txt') + '_readme.txt'
    with open(os.path.join(WorkingFolder, desc), 'w') as readme:
        readme.write('Gate_step = {}'.format(step))

    return outfile