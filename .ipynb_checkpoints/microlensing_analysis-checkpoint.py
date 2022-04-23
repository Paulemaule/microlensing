import os
import h5py
import numpy as np
from matplotlib import pyplot as plt

# Defining the paths to the simulation results

runs_path = 'H:/Bachelorarbeit/Daten/runs/'
if not os.path.exists(runs_path): runs_path = '/work/Tit6/paul.zuern/data/runs/'

# Reading the data

for run in sorted(os.listdir(runs_path), key = lambda x: int(x[4:])):
    print('On ' + run)
    run_path = runs_path + run + '/'
    
    RBAR, TSCALE, N = None, None, None
    with h5py.File(run_path + 'snap.40_0.h5part') as f:
        RBAR = f[list(f.keys())[-1]]['000 Scalars'][2]
        TSCALE  = f[list(f.keys())[-1]]['000 Scalars'][10]
        N = f[list(f.keys())[-1]]['000 Scalars'][4]
        
    print(f'   General Properties: N = {int(N)}, RBAR = {RBAR:.3f}, TSCALE = {TSCALE:.3f}')
    
    for snap in sorted(filter(lambda x: ('snap' in x), os.listdir(runs_path + run)), key = lambda x: int(x.split('.')[1][3:])):
        if int(snap.split('.')[1][3:]) == 20: break
        
        print('   On ' + snap)
        
        with h5py.File(run_path + snap) as f:
            print(f.keys())
            for step in sorted(f.keys(), key = lambda x: int(x[5:])):
                if f[step]['000 Scalars'][0] % 1 != 0: continue
                
                i = np.array(f[step]['032 Name'])
                m = np.array(f[step]['023 M'])
                x = np.array([f[step]['001 X1'], f[step]['002 X2'], f[step]['003 X3']])
                
                sort = np.argsort(i)
                i = i[sort]
                m = m[sort]
                x = x[:,sort]
                
                
    break
    
print('END')