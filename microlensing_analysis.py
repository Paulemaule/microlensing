import os
from time import time
import h5py
import numpy as np
from scipy.spatial import KDTree

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
        
    print(f'  General Properties: N = {int(N)}, RBAR = {RBAR:.3f}, TSCALE = {TSCALE:.3f}')
    
    # The dictionary that will store the found microlensing events
    events = {'time': [], 'id_ffp': [], 'id_star': [], 'x1_ffp': [], 'x2_ffp': [], 'x3_ffp': [], 'x1_star':[], 'x2_star':[], 'x3_star':[]}
    
    time_average = 5
    
    snaps = sorted(filter(lambda x: ('snap' in x), os.listdir(runs_path + run)), key = lambda x: int(x.split('.')[1][3:]))
    for snap in snaps:
        #if int(snap.split('.')[1][3:]) == 20: break
        start_time = time()
        time_remaining = time_average * (int(snaps[-1].split('.')[1][3:]) - int(snap.split('.')[1][3:]))
        print('    Snapshot {} of {} | time remaining: {:.2f} min  '.format(snap.split('.')[1][3:], snaps[-1].split('.')[1][3:], time_remaining))
        
        with h5py.File(run_path + snap) as f:
            for step in sorted(f.keys(), key = lambda x: int(x[5:])):
                # Exclude the timesteps with non integer times
                if f[step]['000 Scalars'][0] % 1 != 0: continue
                
                # Read the date
                i = np.array(f[step]['032 Name'])
                m = np.array(f[step]['023 M'])
                x = np.array([f[step]['001 X1'], f[step]['002 X2'], f[step]['003 X3']])
                
                sort = np.argsort(i)
                i = i[sort]
                m = m[sort]
                x = x[:,sort]
                
                stars = m > 1e-3
                ffps = m < 1e-3
                
                # Choose the axis along which the events are to be observed
                axes = [0,1]
                
                # Construct the KD-Tree
                leafsize = 20
                tree_stars = KDTree(x[axes,:][:,stars].T, leafsize = leafsize)
                tree_ffps = KDTree(x[axes,:][:,ffps].T, leafsize = leafsize)
                neighbours = tree_ffps.query_ball_tree(tree_stars, 0.01)
                
                #print(x[axes,:][:,ffps][:,:5])
                #print([item for sublist in neighbours[:5] for item in sublist])
        
        time_average = 0.5 * time_average + 0.5 * (time() - start_time) / 60
        
    print(('   '.join(['{:13s}',] * len(events.keys())).format(*events.keys())))
    for e in np.arange(len(events['time'])):
        print(('   '.join(['{:13s}',] * len(events.keys())).format(*['{:.6g}'.format(events[key][e]) for key in events.keys()])))
    
    break
    
print('END')