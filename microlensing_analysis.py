import os
from time import time
import h5py
import numpy as np
from scipy.spatial import KDTree

# Defining the paths to the simulation results

runs_path = 'H:/Bachelorarbeit/Daten/runs/'
if not os.path.exists(runs_path): runs_path = '/work/Tit6/paul.zuern/data/runs/'

# Defining the neighbour finding funtion

def find_closest_kdtree(pos, axis, mass):
    mass_threshold = 1e-3
    stars = mass > mass_threshold
    ffps = mass < mass_threshold
    
    axis = [0,1]
    
    leafsize = 20
    tree_stars = KDTree(pos[axis,:][:,stars].T, leafsize = leafsize)
    tree_ffps = KDTree(pos[axis,:][:,ffps].T, leafsize = leafsize)
    neighbours = tree_ffps.query_ball_tree(tree_stars, 0.001)
    
    return neighbours

def find_closest_digitize(pos, axis, mass):
    return None

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
    
    snaps = sorted(filter(lambda x: ('snap' in x), os.listdir(runs_path + run)), key = lambda x: int(x.split('.')[1][3:]))
    for snap_num, snap in enumerate(snaps):
        
        time_read = time()
        print(f'    Snapshot \'{snap}\' of {len(snaps)} - remaining {len(snaps) - snap_num}')
        
        with h5py.File(run_path + snap) as f:
            
            print(f'      Time for File reading: {time() - time_read:.2f} sec')
            time_search = time()
            
            # Sort the step keys and filter those that contain data of integer timesteps
            steps = list(filter(lambda x: f[x]['000 Scalars'][0] % 1 == 0, sorted(f.keys(), key = lambda x: int(x[5:]))))
            for step_num, step in enumerate(steps):
                
                # Read the data
                i = np.array(f[step]['032 Name'])
                m = np.array(f[step]['023 M'])
                x = np.array([f[step]['001 X1'], f[step]['002 X2'], f[step]['003 X3']])
                
                sort = np.argsort(i)
                i = i[sort]
                m = m[sort]
                x = x[:,sort]
                
                # Search for neighbours                
                neighbours = find_closest_kdtree(x, [0,1], m)
                print(neighbours)
                
            print(f'      Time for distance calculation of {len(list(steps))} steps: {time() - time_search:.2f} sec')
        
        total_time = time() - time_read
        print(f'      Total time for snapshot: {total_time:.2f} sec -> estimated remaining: {(total_time * (len(snaps) - snap_num) / 60):.1f} min')
        
    print(('   '.join(['{:13s}',] * len(events.keys())).format(*events.keys())))
    for e in np.arange(len(events['time'])):
        print(('   '.join(['{:13s}',] * len(events.keys())).format(*['{:.6g}'.format(events[key][e]) for key in events.keys()])))
    
    break

print('END')