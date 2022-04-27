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
        
    print(f'  General Properties: N = {int(N)}, RBAR = {RBAR:.3f}, TSCALE = {TSCALE:.3f}')
    
    # The dictionary that will store the found microlensing events
    events = {'time': [], 'id_ffp': [], 'id_star': [], 'x1_ffp': [], 'x2_ffp': [], 'x3_ffp': [], 'x1_star':[], 'x2_star':[], 'x3_star':[]}
    
    for snap in sorted(filter(lambda x: ('snap' in x), os.listdir(runs_path + run)), key = lambda x: int(x.split('.')[1][3:])):
        if int(snap.split('.')[1][3:]) == 20: break
        
        print('    On Snapshot ' + snap.split('.')[1][3:])
        
        with h5py.File(run_path + snap) as f:
            for step in sorted(f.keys(), key = lambda x: int(x[5:])):
                # Exclude the timesteps with non integer times
                if f[step]['000 Scalars'][0] % 1 != 0: continue
                print('      At time {:.0f}'.format(f[step]['000 Scalars'][0]))
                
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
                                
                # Calculate the bins into which stars and ffp will be sorted for a faster microlensing search
                bins1, bins2 = np.linspace(min(x[axes[0]]), max(x[axes[0]]), 100), np.linspace(min(x[axes[1]]), max(x[axes[1]]), 500)
                
                # Sort the stars into the bins
                i_stars, j_stars = np.digitize(x[axes[0]][stars], bins1), np.digitize(x[axes[1]][stars], bins2)
                
                # Loop over all free floating particles
                for ffp_num in np.arange(len(m))[ffps]:
                    # Sort the ffp into a bin
                    i_ffp, j_ffp = np.digitize(x[axes[0]][ffp_num], bins1), np.digitize(x[axes[1]][ffp_num], bins2)
                    
                    # Retrieve all stars from the same and adjacent bins
                    vicinity = [-1,0,1]
                    neighbours = np.logical_and(np.isin(i_stars, i_ffp + vicinity), np.isin(j_stars, j_ffp + vicinity))
                    
                    # Calculate the disance in the observation plane
                    distances = np.linalg.norm(x[axes][:,stars][:,neighbours] - x[axes,ffp_num][:,np.newaxis], axis = 0)
                    
                    # Find the neighbour stars that will be lensed by the free floating planet and saving them as events
                    lensed = distances < 0.01
                    
                    for star_num in np.arange(len(m))[stars][neighbours][lensed]:
                        events['time'].append(f[step]['000 Scalars'][0])
                        events['id_ffp'].append(i[ffp_num])
                        events['id_star'].append(i[star_num])
                        events['x1_ffp'].append(x[0,ffp_num])
                        events['x2_ffp'].append(x[1,ffp_num])
                        events['x3_ffp'].append(x[2,ffp_num])
                        events['x1_star'].append(x[0,star_num])
                        events['x2_star'].append(x[1,star_num])
                        events['x3_star'].append(x[2,star_num])
                        
                    #break
        
    print(('   '.join(['{:13s}',] * len(events.keys())).format(*events.keys())))
    for e in np.arange(len(events['time'])):
        print(('   '.join(['{:13s}',] * len(events.keys())).format(*['{:.6g}'.format(events[key][e]) for key in events.keys()])))
    
    break
    
print('END')