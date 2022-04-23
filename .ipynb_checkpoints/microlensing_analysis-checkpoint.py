import os
import h5py
import numpy as np
from matplotlib import pyplot as plt

# Defining the paths to the simulation results

run_path = 'H:/Bachelorarbeit/Daten/runs/'
if not os.path.exists(run_path): run_path = '/work/Tit6/paul.zuern/data/runs/'

# Reading the data

for run in sorted(os.listdir(runs_path), key = lambda x: int(x[4:])):
    print(run)
    print('TEST')
    
print('END')