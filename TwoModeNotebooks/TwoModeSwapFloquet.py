import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from cycler import cycle

import qutip as qp

import scqubits as scq

import xarray as xr


import scipy.linalg as la
import scipy.optimize as opt
import scipy.fft as fft

import colormaps as cmaps
import matplotlib_inline.backend_inline
matplotlib_inline.backend_inline.set_matplotlib_formats('png')

import itertools as itert

import IPython.display

import sys

import copy

import pickle

import json

import pathlib

import os

import inspect

from datetime import datetime

'''
In the kwargs definition in line 76, the 'shifts' is the array of stark shifts that are sampled. 

When you run this, do it as:

nohup python -u TwoModeSwapFloquet.py > outputs/TwoModeSwapFloquet.out &

The TwoModeFloquetSweep.out can have what ever name you want. It is just the log file for all the print statements. 
The last two will be the stark shift and the drive time. The plot will be saved in the Figs folder in the Mode_3_And_5 directory.
'''

path_to_package = str(pathlib.Path(os.getcwd()).parents[0])+'/Main Package/'
print(path_to_package)
sys.path.insert(0,path_to_package)

main_path = str(pathlib.Path(os.getcwd()).parents[0])+'/'

import importlib

import Transmon_Cavity_Model as TCM

model_name = 'Mode_3_And_5'
save_path = 'Model Saves/'+model_name
save_path = main_path+save_path + '/'

Mode35 = TCM.LoadModel(save_path+model_name+'.json')

print(f'Got Model')


state_g10 = ['g',1,0]
state_g01 = ['g',0,1]

op_name_g10_SWAP = 'SWAP_g10' 

default_freq_SWAP_g10 = Mode35.DefaultFrequency(['g',1,0], ['g',0,1])/2

print(f'Default Frequency: {default_freq_SWAP_g10}')

Odeoptions = {'nsteps': 10000, 'max_step':1, 'atol':1e-6, 'rtol':1e-6}#, 'method':'bdf'}
solver_ops = qp.Options(**Odeoptions)
kwargs = {'shifts':np.linspace(0.01, 0.05, 4), 'show_plot':True, 'options': solver_ops, 'save_plot':True, 'debug':True, 'save_path':save_path+'Figs/' }

SWAP_g10_Stark_Shift = Mode35.GetStarkShift(state_g10, state_g01, 0.5, use_fit = False, kwargs = kwargs, freq_d = default_freq_SWAP_g10)
starkshift = SWAP_g10_Stark_Shift['x'][0]
initial_time = 1/SWAP_g10_Stark_Shift['fun']
print(f'stark shift: {starkshift}')
print(f'Initial Drive Time: {initial_time}')