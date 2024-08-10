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


ops = ['q_g_00', 'half-q_e_00', 'sb_f_0_Mode_3', 'q_e_00', 'sb_f_0_Mode_5']


psi0 = Mode35.get_dressed_state(['g', 0, 0])

print(f'ops:\n{ops}')
print(f'Getting c_ops')
C_and_D_Ops = TCM.GetCollapseAndDephasing(Mode35)
#with open("Data/Mode35_c_ops_save.pkl",  'rb') as file:
#    C_and_D_Ops = pickle.load(file)
c_ops = list(C_and_D_Ops.values())
print(f'Got c_ops')
Odeoptions = {'nsteps': 1000000, 'max_step':None, 'atol':1e-6, 'rtol':1e-6}
res = Mode35.Run_Pulse_Sequence(psi0, ops, spps = 10, Odeoptions = Odeoptions, c_ops=c_ops)


with open("Data/Noon_State_n_1.pkl", 'wb') as file:
    pickle.dump(res, file)