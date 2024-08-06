import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from cycler import cycler

import scipy.optimize as opt

import matplotlib.patches as mpatches
from matplotlib.markers import MarkerStyle
import qutip as qp

import itertools as itert

import scqubits as scq

import xarray as xr

import colormaps as cmaps
import matplotlib_inline.backend_inline
matplotlib_inline.backend_inline.set_matplotlib_formats('png')

import pickle

import inspect

import copy as copy


import pathlib

import os

import IPython.display

import json
from datetime import datetime

Transmon_Labels_LK = {'g':0, 'e':1, 'f':2, 'h':3, 'i':4, 'j':5, 'k':6, 'l':7,'m':8, 'o':9, 'p':10}

Transmon_Labels_NK = {}
for  i,key in enumerate(Transmon_Labels_LK):
    Transmon_Labels_NK[i]=key

marker_dict = {'g':'1', 'e':'2', 'f':'3', 'h':'4', 'i':'+', 'j':'x', 'k':'o', 'l':'v', 'm':"^", 'o':'>', 'p':'<'}

def BraKetLabel(state, dag = False):
    if type(state[0])!= str:
        state[0] = Transmon_Labels_NK[state[0]]
    
    label = state[0]
    for i in range(len(state)-1):
        label += f",{state[i+1]}"
    if not dag: 
        return r"$|$LABEL$\rangle$".replace('LABEL', label)
    if dag:
        return r"$\langle$LABEL$|$".replace('LABEL', label)

def LoadModel(file_name):
    with open(file_name, 'r') as file:
        params = json.load(file)
    
    model = Transmon_Cavity_Model(**params['Main_Config'])
    model.op_drive_params_dict = params['op_drive_params_dict']
    
    if 'Extra_Stuff' in params:
        model.Extra_Stuff = params['Extra_Stuff']

    if 'Drive Sequence' in params:
        model.Drive_Sequences = params['Drive Sequences']

    return model

def CopyModel(old_model_name, new_model_name, save_path, Return = False):
    old_save_path = save_path+old_model_name+'/'
    new_save_path = save_path+new_model_name+'/'

    model = LoadModel(old_save_path+old_model_name+'.json')

    model.save_path = new_save_path
    model.model_name = new_model_name

    model.update_model_params()
    model.Save_Params()

    if Return: return model

def Square_Envelope(t):
    return 1

def Guassian_Envelope(t, sigma, mu=0):
    return np.exp(-(t-mu)**2/(2*sigma**2))

def Linear_Envelope(t, t0, y0, m):
    return m*(t-t0)+y0

def Sine_Squared_Envelope(t, ramp_time, offset=0, phi = 0):
    return np.sin((np.pi/2)*(t-offset)/(ramp_time)+phi)**2

def Gaussian_Ramp_Envelope(t, pulse_time, ramp_time = 20, sigma_factor = 2):
    sigma = ramp_time/sigma_factor

    flat_top_time = pulse_time-2*ramp_time

    if (t<= ramp_time):
        return Guassian_Envelope(t, sigma, ramp_time)
    
    if (t>ramp_time and t<flat_top_time+ramp_time):
        return 1
    
    if (t>= flat_top_time+ramp_time):
        return Guassian_Envelope(t, sigma, mu = flat_top_time+ramp_time)

def Sine_Squared_Ramp_Envelope(t, pulse_time, ramp_time = 20):
    flat_top_time = pulse_time-2*ramp_time

    if (t<= ramp_time):
        return Sine_Squared_Envelope(t, ramp_time)
    
    if (t>ramp_time and t<flat_top_time+ramp_time):
        return 1
    
    if (t>= flat_top_time+ramp_time):
        return Sine_Squared_Envelope(t, ramp_time, offset = flat_top_time+ramp_time, phi = np.pi/2)


global_envelope_dict = {}
global_envelope_dict['Guassian_Ramp'] = Gaussian_Ramp_Envelope
global_envelope_dict['Sine_Squared_Ramp'] = Sine_Squared_Ramp_Envelope
global_envelope_dict['Square'] = Square_Envelope
global_envelope_dict['Guassian'] = Guassian_Envelope
global_envelope_dict['Sine_Squared'] = Sine_Squared_Envelope

def tlist(pulse_time, starttime=0, spps=10):
    ''' - Make a list of times
        - spps = sample points per second (nanoseconds for our purposes)
    '''
    tfull = np.linspace(starttime, starttime+pulse_time, round(pulse_time*spps) + 1)
    
    return tfull

def Calibrate_Ramped_Pulse(model, state1, state2, args, t0, tf, steps, op_name, drive_state_kwargs = {}, make_plot = False, debug = False, save_pulse = True):
    if not isinstance(state1, qp.qobj.Qobj): 
        if type(state1[0])==str:
            state1[0] = model.Transmon_Labels_LK[state1[0]]
        
        psi_i = model.get_dressed_state(state1)
    
    else:
        psi_f = state2

    if not isinstance(state2, qp.qobj.Qobj): 
        if type(state2[0])==str:
            state2[0] = model.Transmon_Labels_LK[state2[0]]
        
        psi_f = model.get_dressed_state(state2)
    
    else:
        psi_f = state2
    
   

    if 'freq_d' not in args:
        args['freq_d'] = model.DefaultFrequency(state1, state2)


    final_prob_list = []
    t_array = np.linspace(t0, tf, steps)

    for i in range(len(t_array)):
        t = t_array[i]
        print(f'Doing step {i+1}/{len(t_array)}, t = {t}')
        args['pulse_time'] = t
        args['Envelope Args']['pulse_time'] = args['pulse_time']
        if debug: print(f'  args: {args}')


        res = model.Drive_State(psi_i, args, **drive_state_kwargs)

        state_prob = qp.expect(psi_f.proj(), res.states)
        
        final_prob_list.append(state_prob[-1])
        print(f'Transition prob: {state_prob[-1]}')
        print(f'-----------------------------------------------------------------\n')
        

    plt.plot(t_array, final_prob_list)

    final_prob_list = np.array(final_prob_list)

    args['pulse_time'] = t_array[final_prob_list == final_prob_list.max()][0]
    args['Envelope Args']['pulse_time'] = args['pulse_time']
    
    if make_plot:
        plt.figure(figsize = (5,3), dpi = 200)
        plt.plot(t_array, final_prob_list, marker = 'x', lw = 0.5)
    
    if save_pulse: model.op_drive_params_dict[op_name] = copy.deepcopy(args)

def Calibrate_Guassian_Pulse(model, state1, state2, args, t0, tf, steps, op_name, drive_state_kwargs = {}, make_plot = False, debug = False, save_pulse = True):
    if type(state1[0])==str:
            state1[0] = model.Transmon_Labels_LK[state1[0]]
        
    if type(state2[0])==str:
        state2[0] = model.Transmon_Labels_LK[state2[0]]

    psi_i = model.get_dressed_state(state1)

    psi_f = state2
    if type(state2) != qp.qobj.Qobj:
        psi_f = model.get_dressed_state(state2)

    if 'freq_d' not in args:
        args['freq_d'] = model.DefaultFrequency(state1, state2)


    final_prob_list = []
    t_array = np.linspace(t0, tf, steps)

    for i in range(len(t_array)):
        t = t_array[i]
        print(f'Doing step {i+1}/{len(t_array)}, t = {t}')
        args['pulse_time'] = t
        args['Envelope Args'] = {'sigma': t/4, 'mu':t/2}

        

        res = model.Drive_State(psi_i, args, **drive_state_kwargs)

        state_prob = qp.expect(psi_f.proj(), res.states)
        
        final_prob_list.append(state_prob[-1])
        print(f'Transition prob: {state_prob[-1]}')
        print(f'-----------------------------------------------------------------\n')
        

    plt.plot(t_array, final_prob_list)

    final_prob_list = np.array(final_prob_list)

    if make_plot:
        plt.figure(figsize = (5,3), dpi = 200)
        plt.plot(t_array, final_prob_list, marker = 'x', lw = 0.5)

    args['pulse_time'] = t_array[final_prob_list == final_prob_list.max()][0]
    
    if save_pulse: model.op_drive_params_dict[op_name] = copy.deepcopy(args)

def Calibrate_Sine_Squared_Pulse(model, state1, state2, args, t0, tf, steps, op_name, drive_state_kwargs = {}, make_plot = False, debug = False, save_pulse = True):
    if type(state1[0])==str:
            state1[0] = model.Transmon_Labels_LK[state1[0]]
        
    if type(state2[0])==str:
        state2[0] = model.Transmon_Labels_LK[state2[0]]

    psi_i = model.get_dressed_state(state1)
    psi_f = model.get_dressed_state(state2)

    if 'freq_d' not in args:
        args['freq_d'] = model.DefaultFrequency(state1, state2)


    final_prob_list = []
    t_array = np.linspace(t0, tf, steps)

    for i in range(len(t_array)):
        t = t_array[i]
        print(f'Doing step {i+1}/{len(t_array)}, t = {t}')
        args['pulse_time'] = t
        args['Envelope Args'] = {'ramp_time': t/2, 'offset':0}

        psi0 = model.get_dressed_state(state1)
        

        res = model.Drive_State(psi0, args, **drive_state_kwargs)

        state_prob = qp.expect(model.get_dressed_state(state2).proj(), res.states)
        
        final_prob_list.append(state_prob[-1])
        print(f'Transition prob: {state_prob[-1]}')
        print(f'-----------------------------------------------------------------\n')
        

    plt.plot(t_array, final_prob_list)

    final_prob_list = np.array(final_prob_list)

    if make_plot:
        plt.figure(figsize = (5,3), dpi = 200)
        plt.plot(t_array, final_prob_list, marker = 'x', lw = 0.5)

    args['pulse_time'] = t_array[final_prob_list == final_prob_list.max()][0]
    
    if save_pulse: model.op_drive_params_dict[op_name] = copy.deepcopy(args)

global_calibrator_dict = {}
global_calibrator_dict['Guassian'] = Calibrate_Guassian_Pulse
global_calibrator_dict['Sine_Squared'] = Calibrate_Sine_Squared_Pulse
global_calibrator_dict['Ramped'] = Calibrate_Ramped_Pulse

def Floquet_Modes_Table(model, times, drive_args, do_check=True):
    H_d_floq =  model.Drive_Hamiltonian(drive_args, floquet = True)
    
    omega = (2*np.pi*(drive_args['freq_d']+drive_args['shift']))
    T = np.abs(2*np.pi/omega)
    args = {'w': omega}

    f_modes_0, f_energies = qp.floquet.floquet_modes(H_d_floq, T = T, args=args)

    f_modes_table = qp.floquet.floquet_modes_table(f_modes_0, f_energies, times, H = H_d_floq, T = T, args = args)

    f_modes_t = []
    for t in times:
        f_modes_t.append(qp.floquet_modes_t_lookup(f_modes_table, t, T))
    
    if do_check:
        f_modes_tf = qp.floquet.floquet_modes_t(f_modes_0, f_energies, times[-1], H_d_floq, T, args)
        temp = []
        for i in range(len(f_modes_tf)):
            temp.append(qp.expect(f_modes_t[-1][i].proj(), f_modes_tf[i]))

        print(f'Average overlap between Table Derived final floquet modes and fully simulated: {np.sum(np.array(temp))/len(temp)}')

    return dict(modes_t = f_modes_t, mode_table = f_modes_table, energies = f_energies, f_modes_0 = f_modes_0)

def Floquet_Modes_Time(model, times, drive_args, do_check):
    H_d_floq = model.Drive_Hamiltonian(drive_args, floquet = True)
    
    omega = (2*np.pi*(drive_args['freq_d']+drive_args['shift']))
    T = np.abs(2*np.pi/omega)
    args = {'w': omega}

    f_modes_0, f_energies = qp.floquet.floquet_modes(H_d_floq, T = T, args=args)

    
    f_modes_t = []
    for i in range(len(times)):
        t = times[i]
        print(f'Doing time step {i}/{len(times)}, t = {t}')
        f_modes_t.append(qp.floquet.floquet_modes_t(f_modes_0, f_energies, t, H_d_floq, T, args))
    
    return dict(modes_t = f_modes_t, energies = f_energies)

def Floquet_Mode_Overlaps(floquet_dat, states, times): 
    coords = dict(quasienergy = floquet_dat['energies'], time = times)
    floquet_da = xr.DataArray(dims = list(coords.keys()), coords = coords, name = 'Overlaps')

    for i in range(floquet_da.coords['time'].shape[0]):
        t = floquet_da.coords['time'].values[i]
        for j in range(floquet_da.coords['quasienergy'].shape[0]):
            energy = floquet_da.coords['quasienergy'].values[j]

            floquet_da.loc[{'quasienergy':energy, 'time':t}] = qp.expect(floquet_dat['modes_t'][i][j].proj(), states[i])

    return floquet_da

def GET_TIME():
    return datetime.now().strftime('%m-%d-%y:%H:%M:%S')

def Probability_From_Psi_List(psi_list, state):
    res = []
    proj = state.proj()
    for psi in psi_list:
        res.append(qp.expect(proj, psi))
    
    return np.concatenate(np.array(res))
    
def TrackStates(states, debug = False, extra_sorts=None, states_to_track = [], states_to_track_keys = []):
    '''states: Array where first index is a slice and the second indexes the states in that slice. Extra Sorts needs to be an array'''

    states = states.copy()
    steps = len(states)
    N = len(states[0])

    if extra_sorts is not None:
        if len(extra_sorts.shape) == 2:
            extra_sorts = np.array([extra_sorts])
    indexes = []
    overlaps = []
    for i in range(1, steps):
        ordered_temp = np.array([None]*N)
        overlaps_temp = np.array([0.0]*N)
        for pair in itert.product(np.arange(0,N), repeat = 2):
            overlap = qp.expect(states[i-1][pair[0]].proj(), states[i][pair[1]])

            if debug: print(f'Pair: {pair}, Overlap: {overlap}\n  Current Temp Overlap: {overlaps_temp[pair[1]]}')
            
            if overlap > overlaps_temp[pair[1]]:
                if debug: print(f'     Updating Overlap')
                overlaps_temp[pair[1]] = overlap
                ordered_temp[pair[1]] = pair[0]
                if debug: print(f'      New Overlap: {overlaps_temp[pair[0]]}, Index: {ordered_temp[pair[1]]}')
            if debug: print(f'\n    Ordered_temp: {ordered_temp}')
            
            if debug:
                if pair[1] == N-1:
                    print(f'------------------------------------------------------------------------\n')
                
        overlaps.append(overlaps_temp.copy())
        if type(states[i]) != np.ndarray:
            states_temp = np.empty(N, dtype = object)
            for j in range(N):
                states_temp[j] = states[i][j]
            
            states[i] = states_temp.copy()
            
        
        idx = ordered_temp.argsort()
        indexes.append(idx)
        if debug: print(f'idx: {idx}')
        
        states[i] = states[i][idx]

        if extra_sorts is not None:
            for j in range(extra_sorts.shape[0]):
                extra_sorts[j][i] = extra_sorts[j][i][idx]

    to_return = {'states':states, 'indexes':indexes, 'overlaps':overlaps}
    
    if extra_sorts is not None:
        if len(extra_sorts) == 1:
            extra_sorts = extra_sorts[0]
        to_return['extra_sorts'] = extra_sorts

    if len(states_to_track) > 0:
        overlaps_track = np.zeros((len(states_to_track), steps, N))
        projectors = []
        best_fits = np.zeros((len(states_to_track), steps), dtype = int)
        for state in states_to_track:
            projectors.append(state.proj())

        keys = {}
        for pair in itert.product(np.arange(len(states_to_track)), np.arange(0,steps), np.arange(0,N)):
            overlap = qp.expect(projectors[pair[0]], states[pair[1]][pair[2]])
            if debug: print(f'Looking at index {pair}, overlap: {overlap}')
            if overlap > overlaps_track[pair[0], pair[1], :].max():
                if debug: print(f'  Updating Best Fit: Previous max: {overlaps_track[pair[0], pair[1], :].max()}, New Max Location: {pair[2]}')
                best_fits[pair[0], pair[1]] = pair[2] 

                if pair[1] == 0:
                    if len(states_to_track_keys)>0: 
                        transmon = states_to_track_keys[pair[0]][0]
                        if type(transmon) != str:
                            transmon = Transmon_Labels_NK[transmon]
                        key = str(transmon)+str(states_to_track_keys[pair[0]][1])
                        if debug: print(f'      Key: {key}, loc: {pair[2]}')
                        keys[key] = pair[2]

            if debug:
                if pair[2] == N-1:
                    print(f'------------------------------------------------------------------------\n')
            overlaps_track[pair[0], pair[1], pair[2]] = overlap
        

        to_return['tracked_state_locs'] = best_fits
        to_return['tracked_state_overlaps'] = overlaps_track
        to_return['tracked_state_keys'] = keys

    return to_return

def Floquet_0_Sweep(model, list_of_drive_args, options = None, debug = False):
    floquet_modes = []
    floquet_modes_energies = []
    for i in range(len(list_of_drive_args)):
        drive_args = list_of_drive_args[i]
        if debug: print(f'Doing args {i+1}/{len(list_of_drive_args)}\n args: {drive_args}')

        H_d_floq = model.Drive_Hamiltonian(drive_args, floquet = True)

        if 'shift' not in drive_args:
            drive_args['shift'] = 0
        omega = (2*np.pi*(drive_args['freq_d']+drive_args['shift']))
        T = np.abs(2*np.pi/omega)
        args = {'w': omega}
        
        res= (qp.floquet.floquet_modes(H_d_floq, T = T, args=args, options=options))

        floquet_modes.append(res[0])
        floquet_modes_energies.append(res[1])

    floquet_modes_array = np.empty((len(floquet_modes), len(floquet_modes[0])), dtype = object)
    floquet_modes_energies_array = np.empty((len(floquet_modes_energies), len(floquet_modes_energies[0])))

    for idx, i in np.ndenumerate(floquet_modes_array):
        floquet_modes_array[idx] = floquet_modes[idx[0]][idx[1]]
        floquet_modes_energies_array[idx] = floquet_modes_energies[idx[0]][idx[1]]

    return dict(modes = floquet_modes_array, energies= floquet_modes_energies_array) 

def Plot_State_Evolution_Single_Mode(model, psi_list, x_axis = None, x_ticks = None, use_all_states = False, emphasize_states = 'All', title = None, x_label = None, legend_args = {}, figsize = (10,5), dpi = 200, markersize = 100, lw = 0.2, debug  =False, plots = 'Probability', emphasis_alpha = 0.1, small_amp_cutoff = 1e-3, relative_phase = 'q0'):
        ''' -model: An Instance of 'Transmon_Cavity_Model'
            -psi_list: List of wavefunctions to plot
            -x_axis: List/array of x values. Must be same size as psi_list
            -x_ticks: list/array of values to place the x_ticks
            -use_all_states: Bool, if True, then it will use all the states in the model, not just the surface N_Transmon and N_Cavity
            -emphasize_states: str/list. If string then it will emphasize that state (or all states if 'All'). Can also be a list of states in the form 'f0' (Transmon Mode then Cavity Mode)
            -title: Plot title
            -x_label: label for x axis. Default is 'Time [ns]'
            -legend_args: kwargs for legend
            plots: Takes strings: 'Probability', 'Amplitude', 'Phase' or 'All'. Can also take a list of any number of those
            emphasis_alpha: the alpha value for all the non-emphasized states
            small_amp_cutoff: this number is min value to show up on the plot. So if there are any points where the probability of a state goes below 1e-3 then the alpha at those points is set to zero. Helps with clarity!
        '''
        
        N_Transmon = model.N_Transmon
        N_Cavity = model.N_Cavity[0]
        if use_all_states:
            N_Transmon = model.transmon_truncated_dim
            N_Cavity = model.resonator_truncated_dim[0]
        
        if plots == 'All':
            plots = ['Probability', 'Amplitude', 'Phase']#, 'Final State Probability']
        if type(plots) == str:
            plots = [plots]

        fig, ax = plt.subplots(figsize = figsize, dpi = dpi, nrows = len(plots), layout = 'compressed')

        if len(plots) == 1:
            ax = np.array([ax], dtype = object)

        photon_numbers = np.arange(0, N_Cavity)
        cmap = cmaps.cet_r_bgyrm

        cmap_for_plot = cmap(photon_numbers/photon_numbers.max())

        norm = colors.BoundaryNorm(np.concatenate([photon_numbers, [N_Cavity]])-0.5, cmap.N)
        cmap_for_cbar = cm.ScalarMappable(cmap = cmap, norm = norm)


        if type(emphasize_states) == str:
            emphasize_states = [emphasize_states]


        y_dict = dict(Probability = {}, Amplitude = {}, Phase = {})
        
        for i in range(N_Transmon):
                y_dict['Probability'][Transmon_Labels_NK[i]] = {}
                y_dict['Amplitude'][Transmon_Labels_NK[i]] = {}
                y_dict['Phase'][Transmon_Labels_NK[i]] = {}
        

        if (psi_list[0].type) == 'ket':
            for pair in itert.product(np.arange(N_Transmon), np.arange(N_Cavity)):
                i = pair[0]
                j = pair[1]
                
                ref_state = model.get_dressed_state([i,j])
                z = np.array([ref_state.dag()*state for state in psi_list]).flatten()
                y_dict['Amplitude'][Transmon_Labels_NK[i]][j] = np.abs(z)
                y_dict['Phase'][Transmon_Labels_NK[i]][j] = np.angle(z)
                y_dict['Probability'][Transmon_Labels_NK[i]][j] = np.abs(z)**2
        
        if (psi_list[0].type) == 'oper':
            for pair in itert.product(np.arange(N_Transmon), np.arange(N_Cavity)):
                i = pair[0]
                j = pair[1]
                
                ref_state_proj = model.get_dressed_state([i,j]).proj()
                z = np.array([qp.expect(ref_state_proj, state) for state in psi_list]).flatten()
                y_dict['Probability'][Transmon_Labels_NK[i]][j] = np.abs(z)
        
        for plot_num in range(len(plots)):
            plot = plots[plot_num]
    
            if debug: print(f'Plotting {plot}')
            for pair in  itert.product(np.arange(N_Transmon), np.arange(N_Cavity)):
                i = pair[0]
                j = pair[1]
                trans_mode = Transmon_Labels_NK[i]
                cavity_mode = j

                label = None
                if j == 0:
                    label = f'|{trans_mode},i>'
                

                y = y_dict[plot][trans_mode][cavity_mode]
                if plot == 'Phase':
                    y = (y-y_dict[plot][relative_phase[0]][int(relative_phase[1])])%(2*np.pi)
                    y[y>=np.pi] += -2*np.pi
                
                x = x_axis
                if x is None:
                    x = np.arange(len(y))                
                
                alpha = emphasis_alpha
                if 'All' in emphasize_states:
                    alpha = 1
                
                if [i,j] in emphasize_states or [trans_mode,j] in emphasize_states or f'{trans_mode}{j}' in emphasize_states:
                    alpha = 1

                ax[plot_num].plot(x, y, color = cmap_for_plot[j], lw = lw, alpha = alpha, zorder = -100)
                
                alpha_array = np.ones_like(x)*alpha

                if small_amp_cutoff != None:
                    alpha_array[np.abs(y_dict['Probability'][trans_mode][cavity_mode])<small_amp_cutoff] = 0

                color = np.array([cmap_for_plot[j]]*len(alpha_array))
                alpha_array[0] = alpha
                #print(color.shape)
                if marker_dict[trans_mode] in ['0', '1', '2', '3', '+', 'x']:
                    ax[plot_num].scatter(x, y, marker=marker_dict[trans_mode], c = color, label=label, s=markersize, alpha = alpha_array)
                
                if marker_dict[trans_mode] in ['^', 'v', '<', '>']:
                    g = ax[plot_num].scatter(x, y, marker=marker_dict[trans_mode], edgecolor = color, facecolor = ['none']*len(color), label=label, s=markersize, alpha = alpha_array)
                    g.set_facecolor('none')
                
            if title == None:
                title = 'System Evolution'
            
            if plot_num == 0:
                ax[plot_num].set_title(title)
            ax[plot_num].set_xlabel(x_label)
            ax[plot_num].set_ylabel(plot)
            if x_ticks is None:
                x_ticks =x
            ax[plot_num].set_xticks(x_ticks)

            if plot == 'Probability':
                ax[plot_num].set_ylim([-0.1, 1.1])
                ax[plot_num].hlines(y = 0.5, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
                ax[plot_num].hlines(y = 1, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)

            if plot == 'Amplitude':
                ax[plot_num].set_ylim([-0.1, 1.1])
                ax[plot_num].hlines(y = np.sqrt(0.5), xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
                ax[plot_num].hlines(y = 1, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
            
            if plot == 'Phase':
                ax[plot_num].set_ylim([-0.25-np.pi, 0.25+np.pi])
                ax[plot_num].hlines(y = -np.pi, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
                ax[plot_num].hlines(y = np.pi, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
                ax[plot_num].set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
                ax[plot_num].set_yticklabels(labels = [r'$-\pi$', r'$-\pi/2$', 0, r'$\pi/2$', r'$\pi$'])


        
        
        ax[0].legend(**legend_args)

        cbar = fig.colorbar(cmap_for_cbar, label = 'Photon Number n', ax = ax.ravel().tolist())
        cbar.ax.set_yticks(photon_numbers)
        cbar.ax.set_yticklabels(photon_numbers)
        plt.show()

def Plot_State_Evolution_Two_Modes(model, psi_list, x_axis = None, x_ticks = None, use_all_states = False, emphasize_states = 'All', title = None, x_label = None, legend_args = {}, figsize = (10,5), dpi = 200, markersize = 100, lw = 0.2, debug  =False, plots = 'Probability', emphasis_alpha = 0.1, small_amp_cutoff = 1e-3, relative_phase = 'q0', marker_lw = 0.1):
        ''' -model: An Instance of 'Transmon_Cavity_Model'
            -psi_list: List of wavefunctions to plot
            -x_axis: List/array of x values. Must be same size as psi_list
            -x_ticks: list/array of values to place the x_ticks
            -use_all_states: Bool, if True, then it will use all the states in the model, not just the surface N_Transmon and N_Cavity
            -emphasize_states: str/list. If string then it will emphasize that state (or all states if 'All'). Can also be a list of states in the form 'f0' (Transmon Mode then Cavity Mode)
            -title: Plot title
            -x_label: label for x axis. Default is 'Time [ns]'
            
            -legend_args: kwargs for legend
            plots: Takes strings: 'Probability', 'Amplitude', 'Phase' or 'All'. Can also take a list of any number of those
            emphasis_alpha: the alpha value for all the non-emphasized states
            small_amp_cutoff: this number is min value to show up on the plot. So if there are any points where the probability of a state goes below 1e-3 then the alpha at those points is set to zero. Helps with clarity!
        '''

        marker_dict = {'g':'v', 'e':'^', 'f':'P', 'h':'X', 'i':'*', 'j':'D', 'k':'8', 'l':'s', 'm':"H", 'o':'h', 'p':'d'}
        
        N_Transmon = model.N_Transmon
        N_Cavity = model.N_Cavity
        if use_all_states:
            N_Transmon = model.transmon_truncated_dim
            N_Cavity = model.resonator_truncated_dim
        
        if plots == 'All':
            plots = ['Probability', 'Amplitude', 'Phase']#, 'Final State Probability']
        if type(plots) == str:
            plots = [plots]

        fig, ax = plt.subplots(figsize = figsize, dpi = dpi, nrows = len(plots), layout = 'compressed')

        if len(plots) == 1:
            ax = np.array([ax], dtype = object)

        photon_numbers = np.arange(0, np.array(N_Cavity).max())
        cmap = cmaps.cet_r_bgyrm

        cmap_for_plot = cmap(photon_numbers/photon_numbers.max())

        norm = colors.BoundaryNorm(np.concatenate([photon_numbers, [photon_numbers[-1]+1]])-0.5, cmap.N)
        cmap_for_cbar = cm.ScalarMappable(cmap = cmap, norm = norm)


        if type(emphasize_states) == str:
            emphasize_states = [emphasize_states]


        y_dict = dict(Probability = {}, Amplitude = {}, Phase = {})
        
        for i in range(N_Transmon):
                y_dict['Probability'][Transmon_Labels_NK[i]] = {}
                y_dict['Amplitude'][Transmon_Labels_NK[i]] = {}
                y_dict['Phase'][Transmon_Labels_NK[i]] = {}
        

        if (psi_list[0].type) == 'ket':
            for pair in itert.product(np.arange(N_Transmon), np.arange(N_Cavity[0]), np.arange(N_Cavity[1])):
                i = pair[0]
                j = pair[1]
                k = pair[2]
                index = str(j)+','+str(k)

                ref_state = model.get_dressed_state([i,j,k])

                z = np.array([ref_state.dag()*state for state in psi_list]).flatten()
                y_dict['Amplitude'][Transmon_Labels_NK[i]][index] = np.abs(z)
                y_dict['Phase'][Transmon_Labels_NK[i]][index] = np.angle(z)
                y_dict['Probability'][Transmon_Labels_NK[i]][index] = np.abs(z)**2
        
        if (psi_list[0].type) == 'oper':
            for pair in itert.product(np.arange(N_Transmon), np.arange(N_Cavity[0]), np.arange(N_Cavity[1])):
                i = pair[0]
                j = pair[1]
                k = pair[2]
                index = str(j)+','+str(k)
                
                ref_state_proj = model.get_dressed_state([i,j,k]).proj()
                z = np.array([qp.expect(ref_state_proj, state) for state in psi_list]).flatten()
                y_dict['Probability'][Transmon_Labels_NK[i]][index] = np.abs(z)
        
        for plot_num in range(len(plots)):
            plot = plots[plot_num]
    
            if debug: print(f'Plotting {plot}')
            for pair in  itert.product(np.arange(N_Transmon), np.arange(N_Cavity[0]), np.arange(N_Cavity[1])):
                i = pair[0]
                j = pair[1]
                k = pair[2]
                trans_mode = Transmon_Labels_NK[i]
                cavity_mode = str(j)+','+str(k)

                y = y_dict[plot][trans_mode][cavity_mode]
                if plot == 'Phase':
                    y = (y-y_dict[plot][relative_phase[0]][int(relative_phase[1])])%(2*np.pi)
                    y[y>=np.pi] += -2*np.pi
                
                x = x_axis
                if x is None:
                    x = np.arange(len(y))                
                
                alpha = emphasis_alpha
                if 'All' in emphasize_states:
                    alpha = 1
                
                if [i,j] in emphasize_states or [trans_mode,j] in emphasize_states or f'{trans_mode}{j}' in emphasize_states:
                    alpha = 1

                f_color_1 = cmap_for_plot[j]
                f_color_2 = cmap_for_plot[k]

                line_color = np.sqrt((f_color_1**2+f_color_2**2)/2)

                ax[plot_num].plot(x, y, color = line_color, lw = lw, alpha = alpha, zorder = -100)
                
                alpha_array = np.ones_like(x)*alpha

                if small_amp_cutoff != None:
                    alpha_array[np.abs(y_dict['Probability'][trans_mode][cavity_mode])<small_amp_cutoff] = 0

                alpha_array[0] = alpha
                #print(color.shape)
                
            
                f_color_1_arr = np.array([f_color_1]*len(alpha_array))
                f_color_2_arr = np.array([f_color_2]*len(alpha_array))

                marker = MarkerStyle(marker_dict[trans_mode], fillstyle="left")
                ax[plot_num].scatter(x, y, marker=marker, c = f_color_1_arr, edgecolor = 'black', label=None, s=markersize, alpha = alpha_array, linewidth = marker_lw)
                
                marker = MarkerStyle(marker_dict[trans_mode], fillstyle="right")
                ax[plot_num].scatter(x, y, marker=marker, c = f_color_2_arr, edgecolor = 'black', label=None, s=markersize, alpha = alpha_array, linewidth = marker_lw)

                label = None
                if j == 0 and k == 0:
                    label = f'|{trans_mode},i,j>'
                    plt.plot([],[], label = label, color = f_color_1, marker = marker_dict[trans_mode], lw = 0)
                
            if title == None:
                title = 'System Evolution'
            
            if plot_num == 0:
                ax[plot_num].set_title(title)
            ax[plot_num].set_xlabel(x_label)
            ax[plot_num].set_ylabel(plot)
            if x_ticks is None:
                x_ticks =x
            ax[plot_num].set_xticks(x_ticks)

            if plot == 'Probability':
                ax[plot_num].set_ylim([-0.1, 1.1])
                ax[plot_num].hlines(y = 0.5, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
                ax[plot_num].hlines(y = 1, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)

            if plot == 'Amplitude':
                ax[plot_num].set_ylim([-0.1, 1.1])
                ax[plot_num].hlines(y = np.sqrt(0.5), xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
                ax[plot_num].hlines(y = 1, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
            
            if plot == 'Phase':
                ax[plot_num].set_ylim([-0.25-np.pi, 0.25+np.pi])
                ax[plot_num].hlines(y = -np.pi, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
                ax[plot_num].hlines(y = np.pi, xmin = x.min(), xmax = x.max(), color = 'black', ls = ':', lw = 0.5)
                ax[plot_num].set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
                ax[plot_num].set_yticklabels(labels = [r'$-\pi$', r'$-\pi/2$', 0, r'$\pi/2$', r'$\pi$'])


        
        
        ax[0].legend(**legend_args)

        cbar = fig.colorbar(cmap_for_cbar, label = 'Photon Number n', ax = ax.ravel().tolist())
        cbar.ax.set_yticks(photon_numbers)
        cbar.ax.set_yticklabels(photon_numbers)
        plt.show()

def FindResonanceFloquet(model, state_i, state_f, freq_d, shifts, epsilon, show_plot = False, fig_args = {}, debug = False, options = None, jump_thresh = 10, save_plot = False, save_path = None):

    list_of_drive_args = []
    for shift in shifts:
        list_of_drive_args.append(dict(freq_d = freq_d, shift = shift, epsilon = epsilon))

    if debug: print('Doing Floquet Sweep\n------------------------------------------------------')
    floq_res = Floquet_0_Sweep(model, list_of_drive_args, options = options, debug = debug)
    if debug: print("------------------------------------------------------\nFloquet Sweep Done")
    states_to_track_names = [state_i, state_f]
    states_to_track = [model.get_dressed_state(state_i), model.get_dressed_state(state_f)]

    if debug: print(f'Doing State Tracking\n  (This does not have any updates, it just takes a while)')
    tracking_res = TrackStates(floq_res['modes'], extra_sorts = floq_res['energies'], states_to_track = states_to_track, states_to_track_keys = states_to_track_names)
    if debug: print(f'State Tracking Done! Doing the next bit')
    x = np.copy(shifts)
    ys = np.zeros((tracking_res['tracked_state_locs'].shape[0], len(x)))
    for pair in itert.product(np.arange(tracking_res['tracked_state_locs'].shape[0]), np.arange(tracking_res['tracked_state_locs'].shape[1])):
        ref_state_idx = pair[0]
        x_step = pair[1]
        matched_state_idx = tracking_res['tracked_state_locs'][ref_state_idx, x_step]
        
        ys[pair] = tracking_res['extra_sorts'][x_step, matched_state_idx]/(np.pi)
        

    difs = np.abs(ys[1]-ys[0]) 
    if debug: print(f'Initial Difs: {difs}')
    for i in range(len(difs)-1):
        if debug: print(f'difs[i+1]: {difs[i+1]}, jump_thresh*difs[i]: {jump_thresh*difs[i]}')
        if difs[i+1]>jump_thresh*difs[i]:
            if debug: print(f'    difs[i+1]: {difs[i+1]}')
            difs[i+1] = np.abs(difs[i+1] + 2*(x[i+1]+freq_d))
            if debug: print(f'    difs[i+1]: {difs[i+1]}')

    difs = np.abs(difs)
    if debug: print(f'Final Difs: {difs}')

    if show_plot:
        if 'figsize' not in fig_args:
            fig_args['figsize'] = [10,5]
        if 'dpi' not in fig_args:
            fig_args['dpi'] = 200
        
        fig, ax1 = plt.subplots(**fig_args)
        
        for i in range(ys.shape[0]):
            if i == 0:
                marker = 'x'
            if i == 1:
                marker = '+'
            ax1.plot(x, ys[i], lw = 1, marker = marker, ms = 7, label = BraKetLabel(states_to_track_names[i]))
        
        
        ax2 = ax1.twinx()
        ax2.plot(x, difs, color = 'black', label = '|difs|', alpha = 0.5, marker = 'x')
        ax2.set_ylabel('|Difs|')

    
    to_fit = lambda x, x_shift, y_shift, slope: slope*np.sqrt((x-x_shift)**2+y_shift**2)
    x_shift = x[np.argmin(difs)]
    y_shift = np.min(difs)
    slope = np.abs((np.max(difs) - np.min(difs))/(x[np.argmax(difs)]-x[np.argmin(difs)]))
    p0 = [x_shift, y_shift, slope]
    
    if debug: print(f'Fit Guesses: {p0}') 

    fit_dat = opt.curve_fit(to_fit, x, difs, p0 = p0) 

    

    if show_plot:
        y = []
        x_fit = np.linspace(x[0], x[-1], 51)
        for i in range(len(x_fit)):
            y.append(to_fit(x_fit[i], fit_dat[0][0], fit_dat[0][1], fit_dat[0][2]))
        
        ax2.plot(x_fit,y, label = 'dif fit', color = 'green', lw = 5, alpha = 0.25)

        x_loc = x[np.argmin(difs)]#-(x[1]-x[0])/2
        ax1.axvline(x_loc, color = 'black', lw = 1, ls = '--', label = f'Cossing at {x_loc}')


    to_minimize = lambda x: to_fit(x, fit_dat[0][0], fit_dat[0][1], fit_dat[0][2])

    min_dat = opt.minimize(to_minimize, x0 = np.min(difs))

    if show_plot:
        min_loc = min_dat['x'][0]
        ax1.axvline(min_loc, color = 'black', lw = 1, label = f"Fittied Shift {min_loc}")

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1+lines2, labels1+labels2)

        ax1.set_xlabel(r'Frequency Shift [GHz]')
        ax1.set_ylabel(r'$\frac{1}{\pi}\omega_Q$ [GHz]')
        ax1.set_title(r'$\text{Fitting Stark Shift: }\epsilon = EPS$ [GHz]'.replace('EPS', str(epsilon)))
        
        if save_plot:
            if save_path == None:
                save_path = model.save_path
            fig.savefig(save_path+'FloquetResonanceSweep_'+GET_TIME()+'.pdf')

        plt.show()


    return min_dat

def StarkShiftFitter(model, state_i, state_f, freq_d, shifts, epsilons, fit_order = 8, func_kwargs = {}, print_progress = True, return_mins = False):
    stark_shift_mins = []
    for i in range(len(epsilons)):
        epsilon = epsilons[i]
        if print_progress: print(f'Doing step {i+1}/{len(epsilons)}')
        stark_shift_mins.append(FindResonanceFloquet(model, state_i, state_f, freq_d, np.copy(shifts), epsilon, show_plot = True))
    
    stark_shifts = []
    for i in range(len(stark_shift_mins)):
        stark_shifts.append(stark_shift_mins[i]['x'][0])

    fit = np.polyfit(epsilons, stark_shifts, deg = 8)

    if return_mins:
        return dict(fit = fit, stark_shift_mins = stark_shift_mins)

    return fit

def GetCollapseAndDephasing(model, T_Kappa_C = 1/(50*1000), T_Kappa_D = 1/(192.5*1000), C_Kappa_C = 1/(1*1000*1000)):
    ''' T_Kappa_C: Kappa for Transmon Collapse
        T_Kappa_D: Kappa for Transmon Dephasing
        C_Kappa_C: Kappa for Cavity Collapse
    '''
    Transmon_Collapse = 0
    Transmon_Dephasing = 0
    Cavity_Collapse = 0
    
    for pair in itert.product(np.arange(model.transmon_truncated_dim-1), np.arange(model.resonator_truncated_dim)):
        j = pair[0]
        n = pair[1]
        Transmon_Collapse += np.sqrt(T_Kappa_C)*np.sqrt(j+1)*model.get_dressed_state([j,n])*model.get_dressed_state([j+1, n]).dag()
        Transmon_Dephasing += np.sqrt(T_Kappa_D)*np.sqrt(j)*model.get_dressed_state([j, n])*model.get_dressed_state([j, n]).dag()
    
    for pair in itert.product(np.arange(model.transmon_truncated_dim), np.arange(model.resonator_truncated_dim-1)):
        j = pair[0]
        n = pair[1]
        Cavity_Collapse += np.sqrt(C_Kappa_C)*np.sqrt(n+1)*model.get_dressed_state([j, n])*model.get_dressed_state([j, n+1]).dag()
    
    return dict(T_C = Transmon_Collapse, T_D = Transmon_Dephasing, C_C = Cavity_Collapse)

class Transmon_Cavity_Model:

    def __init__(self, E_J=27, E_C=0.11, E_osc=5.75, g=0.028, \
                            transmon_ng = 0, transmon_ncut = 60,  \
                            transmon_truncated_dim = None, \
                            resonator_truncated_dim = None, N_Transmon = 4, N_Cavity = 5, model_name = None, save_path = None):
        
        '''All constants with units  in GHz'''
        self.E_J = E_J
        self.E_C = E_C
        self.E_osc = E_osc

        if isinstance(g, (int, np.floating, float)): 
            g = [g]
        
        if isinstance(E_osc, (int, np.floating, float)): 
            E_osc = [E_osc]
        
        if isinstance(N_Cavity, (int, np.floating, float)): 
            N_Cavity = [N_Cavity]*len(g)
        
        if isinstance(resonator_truncated_dim, (int, np.floating, float, type(None))): 
            resonator_truncated_dim = [resonator_truncated_dim]*len(N_Cavity)

        


        self.N_Res = len(E_osc)

        if len(g) != len(E_osc) and len(E_osc)!=len(N_Cavity):
            print('E_osc, g and N_Cavity Need to be the same length. Please try again!')
            return None

        self.g = copy.deepcopy(g)

        self.N_Transmon = N_Transmon
        self.N_Cavity = copy.deepcopy(N_Cavity)

        self.transmon_ng = transmon_ng
        self.transmon_ncut = transmon_ncut
        self.resonator_truncated_dim = copy.deepcopy(resonator_truncated_dim)
        
        for i in range(len(self.resonator_truncated_dim)):
            if self.resonator_truncated_dim[i] == None:
                self.resonator_truncated_dim[i] = N_Cavity[i]+2


        self.transmon_truncated_dim = transmon_truncated_dim
        if self.transmon_truncated_dim == None:
            self.transmon_truncated_dim = self.N_Transmon+2

        self.model_name = model_name
        if self.model_name == None:
            self.model_name = 'model_'+str(np.random.randint(100000))
        
        self.save_path = save_path

        self.Extra_Stuff = {}

        self.model_params = {'E_J':self.E_J, 'E_C':self.E_C, 'E_osc':self.E_osc, 'g':self.g,\
            'N_Transmon':self.N_Transmon, 'N_Cavity':self.N_Cavity, \
                'transmon_ng':self.transmon_ng, 'transmon_ncut':self.transmon_ncut, 'transmon_truncated_dim':self.transmon_truncated_dim, \
                    'resonator_truncated_dim':self.resonator_truncated_dim, 'model_name':self.model_name, 'save_path':self.save_path}

        
        self.n_zpt =  (E_J/32/E_C)**(1/4.)


        self.transmon =  scq.Transmon(EJ = E_J, EC = E_C, ng = transmon_ng,\
                                       ncut = transmon_ncut, truncated_dim=self.transmon_truncated_dim)
        
        self.resonators = []
        for i in range(len(E_osc)):
            self.resonators.append(scq.Oscillator(E_osc=E_osc[i], truncated_dim=self.resonator_truncated_dim[i]))

        self.hilbertspace = scq.HilbertSpace([self.transmon, *self.resonators])

        for i in range(len(self.resonators)):
            self.hilbertspace.add_interaction(g=self.g[i],  op1=self.transmon.n_operator, op2 = self.resonators[i].annihilation_operator, add_hc = True)

        self.hilbertspace.generate_lookup()

        self.Transmon_Labels_LK = Transmon_Labels_LK

        self.Transmon_Labels_NK = Transmon_Labels_NK
        

        self.op_drive_params_dict = {}

        #self.Transmon_Cavity_Abstract = Transmon_Cavity_Abstract(N_Transmon=self.N_Transmon, N_Cavity=self.N_Cavity)

        self.Shift_Search_Dict = {}

        self.Drive_Sequences = {}
    
    def update_model_params(self):
        self.model_params = {'E_J':self.E_J, 'E_C':self.E_C, 'E_osc':self.E_osc, 'g':self.g,\
            'N_Transmon':self.N_Transmon, 'N_Cavity':self.N_Cavity, \
                'transmon_ng':self.transmon_ng, 'transmon_ncut':self.transmon_ncut, 'transmon_truncated_dim':self.transmon_truncated_dim, \
                    'resonator_truncated_dim':self.resonator_truncated_dim, 'model_name':self.model_name, 'save_path':self.save_path}
    
    def get_dressed_state(self, ind):
        if type(ind[0]) == str:
            ind[0] = self.Transmon_Labels_LK[ind[0]]
        j = self.hilbertspace.dressed_index((tuple(ind)))
        return self.hilbertspace.eigensys(self.hilbertspace.dimension)[1][j]
    

    def DefaultFrequency(self, state1, state2):
        if type(state1[0]) == str:
            state1[0] = self.Transmon_Labels_LK[state1[0]]
        
        if type(state2[0]) == str:
            state2[0] = self.Transmon_Labels_LK[state2[0]]
        
        E0 = self.hilbertspace.energy_by_bare_index(tuple(state1))
        E1 = self.hilbertspace.energy_by_bare_index(tuple(state2))
        
        return E1-E0

    def get_proj_op(self, state):
        if type(state[0]) == str:
            state[0] = self.Transmon_Labels_LK[state[0]]
        
        return self.get_dressed_state(state).proj()

    def Drive_Hamiltonian(self, args, identity = False, envelope = None, floquet = False):
        
        if floquet:        
            freq = args['freq_d']+args['shift']
            eps = args['epsilon']
            
            drive_coeff = 2*np.pi*eps
            n_op = scq.identity_wrap(self.transmon.n_operator(), self.transmon, [self.transmon, *self.resonators])
            H_d = [2*np.pi*self.hilbertspace.hamiltonian(), [drive_coeff*n_op, 'sin(w * t)']]
            return H_d

        eps = args['epsilon']
        if isinstance(eps, (int, np.floating, float)): 
            eps = [eps]
        
        freqs = args['freq_d']
        shifts = args['shift']
        if isinstance(freqs, (int, np.floating, float)): 
            freqs = [freqs]
        if isinstance(shifts, (int, np.floating, float)): 
            shifts = [shifts]
        
        freqs = np.array(freqs)
        shifts = np.array(shifts)
        freqs = freqs+shifts

        envelope_args = {}
        if 'Envelope Args' in args:
            envelope_args = copy.deepcopy(args['Envelope Args'])

        if 'Envelope' in args:
            if 'pulse_time' in inspect.getfullargspec(global_envelope_dict[args['Envelope']]).args:
                envelope_args['pulse_time'] = args['pulse_time']
            
        if envelope == None:
            envelope = lambda t: global_envelope_dict[args['Envelope']](t, **envelope_args)
            
        
        if 'phi0' in args:
            phi0 = args['phi0']

        if 'phi0' not in args:
            phi0 = 0
                    
        if isinstance(phi0, (int, np.floating, float)): 
            phi0 = [phi0]

        
    
        drive_coeff = lambda t, args: envelope(t)*2*np.pi*np.sum(eps[i]*np.sin(2*np.pi*(freqs[i])*t + phi0[i]) for i in range(len(eps)))

        n_op = scq.identity_wrap(self.transmon.n_operator(), self.transmon, [self.transmon, *self.resonators])
        H_d = [2*np.pi*self.hilbertspace.hamiltonian(), [n_op, drive_coeff]]

        if identity:
            drive_coeff = lambda t, args: t

            n_op = scq.identity_wrap(np.identity(len(self.transmon.n_operator())), self.transmon, [self.transmon, self.resonator])
            H_d = [2*np.pi*self.hilbertspace.hamiltonian(), [n_op, drive_coeff]] # This 2pi is from some factor of hbar vs h

        return H_d
    
    def Drive_State(self,psi0, args, spps = 100, starttime = 0, t_list = None, Odeoptions = None, progress_bar = True,c_ops = None, identity = False, envelope = None, H_d = None):

        if Odeoptions == None:
            Odeoptions = {}
        
        Odeoptions_Default = {'atol': 1e-12, 'rtol':1e-12, 'max_step': 1/(20*100)}
        
        for key in Odeoptions_Default:
            if key not in Odeoptions:
                Odeoptions[key] = Odeoptions_Default[key]

        
        if H_d == None:
            H_d = self.Drive_Hamiltonian(args, identity = identity, envelope = envelope)
        

        if t_list is None:
            t_list = tlist(args['pulse_time'], starttime = starttime, spps = spps)
    
        output = qp.mesolve(H = H_d,
				  	rho0=psi0,
					tlist=t_list,#e_ops=e_ops,
					args=args,
					options=qp.Odeoptions(**Odeoptions),
					progress_bar=progress_bar,
                    c_ops = c_ops)
        
        return output

    def Run_Pulse_Sequence(self, psi0, op_list, spps = 10, save_history = False, Odeoptions = None, c_ops = None, identity = False):
        results = []
        psi_list = []
        psi_list.append(psi0)

        psi = psi_list[0]
        total_time = np.zeros(0)
        t_list = []
        timestamps = {}
        t = 0
        for i in range(len(op_list)):
            phi = 0
            gate = op_list[i]
            timestamp = f'{i+1}: {gate}'
            
            timestamps[timestamp] = []
            
            args = copy.deepcopy(self.op_drive_params_dict[gate])
            print(f"\n Doing Gate {i+1}/{len(op_list)} ({gate}), Gate Time: {args['pulse_time']}\n-------------------------------------------------------------------")
            res = self.Drive_State(psi, args, spps=spps, Odeoptions=Odeoptions, c_ops=c_ops, identity=identity)
            psi = res.states[-1]
            results.append(res)

            t_list.append(res.times+t)
            timestamps[timestamp].append(t)
            t+= args['pulse_time']
            timestamps[timestamp].append(t)
        
        psi_list = []
        for i in range(len(results)):
            for j in range(len(results[i].states)):
                psi_list.append(results[i].states[j])
        
        res = [np.concatenate(np.array(t_list)), psi_list, op_list, timestamps]

        if  type(save_history)==str:
            with open(save_history+'.pkl', 'wb') as file:
                pickle.dump(res, file)
        return res
    
    def Plot_State_Evolution(self, t_list, psi_list, fig_kwargs = {}, debug = False, timestamps = [], plot_every = 1):
        t_list = copy.deepcopy(t_list)
        psi_list = copy.deepcopy(psi_list)
        if len(timestamps)>0:
            new_psi_list = []
            new_t_list = []
            for timestamp in timestamps:
                temp_t_list = np.array(copy.copy(t_list))
                temp_t_list = temp_t_list-timestamp[0]
                temp_t_list2 = temp_t_list[temp_t_list<0]
                if len(temp_t_list2) == 0:
                    start = 0
                else:
                    start = np.where(temp_t_list == np.max(temp_t_list2))[0][0]

                temp_t_list = np.array(copy.copy(t_list))
                temp_t_list = temp_t_list-timestamp[1]
                temp_t_list2 = temp_t_list[temp_t_list>0]
                if len(temp_t_list2) == 0:
                    stop = len(temp_t_list)-1
                else:
                    stop = np.where(temp_t_list == np.min(temp_t_list2))[0][0]

                for i in range(start, stop+1):
                    new_t_list.append(t_list[i])
                    new_psi_list.append(psi_list[i])

            t_list = np.array(new_t_list)
            psi_list = new_psi_list
        
        t_f = t_list[-1]
        psi_f = psi_list[-1]
        t_list = t_list[::plot_every]
        psi_list = psi_list[::plot_every]

        print(list(fig_kwargs.keys()))
        if t_f not in t_list:
            t_list = np.concatenate([t_list, [t_f]])
            psi_list.append(psi_f)

        if debug:
            print(f'len(t_list): {len(t_list)}, len(psi_list): {len(psi_list)}')

        if 'figsize' not in fig_kwargs:
            fig_kwargs['figsize'] = (10,5)

        if 'dpi' not in fig_kwargs:
            fig_kwargs['dpi'] = 200


        if 'title' not in fig_kwargs:
            fig_kwargs['title'] = 'State Evolution'

        if 'x_label' not in fig_kwargs:
            fig_kwargs['x_label'] = 'time [ns]'

        
        if 'x_axis' not in fig_kwargs:
            fig_kwargs['x_axis'] = t_list
        
        
        if 'x_ticks' not in fig_kwargs:
            fig_kwargs['x_ticks'] = np.linspace(t_list[0], t_list[-1], 11)

        if debug:
            print(f'len(t_list): {len(fig_kwargs["x_axis"])}, len(psi_list): {len(psi_list)}')

        if len(self.N_Cavity) == 1:
            Plot_State_Evolution_Single_Mode(self, psi_list = psi_list, **fig_kwargs)
        
        if len(self.N_Cavity) == 2:
            Plot_State_Evolution_Two_Modes(self, psi_list = psi_list, **fig_kwargs)

    def GetStarkShift(self, state1, state2, epsilon, fit_name = None, use_fit = True, make_fit = False, kwargs = {}, save_new_fit = True, freq_d = None):
        if fit_name == None:
            fit_name = f'{state1[0]}{state1[1]}{state2[0]}{state2[1]}'

        fit_exists = False
        
        if 'Stark_Shift_Fits' not in self.Extra_Stuff: 
            self.Extra_Stuff['Stark_Shift_Fits'] = {}

        if fit_name in self.Extra_Stuff['Stark_Shift_Fits']:
            fit_exists = True
        if use_fit and fit_exists:
            fitpoly = np.poly1d(self.Extra_Stuff['Stark_Shift_Fits'][fit_name])
            return fitpoly(epsilon)

        if use_fit and not fit_exists:
            print('No fit exits for this pulse, please calibrate or generate a fit')
            return None

        if freq_d == None:
            freq_d = self.DefaultFrequency(state1, state2)
        if not make_fit:
            return FindResonanceFloquet(self, state1, state2, freq_d, epsilon = epsilon, **kwargs)
        
        if make_fit:
            if 'return_mins' not in kwargs:
                kwargs['return_mins'] = True
            if not kwargs['return_mins']:
                kwargs['return_mins'] = True
            
            fit_dat = StarkShiftFitter(self, state1, state2, freq_d, **kwargs)

            if save_new_fit:
                self.Extra_Stuff['Stark_Shift_Fits'][fit_name] = fit_dat['fit']
            
            polyfit = np.poly1d(fit_dat['fit'])

            return dict(StarkShift = polyfit(epsilon), fit_dat = fit_dat)
        
        print('Did Nothing')
        return None

    def CalibratePulse(self, pulse_type, state1, state2, args, t0, tf, steps, op_name, kwargs = {}):
        global_calibrator_dict[pulse_type](self, state1, state2, args, t0, tf, steps, op_name, **kwargs)

    def Save_Params(self, filename = None):
        if filename == None:
            if type(self.model_name) == str:
                filename = self.model_name

        if self.save_path == None:
            print('Need model save path')
            return None

        TO_SAVE = {'Main_Config':self.model_params}
        
        TO_SAVE['op_drive_params_dict'] = self.op_drive_params_dict

        TO_SAVE['Drive Sequences'] = self.Drive_Sequences
        
        TO_SAVE['Extra_Stuff'] = self.Extra_Stuff

        

        try:
            with open(self.save_path+filename+'.json', 'w') as file:
                file.write(json.dumps(TO_SAVE))
                
        except FileNotFoundError:
            os.makedirs(self.save_path)

            with open(self.save_path+filename+'.json', 'w') as file:
                file.write(json.dumps(TO_SAVE))

def Make_Pulse_Diagram(model, pulse_sequence, op_colors = 'Default', op_style = 'Default', Transmon_Colors = 'Default', spns = 5,  q_scale_factor = 10, sb_scale_factor = 0.1, base_offset = 1, linewidth = 1, ax = None, figsize = (5,3), dpi = 200, title = None, sb_tick_num = 5, q_tick_num = 5, q_ytick_factor = 1000, sb_ylabel_loc = (0.01, 0.33), q_ylabel_loc = (0.01, 0.7), q_ylabel = r'$\epsilon$ (MHz)', sb_ylabel = r'$\epsilon$ (GHz)', ylim_factor = 1.15, wave_lw_factor = 0.1, wave_alpha = 0.5, num_label_skip = 1, q_yticks = 'Default', sb_yticks = 'Default', wave_factor = 1):
    if ax == None:
        fig, ax = plt.subplots(figsize = figsize, dpi = dpi)
    
    time_stamps = []
    time = 0
    xtick_labels = []
    xtick_locs = []

    list_of_ys = []
    list_of_times = []
    max_y = 0
    for i in range(len(pulse_sequence)):
        pulse = pulse_sequence[i]
        pulse_dat = copy.deepcopy(model.op_drive_params_dict[pulse])
        pulse_length = pulse_dat['pulse_time']
        temp_times = np.linspace(0, pulse_length, int(spns*pulse_length))
        
        y = []
        for t in temp_times:
            y.append(global_envelope_dict[pulse_dat['Envelope']](t, **pulse_dat['Envelope Args']))  
        
        y = np.array(y)*pulse_dat['epsilon']

        if y.max()>max_y:
            max_y = y.max()
        
        list_of_ys.append(y)
        list_of_times.append(temp_times)
    
    max_q = 0
    for i in range(len(pulse_sequence)):
        pulse = pulse_sequence[i]
        pulse_dat = copy.deepcopy(model.op_drive_params_dict[pulse])
        pulse_length = pulse_dat['pulse_time']

        pulse_split = pulse.split('_')
        
        op_style = 'solid'
        if pulse_split[0].split('-')[0] == 'half':
            op_style = ':'
        
        if pulse_split[0].split('-')[-1] == 'q':
            if pulse_split[1] == 'g':
                color = 'blue'
            if pulse_split[1] == 'e':
                color = 'red'
        
        if pulse_split[0].split('-')[-1] == 'sb':
            if pulse_split[-1] == '3':
                color = 'green'
            if pulse_split[-1] == '5':
                color = 'orange'
        


        y = list_of_ys[i]
        temp_times = list_of_times[i]        

        offset = 0
        if pulse.split('_')[0] == 'q' or pulse.split('_')[0] == 'half-q':
            offset = (max_y+base_offset)*sb_scale_factor
            y*= q_scale_factor

            if y.max()>max_q:
                max_q = y.max()
        
        if pulse.split('_')[0] == 'sb' or pulse.split('_')[0] == 'half-sb':
            y*= sb_scale_factor

        
        ax.plot(temp_times+time, y+offset, color = color, ls = op_style, lw = linewidth)
        ax.plot(temp_times+time, -y+offset, color = color, ls = op_style, lw = linewidth)

        freq = pulse_dat['freq_d']+pulse_dat['shift']

        ax.plot(temp_times+time, np.sin(2*np.pi*freq/wave_factor*temp_times)*y+offset, color = color, lw = linewidth*wave_lw_factor, alpha = wave_alpha)


        xtick_labels.append(i)
        xtick_locs.append((temp_times+time).mean())

        time += temp_times[-1]

    if title != None:
        ax.set_title(title)


    ax.axhline(0, ls = ':', lw = 0.5, color = 'black')
    ax.axhline((max_y+base_offset)*sb_scale_factor, ls = ':', lw = 0.5, color = 'black')    

    if q_yticks is 'Default':
        q_ytick_locs = np.linspace((max_y+base_offset)*sb_scale_factor-max_q, (max_y+base_offset)*sb_scale_factor+max_q, q_tick_num)
        q_ytick_labels = np.round(np.linspace(-max_q/q_scale_factor, max_q/q_scale_factor, q_tick_num)*q_ytick_factor,3)

    if sb_yticks is 'Default':
        sb_ytick_locs = np.linspace(-max_y*sb_scale_factor, max_y*sb_scale_factor, sb_tick_num)
        sb_ytick_labels = np.round(np.linspace(-max_y, max_y, sb_tick_num), 3)
    
    if q_yticks is not 'Default':
        q_ytick_locs = np.array(q_yticks)*q_scale_factor + (max_y+base_offset)*sb_scale_factor
        q_ytick_labels = np.array(q_yticks)*q_ytick_factor
        print(q_scale_factor)

    if sb_yticks is not 'Default':
        sb_ytick_locs = np.array(sb_yticks)*sb_scale_factor
        sb_ytick_labels = np.array(sb_yticks)
    

    ytick_locs = np.concatenate([sb_ytick_locs, q_ytick_locs])
    ytick_labels = np.concatenate([sb_ytick_labels, q_ytick_labels])
    ax.set_yticks(ticks = ytick_locs, labels = ytick_labels)
    #ax.set_xticks(ticks=xtick_locs, labels=xtick_labels)
    ax.set_xlabel('Time (ns)')

    ax.annotate(sb_ylabel, sb_ylabel_loc, xycoords = 'subfigure fraction', rotation = 90)

    q_ylabel_x = 0
    q_ylabel_y = 0.4
    ax.annotate(q_ylabel, q_ylabel_loc, xycoords = 'subfigure fraction', rotation = 90)    
    label_eps_scale_function = r'$\epsilon$'
    
    ax.set_ylim(ytick_locs[0]*ylim_factor, ytick_locs[-1]*ylim_factor)

    ax_for_top = ax.twiny()
    ax_for_top.plot(np.linspace(0, time),0*np.linspace(0, time), alpha = 0)

    ax_for_top.tick_params(direction = 'in', pad = -15)
    ax_for_top.set_xticks(ticks=xtick_locs[::num_label_skip], labels=xtick_labels[::num_label_skip])
    ax_for_top.set_xlabel('Pulse #')



class Transmon_Cavity_Abstract:
    def __init__(self, N_Transmon = 4,N_Cavity = [5]):
        self.N_Transmon = N_Transmon
        self.N_Cavity = N_Cavity
        if type(self.N_Cavity) == int:
            self.N_Cavity = [self.N_Cavity]


        self.Transmon_Labels_LK = {'g':0, 'e':1, 'f':2, 'h':3, 'i':4, 'j':5}

        self.Transmon_Labels_NK = {}
        for  i,key in enumerate(self.Transmon_Labels_LK):
            self.Transmon_Labels_NK[i]=key

        self.op_dict = {}

        self.Build_Creation_and_Annihilation_Ops()

        
    def Build_Creation_and_Annihilation_Ops(self):
        temp_list = [qp.destroy(self.N_Transmon)]
        for i in range(len(self.N_Cavity)):
            temp_list.append(qp.identity(self.N_Cavity[i]))
        a_T = qp.tensor(temp_list)

        a_C_list = []
        for i in range(len(self.N_Cavity)):
            temp_list = [qp.identity(self.Transmon)]
            for j in range(len(self.N_Cavity)):
                if i != j:
                    temp_list.append(qp.identity(self.N_Cavity)[j])
                if i == j:
                    temp_list.append(qp.destroy(self.N_Cavity[i]))
            a_C_list.append(qp.tensor(temp_list))


        n_T = a_T.dag()*a_T

        n_C_list = []
        for i in a_C_list:
            n_C_list.append(a_C_list[i].dag()*a_C_list[i])

        self.op_dict['a_T']=a_T
        self.op_dict['a_C']=a_C
        self.op_dict['n_T']=n_T
        self.op_dict['n_C']=n_C
        
        def shift_mat(N):
            mat = np.zeros((N,N))
            for i in  range(N-1):
                mat[i,i+1] = 1
            return mat

        
        temp_list = [qp.Qobj(shift_mat(self.N_Transmon))]
        for i in range(len(self.N_Cavity)):
            temp_list.append(qp.identity(self.N_Cavity[i]))
        a_T_Norm  = qp.tensor(temp_list)

        a_C_Norm_list = []
        for i in range(len(self.N_Cavity)):
            temp_list = [qp.identity(self.Transmon)]
            for j in range(len(self.N_Cavity)):
                if i != j:
                    temp_list.append(qp.identity(self.N_Cavity)[j])
                if i == j:
                    temp_list.append(qp.Qobj(shift_mat(self.N_Cavity[i])))
            a_C_Norm_list.append(qp.tensor(temp_list))
        

        self.op_dict['a_T_Norm']=a_T_Norm
        self.op_dict['a_C_Norm']=a_C_Norm_list

    def create_state(self,i,*cavities):
        if type(i) == str:
            i =  self.Transmon_Labels_LK[i]
        
        temp_list = [qp.fock(self.Transmon, i)]
        for j in range(len(cavities)):
            temp_list.append(qp.fock(self.N_Cavity[j], cavities[j]))
        return qp.tensor(temp_list)
    
    def q_maker(self, l):
        '''l: level which q acts on'''
        if type(l) == str:
            l =  self.Transmon_Labels_LK[l]
        fock_T_2 = qp.fock(self.N_Transmon, l)
        fock_T_3 = qp.fock(self.N_Transmon, l+1)
        
        cavity_ones = []
        for i in range(len(self.N_Cavity)):
            cavity_ones.append(qp.identity(self.N_Cavity[i]))

        proj2 = qp.tensor([fock_T_2*fock_T_2.dag()]+cavity_ones)
        proj3 = qp.tensor([fock_T_3*fock_T_3.dag()]+cavity_ones)

        ONE = qp.tensor([qp.identity(self.N_Transmon)]+cavity_ones)

        op1 = proj3*self.op_dict['a_T_Norm'].dag()*proj2-proj2*self.op_dict['a_T_Norm']*proj3
        op2 = (ONE-proj3)*(ONE-proj2)

        return op1+op2
    
    def half_q_maker(self, l, cavity):
        '''l: level which q acts on
            cavity: level cavity acts on'''
        if type(l) == str:
            l =  self.Transmon_Labels_LK[l]
        if type(cavity)==str:
            cavity = int(cavity)
        
        cavity_ones = []
        for i in range(len(self.N_Cavity)):
            cavity_ones.append(qp.identity(self.N_Cavity[i]))

        cavity = qp.fock(self.N_Cavity, cavity)

        trans_i = qp.fock(self.N_Transmon, l)
        trans_ip1 = qp.fock(self.N_Transmon, l+1)

        ONE = qp.tensor([qp.identity(self.N_Transmon)+cavity_ones])

        trans_op = 1/np.sqrt(2)*(trans_i+trans_ip1)*trans_i.dag() 
        trans_op += 1/np.sqrt(2)*(-trans_i+trans_ip1)*trans_ip1.dag()
        trans_op += -trans_i.proj()-trans_ip1.proj()

        return ONE + qp.tensor([trans_op]+cavity_ones)
    
    def Half_Sideband_Maker(self, l, cavity, mode):
        if type(l)==str:
            l = self.Transmon_Labels_LK[l]
        
        if type(cavity)==str:
            cavity = int(cavity)

        fock1_list = []
        fock2_list = []
        for i in range(len(self.N_Cavity)):
            if i!= mode:
                fockl1_list.append(qp.identity(self.N_Cavity[i]))
                fockl2_list.append(qp.identity(self.N_Cavity[i]))
            if i == mode:
                fockl1_list.append(qp.fock(self.N_Cavity[i], cavity))
                fockl2_list.append(qp.fock(self.N_Cavity[i], cavity+1))
        
        fock1 = qp.tensor([qp.fock(self.N_Transmon, l)]+fock1_list)
        fock2 = qp.tensor([qp.fock(self.N_Transmon, l-2)]+fock2_list)

        proj1 = fock1*fock1.dag()
        proj2 = fock2*fock2.dag()
        


        Id = qp.tensor(qp.identity(self.N_Transmon), qp.identity(self.N_Cavity))
        
        op1 = (proj1+proj2+proj2*self.op_dict['a_C_Norm'][mode].dag()*self.op_dict['a_T_Norm']*self.op_dict['a_T_Norm']*proj1-proj1*self.op_dict['a_C_Norm'][mode]*self.op_dict['a_T_Norm'].dag()*self.op_dict['a_T_Norm'].dag()*proj2)/np.sqrt(2)
        op2 = (Id-proj1)*(Id-proj2)

        return op1+op2
    
    def Make_Phase_Shift(self, i,j):
        if type(i)==str:
            i  =  self.Transmon_Labels_LK[i]
            
        state = qp.tensor(qp.fock(self.N_Transmon, i), qp.fock(self.N_Cavity, j))
        proj = state*state.dag()

        op1 = -proj
        op2 = -proj+qp.tensor(qp.identity(self.N_Transmon), qp.identity(self.N_Cavity))

        return op1+op2   

    def op_maker(self, op_name):
        op_name_list = op_name.split('_')
        print(f"op_name: {op_name}, op_name_list: {op_name_list}")
        if op_name_list[0] == 'q':
            self.op_dict[op_name] = self.q_maker(op_name_list[1])
        
        if op_name_list[0] == 'sb':
            half_op = self.Half_Sideband_Maker(op_name_list[1], op_name_list[2])
            self.op_dict[op_name] = half_op*half_op
        
    
        if op_name_list[0] == 'half-q':
            self.op_dict[op_name] = self.half_q_maker(op_name_list[1], op_name_list[2])
        
        if op_name_list[0] == 'half-sb':
            half_op = self.Half_Sideband_Maker(op_name_list[1], op_name_list[2])
            self.op_dict[op_name] = half_op

    def Make_State_Plot(self, init_state, op_list, op_colors='Default', op_style='Default', Transmon_Colors='Default', mode_sep = 0.1, arrow_radius = 0.1, arrow_rad_growth_exp = 1, dpi = 200, figsize = (7,5), width = 0.5, chi = 0.3, sep = 0.5, save = None, op_num_labels_size = 8, op_label_offset = 0, op_label_bbox_kwargs = dict(boxstyle='Circle,pad=0', facecolor ='white', edgecolor='white', alpha = 0.8), init_color = 'xkcd:minty green', final_color = 'xkcd:grapefruit', init_markers_box_style = dict(boxstyle='Circle,pad=0.25', alpha = 1), final_markers_box_style = dict(boxstyle='Square,pad=0.15', alpha = 1), label_fs = 8, init_final_fs = 8, init_label = r'$\psi_i$', final_label = r'$\psi_f$', threshhold  = 1e-5):
        state_history_list = []
        
        shape = copy.copy(self.N_Cavity)
        shape.insert(0, self.N_Transmon)
        to_iterate = [np.arange(0, i) for i in shape]

        zero_array = np.empty(shape, dtype=object)

        state_array = np.empty_like(zero_array)
        init_array = zero_array.copy()        

        for coord in itert.product(*to_iterate):
            state_array[coord] = self.create_state(*coord)
            zero_array[coord] = 0*state_array[coord]

            proj = state_array[coord]*state_array[coord].dag()
            init_array[coord] = proj*init_state
        
        state_history_list.append(init_array)

        state_proj_history = []
        for n_op in range(len(op_list)):
            op_name = op_list[n_op]

            if op_name not in self.op_dict:
                self.op_maker(op_name)                        
                    
            op = self.op_dict[op_name]

            step_array =  np.empty(shape,  dtype=object)
            
            new_state_array = zero_array.copy()

            proj_list = []
            proj_list_key = []
            
            for coord in itert.product(*to_iterate):
                state = state_history_list[n_op][coord]
        
                if state.norm()>threshhold:
                    new_state =  op*state
                    proj_array =  np.zeros(shape, dtype=complex)
                    
                    for coord in itert.product(*to_iterate):
                        if qp.version.version[0]!='5':
                            proj_array[coord] = (state_array[coord].dag()*new_state)[0][0]
                        if qp.version.version[0]=='5':
                            proj_array[coord] = (state_array[coord].dag()*new_state)
            
                    proj_list.append(proj_array.copy())
                    proj_list_key.append(coord)

            proj_array = np.array(proj_list)
            for coord in itert.product(*to_iterate):
                new_state_array[coord] = state_array[coord]*np.sum(proj_array[:, coord])
            
            state_history_list.append(new_state_array.copy())
            state_proj_history.append((proj_array.copy(), np.array(proj_list_key)))

        fig, ax = plt.subplots(dpi = dpi, layout = 'tight', figsize = figsize)

        connections_counter_dict = {}

        delta = lambda x,y: 1 if x==y else 0

        for mode in range(len(self.N_Cavity)):
            for coord in itert.product(np.arange(self.N_Transomn), np.arange(self.N_Cavity[mode])):
                if coord[1] != 0:
                    offset = mode_sep*mode
                if coord[1] == 0:
                    offset = 0
                ax.add_line(plt.Line2D([coord[0], coord[0]+width], [coord[1]+chi*coord[0]+offset,coord[1]+chi*coord[0]+offset], color = 'black'))
        
        for n in range(len(state_proj_history)):
            for m in range(len(state_proj_history[n][1])):
                transmon_start = state_proj_history[n][1][m][0]
                cavity_start = state_proj_history[n][1][m][1]
                
                new_state = []
                for i in range(state_proj_history[n][0][m].shape[0]):
                    for j in range(state_proj_history[n][0][m].shape[1]):
                        if np.abs(state_proj_history[n][0][m][i,j])>1e-10:
                            new_state.append([i,j])
                
                for i in range(len(new_state)):
                    transmon_end = new_state[i][0]
                    cavity_end = new_state[i][1]

                    shift = 0
                    
                    connection = str([(transmon_start,  cavity_start), (transmon_end, cavity_end)])
                    
                    if connection not in connections_counter_dict:
                        connections_counter_dict[connection] = 0

                    connections_counter_dict[connection]+=1

                    n_con =  connections_counter_dict[connection]            

                    if transmon_end<transmon_start:
                        shift = -width

                    x0 = width*(transmon_start+1)+sep*(transmon_start)+shift
                    y0 = cavity_start+chi*transmon_start

                    x1 = width*(transmon_end)+sep*(transmon_end)-shift
                    y1 = cavity_end+chi*transmon_end

                    dx = x1-x0
                    dy = y1-y0

                    str_split = op_list[n].split('_')
                    
                    ls = op_style[str_split[0]]
                    color =op_colors[str_split[0]][str_split[1]]
                    if (transmon_start==transmon_end) & (cavity_end==cavity_start):

                        x  = width*transmon_start+sep*transmon_start+width/2
                        y = cavity_start+chi*transmon_end
                        w  = 0.05*(n_con+1)
                        l  = 0.15*(n_con+1)

                        x0=x
                        y0 = y+l/2
                        arrow = mpatches.Ellipse((x0,y0), w, l, fill =  False, color = color, ls=ls)

                        dist = np.sqrt(dx**2+dy**2)
                        shift = 0
                        label_pos_x = x0
                        label_pos_y = y0+l/2+0.045
                        
                    else:
                        rad = -arrow_radius*(n_con**arrow_rad_growth_exp)*np.sign(dx)*np.sign(dy)
                        arrow = mpatches.FancyArrowPatch((x0,y0), (x1,y1), arrowstyle='->',ls=ls, connectionstyle='arc3, rad='+str(rad), mutation_scale=10, zorder = 100, color = color)
                        
                        label_pos_x = (x1+x0)/2
                        label_pos_y = (y1+y0)/2+(np.abs(rad))*np.sign(dy)
                        
                    
                    old_artists = ax.get_children()
                    ax.add_patch(arrow)
                    new_artists = ax.get_children()

                    new_artist = [x for x in new_artists if x not in old_artists][0]
                    #print(new_artist.)
                    
                    fontsize = op_num_labels_size
                    ax.annotate(n, [label_pos_x,label_pos_y-fontsize/175+op_label_offset], fontsize =  fontsize, ha ='center', color='black', zorder =  200, bbox=op_label_bbox_kwargs)
                    #ax.annotate(n, [label_pos_x,label_pos_y+op_label_offset*np.sign(dy)], fontsize =  fontsize, ha ='center', color='black', zorder =  200, bbox=op_label_bbox_kwargs)
                    #ax.annotate(n, [0.5,0.5*(1+np.sign(dy))], xycoords = arrow, fontsize =  fontsize, ha ='center', color='black', zorder =  200, bbox=op_label_bbox_kwargs)

                    fs = init_final_fs
                    shiftx =  -fs/1000
                    shifty = -fs/175
                    if n==0:
                        color = 'black'
                        x = width*transmon_start+sep*transmon_start+width/2+shiftx
                        y = cavity_start+transmon_start*chi+shifty
                        bbox_dat = copy.deepcopy(init_markers_box_style)
                        bbox_dat['facecolor']= init_color
                        bbox_dat['edgecolor']= 'black'
                        ax.annotate(init_label, (x,y), color=color, fontsize = fs, ha ='center', bbox = bbox_dat)

                    if n==len(state_proj_history)-1:
                        color = 'black'
                        x = width*transmon_end+sep*transmon_end+width/2+shiftx
                        y = cavity_end+transmon_end*chi+shifty
                        bbox_dat = copy.deepcopy(final_markers_box_style)
                        bbox_dat['facecolor']= final_color
                        bbox_dat['edgecolor']= 'black'
                        ax.annotate(final_label, (x,y), color=color, fontsize = fs, ha='center', bbox = bbox_dat)

            for i in range(self.N_Transmon):
                label=r'$|N, i\rangle$'.replace('N', self.Transmon_Labels_NK[i])
                ax.annotate(label, ((width+sep)*i+width/2, i*chi-0.4), ha =  'center', fontsize = label_fs) 

            for i in range(self.N_Cavity):
                label=r'$|i, N\rangle$'.replace('N', str(i))
                ax.annotate(label, (-0.3, i), ha =  'center', fontsize = label_fs) 


            ax.set_axis_off()
            ax.set_ylim(-0.5,   self.N_Cavity)
            ax.set_xlim(-0.5, self.N_Transmon)  

        if type(save) == str:
            plt.savefig(save)
