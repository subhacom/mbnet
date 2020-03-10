# pn_output.py --- 
# 
# Filename: pn_output.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Wed Jan 10 17:01:37 2018 (-0500)

# Code:

"""Generate array of spike times for 830 PNs. These will be fed to KC
population.

"""
from __future__ import print_function
import numpy as np
import nhpp
from config import Q_, ur, h, logger
import nrnutils as nu
import ephys
from timeit import default_timer as timer
from collections import defaultdict

osc_fn = lambda a, f, t: a * np.sin(2 * np.pi * f * t)   # simple sinusoidal oscillation function

def create_pn_output(params):
    
    """This attempts to create inputs a KC receives from its presynaptic PNs
    References: Mazor and Laurent, 2005; Jortner, et al., 2007.

    Here the assumption is `odor_exc_frac` fraction of all PNs is
    activated by a specific odor, i.e. there firing rate goes up
    (exact formula in code).

    On the other hand, spont_exc_frac of the cells are spontaneously
    active.
    
    npn: total number of PNs

    odor_exc_frac: fraction of PNs that spike in response to odor

    spont_exc_frac: fraction of PNs that spike spontaneously

    stim_rate: baseline firing rate of a PN during odor stimulation

    spont_rate: baseline firing rate of a PN without odor 

    osc_freq: frequency of LFP oscillation


    osc_amp_scale: The amplitude of oscillation as a fraction of
    baseline firing rate (default 0.3 from data from 50000 PNs).
    

    Returns a list of array quantitities containing spike times of
    each PN.

    """
    onset = Q_(params['onset']).to('s').m
    duration = Q_(params['duration']).to('s').m
    tail = Q_(params['tail']).to('s').m
    delta = Q_(params['delta']).to('s').m
    npn = params['npn']
    odor_exc_frac = params['odor_exc_frac']
    spont_exc_frac = params['spont_exc_frac']
    odor_inh_frac = params['odor_inh_frac']
    stim_rate = Q_(params['stim_rate']).to('Hz').m
    spont_rate = Q_(params['spont_rate']).to('Hz').m
    osc_freq = Q_(params['osc_freq']).to('Hz').m
    osc_amp_scale = params['osc_amp_scale']
    
    assert (odor_exc_frac + odor_inh_frac) <= 1

    spont_amp = spont_rate * osc_amp_scale
    stim_amp = stim_rate * osc_amp_scale
    # Rate function for spontaneously active, excited by odor
    rate_fn_spont_odor = lambda t: \
                         stim_rate + stim_amp * np.sin(2 * np.pi * osc_freq * t) \
                         if (onset < t < (onset + duration)) \
                            else spont_rate + spont_amp * np.sin(2 * np.pi * osc_freq * t)
    # Odor inhibited, but spontaneously active
    rate_fn_odor_inh = lambda t: \
                       0.0 if (onset < t < (onset + duration)) \
                       else spont_rate + spont_amp * np.sin(2 * np.pi * osc_freq * t)
    # Rate function for spontaneously active, unaffected by odor
    rate_fn_spont = lambda t: \
                    spont_rate + spont_amp * np.sin(2 * np.pi * osc_freq * t)
    # rate function for spontaneously silent, activated by odor
    rate_fn_odor = lambda t: \
                   stim_rate + stim_amp * np.sin(2 * np.pi * osc_freq * t) \
                   if (onset < t < (onset + duration)) \
                      else 0.0
    n_spont_odor = int(npn * spont_exc_frac * odor_exc_frac + 0.5)   # spontanously active, excited by odor
    n_odor = int(npn * (1 - spont_exc_frac) * odor_exc_frac + 0.5)   # silent, activated by odor
    n_odor_inh = int(npn * spont_exc_frac * odor_inh_frac + 0.5)     # spontaneously active inhibited by odor
    n_spont = int(npn * spont_exc_frac * (1 - odor_exc_frac - odor_inh_frac) + 0.5)    
    # Ignore those that are neither spontaneously active nor activated
    # by a given odor: how many such PNs? 
    pn_output = []
    start = timer()
    for ii in range(n_spont_odor):
        spike_times = nhpp.nhpp_thinning(rate_fn_spont_odor,
                                         onset + duration + tail,
                                         delta)
        pn_output.append(Q_(np.array(spike_times), 's'))
    for ii in range(n_odor):
        spike_times = nhpp.nhpp_thinning(rate_fn_odor,
                                         onset + duration + tail,
                                         delta)
        pn_output.append(Q_(np.array(spike_times), 's'))
    for ii in range(n_spont):
        spike_times = nhpp.nhpp_thinning(rate_fn_spont,
                                         onset + duration + tail,
                                         delta)
        pn_output.append(Q_(np.array(spike_times), 's'))
    for ii in range(n_odor_inh):
        spike_times = nhpp.nhpp_thinning(rate_fn_odor_inh,
                                         onset + duration + tail,
                                         delta)
        pn_output.append(Q_(np.array(spike_times), 's'))
    end = timer()    
    logger.info('Time to create {} PN spiketrains {} s'.format(len(pn_output), end - start))
    silent = npn - len(pn_output)
    for ii in range(silent):
        pn_output.append(Q_(np.array([]), 's'))
    return pn_output


def create_shifting_response(params):
    """Create a response pattern with a shifting window into PN population
    over the LFP cycles.

    The step by step process to arrive at this design is in
    try_pn_output.ipynb.
    
    params should have the following key value pairs:

    onset: stimulus / activity onset 
    duration: duration of stimulus
    offdur: duration of off response
    tail: simulation time after the stimulus is over
    trials(=1): number of trials (only `tail` interval between successive stimuli.
    delta: interval at which the Poisson rate is considered.
    npn(=830): total number of PNs    
    odor_exc_shift(=0.1): fraction of excited PN changing state in each LFP cycle
    stim_rate(=20.0Hz): firing rate of excited PN upon odor, i.e. 1 spike in 50 ms    
    spont_rate(=2.6Hz): spontaneous firing rate of a PN
    osc_freq(=20.0Hz): LFP oscillation frequency    
    start_frac: fraction of cells in each group that start spiking simultaneously at the excitation phase of the group

    Returns: 

    pn_spikes - list of lists: i-th entry is a list spike
    times of i-th PN.

    clusters - dict: list of PNs mapped to their stimulus response
    start time. The entries against off response time include bith EI
    and II PNs.

    response_types  - dict: {
      'EE': indices of PNs that spike throughout odor presentation,
      'EI': indices of PNs that spike during first half of odor presentation,
      'IE': indices of PNs that spike during last half of odor presentation,
      'II': indices of PNs that do not spike during  of odor presentation
      'U': spntaneously active and unresponsive to odor stimulus
    }

    """
    npn = params['npn']
    lfp_bin = 50e-3 # second
    spont_rate = params['spont_rate'].to('Hz').m
    delta = params['delta'].to('s').m
    onset = params['onset'].to('s').m
    dur = params['duration'].to('s').m
    tail = params['tail'].to('s').m
    offdur = params['offdur'].to('s').m
    osc_freq = params['osc_freq'].to('Hz').m
    stim_rate = params['stim_rate'].to('Hz').m 
    osc_amp_scale = params['osc_amp_scale']
    start_frac = params['start_frac']
    # Indices of different kinds of PNs.
    unresponsive = range(0, int(npn/3.0))
    EE = list(range(unresponsive[-1]+1,int(unresponsive[-1]+1+npn/6.0)))
    EI = list(range(EE[-1]+1, int(EE[-1]+1+npn/6.0)))
    IE = list(range(EI[-1]+1, int(EI[-1]+1+npn/6.0)))
    II = list(range(IE[-1]+1, npn))
    response_types = {'EE': EE, 'EI': EI, 'IE': IE, 'II': II, 'U': unresponsive}
    clusters = defaultdict(list)
    pn_spikes = []
    start = timer()
    # Spontaneous activity in all neurons
    spont_rate_fn = lambda t: spont_rate
    for ii in range(params['npn']):
        st1 = nhpp.nhpp_thinning(spont_rate_fn, onset, delta)
        st2 = nhpp.nhpp_thinning(spont_rate_fn, tail, delta) + onset + dur + offdur
        pn_spikes.append(list(st1) + list(st2))

    # Unresponsive PNs - first 1/3rd
    for jj, index in enumerate(unresponsive):
        pn_spikes[index] += list(nhpp.nhpp_thinning(spont_rate_fn, dur + offdur, delta) + onset)
    
    # Excited throughout odor
    
    # The rate function: python has dynamic binding - so ee_rate_fn
    # sees the value of `start` at evaluation time.
    ee_rate_fn = lambda t: stim_rate + osc_fn((dur - t) * stim_rate * osc_amp_scale / dur, osc_freq, t)  \
                 if  t > start else 0
    start = 0
    # After every `pnclus` PNs, increase the start time. 
    pnclus = round(params['shifting_frac'] * len(EE))    
    for jj, index in enumerate(EE):
        pn_spikes[index] += list(nhpp.nhpp_thinning(ee_rate_fn, dur, delta) + onset)
        clusters[onset + start].append(index)
        diff = jj - len(EE) * start_frac
        if (diff > 0) and ( round(diff) % pnclus == 0):
            start += lfp_bin
    # PNs active only in first half of odour presentation
    ei_rate_fn = lambda t: stim_rate + osc_fn((dur - t) * stim_rate * osc_amp_scale / dur, osc_freq, t)  \
                 if  (t > start) and (t < dur/2.0) else 0
    start = 0    
    pnclus = round(params['shifting_frac'] * len(EI))    # How many PNs should share the same starting time?
    for jj, index in enumerate(EI):
        pn_spikes[index] += list(nhpp.nhpp_thinning(ei_rate_fn, dur, delta) + onset)
        clusters[onset + start].append(index)
        diff = jj - len(EI) * start_frac
        if (diff > 0) and (round(diff) % pnclus == 0):
            start += lfp_bin

    # PNs active only in last half of odour presentation
    ie_rate_fn = lambda t: stim_rate + osc_fn((dur - t) * 2.0 * stim_rate * osc_amp_scale / dur, osc_freq, t) \
                 if  (t > start) else 0
    start = dur / 2.0
    pnclus = round(params['shifting_frac'] * len(IE))    # How many PNs should share the same starting time?
    for jj, index in enumerate(IE):
        pn_spikes[index] += list(nhpp.nhpp_thinning(ie_rate_fn, dur, delta) + onset)
        clusters[onset + start].append(index)
        diff = jj - len(IE) * start_frac
        if (diff > 0) and (round(diff) % pnclus == 0):
            start += lfp_bin
    # off response
    off_rate_fn =  lambda t: stim_rate + osc_fn((offdur - t) * stim_rate * osc_amp_scale / offdur, osc_freq, t)  \
                   if  (t > start) else 0
    start = 0
    pnclus = round(params['shifting_frac'] * len(EI))    # How many PNs should share the same starting time?
    for jj, index in enumerate(EI):
        pn_spikes[index] += list(nhpp.nhpp_thinning(off_rate_fn, offdur, delta) + onset + dur)
        clusters[onset + dur + start].append(index)
        diff = jj - len(EI) * start_frac
        if (diff > 0) and (round(diff) % pnclus == 0):
            start += lfp_bin
    start = 0
    pnclus = round(params['shifting_frac'] * len(II))    # How many PNs should share the same starting time?
    for jj, index in enumerate(II):
        pn_spikes[index] += list(nhpp.nhpp_thinning(off_rate_fn, offdur, delta) + onset + dur)
        clusters[onset + dur + start].append(index)
        diff = jj - len(II) * start_frac
        if (diff > 0) and (round(diff) % pnclus == 0):
            start += lfp_bin
    qst = []
    for ii, st in enumerate(pn_spikes):
        st = np.array(st)
        st.sort()
        qst.append(Q_(st, 's'))
    pn_spikes = qst
    end = timer()
    logger.info('Finished PN output creation in {} s'.format(end - start))
    return pn_spikes, clusters, response_types
                                        
               
# 
# pn_output.py ends here
