# change_pn_spikes.py ---
# Author: Subhasis Ray
# Created: Thu Feb 28 19:29:46 2019 (-0500)
# Last-Updated: Tue Mar 10 15:25:47 2020 (-0400)
#           By: Subhasis Ray
# Version: $Id$

# Code:
import sys
sys.path.append('/home/rays3/projects/ggn/common')
import argparse
import shutil
import os
import numpy as np
import h5py as h5
import yaml
import nhpp
# from math import fmod
from collections import defaultdict
from pint import UnitRegistry
from matplotlib import pyplot as plt

ur = UnitRegistry()
Q_ = ur.Quantity


osc_fn = lambda a, f, t: a * np.sin(2 * np.pi * f * t)   # simple sinusoidal oscillation function



def make_parser():
    parser = argparse.ArgumentParser(description="""Update existing model file with new PN spike trains and changes in configuration attributes.""")
    parser.add_argument('-i', dest='infile', type=str, help='Input file', required=True)
    parser.add_argument('-o', dest='outfile', type=str, help='output file', required=True)
    parser.add_argument('--npn', dest='npn', type=int, default=-1, help='number of PNs')
    parser.add_argument('--ig_ggn_gmax', type=str, help='ig-ggn gmax')
    parser.add_argument('--kc_ig_weight', type=float, help='kc-ig synaptic weight')
    parser.add_argument('--pn_duration', type=str, help='pn activity duration')
    parser.add_argument('--duration', type=str, help='stimulus duration')
    parser.add_argument('--onset', type=str, help='stimulus onset time')
    parser.add_argument('--tail', type=str, help='tail of simulation (after PN response stops)')
    parser.add_argument('--interactive', action='store_true')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--responsive', type=float, default=-1.0, help='Fraction of cells responding to stimulus')
    parser.add_argument('--start_frac',  type=float, default=-1.0, help='Fraction PNs starting together')
    parser.add_argument('--shifting_frac',  type=float, default=-1.0, help='Fraction PNs recruited in each time bin')
    parser.add_argument('--ei',  type=float, default=-1.0, help='Fraction of responsive cells with EI response')
    parser.add_argument('--ee',  type=float, default=-1.0, help='Fraction of responsive cells with EE response')
    parser.add_argument('--ie',  type=float, default=-1.0, help='Fraction of responsive cells with IE response')
    return parser

pn_params = """
    onset: 0.5s                   # when the stimulus starts (PN spiking rate increases)
    duration: 1.0s                # duration of stimulus
    tail: 1.0s                    # duration of simulation post stimulus start time
    stabilization_time: 0.2s    # not even spontaneous activity - allow each cell to reach its resting Vm
    npn: 830
    delta: 5ms                  # parameter for inhomogeneous Poisson process to generate the PN spike times
    number: 830                 # number of PNs
    odor_exc_frac: 0.2          # fraction of all PNs excited by odor
    odor_inh_frac: 0.1          # fraction of spontaneously active PNs inhibited by odor - guessed
    spont_exc_frac: 0.77        # fraction of PNs spontaneously active
    stim_rate: 25.0Hz           # baseline firing rate of a PN upon stimulus
    spont_rate: 2.6Hz           # baseline firing rate of a PN during spontaneous activity
    osc_freq: 20.0Hz            # LFP oscillaion frequency
    osc_amp_scale: 0.4          # the oscillation amplitude as a fraction of baseline firing rate
    offdur: 0.5s                  # duration of off response
    shifting: true             # whether to make PN population with shifting activity time
    start_frac: 0.9             # what fraction of PNs start together (shifting activity)
    shifting_frac: 0.05          # what fraction will be newly recruited at each time bin
    pndur: 0.2s
"""

def create_shifting_response(params):
    """Create a response pattern with a shifting window into PN population
    over the LFP cycles.

    The step by step process to arrive at this design is in
    try_pn_output.ipynb.

    params should have the following key value pairs:

    onset: stimulus / activity onset
    duration: duration of stimulus
    offdur: duration of off response
    pndur: duration of a single PN's activity after it starts
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
    spont_rate = Q_(params['spont_rate']).to('Hz').m
    delta = Q_(params['delta']).to('s').m
    onset = Q_(params['onset']).to('s').m
    dur = Q_(params['duration']).to('s').m
    tail = Q_(params['tail']).to('s').m
    offdur = Q_(params['offdur']).to('s').m
    osc_freq = Q_(params['osc_freq']).to('Hz').m
    stim_rate = Q_(params['stim_rate']).to('Hz').m
    osc_amp_scale = params['osc_amp_scale']
    start_frac = params['start_frac']
    # Indices of different kinds of PNs.
    n_responsive = int(round(npn * 2.0 / 3))
    if params['responsive'] >= 0:
        n_responsive = int(round(npn * params['responsive']))
    prev = 0
    nxt = npn - n_responsive
    unresponsive = list(range(prev, nxt))
    ee_frac = params['ee'] if params['ee'] >= 0 else 0.25
    ee_count = int(round(n_responsive * ee_frac))
    prev = nxt
    nxt += ee_count
    EE = list(range(prev, nxt))
    ei_frac = params['ei'] if params['ei'] >= 0 else 0.25
    ei_count = int(round(n_responsive * ei_frac))
    prev = nxt
    nxt += ei_count
    EI = list(range(prev, nxt))
    ie_frac = params['ie'] if params['ie'] >= 0 else 0.25
    ie_count = int(round(n_responsive * ie_frac))
    prev = nxt
    nxt += ie_count
    IE = list(range(prev, nxt))
    prev = nxt
    nxt = npn
    II = list(range(prev, nxt))     # rest are all II
    response_types = {'EE': EE, 'EI': EI, 'IE': IE, 'II': II, 'U': unresponsive}
    clusters = defaultdict(list)
    pn_spikes = []
    # start = timer()
    # Spontaneous activity in all neurons
    spont_rate_fn = lambda t: spont_rate
    for ii in range(params['npn']):
        st1 = nhpp.nhpp_thinning(spont_rate_fn, onset, delta)
        st2 = nhpp.nhpp_thinning(spont_rate_fn, tail, delta) + onset + dur + offdur
        pn_spikes.append(list(st1) + list(st2))

    # Unresponsive PNs - first 1/3rd
    for jj, index in enumerate(unresponsive):
        pn_spikes[index] += list(nhpp.nhpp_thinning(spont_rate_fn, dur + offdur, delta) + onset)

    # Use stop time to remove PNs from already active set
    stop_time = Q_(params['pndur']).to('s').m

    # Excited throughout odor (EE)

    # The rate function: python has dynamic binding - so ee_rate_fn
    # sees the value of `start` at evaluation time.
    start = 0
    # After every `pnclus` PNs, increase the start time.
    pnclus = round(params['shifting_frac'] * len(EE))
    for jj, index in enumerate(EE):
        ee_rate_fn = lambda t: stim_rate + osc_fn((dur - t) * stim_rate * osc_amp_scale / dur, osc_freq, t)  \
                 if  (t > start) and (t < start + stop_time) else 1e-6
        st = list(nhpp.nhpp_thinning(ee_rate_fn, dur, delta) + onset)
        pn_spikes[index] += st
        print(f'{onset + start}s: {index}')
        print('', sorted(st))
        print(pn_spikes[index])
        clusters[onset + start].append(index)
        diff = jj - len(EE) * start_frac
        if (diff > 0) and ( round(diff) % pnclus == 0):
            start += lfp_bin
    # PNs active only in first half of odour presentation (EI)
    start = 0
    ei_rate_fn = lambda t: stim_rate + osc_fn((dur - t ) * stim_rate * osc_amp_scale / dur, osc_freq, t)  \
                 if  (t > start) and (t < start + stop_time) else 1e-6
    pnclus = round(params['shifting_frac'] * len(EI))    # How many PNs should share the same starting time?
    for jj, index in enumerate(EI):
        st = list(nhpp.nhpp_thinning(ei_rate_fn, dur, delta) + onset)
        pn_spikes[index] += st
        clusters[onset + start].append(index)
        diff = jj - len(EI) * start_frac
        if (diff > 0) and (round(diff) % pnclus == 0):
            start += lfp_bin

    # PNs active only in last half of odour presentation (IE)
    start = dur / 2.0
    ie_rate_fn = lambda t: stim_rate + osc_fn((dur - t) * 2.0 * stim_rate * osc_amp_scale / dur, osc_freq, t) \
                 if  (t > start) and (t < start + stop_time) else 1e-6
    pnclus = round(params['shifting_frac'] * len(IE))    # How many PNs should share the same starting time?
    for jj, index in enumerate(IE):
        pn_spikes[index] += list(nhpp.nhpp_thinning(ie_rate_fn, dur, delta) + onset)
        clusters[onset + start].append(index)
        diff = jj - len(IE) * start_frac
        if (diff > 0) and (round(diff) % pnclus == 0):
            start += lfp_bin
    # off response
    start = 0
    off_rate_fn =  lambda t: stim_rate + osc_fn((offdur - t) * stim_rate * osc_amp_scale / offdur, osc_freq, t)  \
                   if  (t > start) and (t < start + stop_time) else 1e-6
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
    # end = timer()
    # logger.info('Finished PN output creation in {} s'.format(end - start))
    # print('Finished PN output creation in {} s'.format(end - start))
    return pn_spikes, clusters, response_types


if __name__ == '__main__':
    parser = make_parser()
    args = parser.parse_args()
    pn_cfg = yaml.load(pn_params)
    if (not args.test) and os.path.exists(args.outfile):
        print('Output file {} exists. Do nothing.'.format(args.outfile))
        sys.exit(0)
    with h5.File(args.infile, 'r') as fd:
        config = yaml.load(fd.attrs['config'])
        if args.npn > 0:
            config['pn']['number'] = args.npn
            pn_cfg['npn'] = args.npn
        if args.tail is not None:
            config['stimulus']['tail'] = args.tail
            pn_cfg['tail'] = args.tail
        if args.ig_ggn_gmax is not None:
            config['ig_ggn_syn']['gmax'] = args.ig_ggn_gmax
        if args.kc_ig_weight is not None:
            config['kc_ig_syn']['weight'] = args.kc_ig_weight
        if args.pn_duration is not None:
            pn_cfg['pndur'] = args.pn_duration
        if args.start_frac > 0:
            pn_cfg['start_frac'] = args.start_frac
        if args.shifting_frac > 0:
            pn_cfg['shifting_frac'] = args.shifting_frac
        if args.duration is not None:
            config['stimulus']['duration'] = args.duration
            pn_cfg['duration'] = args.duration
        if args.onset is not None:
            config['stimulus']['onset'] = args.onset
            pn_cfg['onset'] = args.onset
        pn_cfg['responsive'] = args.responsive
        pn_cfg['ei'] = args.ei
        pn_cfg['ie'] = args.ie
        pn_cfg['ee'] = args.ee
        spikes, clusters, response_types = create_shifting_response(pn_cfg)
    if args.test:
        x = []
        y = []
        for ii, st in enumerate(spikes):
            x.append(st.to('ms').m)
            y.append([ii] * len(st))
        x = np.concatenate(x)
        y = np.concatenate(y)
        fig, ax = plt.subplots(nrows=2, sharex='all')
        ax[0].plot(x, y, '|')
        win = 50.0
        ax[1].hist(x, bins=np.arange(0, x.max() + win / 2.0, win))
        plt.show()
    else:
        shutil.copyfile(args.infile, args.outfile)
        with h5.File(args.outfile, 'r+') as fd:
            fd.attrs['config'] = yaml.dump(config)
            try:
                descr = fd.attrs['description'].decode()
            except KeyError:
                descr = ''
            descr = '{}. Edited config from original and modified PN spike trains.'.format(descr)
            fd.attrs['description'] = descr
            pn_st_grp = fd['/data/event/pn/pn_spiketime']
            for name in list(pn_st_grp.keys()):
                del pn_st_grp[name]
            for ii, st in enumerate(spikes):
                name = 'pn_{}'.format(ii)
                ds = pn_st_grp.create_dataset(name, data=st.to('ms').m)
                ds.attrs['unit'] = 'ms'
                ds.attrs['source'] = name
                ds.attrs['field'] = 'pn_spikes'

        if args.interactive:
            with h5.File(args.outfile, 'r') as fd:
                xlist = []
                ylist = []
                for pn, st in fd['/data/event/pn/pn_spiketime'].items():
                    y = [int(pn.split('_')[-1])] * st.shape[0]
                    x = st.value
                    xlist.append(x)
                    ylist.append(y)
                xarr = np.concatenate(xlist)
                yarr = np.concatenate(ylist)
                t = xarr.max()
                win = 50.0
                fig, axes = plt.subplots(nrows=2, sharex='all')
                axes[0].plot(xarr, yarr, '|')
                axes[1].hist(xarr, np.arange(0, t + win / 2.0, win))
                plt.show()



#
# change_pn_spikes.py ends here
