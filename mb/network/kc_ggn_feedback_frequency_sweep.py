# kc_ggn_feedback_frequency_sweep.py ---
#
# Filename: kc_ggn_feedback_frequency_sweep.py
# Description:
# Author: Subhasis Ray
# Maintainer:
# Created: Fri Apr 27 11:53:11 2018 (-0400)
# Version:
# Package-Requires: ()
# Last-Updated: Fri May  4 18:00:03 2018 (-0400)
#           By: Subhasis Ray
#     Update #: 509
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#


# Code:
"""To understand the frequency limits of the KC->GGN feedback circuit
I make a simplified model of a KC<->GGN circuit with a single KC. The
output of the KC is fed to the GGN via a very strong synapse. It is
also possible to introduce jitter by making many synapses with
different delays.

"""
from __future__ import print_function
import sys
import argparse
import numpy as np
from matplotlib import pyplot as plt
import h5py as h5
from config import Q_, h, logger, timestamp, mypid, myjobid, nrn_version
from timeit import default_timer as timer
import ephys
import nrnutils as nu
import networkx as nx
import neurograph as ng
import nsdf
from pn_kc_ggn_network import make_exp2syn_one_one


def make_kc_ggn_alpha_loop(ggn, kc, ggn_kc_syn_params, kc_ggn_syn_params, jitter=None):
    """Create a network with single KC connected via a single strong
    synapse (equivalent to many synapses from many synchronous KCs)

    jitter, if specified, should be (time_range, nsyn) where nsyn is
    the number of KC->GGN synapses to be created and each will have a
    delay increased by a uniform random time between 0 and time_range
    (the base delay is specified in kc_ggn_syn_params).

    """
    ggn_graph = nu.nrngraph(ggn)
    if jitter is not None:
        print('Jitter:', jitter, type(jitter[0]), type(jitter[1]))
        trange = float(jitter[0])
        nsyn = int(jitter[1])
    else:
        nsyn = 1
        trange = 0.0
    alpha_sec_pos = nu.select_random_terminal_segments_by_sid(ggn_graph,
                                                     ng.custom_types['alphaL'],
                                                     nsyn)    
    alpha_sec, alpha_pos = list(zip(*alpha_sec_pos))
    kc_ggn_syn_info = nu.make_exp2syn_one_one([kc.soma] * nsyn, [0.5] * nsyn,
                                              alpha_sec, alpha_pos,
                                              kc_ggn_syn_params)
    if nsyn > 1:
        delays = np.random.uniform(low=0.0, high=trange, size=nsyn)
        for netcon, delay in zip(kc_ggn_syn_info['netcons'], delays):
            netcon.delay += delay        
    
    ca_sec_pos = nu.select_random_terminal_segments_by_sid(
        ggn_graph, ng.custom_types['LCA'], 1)
    ca_sec, ca_pos = list(zip(*ca_sec_pos))
    ggn_kc_syn_info = nu.make_gradedsyn_one_one(ca_sec, ca_pos, [kc.soma], [0.5],
                                                ggn_kc_syn_params)
    return {'kc_ggn_syn_info': kc_ggn_syn_info,
            'ggn_kc_syn_info': ggn_kc_syn_info,
            'ggn_pre': ca_sec_pos,
            'ggn_post': alpha_sec_pos,
    }


def make_driving_current(sec, pos, amp, freq, tstart, tend):
    inject = h.Izap(sec(pos))
    setattr(inject, 'del', tstart)
    inject.dur = tend - tstart
    inject.f0 = inject.f1 = freq
    inject.amp = amp
    return inject


def init_run_save(group_vec_list, dsname, tstop, attrs):
    """init, run and save the data in fd, creating the dataset `dsname`.

    vec_dict: {field name: Vector recording the field}

    tstop: stop time of simulation

    attrs: dict containing additional attributes for the dataset.
    """
    h.tstop = tend
    logger.info('Init: saving in dataset {}'.format(dsname))
    h.init()
    h.run()
    print('Finished run: saving in dataset {}'.format(dsname))
    for group, vec_dict in group_vec_list:
        data = np.vstack([np.asarray(vec.x) for vec in vec_dict.values()])
        ds = group.create_dataset(dsname, data=data)
        ds.attrs['fields'] = list(vec_dict.keys())
        for k, v in attrs.items():
            ds.attrs[k] = v
        

def run_freq_sweep(group_vec_izap_list, tend=3.1e3,
                   amp=0.1, f0=0.0, f1=2, step=0.2, log=True):
    """Run a frequency sweep from f0 to f1. If log is true, then
       frequencies are taken as expressed in log10 scale. Otherwise
       linear.

    """
    for ii, f in enumerate(np.arange(f0, f1, step)):        
        freq = 10**f if log else f
        logger.info('Frequency {} Hz'.format(freq))
        for item in group_vec_izap_list:
            izap = item[-1]
            izap.f0 = izap.f1 = freq
        group_vec_list = [a[:2] for a in group_vec_izap_list]
        dsname = 'freq_{}'.format(ii)
        init_run_save(group_vec_list, dsname, tend, {'f': freq})
        
        
def run_amp_sweep(group_vec_iclamp_list, tend=3.1e3,
                  amp0=0.01, amp1=1.0, step=0.2):
    for ii, amp in enumerate(np.arange(amp0, amp1, step)):
        logger.info('Amplitude {} nA'.format(amp))
        for item in group_vec_iclamp_list:
            iclamp = item[-1]
            iclamp.amp = amp
        dsname = 'amp_{}'.format(ii)
        group_vec_list = [a[:2] for a in group_vec_iclamp_list]
        init_run_save(group_vec_list, dsname, tend, {'amp': amp})
        

def make_parser():
    parser = argparse.ArgumentParser(description='Simulate KC<->GGN feedback'
                                     ' frequency response')
    parser.add_argument('--kc-file', type=str, dest='kcfile',
                        required=True, help='KC cell template file'
                        ' (.hoc)')
    parser.add_argument('--kc', type=str, dest='kc', required=True,
                        help='KC cell template name in template file')
    # parser.add_argument('--kc-count', type=int, dest='kccount',
    #                     required=True, action='append',
    #                     help='Number of copies of KC instances')
    parser.add_argument('--ggn-file', type=str, dest='ggnfile',
                        required=True, help='GGN cell template file'
                        '(.hoc)')
    parser.add_argument('--ggn', type=str, dest='ggn', required=True,
                        help='GGN cell template name in template '
                        'file')
    parser.add_argument('--nkc', type=int, dest='nkc', default=50000,
                        help='number of KCs (KC->GGN synaptic '
                        'conductance multiplied by this number)')
    parser.add_argument('--out', type=str, dest='outprefix', default='/data/rays3/ggn/kc_ggn_feedback/'
                        'kc_ggn_freq_sweep',
                        help='Output data file prefix')
    parser.add_argument('-j', '--jitter', type=float, nargs=2, help='(trange, nsyn) add '
                        'uniformly distributed jitter with (1) in time'
                        ' range `trange` specified in ms and (2) number of '
                        'synapses `nsyn`.')
    parser.add_argument('--freq', type=float, dest='freq', nargs='+',
                        default=[0], help='The frequency (in Hz) of'
                        ' injected current, or three numbers '
                        'specifying the start, stop and increment')
    parser.add_argument('--flog', action='store_true', dest='flog',
                        help='Frequency sweep in log scale')
    parser.add_argument('--amp', type=float, dest='amp', nargs='+', default=[1e-2, 1e-1, 5e-3], 
                        help='The amplitude (in nA) of injected current, or three numbers'
                        ' specifying the start, stop and increment')
    return parser


if __name__ == '__main__':
    args = make_parser().parse_args()
    logger.info('Command line: {}'.format(' '.join(sys.argv)))
    ggn_kc_syn_params = {
        'vmid': Q_('-40mV').to('mV').m,
        'vslope': Q_('5.0mV').to('mV').m,
        'e': Q_('-80mV').to('mV').m,
        'gbar': Q_('1e-3uS').to('uS').m,
        'tau': Q_('4.0ms').to('ms').m
    }
    # Note: I increased threshold to avoid transmission during low  frequency stimulus
    # where KC can have sustained depolarization.
    kc_ggn_syn_params = {
        'threshold': Q_('-10.0mV').to('mV').m,
        'delay': Q_('1.0ms').to('ms').m,
        'e': Q_('0.0mV').to('mV').m,
        'tau1': Q_('13.333ms').to('ms').m,
        'tau2': Q_('13.333ms').to('ms').m,
        'gmax': Q_('5pS').to('uS').m * args.nkc
    }
    kc = nu.create_cell(args.kc, filename=args.kcfile)
    kc_solo = nu.create_cell(args.kc, filename=args.kcfile)
    ggn = nu.create_cell(args.ggn, filename=args.ggnfile)
    model = make_kc_ggn_alpha_loop(ggn, kc, ggn_kc_syn_params,
                                   kc_ggn_syn_params, args.jitter)

    tstart = 0.5e3
    tend = 3.1e3
    if len(args.freq) > 1:
        inject = make_driving_current(kc.soma, 0.5, args.amp[0], 1.0, tstart, tend)
        inject_solo = make_driving_current(kc_solo.soma, 0.5, args.amp[0], 1.0, tstart, tend)
    else:
        inject = ephys.setup_current_clamp(kc.soma, delay=Q_(tstart, 'ms'),
                                           duration=Q_((tend - tstart), 'ms'),
                                           amplitude=Q_(args.amp[0], 'nA'), pos=0.5)
        inject_solo = ephys.setup_current_clamp(kc_solo.soma, delay=Q_(tstart, 'ms'),
                                                duration=Q_((tend - tstart), 'ms'),
                                                amplitude=Q_(args.amp[0], 'nA'), pos=0.5)
    tvec = h.Vector()
    kc_vvec = h.Vector()
    ggn_alpha_vvec = h.Vector()
    ggn_ca_vvec = h.Vector()
    ivec = h.Vector()
    kc_solo_vvec = h.Vector()
    ivec_solo = h.Vector()
    

    tvec.record(h._ref_t)
    kc_vvec.record(kc.soma(0.5)._ref_v)
    ggn_alpha_vvec.record(model['ggn_post'][0][0](model['ggn_post'][0][1])._ref_v)
    ggn_ca_vvec.record(model['ggn_pre'][0][0](model['ggn_pre'][0][1])._ref_v)
    ivec.record(inject._ref_i)

    kc_solo_vvec.record(kc_solo.soma(0.5)._ref_v)
    ivec_solo.record(inject_solo._ref_i)
    
    fb_dict = {'t': tvec, 'i': ivec, 
                'kc': kc_vvec,
                'alpha': ggn_alpha_vvec,
                'ca': ggn_ca_vvec}
    solo_dict = {'t': tvec, 'i': ivec_solo, 
                'kc': kc_solo_vvec}
    outfile = '{}_UTC{}_PID{}_JID{}.h5'.format(
        args.outprefix,
        timestamp.strftime('%Y_%m_%d__%H_%M_%S'),
        mypid, myjobid)
    logger.info('Starting simulation')
    with h5.File(outfile, 'w') as fd:
        fd.attrs['command_line'] = str(sys.argv)
        fd.attrs['ggn_kc_syn'] = str(ggn_kc_syn_params)
        fd.attrs['kc_ggn_syn'] = str(kc_ggn_syn_params)
        syntype = np.dtype([('sec', h5.special_dtype(vlen=str)),
                            ('pos', np.float32)])
        ds_alpha_syn = fd.create_dataset('ggn_post',
                                         (len(model['ggn_post']),),
                                         dtype=syntype)
        ds_alpha_syn[:] = [(sec.name(), pos) for sec, pos in model['ggn_post']]
        ds_ca_syn = fd.create_dataset('ggn_pre', (len(model['ggn_pre']),),
                                         dtype=syntype)
        ds_ca_syn[:] = [(sec.name(), pos) for sec, pos in model['ggn_pre']]
        grp_feedback = fd.create_group('kc_ggn')
        grp_solo = fd.create_group('kc')
        grp_vec_inject = [(grp_feedback, fb_dict, inject),
                          (grp_solo, solo_dict, inject_solo)]
        if len(args.freq) > 1:
            logger.info('Starting frequency sweep')
            run_freq_sweep(grp_vec_inject, tend=tend,
                           amp=args.amp[0], f0=args.freq[0],
                           f1=args.freq[1], step=args.freq[2],
                           log=args.flog)
        else:
            logger.info('Starting amplitude sweep')
            run_amp_sweep(grp_vec_inject, tend=tend,
                          amp0=args.amp[0], amp1=args.amp[1], step=args.amp[2])
    logger.info('Data saved in {}'.format(outfile))
            
    # fig, axes = plt.subplots(nrows=2, ncols=1, sharex='all')
   # axes[0].legend()
    # axes[1].legend()
    # plt.show()


# 
# kc_ggn_feedback_frequency_sweep.py ends here
