# kc_ggn_nofeedback.py ---
#
# Filename: kc_ggn_nofeedback.py
# Description:
# Author: Subhasis Ray
# Maintainer:
# Created: Fri Apr 27 11:53:11 2018 (-0400)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Jun 11 09:47:43 2019 (-0400)
#           By: Subhasis  Ray
#     Update #: 531
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#


# Code:
"""This is a control case for isolated KC vs KC with feedback
inhibition from GGN. Here KC receives baseline inhibition from GGN,
but its activity has no effect on GGN.

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


def make_ggn_kc_conn(ggn, kc, ggn_kc_syn_params):
    """Create a network with single KC receiving inhibition from an
    isolated GGN.

    """
    ggn_graph = nu.nrngraph(ggn)
    ca_sec_pos = nu.select_random_terminal_segments_by_sid(
        ggn_graph, ng.custom_types['LCA'], 1)
    ca_sec, ca_pos = list(zip(*ca_sec_pos))
    ggn_kc_syn_info = nu.make_gradedsyn_one_one(ca_sec, ca_pos, [kc.soma], [0.5],
                                                ggn_kc_syn_params)
    return {'ggn_kc_syn_info': ggn_kc_syn_info,
            'ggn_pre': ca_sec_pos,            
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
    parser.add_argument('--ggn-file', type=str, dest='ggnfile',
                        required=True, help='GGN cell template file'
                        '(.hoc)')
    parser.add_argument('--ggn', type=str, dest='ggn', required=True,
                        help='GGN cell template name in template '
                        'file')
    parser.add_argument('--out', type=str, dest='outprefix', default='/data/rays3/ggn/kc_ggn_feedback/'
                        'kc_ggn_freq_sweep',
                        help='Output data file prefix')
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
        'gmax': Q_('5pS').to('uS').m
    }
    kc = nu.create_cell(args.kc, filename=args.kcfile)
    kc_solo = nu.create_cell(args.kc, filename=args.kcfile)
    ggn = nu.create_cell(args.ggn, filename=args.ggnfile)
    model = make_ggn_kc_conn(ggn, kc, ggn_kc_syn_params)
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
    ggn_ca_vvec = h.Vector()
    ivec = h.Vector()
    kc_solo_vvec = h.Vector()
    ivec_solo = h.Vector()
    

    tvec.record(h._ref_t)
    kc_vvec.record(kc.soma(0.5)._ref_v)
    ggn_ca_vvec.record(model['ggn_pre'][0][0](model['ggn_pre'][0][1])._ref_v)
    ivec.record(inject._ref_i)

    kc_solo_vvec.record(kc_solo.soma(0.5)._ref_v)
    ivec_solo.record(inject_solo._ref_i)
    
    fb_dict = {'t': tvec, 'i': ivec, 
                'kc': kc_vvec,
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
        syntype = np.dtype([('sec', h5.special_dtype(vlen=str)),
                            ('pos', np.float32)])
        ds_ca_syn = fd.create_dataset('ggn_pre', (len(model['ggn_pre']),),
                                         dtype=syntype)
        ds_ca_syn[:] = [(sec.name(), pos) for sec, pos in model['ggn_pre']]
        grp_feedback = fd.create_group('ggn_kc')
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


## samp[le command line    
## network/kc_ggn_nofeedback.py --kc-file=cell_templates/kc_1_comp.hoc --kc=KC --ggn-file=cell_templates/GGN_20170309_sc.hoc --ggn=GGN_20170309_sc --amp 0.01 0.035 0.001 --out /data/rays3/ggn/kc_ggn_feedback/kc_ggn_amp_sweep_nofeedback
    

# 
# kc_ggn_nofeedback.py ends here
