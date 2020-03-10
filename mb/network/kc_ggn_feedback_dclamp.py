# kc_ggn_feedback_dclamp.py --- 
# Author: Subhasis Ray
# Created: Tue Aug 20 10:58:08 2019 (-0400)
# Last-Updated: Wed Dec 11 17:32:49 2019 (-0500)
#           By: Subhasis Ray
# Version: $Id$

# Code:
"""This script for testing expansion of the dynamic range of a KC due to GGN inhibition.

Instead of running a whole simulation in the full network, we play the
GGN membrane potential back to the KC.

We use GGN Vm from two simulations, one with low PN activity and
another with high PN activity.

"""
from __future__ import print_function
import os
import sys
sys.path += ['D:/subhasis_ggn/model/mb', 'D:/subhasis_ggn/model/mb/network', 'D:/subhasis_ggn/model/nrn']
import argparse
import numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
import h5py as h5
from config import Q_, h, logger, timestamp, mypid, myjobid, nrn_version
from timeit import default_timer as timer
import ephys
import nrnutils as nu
import neurograph as ng
import nsdf


GGN_KC_SYN_PARAMS = {
    'vmid': Q_('-40mV').to('mV').m,
    'vslope': Q_('5.0mV').to('mV').m,
    'e': Q_('-80mV').to('mV').m,
    'gbar': Q_('1e-3uS').to('uS').m,
    'tau': Q_('4.0ms').to('ms').m
}

# keep global reference of created model components so that they are
# not garbage collected when out of scope
model_dict = {}

def make_kc_with_dynaclamp(kc_name, kc_file, inject, tstart, tend, ggn_vm=None):
    """Read KC model from `kc_file`, inject current `inject` nA, apply
    dynamic clamp `ggn_vm`, which should be a 2D array with time (ms)
    in column 0, and voltage (mV) in column 1.

    """
    global model_dict
    kc = nu.create_cell(kc_name, filename=kc_file)
    model_dict[kc] = None
    iclamp = ephys.setup_current_clamp(kc.soma, pos=0.5, delay=Q_(tstart, 'ms'),
                                       duration=Q_((tend - tstart), 'ms'),
                                       amplitude=Q_(inject, 'nA'))
    model_dict[iclamp] = None
    ggn_g_vec = None
    if ggn_vm is not None:
        syn = h.GradedSyn(kc.soma(0.5))
        for attr, value in GGN_KC_SYN_PARAMS.items():
            setattr(syn, attr, value)
        model_dict[syn] = None
        ggn_comp = h.Section('ggn')
        model_dict[ggn_comp] = None
        h.setpointer(ggn_comp(0.5)._ref_v, 'vpre', syn)
        ggn_vm_vec = h.Vector(ggn_vm[:, 1])
        tvec = h.Vector(ggn_vm[:, 0])
        model_dict[tvec] = None
        # vec.play(var_reference, t, continuous) for interpolating        
        ret = ggn_vm_vec.play(ggn_comp(0.5)._ref_v, tvec, 1)
        print('####', ret)
        model_dict[ggn_vm_vec] = None
        ggn_g_vec = h.Vector()
        ggn_g_vec.record(syn._ref_g)
        model_dict[ggn_g_vec] = None
    kc_vm_vec = h.Vector()
    kc_vm_vec.record(kc.soma(0.5)._ref_v)
    model_dict[kc_vm_vec] = None
    print('Built model')
    return (kc_vm_vec, ggn_g_vec)

def make_parser():
    parser = argparse.ArgumentParser(description='Simulate KC with GGN inhibition at multiple current injections')
    parser.add_argument('--kc-file', type=str, dest='kc_file',
                        required=True, help='KC cell template file'
                        ' (.hoc)')
    parser.add_argument('--kc', type=str, dest='kc', required=True,
                        help='KC cell template name in template file')
    parser.add_argument('--ggn-vm-file', type=str, action='append', dest='ggn_vm_file',
                        required=True,
                        help='CSV file with column 0 time in ms, column 1 GGN Vm in mV')
    parser.add_argument('--istart', type=str, dest='istart', help='Starting amplitude of current (with unit)')
    parser.add_argument('--iend', type=str, dest='iend', help='Ending amplitude of current (with unit)')
    parser.add_argument('--di', type=str, dest='di', help='Current increments (with unit)')
    parser.add_argument('--tstart', type=str, dest='tstart', help='Current injection start time (with unit)')
    parser.add_argument('--tend', type=str, dest='tend', help='Current injection end time (with unit)')    
    return parser


def main():
    parser = make_parser()
    args = parser.parse_args()
    logger.info('Command line args: {}'.format(str(sys.argv)))
    print(args.ggn_vm_file)
    # KCs with GGN inhibition
    inhibited_vec = defaultdict(list)
    solo_vec_list = []
    tstart = Q_(args.tstart).to('ms').m
    tend = Q_(args.tend).to('ms').m
    istart = Q_(args.istart).to('nA').m
    iend = Q_(args.iend).to('nA').m
    di = Q_(args.di).to('nA').m
    irange = np.arange(istart, iend + di/2.0, di)
    logger.info('Starting current: {} nA'.format(istart))
    logger.info('End current: {} nA'.format(iend))
    logger.info('Increment: {} nA'.format(di))
    logger.info('current range: {}'.format(irange))
    ggn_vm = {}
    for input_file in args.ggn_vm_file:
        ggn_vm[input_file] = np.loadtxt(input_file)
    for inject in irange:
        for input_file, vm in ggn_vm.items():
            kc_vvec, ggn_gvec = make_kc_with_dynaclamp(args.kc, args.kc_file, inject, tstart, tend, vm)
            inhibited_vec[input_file].append((kc_vvec, ggn_gvec))
        # KC without any inhibition
        kc_vvec, ggn_gvec = make_kc_with_dynaclamp(args.kc, args.kc_file, inject, tstart, tend)
        solo_vec_list.append(kc_vvec)
    tvec = h.Vector()
    tvec.record(h._ref_t)
    h.tstop = tend
    print('Init')    
    h.init()
    print('Run')
    h.run()
    print('Finished simulation')
    fig, ax = plt.subplots(nrows=len(irange)+1, ncols=len(ggn_vm)+1, sharex='all', sharey='all')
    t = np.array(tvec.x)
    solo_data = []
    for ii, vvec in enumerate(solo_vec_list):
        ax[ii+1, 0].plot(tvec, vvec, color='#e66101')
        solo_data.append(np.array(vvec.x))
    combined = np.vstack(solo_data)
    
    prefix = 'UTC' + timestamp.strftime('%Y%m%d_%H%M%S')
    fname = '{}_solo_kc.npz'.format(prefix)
    np.savez(fname,
             t=t,
             vm=combined,
             inject=irange)
    logger.info('Saved solo KC data in {}'.format(fname))
    for jj, input_file in enumerate(args.ggn_vm_file):
        fname = '{}_{}.npz'.format(prefix, os.path.basename(input_file))
        data = []
        kc_vm_list = inhibited_vec[input_file]
        for ii, (vvec, gvec) in enumerate(kc_vm_list):
            data.append(np.array(vvec.x))
            ax[ii+1, jj+1].plot(tvec, vvec, color='#e66101')
            ax[ii+1, 0].set_ylabel('{} pA'.format(irange[ii]*1e3))
            # ax[ii+1, 0].set_ylabel('{} pA'.format(int(np.round(irange[ii]*1e3))))  # to avoid decimal point when integer values 
        # ax[0, jj+1].plot(tvec, gvec)
        # ax[0, jj+1].plot(ggn_vm[input_file][:,0], ggn_vm[input_file][:,1])
        ax[0, jj+1].set_title(input_file)
        combined = np.vstack(data)
        np.savez(fname, combined=combined, irange=irange, ggn_vm=ggn_vm[input_file])
        logger.info('Saved data from dynamic clamp with input from {} in {}'.format(
            input_file, fname))
    for axis in ax.flat:
        axis.set_xlim(250, 1750)
        fig.set_size_inches(210/25.4, 290/25.4)
    fig.tight_layout()
    fig.savefig('{}_KC_dynamic_range_with_ggn_vm.svg'.format(prefix))
    plt.show()
    print('End')

    
if __name__ == '__main__':
    main()
    

# 
# kc_ggn_feedback_dclamp.py ends here
