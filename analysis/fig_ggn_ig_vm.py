# fig_ggn_ig_vm.py --- 
# Author: Subhasis  Ray
# Created: Wed Dec  5 10:58:13 2018 (-0500)
# Last-Updated: Mon Apr 15 15:01:47 2019 (-0400)
#           By: Subhasis  Ray
# Version: $Id$

# Code:

"""Plot reproduction of hyperpolarization of GGN membrane potential """
from __future__ import print_function
import sys
if sys.platform == 'win32':
    sys.path.append('D:/subhasis_ggn/model/analysis')
else:
    sys.path += ['/home/rays3/projects/ggn/analysis',
                 '/home/rays3/projects/ggn/morphutils',
                 '/home/rays3/projects/ggn/nrn',
                 '/home/rays3/apps/matrix-entropy',
                 '/home/rays3/apps/Spike-Contrast/Python3']
import os
import argparse
import h5py as h5
import numpy as np
import time
# import random
import matplotlib
import matplotlib.colors as colors
from matplotlib import pyplot as plt
# import yaml
# import pint
import collections
import operator
import itertools as it
import pandas as pd
import yaml
import network_data_analysis as nda
# from matplotlib.backends.backend_pdf import PdfPages
import pn_kc_ggn_plot_mpl as myplot
import clust_fixed_net_multi_trial as cl

plt.rc('font', size=11)

analysis_dir = '/home/rays3/projects/ggn/analysis'
# datadir = '/data/rays3/ggn/olfactory_network'
datadir = '/data/rays3/ggn/fixed_net'
os.chdir(analysis_dir)

fig, ax = plt.subplots(nrows=3, ncols=2, sharex='all', sharey='row')
# jid = '13081213'   # this reproduces GGN Vm well, PN->KC clustered, same KC cluster receives coactive PNs
jid = '21841553'   # new simulation 
fname = nda.find_h5_file(jid, datadir)
originals = []
with h5.File(fname, 'r') as fd:
    ax[0, 0].set_title(f'{jid}')
    print(fname)
    orig = fd.attrs['original']
    originals.append(orig)
    print(f'original:{orig}')
    with h5.File(orig, 'r') as forig:
        print(yaml.dump(nda.load_config(forig), default_flow_style=False))
    # fig, ax = plt.subplots(nrows=2, ncols=1, sharex='all')
    myplot.plot_ggn_vm(ax[1, 0], fd, fd['/data/uniform/ggn_output/GGN_output_Vm'], 'LCA', 1, color='black', alpha=1.0)
    ax[1, 0].set_ylim(-60, -45)
    ax[1, 0].set_yticks([-55, -50])
    ig_vm = fd['/data/uniform/ig/IG_Vm']
    dt = ig_vm.attrs['dt']
    ig_vm = ig_vm[0, :]
    t = np.arange(len(ig_vm)) * dt
    ax[2, 0].plot(t, ig_vm, color='black')
    ax[2, 0].hlines(y=-60.0, xmin=500, xmax=1500, color='gray', lw=10)
    ax[2, 0].hlines(y=-60.0, xmin=1500, xmax=2000, color='lightgray', lw=10)
    for axis in ax.flat:
        [sp.set_visible(False) for sp in axis.spines.values()]
        axis.tick_params(right=False, top=False)    
    ax[1, 0].tick_params(bottom=False)
    xticks = [200.0, 1000.0, 2000.0, 3000.0]    
    ax[2, 0].set_xticks(xticks)
    ax[2, 0].set_xticklabels([x/1000.0 for x in xticks])
    ax[2, 0].set_xlabel('Time (s)')
    ax[2, 0].set_xlim(200.0, 3000.0)
    ax[2, 0].set_ylabel('Membrane potential (mV)')
    x = []
    y = []
    for pn, st in fd['/data/event/pn/pn_spiketime'].items():
        y.append([int(pn.split('_')[-1])] * st.shape[0])
        x.append(st.value)
    ax[0, 0].set_xlabel('Time (ms)')
    ax[0, 0].set_ylabel('PN #')
    ax[0, 0].plot(np.concatenate(x[::10]), np.concatenate(y[::10]), color='#fdb863', marker='s', ms=3, ls='none')
    ax[0, 0].set_xlim(200.0, 2500.0)

jid = '21841558'   # new simulation
fname = nda.find_h5_file(jid, datadir)
with h5.File(fname, 'r') as fd:
    print(fname)
    orig = fd.attrs['original']
    originals.append(orig)
    print(f'original:{orig}')
    with h5.File(orig, 'r') as forig:
        print(yaml.dump(nda.load_config(forig), default_flow_style=False))
    myplot.plot_ggn_vm(ax[1, 1], fd, fd['/data/uniform/ggn_output/GGN_output_Vm'], 'LCA', 1, color='black', alpha=1.0)
    ig_vm = fd['/data/uniform/ig/IG_Vm']
    dt = ig_vm.attrs['dt']
    ig_vm = ig_vm[0, :]
    t = np.arange(len(ig_vm)) * dt
    ax[2, 1].plot(t, ig_vm, color='black')
    ax[2, 1].hlines(y=-60.0, xmin=500, xmax=1500, color='gray', lw=10)
    ax[2, 1].hlines(y=-60.0, xmin=1500, xmax=2000, color='lightgray', lw=10)
    for axis in ax.flat:
        [sp.set_visible(False) for sp in axis.spines.values()]
        axis.tick_params(right=False, top=False)    
    ax[1, 1].tick_params(bottom=False)
    xticks = [200.0, 1000.0, 2000.0, 3000.0]    
    ax[2, 1].set_xticks(xticks)
    ax[2, 1].set_xticklabels([x/1000.0 for x in xticks])
    ax[2, 1].set_xlabel('Time (s)')
    ax[2, 1].set_xlim(200.0, 3000.0)
    ax[2, 1].set_ylabel('Membrane potential (mV)')
    x = []
    y = []
    for pn, st in fd['/data/event/pn/pn_spiketime'].items():
        y.append([int(pn.split('_')[-1])] * st.shape[0])
        x.append(st.value)
    ax[0, 1].set_title(f'{jid}')
    ax[0, 1].set_xlabel('Time (ms)')
    ax[0, 1].set_ylabel('PN #')
    ax[0, 1].plot(np.concatenate(x[::10]), np.concatenate(y[::10]), color='#fdb863', marker='s', ms=3, ls='none')
    ax[0, 1].set_xlim(200.0, 2500.0)
fig.subplots_adjust(top=0.95, bottom=0.15, left=0.2, right=0.95, hspace=0.1, wspace=0.1)
fig.set_frameon(False)
fig.set_size_inches(9/2.54, 12/2.54)
fig.savefig('pn_ggn_ig.svg'.format(jid))
plt.show()
    # print(nda.load_config(fd))

# fdl = h5.File(originals[0], 'r')
# fdr = h5.File(originals[1], 'r')

# pn_kc_syn_path = '/data/static/pn_kc_synapse/pn_kc_synapse'
# ggn_kc_syn_path = '/data/static/ggn_kc_synapse/ggn_kc_synapse'
# kc_ggn_syn_path = '/data/static/kc_ggn_synapse/kc_ggn_synapse'
# kc_ig_syn_path =  '/data/static/kc_ggn_synapse/kc_ggn_synapse'

# for left, right in zip(fdl[pn_kc_syn_path][:, 0], fdr[pn_kc_syn_path][:, 0]):
#     # print(left)
#     # print(right)
#     assert left == right

# fdl.close()
# fdr.close()
    

# jids = [
#     '13014324',  # 	3 	20 	1 	12028 	-34
#     '13014327',  # 	3 	20 	1 	17923 	-30
#     '13014342',  # 	2.5 	20 	1 	2772 	-43
#     '13014344',  # 	3.5 	20 	0.9 	3870 	-45
#     '13014347',  # 	3.5 	20 	0.9 	3880 	-46
# ]

# for jid in jids:
#     fname = nda.find_h5_file(jid, datadir)
#     with h5.File(fname, 'r') as fd:
#         fig, ax = plt.subplots(nrows=2, ncols=1, sharex='all')
#         myplot.plot_ggn_vm(ax[1, 1], fd, fd['/data/uniform/ggn_output/GGN_output_Vm'], 'LCA', 1, color='black', alpha=1.0)
#         ig_vm = fd['/data/uniform/ig/IG_Vm']
#         dt = ig_vm.attrs['dt']
#         ig_vm = ig_vm[0, :]
#         t = np.arange(len(ig_vm)) * dt
#         ax[1].plot(t, ig_vm, color='black')
#         for axis in ax.flat:
#             [sp.set_visible(False) for sp in axis.spines.values()]
#             axis.tick_params(right=False, top=False)    
#         ax[1, 1].tick_params(bottom=False)
#         xticks = [200.0, 1000.0, 2000.0, 3000.0]
#         ax[1].set_xticks(xticks)
#         ax[1].set_xticklabels([x/1000.0 for x in xticks])
#         ax[1].set_xlabel('Time (s)')
#         ax[1].set_xlim(200.0, 3500.0)
#         fig.set_frameon(False)
#         fig.set_size_inches(9/2.54, 9/2.54)
#         # fig.savefig('ggn_ig_vm_dep_hyp.svg', transparent=True)
#         plt.show()
    
# 
# fig_ggn_ig_vm.py ends here
