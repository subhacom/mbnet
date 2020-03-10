# fig_ggn_hyperpolarization.py --- 
# Author: Subhasis  Ray
# Created: Tue Feb 12 14:53:34 2019 (-0500)
# Last-Updated: Tue Feb 12 15:37:32 2019 (-0500)
#           By: Subhasis  Ray
# Version: $Id$

# Code:
"""Supplementary figures showing GGN hyperpolarization when PN->KC
connection is not clustered"""


from __future__ import print_function
import sys
if sys.platform == 'win32':
    sys.path.append('D:/subhasis_ggn/model/analysis')
else:
    sys.path += ['/home/rays3/projects/ggn/analysis',
                 '/home/rays3/projects/ggn/morphutils',
                 '/home/rays3/projects/ggn/nrn']
import os
import h5py as h5
import numpy as np
import random
from matplotlib import pyplot as plt

import network_data_analysis as nda
import pn_kc_ggn_plot_mpl as myplot


datadir = '/data/rays3/ggn/olfactory_network'

#jid = '20161955'
# jid = '20162331'
jid = '20162116'
fname = nda.find_h5_file(jid, datadir)
fig, ax = plt.subplots()
with h5.File(fname, 'r') as fd:
    spike_counts = [ds.shape[0] for ds in fd['/data/event/kc/kc_spiketime'].values()]
    spiking = len(np.flatnonzero(spike_counts))
    print('spiking kcs:{}'.format(spiking))
    myplot.plot_ggn_vm(
        ax, fd,
        fd['/data/uniform/ggn_basal/GGN_basal_Vm'], 'dend_b', 1)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('GGN Vm (mV)')
    ax.tick_params(right=False, top=False)
    ax.set_xlim(500, 2500)
fig.set_frameon(False)
fig.set_size_inches(20/2.54, 10/2.54)
fig.savefig('supp_4_ggn_ig_hyperpolarization.svg')    



# 50.0,1.4,2.7,20162064
jid = '20162064'
fname = nda.find_h5_file(jid, datadir)
fig, ax = plt.subplots()
with h5.File(fname, 'r') as fd:
    spike_counts = [ds.shape[0] for ds in fd['/data/event/kc/kc_spiketime'].values()]
    spiking = len(np.flatnonzero(spike_counts))
    print('spiking kcs:{}'.format(spiking))
    myplot.plot_ggn_vm(ax, fd,
                       fd['/data/uniform/ggn_basal/GGN_basal_Vm'], 'dend_b', 1)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('GGN Vm (mV)')
    ax.set_xlim(500, 2500)
    ax.tick_params(right=False, top=False)
fig.set_frameon(False)
fig.set_size_inches(20/2.54, 10/2.54)
fig.savefig('supp_5_ggn_ig_de_hyperpolarization.svg')    




# 
# fig_ggn_hyperpolarization.py ends here
