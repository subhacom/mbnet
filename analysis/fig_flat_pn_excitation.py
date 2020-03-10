# fig_flat_pn_excitation.py --- 
# Author: Subhasis Ray
# Created: Wed Dec 19 13:36:58 2018 (-0500)
# Last-Updated: Thu Dec 20 19:31:50 2018 (-0500)
#           By: Subhasis Ray
# Version: $Id$

# Code:
"""This script plots the GGN response when the PN excitation is flat
with a constant active population (not shifting)"""

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
from matplotlib import gridspec
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

analysis_dir = 'D:/subhasis_ggn/model/analysis'
datadir = 'D:/biowulf_stage/olfactory_network'
os.chdir(analysis_dir)

# jid = 9674116  # this one produced spikes in 4473 KCs out of 50k

# jids = [
#     '9674118',
#     '9674119',
#     '9673663',
#     '9674019',
#     '9673664',
#     '9674069',
#     '9674087',
#     '9674120',
#     '9673863']
# ^-These will not do - constant GGN->KC inhibition with 2% weakly inhibited


jid = '9674118'
fname = nda.find_h5_file(jid, datadir)

gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3, 1], hspace=0.05)
fig = plt.figure()
ax0 = fig.add_subplot(gs[0]) 
ax1 = fig.add_subplot(gs[1], sharex=ax0)
axes = [ax0, ax1]
with h5.File(fname, 'r') as fd:
    print('jid: {} spiking KCs: {}'.format(jid, len(nda.get_spiking_kcs(fd))))
    print(yaml.dump(nda.load_config(fd), default_style=''))
    # pn_st = []
    # pn_id = []
    # for pn in fd[nda.pn_st_path].values():
    #     pn_st.append(pn.value)
    #     pn_id.append([int(pn.name.rpartition('_')[-1])] * len(pn))    
    # ax0.plot(np.concatenate(pn_st[::10]), np.concatenate(pn_id[::10]), 'k,')
    kc_x, kc_y = nda.get_event_times(fd[nda.kc_st_path])
    ax0.plot(np.concatenate(kc_x[::10]), np.concatenate(kc_y[::10]), 'k,')
    myplot.plot_ggn_vm(ax1, fd, fd['/data/uniform/ggn_basal/GGN_basal_Vm'], 'dend_b', 1, color='k')
    ax1.set_ylim(-51, -45)
    ax1.set_xlim(700, 2500)
    xticks = np.array([200.0, 1000.0, 2000.0]) + 500 # Here the onset was at 1s, 500 ms past the newer simulations
    ax1.set_xticks(xticks)
    ax1.set_xticklabels((xticks - 500)/1000.0)
    ax1.set_yticks([-50.0, -40.0])
for axis in axes:
    [sp.set_visible(False) for sp in axis.spines.values()]
    axis.tick_params(right=False, top=False)
for axis in axes[:-1]:
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)

fig.subplots_adjust(hspace=0.01)    
fig.set_frameon(False)
fig.set_size_inches(5.5/2.54, 5.5/2.54)
fig.savefig('pn_kc_ggn_noshift_weak_inh.svg')
plt.show()
myplot.plot_kc_spike_count_hist(fname)
plt.show()

# 
# fig_flat_pn_excitation.py ends here
