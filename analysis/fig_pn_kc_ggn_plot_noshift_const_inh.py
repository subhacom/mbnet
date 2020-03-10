# fig_pn_kc_ggn_plot_noshift_const_inh.py --- 
# Author: Subhasis  Ray
# Created: Thu Dec 20 16:19:02 2018 (-0500)
# Last-Updated: Mon Mar 18 18:53:54 2019 (-0400)
#           By: Subhasis  Ray
# Version: $Id$

# Code:
"""Plot an example simulation where the PN odor response is not shifting and GGN->KC gmax is constant."""

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
import yaml
import pint
from collections import defaultdict
import pandas as pd
import network_data_analysis as nda
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import pn_kc_ggn_plot_mpl as myplot
import neurograph as ng
import timeit
from sklearn import cluster, preprocessing

plt.rc('font', size=11)

_ur = pint.UnitRegistry()
Q_ = _ur.Quantity
datadir = 'D:\\biowulf_stage\\olfactory_network\\'
datadir = '/data/rays3/ggn/olfactory_network/'

# jid = '16034794'
jid = '22295165'
fname = nda.find_h5_file(jid, datadir)

gs = gridspec.GridSpec(nrows=3, ncols=1, height_ratios=[2, 2, 1], hspace=0.05)
fig = plt.figure()
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax0)
ax2 = fig.add_subplot(gs[2], sharex=ax0)

axes = [ax0, ax1, ax2]


with h5.File(fname, 'r') as fd:
    print(jid, len(nda.get_spiking_kcs(fd)))
    config = nda.load_config(fd)
    print(yaml.dump(config, default_style=''))
    pn_st = []
    pn_id = []
    for pn in fd[nda.pn_st_path].values():
        pn_st.append(pn.value)
        pn_id.append([int(pn.name.rpartition('_')[-1])] * len(pn))    
    ax0.plot(np.concatenate(pn_st[::10]), np.concatenate(pn_id[::10]), marker='s', ms=1, color='#fdb863', ls='none')
    kc_x, kc_y = nda.get_event_times(fd[nda.kc_st_path])
    ax1.plot(np.concatenate(kc_x[::10]), np.concatenate(kc_y[::10]), color='#e66101', marker='s', ms=1, ls='none')
    myplot.plot_ggn_vm(axes[2], fd, fd['/data/uniform/ggn_basal/GGN_basal_Vm'], 'dend_b', 1, color='#009292', alpha=1.0)
    ax2.set_ylim(-53, -45)
    ax2.set_xlim(200, 2000)
    xticks = np.array([200.0, 1000.0, 2000.0])
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticks/1000.0)
    ax2.set_yticks([-50.0, -45.0])
    ax2.hlines(y=-53.0, xmin=Q_(config['stimulus']['onset']).to('ms').m, xmax=Q_(config['stimulus']['onset']).to('ms').m + Q_(config['stimulus']['duration']).to('ms').m, lw=3, color='lightgray')
for axis in axes:
    [sp.set_visible(False) for sp in axis.spines.values()]
    axis.tick_params(right=False, top=False)
for axis in axes[:-1]:
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
plotted_pns = sum([len(st) >= 0 for st in pn_st[::10]])
plotted_kcs = sum([len(st) >= 0 for st in kc_x[::10]])
print(f'plotted {plotted_pns} PNs,  {plotted_kcs} KCs')

fig.subplots_adjust(hspace=0.01)    
fig.set_frameon(False)
fig.set_size_inches(5.5/2.54, 10.0/2.54)
fig.savefig('pn_kc_ggn_noshift_const_inh.svg')
plt.show()



# 
# fig_pn_kc_ggn_plot_noshift_const_inh.py ends here
