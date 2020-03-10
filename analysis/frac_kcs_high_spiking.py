# frac_kcs_high_spiking.py --- 
# Author: Subhasis  Ray
# Created: Wed Mar 13 16:26:03 2019 (-0400)
# Last-Updated: Wed Mar 13 16:55:15 2019 (-0400)
#           By: Subhasis  Ray
# Version: $Id$

# Code:
"""Check the fraction of KCs that spiked more than 2 times
"""

from __future__ import print_function
import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

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
import itertools as it
from matplotlib import pyplot as plt
import yaml
import pint
from collections import defaultdict
import pandas as pd
import networkx as nx
import pygraphviz as pgv

import network_data_analysis as nda
from matplotlib.backends.backend_pdf import PdfPages
import pn_kc_ggn_plot_mpl as myplot
import neurograph as ng
import timeit
from sklearn import cluster, preprocessing

datadir = '/data/rays3/ggn/olfactory_network/'


siminfo = pd.read_csv(StringIO("""jid,pn_kc_gmax,shifting,lognorm,spiking_kcs
22072442,2.5pS,y,all,6627
22087964,2.5pS,y,all,6497
22087965,2.5pS,y,all,7933
22087966,2.5pS,y,all,6833
22087967,2.5pS,y,all,7591
22087969,2.4pS,y,all,
22087970,2.4pS,y,all,
22087971,2.4pS,y,all,
22087972,2.4pS,y,all,
22087973,2.4pS,y,all,
"""))

for row in siminfo.itertuples():
    filename = nda.find_h5_file(row.jid, datadir)
    # fig, ax = plt.subplots()
    with h5.File(filename, 'r') as fd:
        config = yaml.load(fd.attrs['config'])
        
        # ax0 = fig.add_subplot(rows, cols, nax)
        spiking_kcs = 0
        hyper_kcs = 0
        spike_counts = []
        for kc, st in fd['/data/event/kc/kc_spiketime'].items():
            spike_counts.append(len(st))
            if len(st) > 0:
                spiking_kcs += 1
            if len(st) > 3:
                hyper_kcs += 1
        print('=' * 20)
        print(row)
        print(row.Index, row.jid, config['pn']['shifting'], config['pn_kc_syn']['gmax'], config['pn_kc_syn']['std'],
              config['kc_ggn_alphaL_syn']['gmax'], config['kc_ggn_alphaL_syn']['std'],
              config['ggn_kc_syn']['gmax'], config['ggn_kc_syn']['std'], spiking_kcs, hyper_kcs, hyper_kcs * 1.0 / spiking_kcs)        
        siminfo.loc[row.Index, 'hyper'] = hyper_kcs
        if np.isnan(siminfo.loc[row.Index, 'spiking_kcs']):
            siminfo.loc[row.Index, 'spiking_kcs'] = spiking_kcs
        print('-' * 20)
        bins = np.arange(0, max(spike_counts)+1)
        # ax.hist(spike_counts, bins=bins)
        # ax.axvline(x=bins[-1], color='red')
        # ax.set_title(f'{row.jid} spiking:{spiking_kcs} hyper:{hyper_kcs}')
# plt.show()




# 
# frac_kcs_high_spiking.py ends here
