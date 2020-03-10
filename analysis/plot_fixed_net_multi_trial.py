# fixed_net_multi_trial_plot.py --- 
# Author: Subhasis  Ray
# Created: Thu Sep 27 16:01:45 2018 (-0400)
# Last-Updated: Fri Sep 28 12:35:14 2018 (-0400)
#           By: Subhasis  Ray
# Version: $Id$

# Code:
"""This is a refinement of analysis_20180927.py
- for plotting and comparing multiple trials of same network"""

from __future__ import print_function
import sys
if sys.platform == 'win32':
    sys.path.append('D:/subhasis_ggn/model/analysis')
else:
    sys.path += ['/home/rays3/projects/ggn/analysis', '/home/rays3/projects/ggn/morphutils', '/home/rays3/projects/ggn/nrn', '/home/rays3/apps/matrix-entropy']
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
import pn_kc_ggn_plot_mpl as myplot
import neurograph as ng
import timeit
import time
import sklearn.neighbors as skn

import neighborhood_functions as nfn
import moving_window_filter as mwf
import calculate_profile as calcprof

datadir = '/data/rays3/ggn/fixed_net/'
kde_bandwidth = 20.0


def get_kc_clusters(fd):
    lca_clus_info = fd[nda.lca_kc_cluster_path]['label', 'sec'][:, 0]
    mca_clus_info = fd[nda.mca_kc_cluster_path]['label', 'sec'][:, 0]                
    lca_kcs = pd.DataFrame(data=lca_clus_info)
    mca_kcs = pd.DataFrame(data=mca_clus_info)
    mca_kcs['cluster'] = map('mca_{}'.format, mca_kcs.label)
    lca_kcs['cluster'] = map('lca_{}'.format, lca_kcs.label)
    kc_clusters = pd.concat([lca_kcs, mca_kcs], ignore_index=True)
    kc_clusters.set_index('sec', drop=True, verify_integrity=True)
    return kc_clusters


def plot_fixed_net(fd, kc_clusters, time, ax=None):
    ret = {}
    fig = None
    if ax is None:
        fig, ax = plt.subplots(nrows=5, sharex='all')
        fig.suptitle(os.path.basename(fd.filename))
    pn_st = [(int(k.split('_')[-1]), v.value) for k, v in fd[nda.pn_st_path].items()]
    pn_st = sorted(pn_st, key=lambda x: x[0])
    pn_x = []
    pn_y = []
    for k, st in pn_st:
        pn_x.append(st)
        pn_y.append(k * np.ones_like(st))
    ax[0].plot(np.concatenate(pn_x), np.concatenate(pn_y), ',')
    # It is more efficient to plot all the spike times in one plot command
    xgrid = np.arange(time[0], time[-1], 10.0)
    # Compute PSTH with Gaussian kernel
    pn_psth, _ = np.histogram(np.concatenate(pn_x), xgrid)
    pn_pdf, _ = nda.kde_pdf(np.concatenate(pn_x), kde_bandwidth, xgrid=xgrid)
    ax[1].plot(xgrid[:-1], pn_psth)
    ax[1].twinx().plot(xgrid, pn_pdf)
    kc_st = {v.attrs['source']: v.value for k, v in fd[nda.kc_st_path].items()}
    ggn_basal = fd['/data/uniform/ggn_basal/GGN_basal_Vm'].value
    kc_idx = 0
    all_kc_x = []
    all_kc_y = []
    for clus, grp in kc_clusters.groupby('cluster'):
        kc_x = []
        kc_y = []
        for idx, row in grp.iterrows():
            st = kc_st[row['sec']]
            kc_x.append(st)
            kc_y.append(np.ones_like(st) * kc_idx)
            kc_idx += 1
        all_kc_x.append(kc_x)
        all_kc_y.append(kc_y)
        ax[2].plot(np.concatenate(kc_x), np.concatenate(kc_y), ',')
    kc_x = np.concatenate([np.concatenate(x) for x in all_kc_x])
    kc_psth, _ = np.histogram(kc_x, xgrid) 
    kc_pdf, _ = nda.kde_pdf(kc_x, kde_bandwidth, xgrid=xgrid)
    ax[3].plot(xgrid[:-1], kc_psth)
    ax[3].twinx().plot(xgrid, kc_pdf)
    ax[4].plot(time, np.mean(ggn_basal, axis=0))
    ax[4].set_ylim(-51, -30)
    ret['pn_x'] = pn_x
    ret['pn_y'] = pn_y
    ret['kc_x'] = all_kc_x
    ret['kc_y'] = all_kc_y
    ret['pn_psth'] = pn_psth
    ret['pn_pdf'] = pn_pdf
    ret['kc_psth'] = kc_psth
    ret['kc_pdf'] = kc_pdf
    ret['xgrid'] = xgrid
    ret['ggn_basal'] = ggn_basal
    return fig, ax, ret

def plot_multi_trial(filenames):
    kc_clusters = None
    orig_data = None
    for fname in filenames:
        print(fname)
        try:
            with h5.File(fname, 'r') as fd:
                if kc_clusters is None:
                    with h5.File(fd.attrs['original'], 'r') as fd_orig:
                        kc_clusters = get_kc_clusters(fd_orig)
                        vuniform = fd_orig['/data/uniform/ggn_basal/GGN_basal_Vm']
                        orig_time = np.arange(vuniform.shape[1]) * vuniform.attrs['dt']
                        ofig, oax, orig_data = plot_fixed_net(fd_orig, kc_clusters, orig_time)
                time = fd['/map/time/time'].value
                fig, ax = plt.subplots(nrows=5, sharex='all')
                fig.suptitle(fd.filename)
                ax[1].plot(orig_data['xgrid'][:-1], orig_data['pn_psth'], '--')
                ax[1].twinx().plot(orig_data['xgrid'], orig_data['pn_pdf'], '--')
                ax[3].plot(orig_data['xgrid'][:-1], orig_data['kc_psth'], '--')
                ax[3].twinx().plot(orig_data['xgrid'], orig_data['kc_pdf'], '--')
                ax[4].plot(orig_time, np.mean(orig_data['ggn_basal'], axis=0), '--')
                fig, ax, ret = plot_fixed_net(fd, kc_clusters, time, ax=ax)
        except IOError as e:
            print('Problem with file', fname)
            print(e)


jids_iid_same_odor = [
    '10126599',
    '10126600',
    '10126601',
    '10126602',
    '10126603',
    '10126604',
    '10126605',
    '10126606',
    '10126607',
    '10126609']

filenames = [nda.find_file_by_jid(jid, datadir)[0] for jid in jids_iid_same_odor]            
plot_multi_trial(filenames)
plt.show()

jids_clustered_same_odor = [
    '10126636',
    '10126637',
    '10126641',
    '10126643',
    '10126646',
    '10126648',
    '10126651',
    '10126653',
    '10126655',
    '10126657']

filenames = [nda.find_file_by_jid(jid, datadir)[0] for jid in jids_clustered_same_odor]            
plot_multi_trial(filenames)
plt.show()

jids_iid_diff_odor = [  # these were after JID9932209 - IID PN->KC conn
    '10152882',
    '10152888',
    '10152892',
    '10152898',
    '10152901',
    '10152906',
    '10152909',
    '10152913',
    '10152917',
    '10152920']

jids_clustered_diff_odor = [ # after JID9932198, clustered PN->KC conn
    '10153056',
    '10153059',
    '10153063',
    '10153066',
    '10153069',
    '10153070',
    '10153073',
    '10153075',
    '10153077',
    '10153087']

for jids in zip(jids_iid_same_odor, jids_iid_diff_odor):
    filenames = [nda.find_file_by_jid(jid, datadir)[0] for jid in jids]
    plot_multi_trial(filenames)
    break
for jids in zip(jids_clustered_same_odor,
                jids_clustered_diff_odor):
    filenames = [nda.find_file_by_jid(jid, datadir)[0] for jid in jids]
    plot_multi_trial(filenames)
    break
plt.show()




# 
# fixed_net_multi_trial_plot.py ends here
