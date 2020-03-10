# pn_kc_ggn_plot.py --- 
# 
# Filename: pn_kc_ggn_plot.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Fri Feb 16 13:08:41 2018 (-0500)
# Last-Updated: Mon Mar  4 11:07:59 2019 (-0500)
#           By: Subhasis  Ray
#     Update #: 1108
# 
# Code:
from __future__ import print_function
import warnings
import sys
import os
from timeit import default_timer as timer
import numpy as np
import random
import h5py as h5
import yaml
import pandas as pd
from matplotlib import pyplot as plt
from pint import UnitRegistry

import sys

import neurograph as ng
from sklearn.neighbors import KernelDensity
import network_data_analysis as nda

ur = UnitRegistry()
Q_ = ur.Quantity


def despine(ax, locs='all'):
    """Remove the line marking axing bounder in locations specified in the list `locs`.
    The locations can be 'top', 'bottom', 'left', 'right' or 'all'
    """
    if locs == 'all':
        locs = ['top', 'bottom', 'left', 'right']
    for loc in locs:
        ax.spines[loc].set_visible(False)

def plot_population_psth(ax, spike_trains, ncell, bins, alpha=0.5, rate_sym='b-', cell_sym='r--'):
    """Plot the population PSTH on axis `ax`, where `spike_trains` is a
    list of arrays, each containing the spike times of one cell. ncell
    is the total number of cell and is used for computing rates per
    cell.

    """
    start = timer()
    cell_counts = np.zeros(len(bins) - 1)
    spike_counts = np.zeros(len(bins) - 1)
    spiking_cell_count = 0
    for train in spike_trains:
        if len(train) > 0:
            spiking_cell_count += 1
        hist, bins = np.histogram(train, bins)
        cell_counts += (hist > 0)
        spike_counts += hist
    print('Total number of spiking cells', spiking_cell_count, 'out of', len(spike_trains))
    spike_rate = spike_counts * 1e3 / (ncell * (bins[1] - bins[0]))
    cell_frac = cell_counts * 1e3 / (ncell * (bins[1] - bins[0]))
    if rate_sym is None:
        ax.plot((bins[:-1] + bins[1:]) / 2.0, spike_rate, 
                label='spikes/ncells/binwidth(s)', alpha=alpha)
    else:
        ax.plot((bins[:-1] + bins[1:]) / 2.0, spike_rate, 
                rate_sym, label='spikes/ncells/binwidth(s)', alpha=alpha)
    if cell_sym is None:
        ax.plot((bins[:-1] + bins[1:]) / 2.0, cell_frac, 
            label='cells spiking/ncells/binwidth(s)', alpha=alpha)
    else:        
        ax.plot((bins[:-1] + bins[1:]) / 2.0, cell_frac, cell_sym,
            label='cells spiking/ncells/binwidth(s)', alpha=alpha)
    end = timer()
    print('Plotted PSTH for {} spike trains in {} s'.format(
        len(spike_trains),
        end - start))
    return ax, spike_rate, cell_frac


def plot_population_KDE(ax, spike_trains, xgrid, bandwidth, color='k',
                        maxamp=1.0):
    """Plot kernel density estimation of spike times of all the spike
    trains.

    Scale the PDF so that the maximum amplitude is same as maxamp
    (useful for bringing it to scale with histogram for example).

    """
#    print(spike_trains)
    spike_times = np.concatenate(spike_trains)
    if len(spike_times) == 0:
        warnings.warn('No spikes in spike trains')
        return
    pdf, xgrid = nda.kde_pdf(spike_times, bandwidth=bandwidth, xgrid=xgrid)
    scale = maxamp / max(pdf) if max(pdf) > 0 else 1.0
    ax.plot(xgrid, pdf * scale / len(spike_trains), color=color, label='KDE')
    return pdf, xgrid


def plot_kc_spikes_by_cluster(ax, fd, ca_side, color='k', marker=',', hline=False, subsample_clus=0, subsample_kc=0):
    """Raster plot KC spike times ordered by cluster no.

    ca_side should be 'lca' or 'mca'.

    subsample_clus: number of clusters to plot if subsampling (> 0).

    subsample_kc:  number of kcs in each cluster to plot if subsampling (> 0).
    """
    kc_spikes = {row['source'].decode(): row['data']
                 for row in fd['/map/event/kc/kc_spiketime'].value}

    cluster_info = nda.extract_cluster_info(fd, ca_side)
    ii = 0
    spike_x, spike_y = [], []
    ret = []
    labels = list(set(cluster_info['label']))
    if subsample_clus > 0:
        labels = np.random.choice(labels, size=subsample_clus, replace=False)
    labelgrp = cluster_info.groupby('label')
    for label in labels:
        group = labelgrp.get_group(label)
        grpx, grpy = [], []
        kcs = np.char.decode(group.index.values.flatten().astype('S'))
        if subsample_kc > 0:
            kcs = np.random.choice(kcs, size=subsample_kc, replace=False)
        for kc in kcs:
            st = fd[kc_spikes[kc]].value
            grpx.append(st)
            grpy.append(np.ones(len(st)) + ii)
            ii += 1
        ret.append(ax.plot(np.concatenate(grpx), np.concatenate(grpy),
                           marker=marker, linestyle='none'))
        spike_x += grpx
        spike_y += grpy
        if hline:
            ax.axhline(y=ii-0.5, color='gray', linewidth=1.0)
    return ret, spike_x, spike_y


def plot_kc_spikes(ax, fd, ca_side='both', color='k', marker=','):
    """Raster plot KC spike times for KCs belonging to the specified side
    of calyx ('lca' or 'mca').

    This function does not care about spatial clusters.

    Returns: the line object, list of spike times and list of their
    y-positions.

    """
    if ca_side == 'both':
        nodes = fd[nda.kc_st_path].keys()
    else:
        nodes = nda.get_kc_spike_nodes_by_region(fd, ca_side)
    spike_x, spike_y = [], []
    fname = fd.filename
    try:
        spike_x, spike_y = nda.get_event_times(
            fd['/data/event/kc/kc_spiketime'],
            nodes=nodes)
    except KeyError:
        dirname = os.path.dirname(fname)
        fname = 'kc_spikes_' + os.path.basename(fd.filename)
        with h5.File(os.path.join(dirname, fname)) as kc_file:
            spike_x, spike_y = nda.get_event_times(kc_file, nodes=nodes)
    if len(spike_x) > 0:
        ret = ax.plot(np.concatenate(spike_x), np.concatenate(spike_y),
                      color=color, marker=marker, linestyle='none')
    else:
        ret = None
    return ret, spike_x, spike_y


def plot_kc_vm(ax, fd, region, count, color='k', alpha=0.5):
    """Plot Vm of `count` KCs from sepcified region."""
    match = nda.get_kc_vm_idx_by_region(fd, region)
    if len(match) == 0:
        return [], []
    selected = random.sample(match, min(count, len(match)))
    kc_vm_node = fd['/data/uniform/kc/KC_Vm']
    try:
        t = np.arange(kc_vm_node.shape[1]) * kc_vm_node.attrs['dt']
    except KeyError:
        t = np.arange(kc_vm_node.shape[1])
    for name, idx in selected:
        ax.plot(t, kc_vm_node[idx, :], label=name, color=color, alpha=alpha)
    return selected


def plot_ggn_vm(ax, fd, dataset, region=None, count=5, color='k', alpha=0.5):
    """Plot Vm of GGN from dataset"""
    sec_list = [sec.decode('utf-8')  for sec in dataset.dims[0]['source']]
    match = []
    if region is None:
        match = [(sec, ii) for ii, sec in enumerate(sec_list)]
    else:
        rsid = ng.name_sid[region]
        for ii, sec in enumerate(sec_list):
            sid = nda.ggn_sec_to_sid(sec)
            if sid == rsid:
                match.append((sec, ii))
    # print(match)
    if len(match) == 0:
        return [], []
    selected = random.sample(match, min(len(match), count))
    try:
        t = np.arange(dataset.shape[1]) * dataset.attrs['dt']
    except KeyError:
        t = np.arange(dataset.shape[1])
    for sec, ii in selected:
        ax.plot(t, dataset[ii, :], label=sec, color=color, alpha=alpha)
    return selected
    
    
def plot_spike_rasters(fname, vm_samples=10, psth_bin_width=50.0,
                       kde_bw=50.0, by_cluster=False):
    """The file `fname` has data from pn_kc_ggn simulation. In the early
    ones I did not record the spike times for KCs. binwidths are in
    ms.

    """
    start = timer()
    print('Processing', fname)
    with h5.File(fname, 'r') as fd:
        try:
            print('Description:', fd.attrs['description'])
        except KeyError:
            print('No description available')
        try:
            config = yaml.load(fd.attrs['config'].decode())
        except KeyError:
            try:
                config = yaml.load(fd['model/filecontents/mb/network/config.yaml'][0].decode())
            except KeyError: # possibly fixed network - look for config in template
                try:
                    original = fd.attrs['original'].decode()
                    with h5.File(original, 'r') as origfd:
                        try:
                            config = yaml.load(origfd.attrs['config'].decode())
                        except AttributeError:  # Handling in Python3?
                            config = yaml.load(origfd.attrs['config'])
                except KeyError:
                    print('No config attribute or model file')
                    pass
        # PN spike raster
        pn_st = fd['/data/event/pn/pn_spiketime']
        fig, axes = plt.subplots(nrows=6, ncols=1, sharex=True)
        try:
            if 'calyx' in fd.attrs['description'].decode():
                descr = 'KC->GGN in alphaL + CA'
            else:
                descr = 'KC->GGN in alphaL only'
        except KeyError:
            descr = ''
        
        fig.suptitle('{} {}'.format(os.path.basename(fname), descr))
        ax_pn_spike_raster = axes[0]
        print('Plotting PN spikes')
        ax_pn_spike_raster.set_title('PN spike raster')
        nodes = [int(node.split('_')[-1]) for node in pn_st]
        nodes = ['pn_{}'.format(node) for node in sorted(nodes)]        
        spike_x, spike_y = nda.get_event_times(pn_st, nodes)
        ax_pn_spike_raster.plot(np.concatenate(spike_x), np.concatenate(spike_y), 'k,')
        simtime = Q_(config['stimulus']['onset']).to('ms').m +  \
                  Q_(config['stimulus']['duration']).to('ms').m + \
                  Q_(config['stimulus']['tail']).to('ms').m
        psth_bins = np.arange(0, simtime, psth_bin_width)
        kde_grid = np.linspace(0, simtime, 100.0)
        ax_pn_psth = ax_pn_spike_raster.twinx()
        _, sr, cf = plot_population_psth(ax_pn_psth, spike_x,
                                         config['pn']['number'], psth_bins)
        plot_population_KDE(ax_pn_psth, spike_x, kde_grid, kde_bw, color='y',
                            maxamp=max(sr))
        ax_kc_lca_spike_raster = axes[1]
        print('Plotting KC PSTH in LCA')
        ax_kc_lca_spike_raster.set_title('KC LCA')
        if by_cluster:
            _, lca_spike_x, lca_spike_y = plot_kc_spikes_by_cluster(ax_kc_lca_spike_raster, fd, 'LCA')
        else:
            _, lca_spike_x, lca_spike_y = plot_kc_spikes(ax_kc_lca_spike_raster, fd, 'LCA')
        ax_kc_lca_psth = ax_kc_lca_spike_raster.twinx()
        _, sr, cf = plot_population_psth(ax_kc_lca_psth, lca_spike_x, len(lca_spike_x), psth_bins)
        plot_population_KDE(ax_kc_lca_psth, lca_spike_x, kde_grid, kde_bw, color='y', maxamp=max(sr))
        ax_kc_lca_psth.legend()
        ax_kc_mca_spike_raster = axes[2]
        print('Plotting KC PSTH in MCA')
        ax_kc_mca_spike_raster.set_title('KC MCA')
        if by_cluster:
            _, mca_spike_x, mca_spike_y = plot_kc_spikes_by_cluster(ax_kc_mca_spike_raster, fd, 'MCA')
        else:
            _, mca_spike_x, mca_spike_y = plot_kc_spikes(ax_kc_mca_spike_raster, fd, 'MCA')
        ax_kc_mca_psth = ax_kc_mca_spike_raster.twinx()
        _, sr, cf = plot_population_psth(ax_kc_mca_psth, mca_spike_x, len(mca_spike_x), psth_bins)
        plot_population_KDE(ax_kc_mca_psth, mca_spike_x, kde_grid, kde_bw, color='y', maxamp=max(sr))
        ax_kc_mca_psth.legend()
        # LCA KC Vm
        ax_kc_lca_vm = axes[3]
        ax_kc_lca_vm.set_title('KC LCA')
        plot_kc_vm(ax_kc_lca_vm, fd, 'LCA', vm_samples)
        # MCA KC Vm
        ax_kc_mca_vm = axes[4]
        ax_kc_mca_vm.set_title('KC MCA')
        plot_kc_vm(ax_kc_mca_vm, fd, 'MCA', vm_samples)
        # GGN MCA Vm, GGN LCA Vm
        ggn_vm_plot = axes[5]
        ggn_vm_plot.set_title('GGN Vm')
        ggn_output_vm = fd['/data/uniform/ggn_output/GGN_output_Vm']
        plot_ggn_vm(ggn_vm_plot, fd, ggn_output_vm, 'LCA', vm_samples, color='r')
        plot_ggn_vm(ggn_vm_plot, fd, ggn_output_vm, 'MCA', vm_samples, color='b')
        lca, = ggn_vm_plot.plot([], color='r', label='LCA')
        mca, = ggn_vm_plot.plot([], color='b', label='MCA')
        alpha, = ggn_vm_plot.plot([], color='k', label='alphaL')
        basal, = ggn_vm_plot.plot([], color='g', label='basal')
        ggn_vm_plot.legend(handles=[lca, mca, alpha, basal])
        # GGN alphaL Vm
        ggn_alphaL_vm = fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm']
        plot_ggn_vm(ggn_vm_plot, fd, ggn_alphaL_vm, 'alphaL', vm_samples, color='k')
        try:
            ggn_basal_vm = fd['/data/uniform/ggn_basal/GGN_basal_Vm']
            plot_ggn_vm(ggn_vm_plot, fd, ggn_basal_vm, 'basal', vm_samples, color='g')
        except KeyError:
            warnings.warn('No basal Vm recorded from GGN')
        end = timer()
        print('Time for plotting {}s'.format(end - start))
        return fig, axes



def plot_data_by_jobid(jid, datadir, save=False, by_cluster=False, figdir='figures'):
    flist = os.listdir(datadir)
    match = [fname for fname in flist if 'JID{}'.format(jid) in fname]
    if len(match) > 1:
        print('Two files with same jobid.', match)
    for fname in match:
        if fname.endswith('.h5'):
            break
    fig, ax = plot_spike_rasters(os.path.join(datadir, fname), by_cluster=by_cluster)
    if save:
        figfile = os.path.join(figdir, fname.rpartition('.h5')[0] + '.png')
        fig.savefig(figfile)
        print('Saved figure in', figfile)
    return fig, ax


def plot_kc_spike_count_hist(fname, bins=None, save=False, figdir='figures'):
    """Plot histogram of spike counts recorded as 1D datasets under group
    path in file fname"""
    with h5.File(fname, 'r') as fd:
        kc_st_grp = fd['/data/event/kc/kc_spiketime']
        lca_kcs = nda.get_kc_spike_nodes_by_region(fd, 'LCA')
        lca_spike_count = [len(kc_st_grp[kc]) for kc in lca_kcs] 
        mca_kcs = nda.get_kc_spike_nodes_by_region(fd, 'mca')
        mca_spike_count = [len(kc_st_grp[kc]) for kc in mca_kcs]        
        fig, axes = plt.subplots(nrows=2, ncols=1, sharey='all', sharex='all')
        if bins is None:
            bins = np.arange(max(lca_spike_count + mca_spike_count) + 1)
        hist, bins, patches = axes[0].hist(lca_spike_count, bins=bins)
        axes[0].arrow(max(lca_spike_count), max(hist)/2, 0, -max(hist)/2.0,
                      head_width=0.5, head_length=max(hist)/20.0,
                      length_includes_head=True)
        axes[0].set_title('LCA')
        hist, bins, patches = axes[1].hist(mca_spike_count, bins=bins)
        axes[1].arrow(max(mca_spike_count), max(hist)/2, 0, -max(hist)/2,
                      head_width=0.5, head_length=max(hist)/20.0,
                      length_includes_head=True)
        axes[1].set_title('MCA')
        fname = os.path.basename(fname)
        fig.suptitle(os.path.basename(fname))
        fname = fname.rpartition('.h5')[0]
        if save:
            fname = os.path.join(figdir, fname) + '_kc_spikecount_hist.png'
            fig.tight_layout()
            fig.savefig(fname)
            print('Saved spike count histogram of KCs in', fname)
    return fig, axes


def plot_kc_clusters_with_presynaptic_pn(fpath, region, clusters):
    """Plot the spike rasters for KCs by cluster labels in clusters and
    the spike rasters for presynaptic PNs"""
    lca_clusters = nda.get_kc_clusters(fpath, region)
    print(lca_clusters.shape)
    lca_clusters = pd.DataFrame(data=lca_clusters[:, 0])
    lca_grp = lca_clusters.groupby('label')
    sec_st = {}
    pn_st = {}
    with h5.File(fpath, 'r') as fd:
        sec_st_path = nda.get_kc_event_node_map(fd)
        for kc, path in sec_st_path.items():
            sec_st[kc] = fd[path].value
        pn_kc_syn = fd[nda.pn_kc_syn_path]['pre', 'post']
        for node in fd['/data/event/pn/pn_spiketime'].values():        
            pn_st[node.attrs['source'].decode('utf-8')] = node.value

    keys = np.random.choice(lca_grp.groups.keys(), size=5, replace=False)
    pn_kc_syn = pd.DataFrame(pn_kc_syn[:, 0])            
    kc_pre_pn = pn_kc_syn.groupby('post')

    fig, axes = plt.subplots(nrows=len(keys), ncols=2, sharex='all')
    prev = set()    
    for ii, k in enumerate(keys):
        seclist = lca_grp.get_group(k)['sec']
        cluster_pre = set()
        for jj, sec in enumerate(seclist):
            st = sec_st[sec]
            axes[ii, 1].plot(st, jj * np.ones(len(st)), 'k,')
            current = set(kc_pre_pn.get_group(sec)['pre'])
            cluster_pre.update(current)
            print(len(current.intersection(prev)))
            prev = current
        for jj, pre in enumerate(cluster_pre):
            axes[ii, 0].plot(pn_st[pre], jj * np.ones(len(pn_st[pre])), 'k,')
        print('#################')
    return fig, axes


def compare_data(leftfiles, rightfiles, leftheader, rightheader):
    """Compare two simulations side by side"""
    figs = []
    axeslist = []
    psthaxlist = []
    for left, right in zip(leftfiles, rightfiles):
        fig, axes = plt.subplots(nrows=6, ncols=2, sharey='row')
        psth_axes = []
        
        for ii, fname in enumerate([left, right]):
            fpath = os.path.join(datadir, fname)                         
            with h5.File(fpath, 'r') as fd:
                config = nda.load_config(fd)
                bins = np.arange(0, nda.get_simtime(fd)+0.5, 50.0)
                try:
                    pns = list(fd[nda.pn_st_path].keys())
                except KeyError:
                    print('Could not find PNs in', fname)
                    return figs, axeslist, psthaxlist
                pns = sorted(pns, key=lambda x: int(x.split('_')[-1]))
                pn_st, pn_y = nda.get_event_times(fd[nda.pn_st_path], pns)
                axes[0, ii].plot(np.concatenate(pn_st), np.concatenate(pn_y), ',')
                psth_ax = axes[0, ii].twinx()
                psth_axes.append(psth_ax)
                plot_population_psth(psth_ax, pn_st, config['pn']['number'], bins)
                lines, kc_st, kc_y = plot_kc_spikes_by_cluster(axes[1, ii], fd, 'LCA')
                plot_population_psth(axes[2, ii], kc_st, len(kc_st), bins, rate_sym='b^', cell_sym='rv')
                stiminfo = nda.get_stimtime(fd)
                stimend = stiminfo['onset'] + stiminfo['duration'] + stiminfo['offdur']
                rates = [len(st[(st > stiminfo['onset']) & (st < stimend)]) * 1e3 
                         / (stimend - stiminfo['onset']) for st in kc_st]
                print(rates[:5])
                axes[3, ii].hist(rates, bins=np.arange(21))
                axes[3, ii].set_xlabel('Firing rate')
                plot_kc_vm(axes[4, ii], fd, 'LCA', 5)
                plot_ggn_vm(axes[5, ii], fd,
                                   fd['/data/uniform/ggn_output/GGN_output_Vm'],
                                   'LCA', 5, color='r')
                plot_ggn_vm(axes[5, ii], fd,
                                   fd['/data/uniform/ggn_basal/GGN_basal_Vm'],
                                   'basal', 5, color='g')
                axes[5, ii].set_ylim((-53, -35))
                axes[0, ii].set_title('{}\nFAKE? {}'.format(fname, nda.load_config(fd)['kc']['fake_clusters']))
        time_axes = [axes[ii, jj] for ii in [0, 1, 2, 4, 5] for jj in [0, 1]]
        for ax in time_axes[:-1]:
            ax.set_xticks([])
        axes[0, 0].get_shared_x_axes().join(*time_axes)
        axes[2, 0].get_shared_x_axes().join(*axes[2, :])
        # psth_axes[0].get_shared_y_axes().join(*psth_axes)
        psth_axes[0].autoscale()
        # axes[-1, -1].autoscale()
        fig.text(0.1, 0.95, leftheader, ha='left', va='bottom')
        fig.text(0.6, 0.95, rightheader, ha='left', va='bottom')
        fig.set_size_inches(15, 10)
        # fig.tight_layout()
        figs.append(fig)
        axeslist.append(axes)
        psthaxlist.append(psth_axes)
    return figs, axeslist, psthaxlist


datadir = '/data/rays3/ggn/olfactory_network'
files = [
    'pn_kc_ggn_UTC2018_02_15__14_47_23-PID127486-JID61581731.h5',
    'pn_kc_ggn_UTC2018_02_15__14_47_24-PID24307-JID61581730.h5',
    'pn_kc_ggn_UTC2018_02_16__23_55_29-PID78196-JID61754666.h5',
    'pn_kc_ggn_UTC2018_02_16__23_55_31-PID51661-JID61754888.h5',
    'pn_kc_ggn_UTC2018_02_17__00_01_30-PID113209-JID61756334.h5',
    'pn_kc_ggn_UTC2018_02_17__00_01_30-PID113210-JID61756338.h5'
    ]




import argparse 

def make_parser():
    parser = argparse.ArgumentParser(
        description='Plot data from NSDF simulation file for olfactory network.')
    parser.add_argument('-f', '--file', action='append')
    parser.add_argument('-c', '--cluster', action='store_true')
    parser.add_argument('-b', '--bandwidth', type=float, default=50.0,
                        help='Kernel Density estimation bandwidth')
    parser.add_argument('-s', '--save', action='store_true',
                        help='Save figures in files')
    return parser
                        
        
if __name__ == '__main__':    
    args = make_parser().parse_args()
    for fname in args.file:
        fig, axes = plot_spike_rasters(fname, kde_bw=args.bandwidth, by_cluster=args.cluster)
        fhist, fax = plot_kc_spike_count_hist(fname)
        if args.save:
            factivity = '{}_activity.png'.format(os.path.basename(fname))
            fig.savefig(factivity)
            histfile = '{}_kc_spike_count_hist.png'.format(os.path.basename(fname))
            fhist.savefig(histfile)
            print('Saved figures in {} and {}'.format(factivity, histfile))
    plt.show()
        
# 
# pn_kc_ggn_plot.py ends here
