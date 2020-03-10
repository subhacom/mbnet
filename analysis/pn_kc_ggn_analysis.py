# -*- coding: utf-8 -*-

"""
Created on Tue Feb  6 12:57:22 2018

@author: Subhasis
"""


from __future__ import print_function
import sys
if sys.platform == 'win32':
    sys.path.append('D:/subhasis_ggn/model/analysis')
else:
    sys.path += ['/home/rays3/projects/ggn/analysis', '/home/rays3/projects/ggn/morphutils', '/home/rays3/projects/ggn/nrn']
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

plt.style.use('classic')

_ur = pint.UnitRegistry()
Q_ = _ur.Quantity
datadir = 'D:\\biowulf_stage\\olfactory_network\\'
datadir = '/data/rays3/ggn/olfactory_network/'
spike_count_bins = [0, 1, 2, 5, 10, 20, 40, 80]


def find_file_by_jid(jid, datadir):
    flist = os.listdir(datadir)
    match = [fname for fname in flist if 'JID{}'.format(jid) in fname]
    if len(match) > 1:
        print('Two files with same jobid.', match)
    return [os.path.join(datadir, fname) for fname in match]


def plot_data_by_jobid(jid, datadir, save=False, by_cluster=False, figdir='figures'):
    flist = os.listdir(datadir)
    match = [fname for fname in flist if 'JID{}'.format(jid) in fname]
    if len(match) > 1:
        print('Two files with same jobid.', match)
    fname = match[0]
    fig, ax = myplot.plot_spike_rasters(os.path.join(datadir, fname), by_cluster=by_cluster)
    if save:
        figfile = os.path.join(figdir, fname.rpartition('.h5')[0] + '.png')
        fig.savefig(figfile)
        print('Saved figure in', figfile)
    return fig, ax

def get_config_by_jid(jid, datadir, **kwargs):
    """Example:
        
        get_config_by_jid('64819232', datadir, default_flow_style=False, default_style='' )
        """
    fname = find_file_by_jid(jid, datadir)[0]
    with h5.File(fname, 'r') as fd:
        return yaml.dump(yaml.load(fd['/model/filecontents/mb/network/config.yaml'].value[0]), **kwargs)
        



# #fd = h5.File('pn_kc_ggn_2018_02_06__00_55_42.h5', 'r')
# fd = h5.File('pn_kc_ggn_UTC2018_02_08__23_18_08-PID9812-JID0.h5', 'r')
# kc_vm = fd['/data/uniform/kc/KC_Vm']
# t = np.arange(kc_vm.shape[1]) * kc_vm.attrs['dt']
# fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)
# pn_st = fd['/data/event/pn/pn_spiketime']
# for ii, node in enumerate(pn_st):
#     axes[0].plot(pn_st[node], np.ones(len(pn_st[node]))*ii, 'b|')
# for ii in random.sample(range(kc_vm.shape[0]), 5):
#     axes[1].plot(t, kc_vm[ii,:])
# plt.show()

# ggn_vm = fd['/data/uniform/ggn_output/GGN_output_Vm']
# for ii in random.sample(range(ggn_vm.shape[0]), 5):
#     axes[2].plot(t, ggn_vm[ii,:])
# plt.show()    
# fd.close()

# fd = h5.File('pn_kc_ggn_UTC2018_02_08__23_34_06-PID76972-JID61032510', 'r')
# fd.close()

# ## With KC->GGN synapses both in alpha lobe as well as calyx
# fd = h5.File('D:\\biowulf_stage\\olfactory_network\\pn_kc_ggn_UTC2018_02_08__23_34_06-PID76972-JID61032510.h5', 'r')
# print(fd.attrs['description'])
# kc_vm = fd['/data/uniform/kc/KC_Vm']
# t = np.arange(kc_vm.shape[1]) * kc_vm.attrs['dt']
# fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)
# pn_st = fd['/data/event/pn/pn_spiketime']
# for ii, node in enumerate(pn_st):
#     axes[0].plot(pn_st[node][:]*1000.0, np.ones(len(pn_st[node]))*ii, 'b|')
# for ii in random.sample(range(kc_vm.shape[0]), 5):
#     axes[1].plot(t, kc_vm[ii,:])
# ggn_vm = fd['/data/uniform/ggn_CA_input/GGN_CA_input_Vm']
# for ii in random.sample(range(ggn_vm.shape[0]), 5):
#     axes[2].plot(t, ggn_vm[ii,:], 'g-')
# ggn_vm = fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm']
# for ii in random.sample(range(ggn_vm.shape[0]), 5):
#     axes[2].plot(t, ggn_vm[ii,:], 'b-')
# ggn_vm = fd['/data/uniform/ggn_output/GGN_output_Vm']
# for ii in random.sample(range(ggn_vm.shape[0]), 5):
#     axes[2].plot(t, ggn_vm[ii,:], 'r-')
# plt.show()    
# fd.close()

# fd = h5.File('D:\\biowulf_stage\\olfactory_network\\pn_kc_ggn_UTC2018_02_08__23_34_06-PID76972-JID61032510.h5', 'r')
# kc_vm = fd['/data/uniform/kc/KC_Vm']
# t = np.arange(kc_vm.shape[1]) * kc_vm.attrs['dt']
# fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)
# pn_st = fd['/data/event/pn/pn_spiketime']
# for ii, node in enumerate(pn_st):
#     axes[0].plot(pn_st[node], np.ones(len(pn_st[node]))*ii, 'b|')
# for ii in random.sample(range(kc_vm.shape[0]), 5):
#     axes[1].plot(t, kc_vm[ii,:])
# ggn_vm = fd['/data/uniform/ggn_output/GGN_output_Vm']
# for ii in random.sample(range(ggn_vm.shape[0]), 5):
#     axes[2].plot(t, ggn_vm[ii,:])
# plt.show()    


# fd.close()


def plot_pn_kc_ggn(fname):
    """Plot the PN spike raster, a sample of KC Vms and GGN Vms"""
    with h5.File(fname, 'r') as fd:
        kc_vm = fd['/data/uniform/kc/KC_Vm']
        t = np.arange(kc_vm.shape[1]) * kc_vm.attrs['dt']
        fig, axes = plt.subplots(nrows=5, ncols=1, sharex=True)
        fig.suptitle(fd.attrs['description'])
        config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'][0])
        pn_st = fd['/data/event/pn/pn_spiketime']
        for ii, node in enumerate(pn_st):
            axes[0].plot(pn_st[node], np.ones(len(pn_st[node])) * ii, 'b|')
        axes[0].set_title('PN spike raster')    
        kc_lca = int(kc_vm.shape[0] * config['kc']['lca_frac']+0.5)
        for ii in random.sample(range(kc_lca), 50):
            axes[1].plot(t, kc_vm[ii, :])
        axes[1].set_title('KC LCA')    
        for ii in random.sample(range(kc_vm.shape[0] - kc_lca, kc_vm.shape[0]), 50):
            axes[2].plot(t, kc_vm[ii, :])    
        axes[2].set_title('KC MCA')    
        ggn_out_vm = fd['/data/uniform/ggn_output/GGN_output_Vm']
        for ii in random.sample(range(ggn_out_vm.shape[0]), 5):
            axes[3].plot(t, ggn_out_vm[ii, :])
        axes[3].set_title('GGN output')
        ggn_alphaL_vm = fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm']
        for ii in random.sample(range(ggn_alphaL_vm.shape[0]), 5):
            axes[4].plot(t, ggn_alphaL_vm[ii, :])
        axes[4].set_title('GGN alphaL input')           
        plt.show()    
        
        
        


# plot_pn_kc_ggn('D:\\biowulf_stage\\olfactory_network\\pn_kc_ggn_UTC2018_02_12__23_37_04-PID101025-JID61328421.h5')
# plot_pn_kc_ggn('D:\\biowulf_stage\\olfactory_network\\pn_kc_ggn_UTC2018_02_12__23_37_08-PID129559-JID61328422.h5')

def plot_pn_kc_ggn_spike_raster(fname, threshold=-20.0):
    """Plot PN, KC spike rasters"""
    with h5.File(fname, 'r') as fd:
        fig, axes = plt.subplots(nrows=5, ncols=1, sharex=True)        
        fig.suptitle(fd.attrs['description'])
        kc_vm = fd['/data/uniform/kc/KC_Vm']
        t = np.arange(kc_vm.shape[1]) * kc_vm.attrs['dt']        
        config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'][0])
        pn_st = fd['/data/event/pn/pn_spiketime']
        for ii, node in enumerate(pn_st):
            axes[0].plot(pn_st[node], np.ones(len(pn_st[node]))*ii, 'b|')
        axes[0].set_title('PN spike raster')    
        kc_lca = int(kc_vm.shape[0] * config['kc']['lca_frac']+0.5)
        for ii in range(kc_lca):
            spike_idx = np.flatnonzero((kc_vm[ii, :-1] < threshold) & 
                                       (kc_vm[ii, 1:] > threshold))
            axes[1].plot(t[spike_idx], np.ones(spike_idx.shape) * ii, 'k|')
        axes[1].set_title('KC LCA')    
        for ii in range(kc_vm.shape[0] - kc_lca, kc_vm.shape[0]):
            spike_idx = np.flatnonzero((kc_vm[ii, :-1] < threshold) & 
                                       (kc_vm[ii, 1:] > threshold))
            axes[2].plot(t[spike_idx], np.ones(spike_idx.shape) * ii, 'k|')
        axes[2].set_title('KC MCA')    
        ggn_out_vm = fd['/data/uniform/ggn_output/GGN_output_Vm']
        for ii in random.sample(range(ggn_out_vm.shape[0]), 5):
            axes[3].plot(t, ggn_out_vm[ii,:])
        axes[3].set_title('GGN output')
        ggn_alphaL_vm = fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm']
        for ii in random.sample(range(ggn_alphaL_vm.shape[0]), 5):
            axes[4].plot(t, ggn_alphaL_vm[ii,:])
        axes[4].set_title('GGN alphaL input')           
        plt.show()


        # This is for extracting KC spike times.
def dump_kc_spiketimes():
    """2018-03-07 - making it a function to avoid rerunning by mistake"""
    filenames = [ 'pn_kc_ggn_UTC2018_02_08__23_34_06-PID76972-JID61032510.h5',
                  'pn_kc_ggn_UTC2018_02_08__23_34_08-PID78555-JID61032632.h5',
                  'pn_kc_ggn_UTC2018_02_12__23_37_04-PID101025-JID61328421.h5',
                  'pn_kc_ggn_UTC2018_02_12__23_37_08-PID129559-JID61328422.h5',
                  'pn_kc_ggn_UTC2018_02_13__18_31_37-PID99588-JID61403623.h5',
                  'pn_kc_ggn_UTC2018_02_13__18_31_41-PID41339-JID61403628.h5',]

    threshold = -20
    for fname in filenames:
        fpath = os.path.join(datadir, fname)
        spike_fname = 'kc_spikes_' + fname
        spike_fpath = os.path.join(datadir, spike_fname)
        with h5.File(fpath, 'r') as fd:
            with  h5.File(spike_fpath, 'w') as outfile:
                kc_vm = fd['/data/uniform/kc/KC_Vm']

                for ii in range(kc_vm.shape[0]):
                    st_idx = np.flatnonzero((kc_vm[ii, :-1] < threshold) &
                                            (kc_vm[ii, 1:] >= threshold))
                    st = st_idx * kc_vm.attrs['dt']
                    ds = outfile.create_dataset(str(ii), data=st)
                    ds.attrs['unit'] = 'ms'

                outfile.attrs['original_file'] = fname
            print('Saved', spike_fpath)
    
        
                             
    
# plot_pn_kc_ggn_spike_raster('D:\\biowulf_stage\\olfactory_network\\pn_kc_ggn_UTC2018_02_12__23_37_04-PID101025-JID61328421.h5')
# plot_pn_kc_ggn_spike_raster('D:\\biowulf_stage\\olfactory_network\\pn_kc_ggn_UTC2018_02_12__23_37_08-PID129559-JID61328422.h5')
        
#plot_pn_kc_ggn('D:\\biowulf_stage\\olfactory_network\\pn_kc_ggn_UTC2018_02_13__18_31_41-PID41339-JID61403628.h5')
#plot_pn_kc_ggn('D:\\biowulf_stage\\olfactory_network\\pn_kc_ggn_UTC2018_02_13__18_31_37-PID99588-JID61403623.h5')

def plot_spike_rasters(fname):
    """Plot PN and KC spikes along with GGN"""
    with h5.File(fname, 'r') as fd:
        fig, axes = plt.subplots(nrows=6, ncols=1, sharex=True)            
        fig.suptitle(fd.attrs['description'])
        config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'][0])
        
        # PN spike raster
        
        # LCA KC spike raster, MCA KC spike raster
        # LCA KC Vm
        # MCA KC Vm
        # GGN MCA Vm, GGN LCA Vm
        # GGN alphaL Vm
        kc_vm = fd['/data/uniform/kc/KC_Vm']
        t = np.arange(kc_vm.shape[1]) * kc_vm.attrs['dt']        
        pn_st = fd['/data/event/pn/pn_spiketime']
        for ii, node in enumerate(pn_st):
            axes[0].plot(pn_st[node], np.ones(len(pn_st[node]))*ii, 'b|')
        axes[0].set_title('PN spike raster')    
        kc_lca = int(kc_vm.shape[0] * config['kc']['lca_frac']+0.5)
        for ii in range(kc_lca):
            spike_idx = np.flatnonzero((kc_vm[ii, :-1] < threshold) & 
                                       (kc_vm[ii, 1:] > threshold))
            axes[1].plot(t[spike_idx], np.ones(spike_idx.shape) * ii, 'k|')
        axes[1].set_title('KC LCA')    
        for ii in range(kc_vm.shape[0] - kc_lca, kc_vm.shape[0]):
            spike_idx = np.flatnonzero((kc_vm[ii, :-1] < threshold) & 
                                       (kc_vm[ii, 1:] > threshold))
            axes[2].plot(t[spike_idx], np.ones(spike_idx.shape) * ii, 'k|')
        axes[2].set_title('KC MCA')    
        ggn_out_vm = fd['/data/uniform/ggn_output/GGN_output_Vm']
        for ii in random.sample(range(ggn_out_vm.shape[0]), 5):
            axes[3].plot(t, ggn_out_vm[ii,:])
        axes[3].set_title('GGN output')
        ggn_alphaL_vm = fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm']
        for ii in random.sample(range(ggn_alphaL_vm.shape[0]), 5):
            axes[4].plot(t, ggn_alphaL_vm[ii,:])
        axes[4].set_title('GGN alphaL input')           
        plt.show()    
        
# Check the Fourier transform of spontaneous acticity
## 

# Check spatial structure

# fcalyx = 'corr_net_UTC2018_02_27__18_28_51-PID89172-JID62378871.h5'
# fd = h5.File(os.path.join(datadir, fcalyx), 'r')
# fd = h5.File('D:/subhasis_ggn/model/mb/corr_net_UTC2018_03_07__00_44_53-PID22896-JID0.h5', 'r')
# fd.close()

fd = h5.File('D:/subhasis_ggn/model/mb/corr_net_UTC2018_03_07__18_57_35-PID4488-JID0.h5', 'r')
lca_labels = fd['/data/static/lca_cluster_labels/lca_cluster_labels']
print(len(lca_labels))
for entry in lca_labels:
    print(entry)
    break
kc_ggn_ca_syn = fd['/data/static/kc_ggn_CA_synapse/kc_ggn_CA_synapse']
ggn_kc_syn =  fd['/data/static/ggn_kc_synapse/ggn_kc_synapse']
print(ggn_kc_syn[0])
print(kc_ggn_ca_syn[0])
print(lca_labels[0])
sec_name = lca_labels[0]['sec']
for ii, sec in enumerate(ggn_kc_syn['pre']):
    if sec == sec_name:
        print(ggn_kc_syn[ii])
for ii, sec in enumerate(kc_ggn_ca_syn['post']):
    if sec == sec_name:
        print(kc_ggn_ca_syn[ii])
matching_pre_kc = kc_ggn_ca_syn[kc_ggn_ca_syn['post'] == sec_name]
matching_post_kc = ggn_kc_syn[ggn_kc_syn['pre'] == sec_name]

## This section is to confirm that pre and post synaptic KC
## populations corresponding to a single cluster label are identical
labels = sorted(list(set(lca_labels['label'].flatten())))
clus_kc = {}
for label in labels:
    label_info = lca_labels[lca_labels['label'] == label]
    idx = np.flatnonzero(lca_labels['label'] == label)
    kc_by_idx = set(['KC[{}].soma'.format(ii) for ii in idx])    
    print(label, len(label_info))
    print(kc_by_idx)
    label_info = set([(row['sec'], row['pos']) for row in label_info])
    print(label, len(label_info))
    kcs = []
    matching_pre_syn = []
    matching_post_syn = []
    for row in label_info:
        print('section info', row)
        matching_pre_syn += list(kc_ggn_ca_syn[(kc_ggn_ca_syn['post'] == row[0]) &
                                               (kc_ggn_ca_syn['postpos'] == row[1])].flatten())
        matching_post_syn += list(ggn_kc_syn[(ggn_kc_syn['pre'] == row[0]) &
                                      (ggn_kc_syn['prepos'] == row[1])].flatten())
    pre_kc = set(row[0] for row in matching_pre_syn)
    post_kc = set(row[1] for row in matching_post_syn)
    print('Number of mismatch in pre and post synaptic KCs for this label', len(pre_kc - post_kc))
    assert pre_kc == post_kc
    clus_kc[label] = pre_kc
    tmp = set([kc.decode() for kc in pre_kc])
    print( tmp - kc_by_idx)    
    break
##

### Check if GGN->KC synapses have 1-1 correspondence with cluster
### labels.
for label_info, syn_info in zip(lca_labels, ggn_kc_syn):
    print(label_info['label'], label_info['sec'], syn_info['pre'])




labeldf = nda.extract_cluster_info(fd, 'lca')
kc_ggn_syn = nda.get_kc_ggn_ca_syn_info(fd)
ggn_kc_syn = nda.get_ggn_kc_syn_info(fd)

print(labeldf)
ggn_kc_syn_label = labeldf.join(ggn_kc_syn)
assert np.all(ggn_kc_syn_label['ggn_sec'] == ggn_kc_syn_label['pre'])
kc_ggn_syn_label = labeldf.join(kc_ggn_syn)
ggn_kc_grp = ggn_kc_syn_label.groupby('label')
for label, group in kc_ggn_syn_label.groupby('label'):
    ggn_kc = ggn_kc_grp.get_group(label)
    assert set(group.index) == set(ggn_kc.index)
    assert set(group.post) == set(ggn_kc.pre)
fd.close()
#!! Verified - Subhasis 2018-03-07 18:15 UTC-5:00

#------------------------------------------------------------------
## These were tests for NSDF - static data is written as array of
## arrays - useful for mapping second dimension to other fields, like
## compartment no. etc
#
# test = h5.File('test.h5', 'w')
# a = np.array(('a', 1), dtype=[('name', h5.special_dtype(vlen=str)), ('age', int)])
# b = np.array(('alpha', 2), dtype=[('name', h5.special_dtype(vlen=str)), ('age', int)])
# data = np.vstack((a, b))
# data = np.array([('a', 1), ('alpha', 2)], dtype=[('name', h5.special_dtype(vlen=str)), ('age', int)])
# ds = test.create_dataset('test', data=data)
# test.close()
# test2 = h5.File('test.h5', 'r')
# print(test2['test'])
# test2.close()
#------------------------------------------------------------------

#########################
# Subhasis 2018-03-07 18:18 UTC-5:00
#
# Now try plotting the KC spike trains properly
datadir = 'D:\\biowulf_stage\\olfactory_network\\'
fcalyx = 'corr_net_UTC2018_02_27__18_28_51-PID89172-JID62378871.h5'
fd = h5.File(os.path.join(datadir, fcalyx), 'r')
kc_spikes = {row['source'].decode(): row['data']
             for row in fd['/map/event/kc/kc_spiketime'].value}

cluster_info = nda.extract_cluster_info(fd, 'lca')


fig, ax = plt.subplots(nrows=1, ncols=1)
ii = 0
spike_x, spike_y = [], []
for label, group in cluster_info.groupby('label'):
    print('Label:', label)
    for kc in np.char.decode(group.index.values.flatten().astype('S')):
        st = fd[kc_spikes[kc]].value
        spike_x.append(st)
        spike_y.append(np.ones(len(st)) * ii)
        ii += 1
    ax.axhline(y=ii-0.5, color='red')
ax.plot(np.concatenate(spike_x), np.concatenate(spike_y), 'k,')    
plt.show()
fd.close()

# Final code after trial above was transferred to network_data_analysis.py
datadir = 'D:\\biowulf_stage\\olfactory_network\\'
fname_calyx = 'corr_net_UTC2018_02_27__18_28_51-PID89172-JID62378871.h5'
fd_calyx = h5.File(os.path.join(datadir, fname_calyx), 'r')
fname_alpha = 'corr_net_UTC2018_02_27__18_28_50-PID86391-JID62378874.h5'
fd_alpha = h5.File(os.path.join(datadir, fname_alpha), 'r')
fig_alpha, axes = plt.subplots(nrows=2, ncols=1, sharex=True)
ax_alpha = axes[0]
ax_calyx = axes[1]
nda.plot_kc_spikes_by_cluster(ax_alpha, fd_alpha)
ax_alpha.set_title(fd_alpha.attrs['description'])
nda.plot_kc_spikes_by_cluster(ax_calyx, fd_calyx)
ax_calyx.set_title(fd_calyx.attrs['description'])
plt.show()
fd.close()

######################################################
## Check KC no if it was correct to go by node number
######################################################
# - no - it wasn't - code in pn_kc_ggn_plot_mpl.py
#--------------

fname = 'mb_net_UTC2018_03_13__02_15_25-PID74836-JID63738870.h5'
fda = h5.File(os.path.join(datadir, fname), 'r')
kcs = nda.get_kcs_by_region(fda, 'mca')
idx = nda.get_kc_vm_idx_by_region(fda, 'mca')
fig, ax = plt.subplots(nrows=1, ncols=1)
myplot.plot_kc_spikes_by_cluster(ax, fda, 'lca')
ax.set_title(fda.attrs['description'])
plt.show()


fda.close()

# Wed Mar 14 14:26:34 EDT 2018

# This file has list of simulation outputs by jobid.
#
#cluspn: whether PN->KC connections are shared by KCs in same
# cluster.
#
# ca: kc->ggn conn in CA
# clusca: kc->ggn conn in CA are within same clustered
files = [
    'mb_net_UTC2018_03_13__02_15_25-PID74836-JID63738870.h5',
    'mb_net_UTC2018_03_13__02_35_05-PID4200-JID63739839.h5',
    'mb_net_UTC2018_03_13__02_35_08-PID45474-JID63739862.h5',
    'mb_net_UTC2018_03_13__02_54_09-PID5464-JID63741201.h5',
    'mb_net_UTC2018_03_13__02_15_26-PID10386-JID63738868.h5',
    'mb_net_UTC2018_03_13__02_35_06-PID77583-JID63739829.h5',
    'mb_net_UTC2018_03_13__02_35_08-PID45475-JID63739868.h5',
    # 'mb_net_UTC2018_03_13__15_18_59-PID79964-JID63776172.h5',
    'mb_net_UTC2018_03_13__02_25_11-PID103579-JID63739341.h5',
    'mb_net_UTC2018_03_13__02_35_07-PID67280-JID63739844.h5',
    'mb_net_UTC2018_03_13__02_54_06-PID136087-JID63741200.h5',
#    'mb_net_UTC2018_03_13__15_19_00-PID47752-JID63776173.h5'
]

fname_parts = pd.DataFrame(index=range(len(files)),
                           columns=['filename',
                                    'timestamp',
                                    'pid', 'jid'])
tmp = pd.Series(files).str.split('-')
ts , pid, jid = zip(*tmp)
print(ts)
ts = [x.partition('UTC')[-1] for x in ts]
print(ts)
print(pid)
fname_parts.filename = files
fname_parts.timestamp = ts
fname_parts.pid = [x[3:] for x in pid]
fname_parts.jid = [x[3:-3] for x in jid]
print(fname_parts.jid)
finfo = pd.read_csv('simlist_201803.csv')
finfo.jobid = finfo.jobid.apply(str)
print(finfo[finfo.RA == 70.0])
print(finfo.jobid)
print(fname_parts.jid)
combo = pd.merge(finfo, fname_parts, left_on=['jobid'], right_on=['jid'])
for fname in combo.filename:
    try:
        with h5.File(os.path.join(datadir, fname), 'r') as fd:
            print('##', fname)
            config = yaml.load(fd.attrs['config'])
            print('pn_kc_clusterd?', config['pn_kc_syn']['clustered'])
            print('kc->ggn ca present?', config['kc_ggn_CA_syn']['present'])
            print('kc->ggn ca clustered?', config['kc_ggn_CA_syn']['clustered'])
            print('kc->ggn ca regional?', config['kc_ggn_CA_syn']['regional'])
    except OSError:
        print('Cannot load', fname)
    continue
    
#     try:
#         fig, axes = myplot.plot_spike_rasters(os.path.join(datadir, fname), by_cluster=True)
#     except OSError:
#         print('Could not load', fname)
#         pass
# plt.show()

            
fname = combo[combo.jid == '63738868'].iloc[0].filename
fpath = os.path.join(datadir, fname)
nda.get_shared_pn_frac(fname, verbose=True)             # Verified a little more than 80% shared (81-84%).
#######################################################
# Thu Mar 15 16:06:15 EDT 2018
# Moved file selection to select_good_files function
######################################################
finfo = nda.select_good_files(datadir, 'simlist_201803.csv')
# fname = finfo[finfo.jid == '63879957'].iloc[0].fname
# fpath = os.path.join(datadir, fname)
# print(fpath)
# fd = h5.File(fpath, 'r')
# cfg = yaml.load(fd.attrs['config'])
# print(cfg['kc_ggn_CA_syn'])
# fname = finfo[finfo.jid == '63879957'].iloc[0].fname
# fpath = os.path.join(datadir, fname)
# print(fpath)
# fd = h5.File(fpath, 'r')
# cfg = yaml.load(fd.attrs['config'])
# print(cfg['kc_ggn_CA_syn'])
# fd.close()

# Get a summary about ca-clustering config
for fname in finfo.fname:
    fpath = os.path.join(datadir, fname)
    with h5.File(fpath, 'r') as fd:
        print(fpath.rpartition('JID')[-1])
        print(cfg['kc_ggn_CA_syn'])

for fname in finfo[finfo.clussize == 100].fname:
    fpath = os.path.join(datadir, fname)
    fig, axes = myplot.plot_spike_rasters(fpath, kde_bw=200, by_cluster=True)
plt.show()    

for fname in finfo[finfo.clussize == 1000].fname:
    fpath = os.path.join(datadir, fname)
    fig, axes = myplot.plot_spike_rasters(fpath, kde_bw=200, by_cluster=True)
plt.show()    

fig, axes = plt.subplots(nrows=int((finfo.shape[0] + 0.5) // 2),
                         ncols=2, sharex=True, sharey=True)
Em_ggn = -51.0
for ii in range(finfo.shape[0]):
    row = finfo.iloc[ii]
    ax = axes[ii // 2, ii % 2]
    ax.set_title(str(row.jid))
    fpath = os.path.join(datadir, row.fname)
    print('#', row.jobid, '#')
    with h5.File(fpath, 'r') as fd:
        ds = fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm']
        s_alphaL = myplot.plot_ggn_vm(ax, fd, ds, 'alphaL', 10,
                                      color='k')
        vmax_alphaL = np.max(ds[:, :])
        print('alphaL Vmax =', vmax_alphaL - Em_ggn)
        ds = fd['/data/uniform/ggn_output/GGN_output_Vm']
        s_lca = myplot.plot_ggn_vm(ax, fd, ds, 'lca', 10, color='r')
        vmax_ca = np.max(ds[:, :])
        deflection_ratio = (vmax_ca - Em_ggn) / (vmax_alphaL - Em_ggn)
        print('CA Vmax =', vmax_ca - Em_ggn)
        print('Ratio =', deflection_ratio)
        print(row)
plt.show()    
    

#################################################    
# 2018-03-16 testing old data from SfN meeting
#################################################
def plot_RM_RA_param_sweep_deflection(fname):
    with h5.File(fname, 'r') as fd:
        deflection_map = defaultdict(dict)
        fig, ax = plt.subplots(nrows=1, ncols=1)
        for simname in fd:
            sgrp = fd[simname]
            settings = dict(sgrp.attrs)
            syndata = sgrp['syninfo'].value
            syninfo = dict(sgrp['syninfo'].attrs)
            vmnodes = {sgrp[name].attrs['section']: sgrp[name]
                       for name in sgrp if name.startswith('v_')}
            avg_input = np.mean([vmnodes[synsec].value for synsec in syndata], axis=0)
            lca_nodes = [node for node in vmnodes if 'dend_{}'.format(ng.custom_types['LCA']) in node.decode()]
            lca_sample = random.sample(lca_nodes, 5)
            avg_output = np.mean([vmnodes[sec].value for sec in lca_sample], axis=0)
            input_defl = max(avg_input + 51.0)
            output_defl = max(avg_output + 51.0)
            print(input_defl, output_defl)
            ratio = output_defl / input_defl
            print(ratio)
            ax.plot(avg_output, 'b-', alpha=0.3)
            ax.plot(avg_input, 'r-', alpha=0.3)
            deflection_map[settings['RA']][settings['RM']] = (ratio, input_defl)
        for RA, RM_map in deflection_map.items():
            RM_list = np.array(sorted(list(RM_map.keys())))
            deflections = [RM_map[RM][0] for RM in RM_list]
#            ax.plot(RM_list/1e3, deflections, label=str(RA))
#        ax.legend()    
        plt.show()
        
print('Checking old data from parameter sweep for synchronous input.')
fname = 'D:/subhasis_ggn/model/nrn/data/B_Vm_multi_syn_20170510_171812.h5'
plot_RM_RA_param_sweep_deflection(fname)

# Conclusion: this was inputs at 10 random alpha lobe sections -
# output voltage attenuates to ~20%

# Ran another simulation with new GGN trace:
# python localized_input_output_passive_sweep.py -t 13.33 -n 10 -c GGN_20170309_sc -f ..\mb\cell_templates\GGN_20170309_sc.hoc --RAmin=10 --RAmax=150 --RAdivs=5 -i 8 --gk=0 --gca=0    
#'D:/data/rays3/ggn/single_syn_input/B_Vm_multi_syn_20180316_175240.h5'        
# - this plot does not work because just one RM
# fname = 'D:/data/rays3/ggn/single_syn_input/B_Vm_multi_syn_20180316_181801.h5' 
# - this had too small deflection
#fname = '/data/rays3/ggn/single_syn_input/B_Vm_multi_syn_20180316_182447.h5'
# gsyn increased to 200 nS - RM was too small - should have been specified in ohm-cm2
fname = 'D:/data/rays3/ggn/single_syn_input/B_Vm_multi_syn_20180316_184656.h5'
plot_RM_RA_param_sweep_deflection(fname)

####################################################
# Getting back to network model data
####################################################
# Plot KC spikes
for ii in range(finfo.shape[0]):
    row = finfo.iloc[ii]
    fig = plt.figure(row.jid)
    ax = fig.add_subplot(111)
    fpath = os.path.join(datadir, row.fname)
    with h5.File(fpath, 'r') as fd:
        myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca')
plt.show()

#
# Mon Mar 19 18:29:45 EDT 2018
#
# With reduced GGN->KC conductance, KCs spike way too much.
fname = 'd:\\biowulf_stage\\olfactory_network\\mb_net_UTC2018_03_17__00_07_26-PID34884-JID64085794.h5'
fd = h5.File(fname, 'r')
kc_vm = fd['/data/uniform/kc/KC_Vm']
indices = random.sample(range(kc_vm.shape[0]), 3)
for ii in indices:
    plt.plot(kc_vm[ii,:], alpha=0.5)
plt.show()

fname = 'D:/data/rays3/ggn/single_syn_input/B_Vm_multi_syn_20180319_182619.h5'
plot_RM_RA_param_sweep_deflection(fname)

####################################################
# Getting back to network model data
####################################################
# Plot KC spikes
for ii in range(finfo.shape[0]):
    row = finfo.iloc[ii]
    fig = plt.figure(row.jid)
    ax = fig.add_subplot(111)
    fpath = os.path.join(datadir, row.fname)
    with h5.File(fpath, 'r') as fd:
        myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca')
plt.show()


###############################################################################
## 2018-04-04
## checking the GGN Vm in simulations with low GGN->KC conductance
###############################################################################
plt.rcParams['figure.figsize'] = [10, 8]
plt.rcParams['figure.titlesize'] = 'small'
plt.rcParams['figure.frameon'] = False
plt.rcParams['figure.subplot.hspace'] = 0
plt.rcParams['figure.subplot.wspace'] = 0
plt.rcParams['font.size'] = 8
plt.rcParams['lines.linewidth'] = 1              
plt.rcParams['lines.markersize'] = 3.0
            

jid = '64085710'
flist = os.listdir(datadir)
match = [fname for fname in flist if 'JID{}'.format(jid) in fname]
if len(match) > 1:
    print('Two files with same jobid.', match)
fname = match[0]

fig_noca, ax_noca = myplot.plot_spike_rasters(os.path.join(datadir, fname))
ts = fname.partition('-')[0].partition('UTC')[-1]
fig_noca.savefig('figures/utc_{}_jid_{}.png'.format(ts, jid))

jid = '64085711'
match = [fname for fname in flist if 'JID{}'.format(jid) in fname]
if len(match) > 1:
    print('Two files with same jobid.', match)
fname = match[0]

fig_ca, ax_ca = myplot.plot_spike_rasters(os.path.join(datadir, fname))
ts = fname.partition('-')[0].partition('UTC')[-1]
fig_ca.savefig('figures/utc_{}_jid_{}.png'.format(ts, jid))

#with PdfPages('figures/ggn_vm.pdf') as pdf:
#    meta = {'Title': 'PN->KC<->GGN network simulation',
#                        'Author': 'Subhasis Ray'}
#    pdf.attach_note('jobid: 64085710, gmax(GGN->KC): 100 pS, KC<->GGN cluster size: 1000, no clustering in PN->KC')
#    pdf.savefig(fig_noca)
#    pdf.attach_note('jobid: 64085711, gmax(GGN->KC): 100 pS, KC<->GGN cluster size: 1000, no clustering in PN->KC')
#    pdf.savefig(fig_ca)

# This turns out to be too heavy to display as PDF        

jid = '64085794'
flist = os.listdir(datadir)
match = [fname for fname in flist if 'JID{}'.format(jid) in fname]
if len(match) > 1:
    print('Two files with same jobid.', match)
fname = match[0]

fig_ca, ax_ca = myplot.plot_spike_rasters(os.path.join(datadir, fname))

jid = '64085793'
flist = os.listdir(datadir)
match = [fname for fname in flist if 'JID{}'.format(jid) in fname]
if len(match) > 1:
    print('Two files with same jobid.', match)
fname = match[0]

fig_noca, ax_noca = myplot.plot_spike_rasters(os.path.join(datadir, fname))




fig, ax = plot_data_by_jobid('64210017', datadir, save=True, by_cluster=True)
fname = find_file_by_jid('64210017', datadir)[0]
with h5.File(fname, 'r') as fd:
    print(yaml.dump(yaml.load(fd['/model/filecontents/mb/network/config.yaml'].value[0])))
    
fig, ax = plot_data_by_jobid('64210018', datadir, save=True, by_cluster=True)
fig, ax = plot_data_by_jobid('64211960', datadir, save=True, by_cluster=True)    
fig, ax = plot_data_by_jobid('64440225', datadir, save=True, by_cluster=True) 

fname = find_file_by_jid('64440225', datadir)[0]
with h5.File(fname, 'r') as fd:
    print(yaml.dump(yaml.load(fd['/model/filecontents/mb/network/config.yaml'].value[0])))

print(get_config_by_jid('64819232', datadir, default_flow_style=False, default_style='' ))
fig, ax = plot_data_by_jobid('64819232', datadir, save=True, by_cluster=True)


# 2018-04-06 plot individual KCs and see if it matches Perez-Orive 2004
# This simulation had reduced GGN->KC inhibitory conductance
kcs = []
fname = find_file_by_jid('64085794', datadir)[0]
with h5.File(fname, 'r') as fd:
    config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'][0])
    kcs = [node for node in fd['/data/event/kc/kc_spiketime']]
    fig, axes = plt.subplots(nrows=5, ncols=1)
    for ax in axes:
        st = []
        while len(st) == 0:
            sample = random.choice(kcs)
            st = fd['/data/event/kc/kc_spiketime'][sample][:]
            print(sample, st)
            
        simtime = Q_(config['stimulus']['onset']).to('ms').m +  \
                  Q_(config['stimulus']['duration']).to('ms').m + \
                  Q_(config['stimulus']['tail']).to('ms').m
        kde_grid = np.linspace(0, simtime, 100.0)
        kde_bw = 50.0
        myplot.plot_population_KDE(ax, [st], kde_grid, kde_bw, color='y')
## No, instead of different KCs spiking at different narrow time windows, here each kc spikes throughout the stimulus.        


## 2018-04-11
## Check the PN popolations spiking in 50 ms bins
fname = find_file_by_jid('64085794', datadir)[0]
pn_spikes = {}
with h5.File(fname, 'r') as fd:
    for pn in fd['/data/event/pn/pn_spiketime']:
        pn_spikes[pn] = fd['/data/event/pn/pn_spiketime'][pn].value

plt.figure()                 
for ii, st in enumerate(pn_spikes.values()):
    plt.plot(st, np.ones(st.shape)*ii+1, '|')
    
binwidth = 50.0
bins = np.arange(0, 5000.1, binwidth)
plt.hist(np.concatenate(list(pn_spikes.values())), bins=bins)

spiking_cells = defaultdict(set)
for jj, (left, right) in enumerate(zip(bins[:-1], bins[1:])):
    for pn, st in pn_spikes.items():
        if np.any((st > left) & (st < right)):
            spiking_cells[jj].add(pn)
for ii in range(1, len(spiking_cells)):
    print(ii-1, '->', ii, ':', len(spiking_cells[ii-1]), ' shares', len(spiking_cells[ii].intersection(spiking_cells[ii-1])))
pops = [len(sc) for sc in spiking_cells.values()]                 
min(pops)
max(pops)


###### Checking the KC spike times again
## 2018-04-11 Reconfirmed that spikes in individual KCs are spread out, 
## and not in narrow time bands like experiment

kcs = []
fname = find_file_by_jid('64085794', datadir)[0]
with h5.File(fname, 'r') as fd:
    config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'][0])
    kcs = [node for node in fd['/data/event/kc/kc_spiketime']]
    fig, axes = plt.subplots(nrows=1, ncols=1)
    for ii in range(20):
        st = []
        while len(st) == 0:
            sample = random.choice(kcs)
            st = fd['/data/event/kc/kc_spiketime'][sample][:]
            if len(st) == 0:
                print('No spikes in', sample)
        axes.plot(st, np.ones(st.shape) * ii + 1, '|')
plt.show()
### 2018-05-16 Check simulations with shifting PN activity
## 
fname = find_file_by_jid('976212', datadir)[0]
fd = h5.File(fname, 'r')
for name, st in fd['/data/event/kc/kc_spiketime'].items():
    if len(st) > 0:
        print(name)
## No spike in KCs
fig, axes = plt.subplots(nrows=3, ncols=1, sharex='all')
ii = 0
for name, st in fd['/data/event/pn/pn_spiketime'].items():
    axes[0].plot(st, np.ones(st.shape) * ii, '|')
    ii += 1
ggn_basal = fd['/data/uniform/ggn_basal/GGN_basal_Vm']
rows = np.random.choice(range(ggn_basal.shape[0]), 3, replace=False)
for row in rows:
    vm = ggn_basal[row, :]
    axes[1].plot(np.arange(len(vm)) * ggn_basal.attrs['dt'], vm)

ggn_out = fd['/data/uniform/ggn_output/GGN_output_Vm']
rows = np.random.choice(range(ggn_out.shape[0]), 10, replace=False)
for row in rows:
    vm = ggn_out[row, :]
    axes[2].plot(np.arange(len(vm)) * ggn_out.attrs['dt'], vm)
plt.show()
fd.close()
## Significant attenuation in GGN calyx
fig, ax = plot_data_by_jobid('976212', datadir, save=True, by_cluster=True)
plt.show()

##
fig, ax = plot_data_by_jobid('863985', datadir, save=True, by_cluster=True)
plt.show()

fig, ax = plot_data_by_jobid('1079725', datadir, save=True, by_cluster=True)
plt.show()

fig, ax = plot_data_by_jobid('1079730', datadir, save=True, by_cluster=True)
plt.show()

fig, ax = plot_data_by_jobid('1080421', datadir, save=True, by_cluster=True)
plt.show()

fig, ax = plot_data_by_jobid('1101504', datadir, save=True, by_cluster=True)
plt.show()

# Comparing with Perez-Orive 2002 - reducing GGN conductance increases
# KC spiking? (Although they say KC firind rate did not increase
# significantly, fig 6 shows increase)
fig, ax = plot_data_by_jobid('64085710', datadir, save=True, by_cluster=True)
plt.show()

fig, ax = plot_data_by_jobid('64941216', datadir, save=True, by_cluster=True)
plt.show()

# Check if this clustering is real -

## Randomize the clustered KCs and plot
fname = find_file_by_jid('64941216', datadir)
fd = h5.File(fname[0], 'r')
kc_st = {k: v.value for k, v in fd['/data/event/kc/kc_spiketime'].items()}
fig, ax = plt.subplots(nrows=1, ncols=1)
import random
for ii, kc in enumerate(random.sample(list(kc_st.keys()), len(kc_st))):
    ax.plot(kc_st[kc], ii*np.ones(len(kc_st[kc])), 'k,')
plt.show()
fd.close()
## Looks just like unclustered case.

## can we cluster the random connection case to obtain such clusters?
import pyspike as spk  # PySpike implements some methods by Thomas Kreutz
fname = find_file_by_jid('63879737', datadir)[0]
kc_st = {}
with h5.File(fname, 'r') as fd:
    for kc, st in fd['/data/event/kc/kc_spiketime'].items():
        kc_st[kc] = st.value
profile = spk.spike_profile(kc_st['1000'], kc_st['1003'])
print("SPIKE distance: %.8f" % profile.avrg())
x, y = profile.get_plottable_data()
plt.plot(x, y, '--k')
plt.show()

########################################
## Section START
########################################
from cluster_spiketrains import spiketrains

kc_st2 = spiketrains(kc_st, left=1e3, right=4e3)  ## This turned out to be too much for memory (50000)^2 floats.
spike_counts = kc_st2.spike_count_vec(50.0)
km = kc_st2.cluster_kmeans(100, 50.0)

label_map = defaultdict(list)
for ii, label in enumerate(km.labels_):
    label_map[label].append(ii)

units = kc_st2.units()
ii = 0
spiketimes = []
ypos = []
for label in sorted(list(label_map.keys())):
    for index in label_map[label]:
        st = kc_st2.unit_st[units[index]].spikes
        spiketimes += (list(st))
        ypos += [ii] * len(st)
        ii += 1
plt.plot(spiketimes, ypos, 'k,')
plt.show()
# clustered = kc_st2.cluster_spikedistance()
## Does the KC response expand, as in Perez-Orive 2002 and Lin et al 2014 Nature Neurosci.?
fig, ax = plot_data_by_jobid('63879737', datadir)
plt.show()
# Problem is how can we be sure this measurement is good? The
# spikedistance measure could not be used for size. So trying a
# subsample.
kc3 = np.random.randint(low=0, high=len(kc_st), size=5000)
kc3_st = spiketrains({kc: kc_st[str(kc)] for kc in kc3}, left=1e3, right=4e3)
dist = kc3_st.compute_spike_distance()
## It is also worth trying to cluster KCs by their number of shared
## PNs as a measure of nearness.
## This I did using mcl 
mcl_fname = 'UTC2018_03_14__14_17_08_PID83321_JID63879737.clusters'
clusters = {}
with open(os.path.join(datadir, mcl_fname), 'r') as fd:
    for line in fd:
        print(line)
        print('-----------')
pn_kc_adj = get_pn_kc_adjacency(fname)
########################################
## Section END
########################################

## compare low 
fig, ax = plot_data_by_jobid('64085793', datadir)
plt.show()

## 1585119 is with clustered PN and stimulus start=1s end=2s, GGN is fed with recorded data
jid = '1585119'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
## jid 1585024 is without clustering in PNs
jid = '1585024'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
plt.show()

## These were without clustering KC<->GGN connection, nor were PN->KC clustered
jid = '1576735'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)

jid = '1576726'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
plt.show()

## These were without clustering KC<->GGN connection, but PN->KC clustered
jid = '1576207'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
fig.tight_layout()
jid = '1576217'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
fig.tight_layout()
plt.show()

##
jid = '1298124'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)

jid = '1298126'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
plt.show()

##
jid = '1101503'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
plt.show()

jid = '1585024'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
plt.show()

## with PN->KC gmax=4.5pS and GGN->KC gmax=0.1pS
jid = '1308282'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
plt.show()

jid = '1308286'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir)
plt.show()

nodes = list(range(1000)) + list(range(16000, 17000))
kc_nodes = [str(ii) for ii in nodes]

f1 = nda.find_file_by_jid('1308282', datadir)[0]
with h5.File(f1, 'r') as fd:
    x, y = nda.get_event_times(fd['/data/event/kc/kc_spiketime'], kc_nodes)
    ax = plt.subplot(211)
    plt.plot(np.concatenate(x), np.concatenate(y), ',')

f2 = nda.find_file_by_jid('1308286', datadir)[0]
with h5.File(f2, 'r') as fd:
    x, y = nda.get_event_times(fd['/data/event/kc/kc_spiketime'], kc_nodes)
    plt.subplot(212, sharex=ax, sharey=ax)
    plt.plot(np.concatenate(x), np.concatenate(y), ',')
plt.show()

########################
## Further analyzing 1576207
jid = '1576207'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fname = nda.find_file_by_jid(jid, datadir)[0]
fd = h5.File(fname, 'r')
kc_st = {}
for kc, st in fd['/data/event/kc/kc_spiketime'].items():
    kc_st[kc] = st.value

tmp = kc_st
kc_st = {}
for k, v in tmp.items():
    kc_st[int(k)] = v

x = []
y = []
for ii in range(500):
    st = kc_st[ii]
    if len(st) == 0:
        continue
    x += list(st)
    y += [ii] * st.shape[0]

plt.plot(x, y, '|')
plt.show()

## New data: PN->KC 2pS, GGN->KC 0.1nS,  KC cluster size=1000
ipdb.set_trace()
jid = '1931700'  # KC->GGN CA
fpath = nda.find_file_by_jid(jid, datadir)[0]
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
figh, axh = myplot.plot_kc_spike_count_hist(fpath)
fname = os.path.basename(fpath)
figh.tight_layout()
figh.savefig('C:/Users/Subhasis/Documents/insect_brain/NOTES_DATA/{}_kc_spike_count_hist.png'.format(fname.rpartition('.h5')[0]))
jid = '1931706'  # KC->GGN xCA
fpath = nda.find_file_by_jid(jid, datadir)[0]
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
figh, axh = myplot.plot_kc_spike_count_hist(fpath)
fname = os.path.basename(fpath)
figh.tight_layout()
figh.savefig('C:/Users/Subhasis/Documents/insect_brain/NOTES_DATA/{}_kc_spike_count_hist.png'.format(fname.rpartition('.h5')[0]))
plt.show()

#### hxa: GGN Vm recorded from hxa used for dynamic clamp KC cluster size=1000 
jid = '2009430'
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
fig, ax = plot_data_by_jobid(jid, datadir, save=True)

# What is the number of spikes from a KC?
fname = nda.find_file_by_jid(jid, datadir)[0]
with h5.File(fname, 'r') as fd:
    kc_st_grp = fd['/data/event/kc/kc_spiketime']
    lca_kcs = nda.get_kc_spike_nodes_by_region(fd, 'lca')
    lca_spike_count = [len(kc_st_grp[kc]) for kc in lca_kcs] 
    mca_kcs = nda.get_kc_spike_nodes_by_region(fd, 'mca')
    mca_spike_count = [len(kc_st_grp[kc]) for kc in mca_kcs]
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex='all', sharey='all')
    axes[0].hist(lca_spike_count)
    axes[1].hist(mca_spike_count)
plt.show()    
        
jid = '2007844'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)

jid = '1992040'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
plt.show()

# X X X X These have small number of PNs per KC - not realistic. X X X X
## Low fan in for KCs
jid = '2230234'   # ~40 PN/KC
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
plt.show()
## No spiking
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
# X X X X  X X X X  X X X X  X X X X 


## New simulations - comparing alpha lobe KC input vs calyx input
jid = '2397844'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)

jid = '2398081'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()

# With GGN->KC gmax=0.2nS
jid = '2466867'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

jid = '2466468'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

### This is just for comparing depolarization in older simulations
jid = '64941197'
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
plt.show()

## Old simulations - dump plots
jid = '61403623'
datadir = 'F:/biowulf_ggn/olfactory_network'
fpath = nda.find_file_by_jid(jid, datadir)

fd = h5.File(fpath[2], 'r')
print(list(fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm'].attrs.keys()))
fd.close()
with h5.File(fpath[2], 'r') as fd:
    fig, ax = plt.subplots(nrows=3, ncols=1, sharex='all')
    pn_x, pn_y = nda.get_event_times(fd['/data/event/pn/pn_spiketime'])
    ax[0].plot(np.concatenate(pn_x), np.concatenate(pn_y), 'k,')
    ax[0].set_title('PN spikes')
    with h5.File(fpath[0], 'r') as kcf:
        kc_x, kc_y = nda.get_event_times(kcf)
        ax[1].plot(np.concatenate(kc_x), np.concatenate(kc_y), 'k,')
    ax[1].set_title('KC spikes')
    ggn_alpha_Vm = fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm'].value
    dt = fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm'].attrs['dt']
    index = np.random.randint(0, ggn_alpha_Vm.shape[0], 5)
    for ii in index:
        ax[2].plot(np.arange(ggn_alpha_Vm.shape[1]) * dt, ggn_alpha_Vm[ii, :], 'r', label='alphaL')
    ggn_out_Vm = fd['/data/uniform/ggn_output/GGN_output_Vm'].value
    index = np.random.randint(0, ggn_out_Vm.shape[0], 5)
    for ii in index:
        ax[2].plot(np.arange(ggn_out_Vm.shape[1]) * dt, ggn_out_Vm[ii, :], 'b', label='calyx')
    ax[2].legend()
    ax[2].set_title('GGN Vm')
    fig.suptitle(os.path.basename(fpath[2]))
    fig.savefig('figures/{}'.format(os.path.basename(fpath[2]).replace('.h5', '.png')))

plt.show()

jid = '2706143'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()


jid = '2706144'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()

jid = '2710326'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()


jid = '2710908'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()

jid = '2710961'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()

jid = '2714059'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

jid = '2714061'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

jid = '2714083'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

jid = '2715134'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)


jid = '2716893'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

jid = '2716842'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)


jid = '2716857'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

jid = '2718859'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

# jid = '2718860'  ## failed 

jid = '2718861'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

jid = '2718862'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)

jid = '2718863'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)


jid = '3047095'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)


jid = '3047110'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins, save=True)
plt.show()
config = nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style='' )
print(config)
## SCRATCH: For developing of kc_spiking_correlates.
# jid = '3047110'
# fpath = nda.find_file_by_jid(jid, datadir)
# print(type(fpath[0]))
# fd = h5.File(fpath[0], 'r')
# x = fd['/data/static/ggn_kc_synapse/ggn_kc_synapse']['post', 'gbar']
# ggn_kc = pd.DataFrame({'ginh': x['gbar'][:, 0], 'kc': x['post'][:, 0]})
# y = fd['/data/static/pn_kc_synapse/pn_kc_synapse']['post', 'gmax'][:, 0]
# pn_kc = pd.DataFrame({'gexc': y['gmax'], 'kc': y['post']})
# # pn_kc_gmax = pd.DataFrame(fd['/data/static/pn_kc_synapse/pn_kc_synapse']['post', 'gmax'][:, 0])
# pn_kc_tot = pn_kc.groupby('kc').sum().reset_index()
# combined = pd.merge(ggn_kc, pn_kc_tot, on='kc')

# conf = yaml.load(fd.attrs['config'].decode())
# print(conf['stimulus'])
# print(conf['pn']['offdur'])
# kc_st = nda.get_kc_st(fd)
# fd.close()
jid = '3047110'
fpath = nda.find_file_by_jid(jid, datadir)
# print(type(fpath[0]))
# fd = h5.File(fpath[0], 'r')
# pn_kc_syn_path =  '/data/static/pn_kc_synapse/pn_kc_synapse'
# kc_st_path = '/data/event/kc/kc_spiketime'
# ggn_kc_syn_path = '/data/static/ggn_kc_synapse/ggn_kc_synapse'
# ggn_kc_gbar = pd.DataFrame(fd[ggn_kc_syn_path]['post', 'gbar'][:, 0])
# pn_kc_gmax = pd.DataFrame(fd[pn_kc_syn_path]['post', 'gmax'][:, 0])
# pn_kc_gmax_by_kc = pn_kc_gmax.groupby('post').sum().reset_index()
# conf = yaml.load(fd.attrs['config'].decode())
# delt = Q_(conf['stimulus']['duration']).to('ms').m
# try:
#     delt += Q_(conf['pn']['offdur']).to('ms').m
# except KeyError:
#     pass
# print('PN excited for', delt, 'ms')
# kc_spike_count = [(kc, len(st)) for kc, st in nda.get_kc_st(fd)]
# kc_spike_count = pd.DataFrame(data=kc_spike_count, columns=['kc', 'spikes'])
# kc_spike_rate = pd.DataFrame({'kc': kc_spike_count['kc'],
#                               'spikerate': kc_spike_count['spikes'] / delt})
# pn_ggn = pd.merge(pn_kc_gmax_by_kc, ggn_kc_gbar, on='post')
# pn_ggn_kc_rate = pd.merge(pn_ggn, kc_spike_rate, left_on='post', right_on='kc')

jid = '3047110'
fpath = nda.find_file_by_jid(jid, datadir)
correlates = nda.kc_spiking_correlates(fpath[0])
correlates = correlates.rename(columns={'gbar': 'ginh', 'gmax': 'gexc'})

fig, axes = plt.subplots(nrows=2, ncols=1, sharey='all')
axes[0].scatter(correlates['gexc']*1e3, correlates['spikerate'])
# axes[0].legend()
axes[1].scatter(correlates['ginh']*1e3, correlates['spikerate'])
axes[1].set_xlim(-0.5, 2)
# axes[1].legend()
axes[0].set_xlabel('Maximum total excitatory conductance (nS)')
axes[0].set_title('spikerate vs g_exc')
axes[1].set_ylabel('KC spike rate during PN activity (/s)')
axes[1].set_xlabel('Maximum inihibitory conductance (nS)')
axes[1].set_title('spikerate vs g_inh')
# despine(axes[0])
# despine(axes[1])
fig.set_size_inches(4, 6)
fig.tight_layout()
fig.savefig('spikerate_exc_inh.png', transparent=True)
plt.show()

fig, axes = plt.subplots()
axes.scatter(correlates['gexc']*1e3, correlates['ginh']*1e3, c=correlates['spikerate'], s=3, alpha=0.3)
plt.show()

jid = '3251240'
fig, ax = plot_data_by_jobid(jid, datadir, save=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins, save=True)
plt.show()

## Fri Jun 15 12:24:25 EDT 2018
## Back check for labmeet
jid = '62378871' # Alpha + CA
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()

jid = '62378874'  # alpha only


jid = '3621292'  # alphaL, clustered PN, shifting PN
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)

jid = '3620975'   # alphaL + CA, clustered PN, shifting PN
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)


jid = '3620994'   # alphaL + CA, clustered PN 
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()

print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))


jid = '3621315'   # alphaL, unclustered PN 
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths) 
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()

jid = '3621387'     # alphaL, unclustered PN, no regional separation between KCs.
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()


jid = '3621019'   # alphaL + CA, PN not shifting
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()

jid = '3621435'
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()


###############################
# Thu Jun 21 18:02:59 EDT 2018
###############################

## Clustered PN
jid = '3690512'   # RA=70 ohm-cm
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins, save=True)
plt.show()

jid = '3690513'   # RA=70 ohm-cm
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins, save=True)
plt.show()


jid = '3690520'   # RA=150 ohm-cm
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins, save=True)
plt.show()

#########
### Display just a few clusters - so we can recognize the difference
### in spike times
lca_clusters = nda.get_kc_clusters(fpath, 'lca')
print(lca_clusters.shape)
lca_clusters = pd.DataFrame(data=lca_clusters[:, 0])
lca_grp = lca_clusters.groupby('label')
sec_st = {}
with h5.File(fpath, 'r') as fd:
    sec_st_path = nda.get_kc_event_node_map(fd)
    for kc, path in sec_st_path.items():
        sec_st[kc] = fd[path].value
        
keys = np.random.choice(lca_grp.groups.keys(), size=5, replace=False)

fig, axes = plt.subplots(nrows=len(keys), sharex='all')
for k, ax in zip(keys, axes.flatten()):
    seclist = lca_grp.get_group(k)['sec']
    for ii, sec in enumerate(seclist):
        st = sec_st[sec]
        ax.plot(st, ii * np.ones(len(st)), 'k,')
plt.show()

#########
# Plot presynaptic PN spike raster for each cluster
pn_st = {}
with h5.File(fpath, 'r') as fd:
    pn_kc_syn = fd[nda.pn_kc_syn_path]['pre', 'post']
    for node in fd['/data/event/pn/pn_spiketime'].values():        
        pn_st[node.attrs['source'].decode('utf-8')] = node.value
pn_kc_syn = pd.DataFrame(pn_kc_syn[:, 0])            
kc_pre_pn = pn_kc_syn.groupby('post')

fig, axes = plt.subplots(nrows=len(keys), ncols=2, sharex='all')
prev = set()    
for ii, k in enumerate(keys):
    seclist = lca_grp.get_group(k)['sec']
    cluster_pre = set()
    for jj, sec in enumerate(seclist):
        st = sec_st[sec]
        axes[ii, 1].plot(st, jj * np.ones(len(st)), 'k|')
        current = set(kc_pre_pn.get_group(sec)['pre'])
        cluster_pre.update(current)
        print(len(current.intersection(prev)))
        prev = current
    for jj, pre in enumerate(cluster_pre):
        axes[ii, 0].plot(pn_st[pre], jj * np.ones(len(pn_st[pre])), 'k,')
    print('#################')
# for ax in axes.flat:
#     ax.set_xlim(1000, 2000)
plt.show()


## Unclustered PN, regional connection

jid = '3690515'   # RA=70 ohm-cm
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins, save=True)
plt.show()

fig, axes = myplot.plot_kc_clusters_with_presynaptic_pn(fpath, 'lca', [1, 5, 10, 20])
plt.show()


###########################
## Fri Jun 22 12:00:49 EDT 2018
###########################
## This is for checking clustering of KC activity
jid = '3805935'
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)
fpath = fpaths[0]
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins, save=True)
plt.show()
fig, axes = myplot.plot_kc_clusters_with_presynaptic_pn(fpath, 'lca', [1, 5, 10, 20])
for ax in axes.flat:
    du.despine(ax)
for ax in axes.flatten()[:-1]:
    ax.xaxis.set_visible(False)
fig.set_size_inches(6, 3)
plt.show()



##################################
## Just to confirm that the old simulations did show clustering of KC activity
jid = '62378874'
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)

fig, ax = plt.subplots()
fd = h5.File(fpaths[0], 'r')
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()

jid = '62572477'
fpaths = nda.find_file_by_jid(jid, datadir)
print(fpaths)

fig, ax = plt.subplots()
fd = h5.File(fpaths[0], 'r')
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca')
plt.show()

################
# Mon Jun 25 19:51:16 EDT 2018
# fd = h5.File('d:/tmp/mb_net_UTC2018_06_26__00_38_05-PID11196-JID0.h5', 'r')

# mb_net_UTC2018_06_26__21_06_25-PID10316-JID0.h5 has the following
# config and IPSPs are about 2mV

# ig:
#     filename: D:/subhasis_ggn/model/mb/dclamp_input/ig_spikerate.npy

# ig_ggn_syn:
#     threshold: -20.0mV
#     delay: 1.0ms
#     e: -80.0mV
#     tau1: 1ms
#     tau2: 1ms
#     gmax: 500nS
#     count: 1

fd = h5.File('d:/tmp/mb_net_UTC2018_06_26__21_06_25-PID10316-JID0.h5', 'r')
ggn_vm = fd['/data/uniform/ggn_basal/GGN_basal_Vm']
t = np.arange(ggn_vm.shape[1]) * ggn_vm.attrs['dt']
plt.plot(t, ggn_vm[1, :])
plt.show()
fd.close()

######################
# Wed Jun 27 16:05:12 EDT 2018
# Analysis of simulations with IG IPSPs
jid = '4069523'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
# The spikes seem somewhat clustered. What's going on? KCs are
# clustered by presynaptic GGN segment. This alone is sufficient for
# clustering?
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4069516'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
# ^ although same settings as '4069523' this has more KCs with 2-5 spikes.
fd = h5.File(fpath, 'r')
# print(list(fd['/data/static'].keys()))
# print(list(fd['/data/static/lca_cluster_labels'].keys()))
# labels = fd['/data/static/lca_cluster_labels/lca_cluster_labels']
# print(labels[0])
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()
# ^_ it does not show KC spike clusters.
# v- same as above, but RA=70 ohm cm instead of 150
jid = '4069450'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()
# 4069463 is repeat of above
jid = '4069463'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

# jid = '4066145' has kc->ggn in ca RA=150 ohm-cm
jid = '4066145'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

# jid = '4066141' repeat of above
jid = '4066141'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

#
jid = '4066113'  # alphaL + CA, RA=70 ohm-cm
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

# jid = 4066108
jid = '4066108'  # repeat of 4066113
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4069505'  # alphaL only RA=150 regional
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4069428' # RA = 70, rest as above
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4066103' # RA=70, regional alphaL+CA 
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4066136' # RA=150, regional, alphaL+CA
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4069497'  # RA=150, clustered PN, alphaL
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4066132'  # RA=150, clustered PN, alphaL+CA
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4069497'   # RA=150, clustered PN, alphaL only
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4066098' # regional, clustered pn, alphaL+CA RA=70
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4066084'  # all above + shifting_pn
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4069374'  # no kc->ggn in CA, all else same as above
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

jid = '4069363'  # same as above
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
fd = h5.File(fpath, 'r')
fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
plt.show()
fd.close()

##################################
# Thu Jun 28 11:33:05 EDT 2018
# Check why MCA activity does not have clusters
jid = '4069363'  # same as above
fpath = nda.find_file_by_jid(jid, datadir)[0]
fd = h5.File(fpath, 'r')
mca_cluster_path = '/data/static/mca_cluster_labels/mca_cluster_labels'
ggn_kc_syn_path = '/data/static/ggn_kc_synapse/ggn_kc_synapse'
pn_kc_syn_path = '/data/static/pn_kc_synapse/pn_kc_synapse'
mca_clusters = pd.DataFrame(data=fd[mca_cluster_path]['label', 'sec'][:,0])
ggn_kc_syn = pd.DataFrame(fd[ggn_kc_syn_path]['pre', 'post'][:, 0])
mca_kc = pd.merge(mca_clusters, ggn_kc_syn, left_on='sec',
                  right_on='pre')
pn_kc_syn = pd.DataFrame(data=fd[pn_kc_syn_path]['pre', 'post'][:, 0])
mca_pns = pd.merge(mca_kc, pn_kc_syn, left_on='post', right_on='post')
mca_kc_labels = pd.DataFrame(fd['/data/static/mca_kc_cluster_labels/mca_kc_cluster_labels']['label', 'sec'][:, 0])
print(pn_kc_syn.shape)
print(mca_pns.shape)
fig, ax = plt.subplots()
ii = 0
for pn in mca_pns['pre_y'].unique():
    print(pn)
    pn_st = fd[nda.pn_st_path][pn.decode('utf-8')]
    ax.plot(pn_st, np.ones(pn_st.shape[0]) * ii, '|', ms=2)
    ii += 1
plt.show()
fig, ax = plt.subplots()
ii = 0
for clus, grp in mca_pns.groupby('label'):
    for pn in grp['pre_y'].unique():
        pn_st = fd[nda.pn_st_path][pn.decode('utf-8')]
        ax.plot(pn_st, np.ones(pn_st.shape[0]) * ii, '|', ms=2)
        ii += 1
    break
plt.show()        


lca_cluster_path = '/data/static/lca_cluster_labels/lca_cluster_labels'
lca_clusters = pd.DataFrame(data=fd[lca_cluster_path]['label', 'sec'][:,0])
lca_kc = pd.merge(lca_clusters, ggn_kc_syn, left_on='sec',
                  right_on='pre')
lca_pns = pd.merge(lca_kc, pn_kc_syn, left_on='post', right_on='post')
lca_kc_labels = pd.DataFrame(fd['/data/static/lca_kc_cluster_labels/lca_kc_cluster_labels']['label', 'sec'][:, 0])
print(pn_kc_syn.shape)
print(lca_pns.shape)
fig, ax = plt.subplots()
ii = 0
for pn in lca_pns['pre_y'].unique():
    print(pn)
    pn_st = fd[nda.pn_st_path][pn.decode('utf-8')]
    ax.plot(pn_st, np.ones(pn_st.shape[0]) * ii, '|', ms=2)
    ii += 1
plt.show()
fig, ax = plt.subplots()
ii = 0
for clus, grp in lca_pns.groupby('label'):
    for pn in grp['pre_y'].unique():
        pn_st = fd[nda.pn_st_path][pn.decode('utf-8')]
        ax.plot(pn_st, np.ones(pn_st.shape[0]) * ii, '|', ms=2)
        ii += 1
    break
plt.show()        

fig, ax = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax, fd, 'lca', hline=True)
fig2, ax2 = plt.subplots()
myplot.plot_kc_spikes_by_cluster(ax2, fd, 'mca', hline=True)
plt.show()

fd.close()

#######################################
# Fri Jun 29 11:41:09 EDT 2018
#######################################
# Analyzing simulations with IG, alphaL only, where unresponsive PNs
# spike through stimulus.

jid = '4242763'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()

fd = h5.File(fpath, 'r')

pn_kc_syn = pd.DataFrame(fd[nda.pn_kc_syn_path]['pre', 'post'][:,0])
mca_kc_cluster = defaultdict(list)
for label, kc in fd[nda.mca_kc_cluster_path]['label', 'sec'][:,0]:
    mca_kc_cluster[label].append(kc)
print(mca_kc_cluster[0][0])
print(pn_kc_syn)

pre0 = pn_kc_syn.loc[pn_kc_syn['post'] == mca_kc_cluster[0][0]]['pre']
pre1 = pn_kc_syn.loc[pn_kc_syn['post'] == mca_kc_cluster[0][1]]['pre']
print(len(set(pre0).intersection(set(pre1))))

pre2 = pn_kc_syn.loc[pn_kc_syn['post'] == mca_kc_cluster[1][0]]['pre']
print(len(set(pre0).intersection(set(pre2))))

lca_kc_cluster = defaultdict(list)
for label, kc in fd[nda.lca_kc_cluster_path]['label', 'sec'][:,0]:
    lca_kc_cluster[label].append(kc)

pre0 = pn_kc_syn.loc[pn_kc_syn['post'] == lca_kc_cluster[0][0]]['pre']
pre1 = pn_kc_syn.loc[pn_kc_syn['post'] == lca_kc_cluster[0][1]]['pre']
print(len(set(pre0).intersection(set(pre1))))

pre2 = pn_kc_syn.loc[pn_kc_syn['post'] == lca_kc_cluster[1][0]]['pre']
print(len(set(pre0).intersection(set(pre2))))

fd.close()

## The above had all 415 PNs shared in the same cluster, which does
## not match expectation of about 80% PNs shared between KCs in same
## cluster.
## - later realized that I was not following the shared fraction design for the shifting PNs.
jid = '4069363'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()

fd = h5.File(fpath, 'r')

pn_kc_syn = pd.DataFrame(fd[nda.pn_kc_syn_path]['pre', 'post'][:,0])
mca_kc_cluster = defaultdict(list)
for label, kc in fd[nda.mca_kc_cluster_path]['label', 'sec'][:,0]:
    mca_kc_cluster[label].append(kc)
print(mca_kc_cluster[0][0])
print(pn_kc_syn)

pre0 = pn_kc_syn.loc[pn_kc_syn['post'] == mca_kc_cluster[0][0]]['pre']
pre1 = pn_kc_syn.loc[pn_kc_syn['post'] == mca_kc_cluster[0][1]]['pre']
print(len(set(pre0).intersection(set(pre1))))

pre2 = pn_kc_syn.loc[pn_kc_syn['post'] == mca_kc_cluster[1][0]]['pre']
print(len(set(pre0).intersection(set(pre2))))

lca_kc_cluster = defaultdict(list)
for label, kc in fd[nda.lca_kc_cluster_path]['label', 'sec'][:,0]:
    lca_kc_cluster[label].append(kc)

pre0 = pn_kc_syn.loc[pn_kc_syn['post'] == lca_kc_cluster[0][0]]['pre']
pre1 = pn_kc_syn.loc[pn_kc_syn['post'] == lca_kc_cluster[0][1]]['pre']
print(len(set(pre0).intersection(set(pre1))))

pre2 = pn_kc_syn.loc[pn_kc_syn['post'] == lca_kc_cluster[1][0]]['pre']
print(len(set(pre0).intersection(set(pre2))))

fd.close()

# This showed same issue for MCA, but not for LCA. There was no shifting_pn, just clustered_pn
jid = '4242796'

fpath = nda.find_file_by_jid(jid, datadir)[0]
print(fpath)
print(nda.get_shared_pn_frac(fpath, region='lca', verbose=True))
print(nda.get_shared_pn_frac(fpath, region='mca', verbose=True)) # both have the correct fraction of shared pns

# Check shared PNs between clusters
fd = h5.File(fpath, 'r') 
mca_kc_clus = pd.DataFrame(data=fd[nda.mca_kc_cluster_path]['label', 'sec'][:,0])
pn_kc_syn = pd.DataFrame(data=fd[nda.pn_kc_syn_path]['pre', 'post'][:,0])
    
kcs_0 = mca_kc_clus.loc[mca_kc_clus['label'] == 0]['sec']
kcs_1 = mca_kc_clus.loc[mca_kc_clus['label'] == 1]['sec']
kc0 = np.random.choice(kcs_0, size=1)[0]
kc1 = np.random.choice(kcs_1, size=1)[0]
print(kc0, kc1)
pre0 = pn_kc_syn.loc[pn_kc_syn['post'] == kc0]['pre']
pre1 = pn_kc_syn.loc[pn_kc_syn['post'] == kc1]['pre']
print(len(set(pre0).intersection(pre1)))
intersect = []
for kc1 in kcs_1:
    pre1 = pn_kc_syn.loc[pn_kc_syn['post'] == kc1]['pre']
    intersect.append(len(set(pre0).intersection(pre1)))
    # this consistently showed ~200 shared PNs between clusters
print(max(intersect), min(intersect))    
kc_st_map = fd['/map/event/kc/kc_spiketime']
kc_st_map = pd.DataFrame(data=kc_st_map.value)
fig, ax = plt.subplots()
ii = 0
for label, clus in mca_kc_clus.groupby('label'):
    for index, kc in clus['sec'].iteritems():
        st_ds = kc_st_map.loc[kc_st_map['source'] == kc]['data'].values[0]
        spikes = fd[st_ds].value
        ax.plot(spikes, np.ones(len(spikes)) * ii, 'k,')
        ii += 1
    ax.axhline(ii, color='red', alpha=0.5)
plt.show()    
fd.close()        

## This was a dry run to test connectivity issues
jid = '4299608'
fpath = os.path.join(datadir,
                     'mb_net_UTC2018_06_29__16_56_55-PID40060-JID4299608.h5')
## no shifting pn
jid = '4242796'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()

jid = '4242816'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

########################
# Sat Jun 30 18:44:44 EDT 2018
########################
jid = '4348978'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4349031'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4349034'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
# No spiking

jid = '4349039'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
plt.show()
fpath = nda.find_file_by_jid(jid, datadir)[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
# No spiking

###########################################
# Mon Jul  2 12:32:22 EDT 2018
###########################################
jid = '4448145'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448146'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448147'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448148'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448149'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448150'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448151'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

fig, ax = plt.subplots()
with h5.File(fpath, 'r') as fd:
    x, y = nda.get_event_times(fd[nda.kc_st_path])
    ax.plot(np.concatenate(x), np.concatenate(y), 'k|', ms=3)
    # for ii, ds in enumerate(fd[nda.kc_st_path].values()):
    #     ax.plot(ds.value, np.ones(ds.shape[-1]) + ii, 'k|', ms=3)
plt.show()

fig, ax = plt.subplots()
with h5.File(fpath, 'r') as fd:
    clus = nda.extract_cluster_info(fd, 'mca')
    print(clus)
    ret, x, y = myplot.plot_kc_spikes_by_cluster(ax, fd, 'mca')
plt.show()

print(sum(1 for st in x if len(st) > 0))

jid = '4448152'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))


jid = '4448153'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448154'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448155'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))


jid = '4448156'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448157'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448158'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448159'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448160'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448161'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448162'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448163'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4448164'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

# This was from before changing the options in pn_kc_ggn_network.py
jid = '4439571'  # RA=70 ohm-cm
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4439572'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4439578'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

# Compared 4242763 which spiked a lot with 4439571 which had no
# spike. The former had PN->KC gmax=4.5 pS and GGN->KC gmax = 1nS. The
# latter, 2pS and 0.5pS. Thus PN->KC conductance dominates.
print(nda.get_config_by_jid('4242763', datadir, default_flow_style=False, default_style=''))
print(nda.get_config_by_jid('4439571', datadir, default_flow_style=False, default_style=''))

jid = '4349040'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))



jid = '4242907'  # older simulation - checking if too many KCs have
                 # high firing rate.
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
plt.show()
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

jid = '4242928'  # older simulation - checking if too many KCs have
                 # high firing rate.this is without MCA/LCA separation
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))

plt.show()

jobs = ['4242867', '4242888', '4242907', '4242928']
st_dict = {}
fig, axes = plt.subplots(nrows=len(jobs), ncols=2, sharex='col')
for ii, jid in enumerate(jobs):
    fpath = nda.find_file_by_jid(jid, datadir)
    print(fpath)
    fpath = fpath[0]
    with h5.File(fpath, 'r') as fd:
        lines, st, y = myplot.plot_kc_spikes_by_cluster(axes[ii, 0], fd, 'lca', hline=True)
        st_dict[jid] = st
        sc = [len(spikes) for spikes in st]
        patches, hist, bins = axes[ii, 1].hist(sc, spike_count_bins)
        axes[ii, 1].set_title(fpath)
        print('Plotted', fpath)
fig.set_size_inches(20, 10)
outfile = PdfPages('comparison_of_conn_schemes.pdf')
outfile.savefig(fig)
plt.show()        


###########################################
# Tue Jul  3 12:01:38 EDT 2018
###########################################
# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --pn_kc_clustered_pre --kc_ggn_ca_gmax=2.5pS  --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered | 4529319 |
jid = '4529319'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                        | 4529321 |
jid = '4529321'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

# sbatch slurm/run_mb_net.sh --shifting_pn --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                                          | 4529322 |
jid = '4529322'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --pn_kc_clustered_pre --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional          | 4529326 |
jid = '4529326'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

# Comparing '4529327' with earlier replicate '4448148'
jobs = ['4529326', '4448148']
st_dict = {}
fig, axes = plt.subplots(nrows=1, ncols=2, sharex='all', sharey='all')
for ii, jid in enumerate(jobs):
    fpath = nda.find_file_by_jid(jid, datadir)
    print(fpath)
    fpath = fpath[0]
    with h5.File(fpath, 'r') as fd:
        lines, st, y = myplot.plot_kc_spikes_by_cluster(axes[ii], fd, 'lca', hline=True)
        axes[ii].set_title(fpath)
        print('Plotted', fpath)
fig.set_size_inches(20, 10)
# outfile = PdfPages('comparison_of_conn_schemes.pdf')
# outfile.savefig(fig)
plt.show()        


# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                | 4529327 |
jid = '4529327'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

# sbatch slurm/run_mb_net.sh --shifting_pn --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                                  | 4529333 |
jid = '4529333'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --pn_kc_clustered_pre --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                     | 4529337 |
jid = '4529337'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                           | 4529340 |
jid = '4529340'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

# sbatch slurm/run_mb_net.sh --shifting_pn --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                                             | 4529349 |
jid = '4529349'
fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
fpath = nda.find_file_by_jid(jid, datadir)
print(fpath)
fpath = fpath[0]
fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

# sbatch slurm/run_mb_net.sh --shifting_pn  --pn_kc_clustered --pn_kc_clustered_pre   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                  | 4529350 |
# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                         | 4529351 |
# sbatch slurm/run_mb_net.sh --shifting_pn  --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                                            | 4529352 |
# sbatch slurm/run_mb_net.sh --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                                      | 4529354 |
# sbatch slurm/run_mb_net.sh --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                                                        | 4529440 |
# sbatch slurm/run_mb_net.sh --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                              | 4529441 |
# sbatch slurm/run_mb_net.sh --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                                                | 4529442 |
# sbatch slurm/run_mb_net.sh --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                                         | 4529443 |
# sbatch slurm/run_mb_net.sh --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                                                           | 4529445 |
# sbatch slurm/run_mb_net.sh --pn_kc_clustered   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                                       | 4529446 |
# sbatch slurm/run_mb_net.sh  --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                                                          | 4529447 |
for jid in ['4529350', '4529351', '4529352', '4529354', '4529440', '4529441']
    fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
    fpath = nda.find_file_by_jid(jid, datadir)
    print(fpath)
    fpath = fpath[0]
    fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
    # print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()
# '4529442' cancelled by system for timeout
for jid in ['4529443', '4529445', '4529446', '4529447']:
    fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
    fpath = nda.find_file_by_jid(jid, datadir)
    print(fpath)
    fpath = fpath[0]
    fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
    # print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
plt.show()

for jid in ['4448146', '4448147', '4448149', '4448150', '4448152', '4448153',
            '4448155', '4448156', '4448157', '4448158', '4448159', '4448160',
            '4448161', '4448162', '4448163', '4448164']:
    fpath = nda.find_file_by_jid(jid, datadir)[0]
    try:
        with h5.File(fpath, 'r') as fd:
            try:
                config = yaml.load(fd.attrs['config'].decode())
            except KeyError:
                config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'].value[0])
        with open('{}.yaml'.format(fpath), 'w') as outfile:
            yaml.dump(config, outfile, default_flow_style=False, default_style='')
    except IOError:
        print('JID', jid, 'Could not open file', fpath)


jids = ['4448146', '4448147', '4448149', '4448150', '4448152',
        '4448153', '4448155', '4448156', '4448157', '4448158',
        '4448159', '4448160', '4448161', '4448162', '4448163',
        '4448164', '4529321', '4529322', '4529327', '4529333',
        '4529340', '4529349', '4529351', '4529352', '4529354',
        '4529440', '4529441', '4529442', '4529443', '4529445',
        '4529446', '4529447']

for jid in jids:
    fpaths = nda.find_file_by_jid(jid, datadir)
    for fpath in fpaths:
        if fpath.endswith('h5'):
            break
    try:
        with h5.File(fpath, 'r') as fd:
            try:
                config = yaml.load(fd.attrs['config'].decode())
            except KeyError:
                config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'].value[0])
        with open('{}.yaml'.format(fpath), 'w') as outfile:
            yaml.dump(config, outfile, default_flow_style=False, default_style='')
            print('saved config for {} in {}'.format(jid, outfile.name))            
    except IOError:
        print('JID', jid, 'Could not open file', fpath)
    finally:
        try:
            os.remove(fpath)
            print('Deleted', fpath)
        except IOError:
            print('Could not remove', fpath)




######################################################
# These had 4pS PN->KC gmax and 1 nS GGN->KC gmax


# |-------------------------------------------------------------------------------------------------------------------------------------------------------+---------|
# | sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --pn_kc_clustered_pre --kc_ggn_ca_gmax=2.5pS  --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered | 4553621 |
# | sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                        | 4553622 |
# | sbatch slurm/run_mb_net.sh --shifting_pn --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                                          | 4553623 |
# | sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --pn_kc_clustered_pre --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional          | 4553624 |
# | sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                | 4553625 |
# | sbatch slurm/run_mb_net.sh --shifting_pn --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                                  | 4553626 |
# | sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --pn_kc_clustered_pre --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                     | 4553627 |
# | sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                           | 4553628 |
# | sbatch slurm/run_mb_net.sh --shifting_pn --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                                             | 4553629 |
# | sbatch slurm/run_mb_net.sh --shifting_pn  --pn_kc_clustered --pn_kc_clustered_pre   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                  | 4553630 |
# | sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                         | 4553631 |
# | sbatch slurm/run_mb_net.sh --shifting_pn  --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                                            | 4553632 |
# | sbatch slurm/run_mb_net.sh --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                                      | 4553633 |
# | sbatch slurm/run_mb_net.sh --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                                                        | 4553635 |
# | sbatch slurm/run_mb_net.sh --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                              | 4553636 |
# | sbatch slurm/run_mb_net.sh --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                                                | 4553638 |
# | sbatch slurm/run_mb_net.sh --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                                         | 4553639 |
# | sbatch slurm/run_mb_net.sh --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                                                           | 4553640 |
# | sbatch slurm/run_mb_net.sh --pn_kc_clustered   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                                       | 4553641 |
# | sbatch slurm/run_mb_net.sh  --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                                                          | 4553642 |

jidlist = ['4553621', '4553622', '4553623', '4553624', '4553625',
           '4553626', '4553627', '4553628', '4553629', '4553630',
           '4553631', '4553632', '4553633', '4553635', '4553636',
           '4553638', '4553639', '4553640', '4553641', '4553642']


fig_ax_list = []
for jid in jidlist[:len(jidlist)//2]:
    fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
    fpath = nda.find_file_by_jid(jid, datadir)
    print(fpath)
    fpath = fpath[0]
    fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
    print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
    fig_ax_list += [(fig, ax), (fh, axh)]
plt.show()

fig_ax_list = []
for jid in jidlist[len(jidlist)//2:]:
    fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
    fpath = nda.find_file_by_jid(jid, datadir)
    print(fpath)
    fpath = fpath[0]
    fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
    print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
    fig_ax_list += [(fig, ax), (fh, axh)]
plt.show()


#
# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --pn_kc_clustered_pre --kc_ggn_ca_gmax=2.5pS  --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered | 4557692 |
# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                        | 4557694 |
# sbatch slurm/run_mb_net.sh --shifting_pn --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                                          | 4557696 |
# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --pn_kc_clustered_pre --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional          | 4557697 |
# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                | 4557699 |
# sbatch slurm/run_mb_net.sh --shifting_pn --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                                  | 4557701 |
# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --pn_kc_clustered_pre --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                     | 4557703 |
# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                           | 4557705 |
# sbatch slurm/run_mb_net.sh --shifting_pn --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                                             | 4557707 |
# sbatch slurm/run_mb_net.sh --shifting_pn  --pn_kc_clustered --pn_kc_clustered_pre   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                  | 4557709 |
# sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                         | 4557711 |
# sbatch slurm/run_mb_net.sh --shifting_pn  --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                                            | 4557713 |
# sbatch slurm/run_mb_net.sh --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                                      | 4557714 |
# sbatch slurm/run_mb_net.sh --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --kc_ggn_clustered                                                        | 4557716 |
# sbatch slurm/run_mb_net.sh --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                              | 4557718 |
# sbatch slurm/run_mb_net.sh --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS --regional                                                                | 4557720 |
# sbatch slurm/run_mb_net.sh --pn_kc_clustered --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                                         | 4557721 |
# sbatch slurm/run_mb_net.sh --kc_ggn_ca_gmax=2.5pS --kc_ggn_alpha_gmax=2.5pS                                                                           | 4557723 |
# sbatch slurm/run_mb_net.sh --pn_kc_clustered   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                                       | 4557725 |
# sbatch slurm/run_mb_net.sh  --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=5.0pS                                                                          | 4557727 |

jidlist = ['4557692', '4557694', '4557696', '4557697', '4557699', '4557701',
        '4557703', '4557705', '4557707', '4557709', '4557711', '4557713',
        '4557714', '4557716', '4557718', '4557720', '4557721', '4557723',
        '4557725', '4557727']


fig_ax_list = []
for jid in jidlist:
    fig, ax = myplot.plot_data_by_jobid(jid, datadir, by_cluster=True)
    # fpath = nda.find_file_by_jid(jid, datadir)
    # print(fpath)
    # fpath = fpath[0]
    # fh, axh = myplot.plot_kc_spike_count_hist(fpath, bins=spike_count_bins)
    # print(nda.get_config_by_jid(jid, datadir, default_flow_style=False, default_style=''))
    fig_ax_list += [(fig, ax), (fh, axh)]
plt.show()

jidlist = ['4557694', '4557696', '4557699', '4557701', '4557705', '4557707',
           '4557711', '4557713', '4557714', '4557716', '4557718', '4557720',
           '4557721', '4557723', '4557725', '4557727']



fpaths = [nda.find_file_by_jid(jid, datadir)[0] for jid in jidlist]
for fpath in fpaths:
    assert fpath.endswith('.h5')
    print(fpath)
# Delete the files with no KC spikes after saving config
nda.dump_config_and_delete_datafile(fpaths)

#####################
# Try clustering a subset of KCs based on PNs
fpath = nda.find_file_by_jid('4553621', datadir)
print(fpath)
fpath = fpath[0]
fd = h5.File(fpath, 'r')
kc_st = pd.DataFrame(data=fd['/map/event/kc/kc_spiketime'].value)
keep = np.random.choice(np.arange(len(kc_st)), size=100, replace=False)
kc_st = kc_st.iloc[keep].copy()
print(kc_st)
kc_st.drop_duplicates(subset=['source'], inplace=True)
syninfo = pd.DataFrame(data=fd['/data/static/pn_kc_synapse/pn_kc_synapse'].value[:,0])
presyn_list = []
ii = 0
for index, row in kc_st.iterrows():
    presyn = syninfo[syninfo['post'] == row['source']]['pre']
    presyn_list.append(set(presyn))
    ii += 1
    print('Done', ii, 'of', len(kc_st))
    
distance_matrix = np.ones((len(kc_st), len(kc_st)), dtype=float)
for ii in range(len(kc_st)):
    print(ii)
    for jj in range(ii):
        print('-', jj)
        distance_matrix[ii, jj] = distance_matrix[jj, ii] = len(presyn_list[ii].intersection(presyn_list[jj]))
distance_matrix = 1.0 / distance_matrix
distance_matrix /= distance_matrix.max()
np.fill_diagonal(distance_matrix, 0.0)

from sklearn.cluster import dbscan
core, labels = dbscan(distance_matrix, metric='precomputed')
clus  = defaultdict(list)
for ii, label in enumerate(labels):
    clus[label].append(ii)
# HDBSCAN has import errors both from pip and from conda
fig, ax = plt.subplots()
ii = 0
for label, idxlist in clus.items():
    for idx in idxlist:
        st  = fd[kc_st.iloc[idx]['data']].value
        ax.plot(st, np.ones(len(st)) + ii, 'k|', ms=3)
        ii += 1
plt.show()                 

clus = nda.cluster_kcs_by_pn(fd, 100)
fig, ax = plt.subplots()
ii = 0
for label, rows in clus.items():
    for row in rows:
        st = fd[row['data']].value
        ax.plot(st, np.zeros(len(st)) + ii, 'k|', ms=3)
        
fd.close()

#####################################
# Thu Jul  5 13:27:41 EDT 2018
#####################################

jid = '4553638'
fpath = nda.find_file_by_jid('4553621', datadir)
print(fpath)
fpath = fpath[0]
with h5.File(fpath, 'r') as fd:
    clus = nda.cluster_kcs_by_pn(fd, 100)
    print(clus)
    fig, ax = plt.subplots()
    ii = 0
    for label, dslist in clus.items():
        for row in dslist:
            st = fd[row['data']].value
            ax.plot(st, ii + np.ones(len(st)))
            ii += 1
plt.show()            
    
################################################################################
# Thu Jul 12 10:47:33 EDT 2018
################################################################################
def get_gmax(directory):
    """Go through all HDF5 files in directory and try to extract the gmax
    for the synapses from config"""
    ret = []
    with fname in os.listdir(directory):
        if not fname.endswith('.h5'):
            continue
        with h5.File(os.path.join(directory, fname), 'r') as fd:
            config = None
            try:
                config = yaml.load(fd.attrs['config'].decode())
            except KeyError:
                config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'].value[0])
            ret.append(fname, {'PN->KC': config['pn_kc_synapse']['gmax'],
                    'GGN->KC': config['ggn_kc_synapse']['gmax'],
                    'KC->GGN (CA)': config['kc_ggn_CA_synapse']['gmax'],
                    'KC->GGN (alphaL)': cofig['kc_ggn_alphaL_synapse']['gmax']
            })
    return ret
            
