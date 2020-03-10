# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 13:55:10 2019

@author: Subhasis
"""


import sys
import os
import shutil
import numpy as np
import h5py as h5
import pandas as pd
import yaml
import network_data_analysis as nda


"""Find the highest spiking KCs so we can disconnect them and repeat the 
network model simulation"""

# jids = [
#     22072442,
#     22087964,
#     22087965,
#     22087966,
#     22087967,
#     22087970,
#     22087971,
#     22087972,
#     22087973,
#     ]
# datadir = 'Y:/Subhasis/ggn_model_data/olfactory_network'
datadir = '/data/rays3/ggn/fixed_net'
templatedir = '/data/rays3/ggn/fixed_net_templates'

jids = [
    30184869, 
    30184873,
    30184874,
    30184876,
    30184878,
    30184880,
    30184882,
    30184884,
    30184886,
    30184888,
]


for jid in jids:
    fpath = nda.find_h5_file(jid, datadir)
    forig_path = fpath
    kc_spike_count = []
    with h5.File(fpath, 'r') as fd:
        for kc, spikes in fd[nda.kc_st_path].items():
            kc_spike_count.append((spikes.attrs['source'], len(spikes)))
        try:
            forig_path = fd.attrs['original']
        except KeyError:
            forig_path = fpath
        syn_orig = pd.DataFrame(data=fd[nda.pn_kc_syn_path][:, 0])
        
    kc_spikes = pd.DataFrame(kc_spike_count, columns=['kc', 'spikes'])
    #kc_spikes['spikes'].plot('hist')
    print('# kcs > 5 spikes:', len(kc_spikes[kc_spikes['spikes'] > 5]))
    print('# kcs > 10 spikes:', len(kc_spikes[kc_spikes['spikes'] > 10]))
    print('# spiking kcs:', len(kc_spikes[kc_spikes['spikes'] > 0]))
    orig_0 = set(np.where(syn_orig['gmax'] == 0)[0])
    limits = [5]
    for lim in limits:
        over = kc_spikes[kc_spikes['spikes'] > lim]['kc']        
        print(f'{len(over)} kcs spiked more than {lim} spikes')
        over_syn = pd.merge(syn_orig, over, left_on='post', right_on='kc')
        print('Synapses to > lim kcs', len(over_syn))
        print('Synapse to these kcs set to 0?', len(np.flatnonzero(over_syn['gmax'].values == 0)))
        fname = os.path.basename(fpath).rpartition('.')[0]
        out_fname = f'{orig_fname}_kc{lim}.h5'
        outfile = os.path.join(templatedir, out_fname)
        shutil.copyfile(forig_path, outfile)
        print('Copying data from {} as {} to update using {} KC spiking'.format(forig_path, outfile, fpath))
        with h5.File(outfile, 'a') as ofd:
            syndf = ofd[nda.pn_kc_syn_path]
            # Set the conductances to each KC spiking more than 5 spikes to 0
            changed_syn_count = 0
            for ii, kc in over.iteritems():
                idx = np.where(syn_orig['post'] == kc)[0]
                print('Common 0 syn:', len(set(idx).intersection(orig_0)))
                changed_syn_count += len(idx)
                syndf[idx, 0, 'gmax'] = 0.0
            print('Modified: synapses set to 0 conductance:', changed_syn_count)
        ofd.close()
        print('Original: synapses with 0 condictance:', len(syn_orig[syn_orig['gmax'] == 0.0]))
        
        with h5.File(outfile, 'r') as o2:
            syn_new = pd.DataFrame(data=o2[nda.pn_kc_syn_path][:, 0])
            print('# synapses in updated file', len(syn_new))
            print('# shape of synapse data in updated file', syn_new.shape)
            print('# synapses with 0 conductance', len(syn_new[syn_new['gmax'] == 0.0]))
            assert (len(syn_new[syn_new['gmax'] == 0.0]) - len(syn_orig[syn_orig['gmax'] == 0.0])) == changed_syn_count


# fd.close()

# # Test with some concrete samples
# forig_path = '/data/rays3/ggn/olfactory_network/mb_net_UTC2019_03_09__18_28_19-PID22056-JID22087969.h5'
# f_kc5_template_path = '/data/rays3/ggn/fixed_net_templates/mb_net_UTC2019_03_09__18_28_19-PID22056-JID22087969_kc5.h5'
# f_kc5_path = '/data/rays3/ggn/fixed_net/fixed_net_UTC2019_06_26__20_46_06-PID43476-JID30184880.h5'

# forig = h5.File(forig_path, 'r')
# f_tmp = h5.File(f_kc5_template_path, 'r')
# f_kc5 = h5.File(f_kc5_path, 'r')
# print(f_kc5.attrs['original'])

# orig_spike_count = pd.DataFrame([(st.attrs['source'], len(st)) for st in forig[nda.kc_st_path].values()], columns=['kc', 'spikes'])
# print(orig_spike_count.head())

# tmp_spike_count = pd.DataFrame([(st.attrs['source'], len(st)) for st in f_tmp[nda.kc_st_path].values()], columns=['kc', 'spikes'])
# print(tmp_spike_count.head())

# kc5_spike_count = pd.DataFrame([(st.attrs['source'], len(st)) for st in f_kc5[nda.kc_st_path].values()], columns=['kc', 'spikes'])
# print(kc5_spike_count.head())

# syn_orig = pd.DataFrame(data=forig[nda.pn_kc_syn_path][:, 0])
# print(syn_orig.head())

# syn_tmp = pd.DataFrame(data=f_tmp[nda.pn_kc_syn_path][:, 0])
# print(syn_tmp.head())

# syn_kc5 = pd.DataFrame(data=f_kc5[nda.pn_kc_syn_path][:, 0])
# print(syn_kc5.head())

# high_orig = orig_spike_count[orig_spike_count.spikes > 5]
# high_kc5 = kc5_spike_count[kc5_spike_count.spikes > 5]

# print('Original > 5:', len(high_orig))
# print('KC5 > 5:', len(high_kc5))

# syn_high_orig = pd.merge(high_orig, syn_orig, left_on='kc', right_on='post')
# print('Synapses to high KCs orig', len(syn_high_orig))

# syn_high_kc5 = pd.merge(high_kc5, syn_kc5, left_on='kc', right_on='post')
# print('Synapses to high KCs KC5', len(syn_high_kc5))
# print(len(pd.merge(syn_kc5, high_kc5, right_on='kc', left_on='post')))

# tmp_syn_0 = syn_tmp[syn_tmp.gmax == 0]
# print('0 syn in template', len(tmp_syn_0))

# kc5_syn_0 = syn_kc5[syn_kc5.gmax == 0]
# print('0 syn in kc5', len(kc5_syn_0))

# intersec = pd.merge(kc5_syn_0, syn_high_kc5, left_on=['pre', 'post'], right_on=['pre', 'post'])
# print('Common synapses', len(intersec))

# for fd in [forig, f_tmp, f_kc5]:
#     fd.close()
