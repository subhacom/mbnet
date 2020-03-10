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
import argparse
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

DATADIR = '/data/rays3/ggn/fixed_net'
TEMPLATEDIR = '/data/rays3/ggn/fixed_net_templates'

def make_parser():
    parser = argparse.ArgumentParser(description='Disconnect high firing KCs from data file to create a template network')
    parser.add_argument('--limit', type=int, default=5, help='Upper limit of allowed KC spike count')
    parser.add_argument('--sdir', type=str, help='source directory')
    parser.add_argument('--tdir', type=str, help='target directory')
    parser.add_argument('--jid', type=str, help='JID of source dataset')
    return parser


def remove_high_firing_kcs(jid, limit, sdir, tdir):
    """Remove KCs firing more than `limit` from dataset of `jid`. `sdir`
    points to directory containing the data file, `tdir` is where the
    output file will be written. If source file is called `x.h5`,
    output file will be `x_kc{limit}.h5`

    """
    fpath = nda.find_h5_file(jid, sdir)
    kc_spike_count = []
    with h5.File(fpath, 'r') as fd:
        for kc, spikes in fd[nda.kc_st_path].items():
            kc_spike_count.append((spikes.attrs['source'], len(spikes)))
        try:
            forig_path = fd.attrs['original']
            print(jid, ': original template:', forig_path)
        except KeyError:
            forig_path = fpath
        syn = fd[nda.pn_kc_syn_path]
        if len(syn.shape) == 2:
            syn_orig = pd.DataFrame(data=syn[:, 0])
        else:
            syn_orig = pd.DataFrame(data=syn[:])
    kc_spikes = pd.DataFrame(kc_spike_count, columns=['kc', 'spikes'])
    print('# kcs > 5 spikes:', len(kc_spikes[kc_spikes['spikes'] > 5]))
    print('# kcs > 10 spikes:', len(kc_spikes[kc_spikes['spikes'] > 10]))
    print('# spiking kcs:', len(kc_spikes[kc_spikes['spikes'] > 0]))
    orig_0 = set(np.where(syn_orig['gmax'] == 0)[0])
    over = kc_spikes[kc_spikes['spikes'] > limit]
    print('{} kcs spiked more than {} spikes'.format(len(over), limit))
    over_syn = pd.merge(syn_orig, over, left_on='post', right_on='kc')
    print('Synapses to > limit kcs', len(over_syn))
    print('Synapse to these kcs set to 0?', len(np.flatnonzero(over_syn['gmax'].values == 0)))
    # This is tricky - a simulation based on a template used to
    # generate a datafile with external reference to the synapse
    # datasets in the template. So attempt to generate another
    # template by updating synapses in the produced datafile
    # referred back to the original template.
    fname = os.path.basename(fpath).rpartition('.')[0]
    out_fname = '{}_kc{}.h5'.format(fname, limit)
    outfile = os.path.join(tdir, out_fname)
    if os.path.exists(outfile):
        print(f'File already exists: {outfile}')
        return 0    
    print('Copying data from {} as {} to update using {} KC spiking'.format(forig_path, outfile, fpath))
    shutil.copyfile(forig_path, outfile)
    print('Disabling PN synapses to KCs firing > {} spikes'.format(limit))
    changed_syn_count = 0
    with h5.File(outfile, 'a') as ofd:
        syndf = ofd[nda.pn_kc_syn_path]
        # Set the conductances to each KC spiking more than 5 spikes to 0
        for row in over.itertuples():
            idx = np.where(syn_orig['post'] == row.kc)[0]
            # print('Common 0 syn:', len(set(idx).intersection(orig_0)))
            changed_syn_count += len(idx)
            syndf[idx, 0, 'gmax'] = 0.0
        print('Modified: synapses set to 0 conductance:', changed_syn_count)
    ofd.close()
    print('Original: synapses with 0 conductance:', len(syn_orig[syn_orig['gmax'] == 0.0]))

    with h5.File(outfile, 'r') as o2:
        syn_new = pd.DataFrame(data=o2[nda.pn_kc_syn_path][:, 0])
        print('# synapses in updated file', len(syn_new))
        print('# shape of synapse data in updated file', syn_new.shape)
        print('# synapses with 0 conductance', len(syn_new[syn_new['gmax'] == 0.0]))
        assert (len(syn_new[syn_new['gmax'] == 0.0]) - len(syn_orig[syn_orig['gmax'] == 0.0])) == changed_syn_count
        print('Finished checking updated file')
    return changed_syn_count


if __name__ == '__main__':
    parser = make_parser()
    args = parser.parse_args()
    changed_syn_count = remove_high_firing_kcs(args.jid, args.limit, args.sdir, args.tdir)
    if changed_syn_count == 0:
        sys.exit(1)

    
    
