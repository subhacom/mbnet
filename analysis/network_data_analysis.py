# network_data_analysis.py --- 
# 
# Filename: network_data_analysis.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Wed Mar  7 17:35:29 2018 (-0500)
# Version: 
# Package-Requires: ()
# Last-Updated: Wed Jun 26 10:17:52 2019 (-0400)
#           By: Subhasis  Ray
#     Update #: 532
# URL: 
# Doc URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Code:
from __future__ import print_function
import sys
sys.path.append('/home/rays3/projects/ggn/morphutils')
import os
import h5py as h5
import numpy as np
import random
from collections import defaultdict
import yaml
import pint
import pandas as pd
from timeit import default_timer as timer
import neurograph as ng
import sklearn.cluster as skc
import sklearn.neighbors as skn
import itertools as it
import pyspike as spk
import warnings


_ur = pint.UnitRegistry()
Q_ = _ur.Quantity

ggn_kc_syn_path = '/data/static/ggn_kc_synapse/ggn_kc_synapse'
pn_kc_syn_path =  '/data/static/pn_kc_synapse/pn_kc_synapse'
kc_st_path = '/data/event/kc/kc_spiketime'
pn_st_path = '/data/event/pn/pn_spiketime'
mca_cluster_path = '/data/static/mca_cluster_labels/mca_cluster_labels'
mca_kc_cluster_path = '/data/static/mca_kc_cluster_labels/mca_kc_cluster_labels'
lca_cluster_path = '/data/static/lca_cluster_labels/lca_cluster_labels'
lca_kc_cluster_path = '/data/static/lca_kc_cluster_labels/lca_kc_cluster_labels'
ggn_alphaL_vm_path = '/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm'
ggn_basal_vm_path = '/data/uniform/ggn_basal/GGN_basal_Vm'
ggn_output_vm_path = '/data/uniform/ggn_output/GGN_output_Vm'

def ggn_sec_to_sid(name):
    sec_name = name.rpartition('.')[-1].partition('[')[0]
    if sec_name.startswith('soma'):
        return 1
    elif sec_name.startswith('dend_'):
        return int(sec_name.partition('_')[-1])
    elif sec_name == 'dend':
        return 3


def get_ggn_vm(fd, region, n=None):
    """Get the Vm at `n` random GGN sections from specified region"""
    if region == 'basal':
        ds = fd[ggn_basal_vm_path]
    elif region.startswith('alpha'):
        ds = fd[ggn_alphaL_vm_path]
    elif (region == 'output') or (region == 'calyx'):
        ds = fd[ggn_output_vm_path]
    t = np.arange(0, ds.shape[1]) * ds.attrs['dt']
    if (n is None) or (n < ds.shape[0]):
        return ds[()], t
    idx = np.random.randint(0, ds.shape[0], n)
    return ds[idx, :], t


def get_ggn_kc_syn_info(fd):
    """Get the GGN->KC syninfo with KC as the index. There is one synapse for each KC, so this makes sense.
    Also, the KCs 1-1 correspond to the clustering data."""    
    ggn_kc_syn =  fd['/data/static/ggn_kc_synapse/ggn_kc_synapse'].value.flatten()
    presec = np.char.decode(ggn_kc_syn['pre'].astype('S'))
    postsec = np.char.decode(ggn_kc_syn['post'].astype('S'))
    syndf = pd.DataFrame(data={'pre': presec, 'post': postsec, 'pos': ggn_kc_syn['prepos'][0]})
    syndf = syndf.set_index('post', drop=False, verify_integrity=True)
    return syndf


def get_kc_ggn_ca_syn_info(fd):
    """Get the KC->GGN synapse info with KC as the index, since there is
    just one such synapse for each KC.
    
    Returns: pandas DataFrame with (pre, post, pos) where pre is the
    presynaptic section name, post is postsynaptic section name, pos
    is postsynaptic segment position (it varies only on the GGN side,
    all KCs connect at soma(0.5)

    """
    kc_ggn_syn = fd['/data/static/kc_ggn_CA_synapse/kc_ggn_CA_synapse'].value.flatten()
    presec = np.char.decode(kc_ggn_syn['pre'].astype('S'))
    postsec = np.char.decode(kc_ggn_syn['post'].astype('S'))
    syndf = pd.DataFrame(data={'pre': presec, 'post': postsec, 'pos': kc_ggn_syn['postpos']})
    syndf = syndf.set_index('pre', drop=True, verify_integrity=True)
    return syndf
    

def extract_cluster_info(fd, ca_side):
    """For simulations prior to 2018-03-07 the cluster labels for KCs were
    not stored separately. Nor did NSDF maintain order of data
    addition. Thus the only info associating KCs with corresponding
    GGN segments and cluster is the cluster labels for the GGN
    segments (section + position) and the GGN->KC synapse
    corresponding to that segment.

    The ordering is maintained in a way because the source dimension scale in
    /map for this dataset has lca_pos_{ii} for integers ii
    corresponding to cluster label indices.

    Mon Jul 2 14:20:03 EDT 2018: realized that this does not work for
    MCA cluster, which start after LCA clusters have been assigned
    KCs. So modified to try {side}_kc_cluster_labels first, then the
    {side}_cluster_labels.


    Returns: DataFrame with (cluster_id, ggn_sec, ggn_sec_pos, x, y,
    z, kc)

    """
    path = '/data/static/{}_kc_cluster_labels/{}_kc_cluster_labels'.format(ca_side, ca_side)
    try:
        labeldf = pd.DataFrame(fd[path].value[:,0])
        labeldf.rename(columns={'sec': 'kc'}, inplace=True)
        return labeldf.set_index('kc', drop=True, verify_integrity=True)
    except KeyError:
        print('Older format, no dataset {}'.format(path))
    if ca_side == 'lca':
        start = 0
    else:
        start = fd['/data/static/lca_cluster_labels/lca_cluster_labels'].shape[-1]
        
    labels = fd['/data/static/{}_cluster_labels/{}_cluster_labels'.format(ca_side, ca_side)]
    sources = labels.dims[0]['source']
    labels = labels.value.flatten()
    # The sources are of the form "lca_pos_{ii}" - extract ii
    source_idx = [int(source.decode().rpartition('_')[-1]) for source in sources]
    label_kc = ['KC[{}].soma'.format(start+ii) for ii in source_idx]
    label_data = [(row['label'], 
                   row['pos'], row['x'], row['y'], row['z'], kc)
                  for row, kc in zip(labels, label_kc)]
    labeldf = pd.DataFrame(data=label_data, columns=('label',
                                                     'pos', 'x',
                                                     'y', 'z', 'kc'))
    labeldf = labeldf.set_index('kc', drop=True, verify_integrity=True)
    return labeldf


def get_kc_st(fd):
    """Get a list of KCs and their spike times"""
    return [(node.attrs['source'].decode('utf-8'), node.value)  \
            for node in fd[kc_st_path].values()]        


def get_event_times(group, nodes=None):
    """Get the spike times under group and the index of the spike train for raster plots.

    group: the HDF5 group under which the spike times are stored.

    nodes: list of str specifying the dataset names under `group`.

    Returns: (spike_x, spike_y) where spike_x is a list of arrays
    containing spike times and spike_y is a list of arrays of the same
    size containing the index of that spike train.

    """
    start = timer()
    spike_x = []
    spike_y = []
    if nodes is None:
        nodes = group
    for ii, node in enumerate(nodes):
        st = group[node][:]
        spike_x.append(st)
        sy = np.zeros(st.shape)
        sy[:] = ii
        spike_y.append(sy)
    end = timer()
    print('get_event_times: {}s'.format(end - start))
    return spike_x, spike_y


def get_kc_event_node_map(fd):
    """Return a dict of kc section<-> spike train data node name"""
    grp = fd[kc_st_path]
    ret = {}
    for node in grp.values():
        ret[node.attrs['source']] = node.name
    return ret

def get_kcs_by_region(fd, region):
    """Returns a set containing the section names of kcs in the region
    ('lca' or 'mca' or 'all').

    fd (h5py.File) is the open file handle.

    """
    ggn_kc_syn = get_ggn_kc_syn_info(fd)
    # The dendrites are named dend_{n}[iii] where n is the numerical
    # type in the SWC file of the GGN morphology
    if region == 'all':
        grp = fd[kc_st_path]
        matching_kcs = [grp[node].attrs['source'] for node in grp]
    else:
        dend_pattern = 'dend_{}'.format(ng.custom_types[region.upper()])
        matching_kcs = ggn_kc_syn.post[ggn_kc_syn.pre.str.find(dend_pattern) >= 0]
    return set(matching_kcs)


def get_kc_spike_nodes_by_region(fd, region):
    """Select dataset names for KC spiketimes where the KC belongs to the
    specified region ('lca' or 'mca' or 'all') of the calyx"""
    start = timer()
    try:
        kc_st = pd.DataFrame(data=fd['/map/event/kc/kc_spiketime'].value)
        if region == 'all':
            kc_spike_nodes = [fd[node].name for index, node in kc_st['data'].iterrows()]
        else:
            kcs = pd.DataFrame(
                data=fd['/data/static/{}_cluster_labels/{}_cluster_labels'.format(
                    region, region)].value['sec'])
            kc_st = pd.merge(kcs, kc_st, left_on='sec', right_on='source')
            kc_spike_nodes = [fd[node].name for index, node in kc_st['data'].iterrows()]            
    except KeyError:
        # Find KCs that are postsynaptic to a GGN dendrite corresponding to region
        matching_kcs = get_kcs_by_region(fd, region)
        ## Now retrieve the map between KC name and spike train dataset
        kc_grp = fd['/data/event/kc/kc_spiketime']
        kc_spike_nodes = [node for node in kc_grp
                          if kc_grp[node].attrs['source'].decode()
                          in matching_kcs]
    end = timer()
    print('get_kc_spike_nodes_by_region: {} s'.format(end - start))
    return kc_spike_nodes


def get_kc_vm_idx_by_region(fd, region):
    """Obtain the row indices in 2D Vm array for KCs in `region`. Returns
    a list containing (KC name, row index)"""
    kc_vm_node = fd['/data/uniform/kc/KC_Vm']
    kcs = get_kcs_by_region(fd, region)    
    kc_vm_name_index = {name.decode(): ii for ii, name in enumerate(kc_vm_node.dims[0]['source'])}
    match = kcs.intersection(kc_vm_name_index.keys())
    print(region, list(kcs)[0], list(kc_vm_name_index.keys())[0], len(match))
    return [(name, kc_vm_name_index[name]) for name in match]


def get_shared_pn_frac(filename, region='lca', verbose=False):
    frac = {}
    with h5.File(filename, 'r') as fd:
        config = yaml.load(fd.attrs['config'].decode())
        shared_pn_frac = config['kc']['shared_pn_frac']
        kc_clus = fd['/data/static/{}_kc_cluster_labels/{}_kc_cluster_labels'.format(region, region)]
        kc_clus = pd.DataFrame(data=kc_clus.value.flatten())
        pn_kc_syn = fd['/data/static/pn_kc_synapse/pn_kc_synapse']
        pn_kc_syn = pd.DataFrame(data=pn_kc_syn.value.flatten())
        for label, clus in kc_clus[['label', 'sec']].groupby('label'):
            kcs = random.sample(list(clus['sec']), 2)
            pn1 = set(pn_kc_syn.pre[pn_kc_syn.post == kcs[0]])
            pn2 = set(pn_kc_syn.pre[pn_kc_syn.post == kcs[1]])
            shared = pn1.intersection(pn2)
            frac[label] = 1.0 * len(shared) / len(pn1)
            if verbose:
                print('cluster', label, ': shared pn frac:', frac[label])
        cfrac = config['kc']['shared_pn_frac']
        msdiff = (np.array(frac.values()) - cfrac)**2 / len(frac)
        if verbose:
            print('mean square difference with config', msdiff)
    return frac


def get_files_by_jid(datadir):
    """Return a dataframe with filename as fname and jobid as jid
    columns"""
    jid_fmap = {}
    for fname in os.listdir(datadir):
        # File names are 'xxx-JID{jobid}.h5'
        if fname.endswith('h5'):
            jid = fname.rpartition('-JID')[-1][:-3]
            jid_fmap[fname] = jid
    return pd.DataFrame(data={'fname': list(jid_fmap.keys()),
                              'jid': list(jid_fmap.values())})

def find_file_by_jid(jid, datadir, verbose=False):
    flist = os.listdir(datadir)
    match = [fname for fname in flist if 'JID{}'.format(jid) in fname]
    if verbose and (len(match) > 1):
        print('{} files with same jobid.', len(match), match)
    return [os.path.join(datadir, fname) for fname in match]


def find_h5_file(jid, datadir):
    flist = find_file_by_jid(jid, datadir)
    match = [fname for fname in flist if fname.endswith('.h5')]
    if len(match) > 1:
        raise Exception('Two HDF5 files with same JID: {}'.format(jid))
    elif len(match) == 0:
        raise Exception('No match for jid {} in {}'.format(jid, datadir))
    return match[0]


def find_spike_distance_file(jid, datadir):
    dfpath = find_h5_file(jid, datadir)
    sfpath = '{}.spike_distance.npz'.format(dfpath.partition('.h5')[0])
    if not os.path.exists(sfpath):
        raise Exception('No spike distance file for {}'.format(dfpath))
    return sfpath


def get_config_by_jid(jid, datadir, **kwargs):
    """Example:
        
        get_config_by_jid('64819232', datadir, default_flow_style=False, default_style='' )
        """
    fname = find_file_by_jid(jid, datadir)[0]
    with h5.File(fname, 'r') as fd:
        try:
            return yaml.dump(yaml.load(fd.attrs['config'].decode()), **kwargs)
        except KeyError:
            return yaml.dump(yaml.load(fd['/model/filecontents/mb/network/config.yaml'].value[0]), **kwargs)


def select_good_files(datadir, simlistfile):
    finfo = pd.read_csv(simlistfile)
    finfo.jobid = finfo.jobid.apply(str)
    flist = get_files_by_jid(datadir)
    combo = pd.merge(finfo, flist, left_on=['jobid'], right_on=['jid'])
    return combo


def get_kc_gmax_inh(fd):
    """Extract inhibitory conductance on KCs"""
    return fd[ggn_kc_syn_path]['post', 'gbar']


def get_kc_gmax_exc(fd):
    """Extract excitatory conductances on KCs. No attempt is made to
    separate out inputs from different PNs

    """
    return fd[pn_kc_syn_path]['pre', 'gmax']


def kc_spiking_correlates(fname):
    """Correlate the spiking rate in KC during the PN activity window and
    with the inhibitory and excitatory synaptic conductances"""
    with h5.File(fname, 'r') as fd:
        ggn_kc_gbar = pd.DataFrame(fd[ggn_kc_syn_path]['post', 'gbar'][:, 0])
        pn_kc_gmax = pd.DataFrame(fd[pn_kc_syn_path]['post', 'gmax'][:, 0])
        pn_kc_gmax_by_kc = pn_kc_gmax.groupby('post').sum().reset_index()
        conf = yaml.load(fd.attrs['config'].decode())
        delt = Q_(conf['stimulus']['duration']).to('ms').m
        try:
            delt += Q_(conf['pn']['offdur']).to('ms').m
        except KeyError:
            pass
        print('PN excited for', delt, 'ms')
        kc_spike_count = [(kc, len(st)) for kc, st in get_kc_st(fd)]
        kc_spike_count = pd.DataFrame(data=kc_spike_count, columns=['kc', 'spikes'])
        kc_spike_rate = pd.DataFrame({'kc': kc_spike_count['kc'],
                                      'spikerate': kc_spike_count['spikes'] / delt})
        pn_ggn = pd.merge(pn_kc_gmax_by_kc, ggn_kc_gbar, on='post')
        pn_ggn['post'] = pn_ggn['post'].str.decode('utf-8')
        pn_ggn_kc_rate = pd.merge(pn_ggn, kc_spike_rate, left_on='post', right_on='kc')
        pn_ggn_kc_rate.drop('post', axis=1, inplace=True)
        print(pn_ggn_kc_rate)        
        return pn_ggn_kc_rate.rename(columns={'gbar': 'ginh', 'gmax': 'gexc'})


def dump_config_and_delete_datafile(filenames):
    for fname in filenames:
        try:
            with h5.File(fname, 'r') as fd:
                try:
                    config = yaml.load(fd.attrs['config'].decode())
                except KeyError:
                    config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'].value[0])
            with open('{}.yaml'.format(fname), 'w') as outfile:
                yaml.dump(config, outfile, default_flow_style=False, default_style='')
                print('saved config for {} in {}'.format(fname, outfile.name))            
                try:
                    os.remove(fname)
                    print('Deleted', fname)
                except IOError:
                    print('Could not remove', fname)
        except IOError:
            print('Could not open file', fname)

            
def get_kc_similarity(fd, nkc):
    tstart = timer()
    kc_st_ds = pd.DataFrame(data=fd['/map/event/kc/kc_spiketime'].value)
    keep = np.random.choice(np.arange(len(kc_st_ds)), size=nkc, replace=False)
    kc_st_ds = kc_st_ds.iloc[keep].copy()
    syninfo = pd.DataFrame(fd['/data/static/pn_kc_synapse/pn_kc_synapse'].value[:, 0])
    tend = timer()
    print('Read syninfo in {} s'.format(tend - tstart))
    sys.stdout.flush()
    presyn_list = []
    ii = 0
    tstart = timer()
    kcgroups = syninfo.groupby('post')
    for index, row in kc_st_ds.iterrows():
        start = timer()
        presyn = kcgroups.get_group(row['source'])['pre']
        presyn_list.append(set(presyn))
        ii += 1
        end = timer()
        # print('Done {} of {} in {} s'.format(ii, len(kc_st_ds), end - start))
        # sys.stdout.flush()
    tend = timer()
    print('Finished finding presynaptic PNs for {} KCs in {} s'.format(nkc, tend - tstart))
    sys.stdout.flush()
    tstart = timer()
    similarity_matrix = np.zeros((len(kc_st_ds), len(kc_st_ds)), dtype=float)
    for ii in range(len(kc_st_ds)):
        for jj in range(ii+1):
            similarity_matrix[ii, jj] = similarity_matrix[jj, ii] = \
                len(presyn_list[ii].intersection(presyn_list[jj]))
    return similarity_matrix


def cluster_kcs_by_pn(fd, nkc, pn_per_kc=415, eps=0.5/415, min_samples=5):
    """Try to cluster KCs by how many PNs they share.

    Use DBSCAN algorithm from scikits-learn.
    """
    tstart = timer()
    similarity_matrix = get_kc_similarity(fd, nkc)
    distance_matrix = 1 / similarity_matrix - 1.0 / pn_per_kc
    tend = timer()
    print('Computing distance matrix: Simple loop {} s'.format(tend - tstart))
    tstart = timer()
    core, labels = skc.dbscan(distance_matrix, metric='precomputed', eps=eps, min_samples=min_samples)
    tend = timer()
    print('DBSCAN done in {} s.'.format(tend - tstart))
    sys.stdout.flush()
    clus = defaultdict(list)
    kc_st_ds = pd.DataFrame(data=fd['/map/event/kc/kc_spiketime'].value)
    for ii, label in enumerate(labels):
        clus[label].append(kc_st_ds.iloc[ii])
    return clus, distance_matrix


def load_config(fd):
    """Load and return the config from open file fd as an YAML structure"""
    config = None
    try:
        config = yaml.load(fd.attrs['config'].decode())
    except KeyError:
        config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'].value[0])
    return config


def get_simtime(fd):
    config = load_config(fd)
    simtime = Q_(config['stimulus']['onset']).to('ms').m +  \
        Q_(config['stimulus']['duration']).to('ms').m + \
        Q_(config['stimulus']['tail']).to('ms').m
    return simtime


def get_stimtime(fd):
    """Returns the timing parameters (in ms) for stimulus,

    onset: time when stimulus started,
    duration: interval for which stimulus was on,
    offdur: duration of OFF response"""
    config = load_config(fd)
    return {'onset': Q_(config['stimulus']['onset']).to('ms').m,
            'duration': Q_(config['stimulus']['duration']).to('ms').m,
            'offdur': Q_(config['pn']['offdur']).to('ms').m}
    


def kde_pdf(spike_times, bandwidth=50.0, xgrid=None, kernel='gaussian'):
    """Compute the probability density function using KernelDensity
    estimation with specified bandwidth and evaluated at positions in
    xgrid.

    Return (pdf, xgrid)

    """
    if len(spike_times) == 0:
        warnings.warn('No spikes in spike trains')
        return (None, None)
    kde = skn.KernelDensity(kernel=kernel, bandwidth=bandwidth)
    kde.fit(spike_times[:, np.newaxis])
    if xgrid is None:
        xgrid = np.arange(min(spike_times), max(spike_times), bandwidth/2.0)
    log_pdf = kde.score_samples(xgrid[:, np.newaxis])    
    pdf = np.exp(log_pdf)
    return pdf, xgrid


def get_kc_clusters(fd):
    myfd = fd
    try:
        lca_clus_info = myfd[lca_kc_cluster_path]['label', 'sec'][:, 0]
    except KeyError:
        myfd = h5.File(fd.attrs['original'], 'r')
        lca_clus_info = myfd[lca_kc_cluster_path]['label', 'sec'][:, 0]
    mca_clus_info = myfd[mca_kc_cluster_path]['label', 'sec'][:, 0]
    if myfd != fd:
        myfd.close()
    lca_kcs = pd.DataFrame(data=lca_clus_info)
    mca_kcs = pd.DataFrame(data=mca_clus_info)
    mca_kcs['cluster'] = ['mca_{}'.format(label) for label in mca_kcs.label]
    lca_kcs['cluster'] = ['lca_{}'.format(label) for label in lca_kcs.label]
    kc_clusters = pd.concat([lca_kcs, mca_kcs], ignore_index=True)
    kc_clusters.set_index('sec', drop=True, verify_integrity=True)
    return kc_clusters


def get_spiking_kcs(fd):
    return [node.attrs['source'].decode('utf-8') for node in fd[kc_st_path].values() if len(node) > 0]
    

def cluster_spike_trains(spike_trains, tstart, tend, interval=None, eps=0.01, measure='SPIKE_distance'):
    """Cluster a list of spike trains by SPIKE distance measure.

    tstart and tend are the start and end time of recording (used for creating SpikeTrain object from PySpike.

    interval is (t0, t1) the time interval during which to compute
    SPIKE_distance. If None, (tstart, tend) is used.

    eps: epsilon parameter for DBSCAN algorithm. This is the maximum
    distance between two samples for them to be considered in the same
    neighborhood.

    All spike trains should be nonempty.

    Return (cluster_info, spike_distance_matrix)
    where
    cluster_info: the result of applying DBSCAN. This gives 
    spike_distance_matrix: an NxN distance matrix for N spike trains.

    """
    if interval is None:
        interval = (tstart, tend)
    st_list = [spk.SpikeTrain(st, (tstart, tend)) for st in spike_trains]
    print('Interval:', interval)
    print('tstart: {} tend: {}'.format(tstart, tend))
    print('Number of spike trains:', len(st_list))
    
    dist = spk.spike_distance_matrix(st_list, interval=interval)
    clus = skc.DBSCAN(eps=eps, metric='precomputed').fit(dist)
    return clus, dist
    

# 
# network_data_analysis.py ends here
