#!/usr/bin/env python
# fixed_network_changing_stim.py ---
# Author: Subhasis  Ray
# Created: Wed Sep 19 17:28:00 2018 (-0400)
# Last-Updated: Tue Oct  1 17:35:39 2019 (-0400)
#           By: Subhasis  Ray
# Version: $Id$

# Code:
"""This script is for simulating a reduced and fixed Mushroom body
network with connectivity info specified in input files.

It should take the base spike trains for PNs from input file and add some
random noise to represent noisy signal. It can be used for simulating
multiple trials of the same stimulus. The identity of the
activated/inhibited PNs will represent the odor identity in a manner
similar to Assisi, et al, 2007.

"""
import sys
import os
from datetime import datetime
import argparse
import numpy as np
import h5py as h5
import yaml
from timeit import default_timer as timer
from neuron import h
import nsdf

import config as cfg
import nrnutils as nu
import neurograph as ng
import ephys
from ig import IGIzhi

h.load_file('stdrun.hoc')

ggn_file = os.path.join(os.environ['HOME'], 'projects', 'ggn', 'mb',
                        'cell_templates', 'GGN_20170309_sc.hoc')
kc_file = os.path.join(os.environ['HOME'], 'projects', 'ggn', 'mb',
                       'cell_templates', 'kc_1_comp.hoc')
kc_name = 'KC'
ggn_name = 'GGN_20170309_sc'

# Load the cell templates
if not hasattr(h, ggn_name):
    h.xopen(ggn_file)
if not hasattr(h, kc_name):
    h.xopen(kc_file)

kc_st_path = '/data/event/kc/kc_spiketime'
pn_st_path = '/data/event/pn/pn_spiketime'
pn_kc_syn_path = '/data/static/pn_kc_synapse/pn_kc_synapse'
ggn_kc_syn_path = '/data/static/ggn_kc_synapse/ggn_kc_synapse'
kc_ggn_alphaL_syn_path = '/data/static/kc_ggn_alphaL_synapse/kc_ggn_alphaL_synapse'
kc_model_path = '/model/modeltree/olf/kc'
kc_ig_syn_path = '/data/static/kc_ig_synapse/kc_ig_synapse'

model_dict = {}   # global to keep refs of created objects
kc_name_sec_dict = {}   # global to keep refs kc sections by name


def dither_spiketrain(st_list, cell_shift=0, dither=5.0):
    """Dither the spike trains with uniform random time between 0 and
    `dither` (so each spike is shifted randomly in either direction by
    0 to dither time). In addition the spike trains are assigned to a
    the cells in cell_list shifted by cell_shift.

    st_list: list of spike trains for.

    Return modified list of spike trains, dithered and rotated.

    """
    if (cell_shift == 0) and (dither == 0):
        return st_list
    altered = []
    for ii, st in enumerate(st_list):
        if dither > 0:
            newst = np.array(st) + (np.random.sample(len(st)) - 0.5) * 2 * dither
            if np.any(newst < 0):
                cfg.logger.info('Removing negative spike time from PN {}'.format(ii))
            newst = np.sort(newst[newst > 0])  # This is important for NetCon
        else:
            newst = st
        altered.append(newst)
    return altered[cell_shift:] + altered[:cell_shift]


def create_ggn():
    """Create a singleton GGN"""
    global model_dict
    cfg.logger.info('Started creating GGN')
    for obj in dir():
        ## Weird, but NEURON hoc.HocObject is not a type object though type(h)
        ## prints hoc.HocObject
        if type(obj) == type(h):  
            if obj.hname().startswith(ggn_name):
                return obj
    cfg.logger.info('Finished creating GGN')
    ggn = eval('h.{}()'.format(ggn_name))
    model_dict['ggn'] = ggn
    return ggn


def load_pn_spikes(fd, pn_st_path=pn_st_path):
    """Read PN spike times from HDF5 file descriptor fd.

    Returns a list of tuples mapping PN names to spike trains sorted
    in order of PN index (assuming PN names are of PN_{index} form).

    """
    pn_st = {}
    for pn, st in fd[pn_st_path].items():
        pn_st[pn] = st.value
    return sorted(pn_st.items(), key=lambda x: int(x[0].split('_')[-1]))        


def create_pn_output(spike_trains):
    """Create stimvec and vecstim for proxy of PN spikes"""
    global model_dict
    cfg.logger.info('Started creating PN output vecstims and stimvecs')
    stimvecs = []
    vecstims = []
    for st in spike_trains:
        stimvecs.append(h.Vector(st))
        vecstims.append(h.VecStim())
        vecstims[-1].play(stimvecs[-1])
    model_dict['vecstim'] = vecstims
    model_dict['stimvec'] = stimvecs
    cfg.logger.info('Finished creating PN output vecstims and stimvecs')
    return stimvecs, vecstims


def create_kcs(fd, cellname='KC', config={}):
    """Create KCs from data in fd. Returns a list of hoc KC objects."""
    global model_dict
    cfg.logger.info('Started KC creation')
    try:
        nkc = len(fd[kc_model_path])
    except KeyError:  # some generated datafiles do not store full model info
        nkc = config['number']
        
    cellclass = eval('h.{}'.format(cellname))
    cfg.logger.info('Finished KC creation')
    kcs = [cellclass() for ii in range(nkc)]
    model_dict['kc'] = kcs
    return kcs


def create_pn_kc_conn(fd, pn_vecstim_dict, kc_name_sec_dict,
                      path=pn_kc_syn_path, config=None):
    """Create PN->KC connection from dataset at specified path

    pn_vecstim_list should be a dict from pn name to corresponding
    vecstim

    kc_name_sec_dict should be a dict from kc soma section names to
    the sections themselves

    returns a dict containing kc names mapped to their synapses and a
    list of all netcons from PNs to KCs.

    """
    global model_dict
    cfg.logger.info('Started PN->KC connections')
    if config is None:
        config = {}
    syn_dict = {}
    netcons = []
    ncdict = {}
    syn_data = fd[path]
    # This is to take care of discrepancy between NSDF saving synapses
    # as dataset of dim [nx1] vs my kc removal code saving them as dim
    # [n]
    if len(syn_data.shape) == 2:
        syn_data = syn_data[:, 0]
    for ii, row in enumerate(syn_data):
        if row['post'] not in syn_dict:
            target = kc_name_sec_dict[row['post']]
            syn = h.Exp2Syn(target(row['postpos']))
            syn.e = row['e']
            syn.tau1 = row['tau1']
            syn.tau2 = row['tau2']
            syn_dict[row['post']] = syn
        else:
            syn = syn_dict[row['post']]
        vecstim = pn_vecstim_dict[row['pre']]
        if 'delay' in row.dtype.names:
            delay = row['delay']
        else:
            delay = config.get('delay', '0.0ms')
            delay = cfg.Q_(delay).to('ms').m
        thresh = -20.0            
        nc = h.NetCon(vecstim, syn, thresh, delay, row['gmax'])
        netcons.append((row['pre'], row['post'], nc))
        assert (vecstim, syn) not in ncdict
        ncdict[(vecstim, syn)] = nc
        # print(ii, nc, nc.pre(), nc.syn())
        # sys.stdout.flush()
    model_dict['pn_kc_conn'] = (netcons, syn_dict)
    cfg.logger.info('Finished PN->KC connections')
    return netcons, syn_dict


def create_ggn_kc_conn(fd, kc_name_sec_dict, ggn_name_sec_dict,
                       path=ggn_kc_syn_path, config=None):
    """Create graded synapses from GGN sections to KC based on synaptic
    connection data in fd at path `ggn_kc_syn_path`."""
    global model_dict
    cfg.logger.info('Started GGN->KC connection')
    model_dict['ggn_kc_conn'] = []
    syn_dict = {}
    syn_data = fd[path]
    if len(syn_data.shape) == 2:
        syn_data = syn_data[:, 0]
    for row in syn_data:
        kc = kc_name_sec_dict[row['post']]
        syn = h.GradedSyn(kc(row['postpos']))
        syn.vmid = row['vmid']
        syn.vslope = row['vslope']
        syn.e = row['e']
        syn.gbar = row['gbar']
        syn.tau = row['tau']
        h.setpointer(ggn_name_sec_dict[row['pre']](row['prepos'])._ref_v,
                     'vpre', syn)
        model_dict['ggn_kc_conn'].append({'pre': row['pre'],
                                          'prepos': row['prepos'],
                                          'post': row['post'],
                                          'postpos': row['postpos'],                                          
                                          'vmid': syn.vmid,
                                          'vslope': syn.vslope,
                                          'e': syn.e,
                                          'gbar': syn.gbar,
                                          'tau': syn.tau})
        syn_dict[row['post']] = syn
        
    model_dict['ggn_kc_syn'] = syn_dict
    cfg.logger.info('Finished GGN->KC connection')
    return syn_dict


def create_kc_ggn_conn(fd, kc_name_sec_dict, ggn_name_sec_dict,
                       path=kc_ggn_alphaL_syn_path, config=None):
    """Created excitatory synapses from KCs to GGN segments"""
    global model_dict
    cfg.logger.info('Started KC->GGN connection')
    if config is None:
        config = {}
    model_dict['kc_ggn_syn'] = []
    model_dict['kc_ggn_nc'] = []
    model_dict['kc_ggn_conn'] = []
    syn_dict = {}
    nc_dict = {}
    syn_data = fd[path]
    if len(syn_data.shape) == 2:
        syn_data = syn_data[:, 0]
    # for row_arr in fd[kc_ggn_alphaL_syn_path]:
        # row = row_arr[0]
    for row in syn_data:
        postsec = ggn_name_sec_dict[row['post']]
        syn = h.Exp2Syn(postsec(row['postpos']))
        syn.e = row['e']
        syn.tau1 = row['tau1']
        syn.tau2 = row['tau2']
        presec = kc_name_sec_dict[row['pre']]
        if 'delay' in row.dtype.names:
            delay = row['delay']
        else:
            delay = config.get('delay', '1.0ms')
            delay = cfg.Q_(delay).to('ms').m
        nc = h.NetCon(presec(row['prepos'])._ref_v, syn, -20.0,
                      delay, row['gmax'], sec=presec)
        syn_dict[row['pre']] = syn
        nc_dict[row['pre']] = nc
        model_dict['kc_ggn_syn'].append(syn)
        model_dict['kc_ggn_nc'].append(nc)
        model_dict['kc_ggn_conn'].append({'pre': row['pre'], 'prepos':
                                          row['prepos'], 'post':
                                          row['post'], 'postpos':
                                          row['postpos'], 'nc': nc,
                                          'syn': syn})
    cfg.logger.info('Finished KC->GGN connection')
    return syn_dict, nc_dict


def setup_ig(fd,  ggn_name_sec_dict, config=None):
    """Skipping model structure and just create IG from config"""
    global model_dict
    if config is None:
        config = {}
    if 'ig' not in config:
        return    
    ig = IGIzhi(inject=config['ig']['inject'],
                gbarGABA=cfg.Q_(config['ggn_ig_syn']['gmax']).to('uS').m)
    ig.vmidGraded = cfg.Q_(config['ggn_ig_syn']['vmid']).to('mV').m
    ig.vslopeGraded = cfg.Q_(config['ggn_ig_syn']['vslope']).to('mV').m
    ggn =  model_dict['ggn']
    # presec = ggn_name_sec_dict[config['ggn_ig_syn']['source']]
    presec = eval('ggn.{}'.format(config['ggn_ig_syn']['source']))
    h.setpointer(presec(0.5)._ref_v, 'vpreGraded', ig.izhi)
    cfg.logger.info('Created IG with injection {}'.format(config['ig']['inject']))
    model_dict['ig'] = ig
    # post_sec = ggn_name_sec_dict[config['ig_ggn_syn']['target']]
    post_sec = eval('ggn.{}'.format(config['ig_ggn_syn']['target']))
    make_ig_ggn_syn(ig, config['ig_ggn_syn'], post_sec)
    if kc_ig_syn_path in fd:
        connect_kcs_to_ig(fd, config['kc_ig_syn'])
    else:
        connect_kcs_to_ig(model_dict['kc'],  config['kc_ig_syn'])
    

def make_ig_ggn_syn(ig, params, sec):
    # This is ugly but only once
    global model_dict
    ig_syn = h.Exp2Syn(sec(0.5))
    ig_syn.tau1 = cfg.Q_(params['tau1']).to('ms').m
    ig_syn.tau2 = cfg.Q_(params['tau2']).to('ms').m
    ig_syn.e = cfg.Q_(params['e']).to('mV').m
    ig_nc = h.NetCon(ig.izhi._ref_V, ig_syn,
                          cfg.Q_(params['threshold']).to('mV').m,
                          cfg.Q_(params['delay']).to('ms').m,
                          cfg.Q_(params['gmax']).to('uS').m,
                          sec=ig.sec)
    model_dict['ig_ggn_syn'] = ig_syn
    model_dict['ig_ggn_nc'] = ig_nc
    cfg.logger.info('Completed creation of synapse from IG to GGN')
    cfg.logger.info('tau1: {}, tau2: {}, e: {}, thresh: {}, delay: {}, gmax: {}'.format(ig_syn.tau1,
                                                                                        ig_syn.tau2,
                                                                                        ig_syn.e,
                                                                                        ig_nc.threshold,
                                                                                        ig_nc.delay,
                                                                                        ig_nc.weight[0]))
                                                                 

def connect_kcs_to_ig(kcs=None, params=None):
    """Connect single compartmental KCs to IG"""
    global model_dict
    if isinstance(kcs, h5.File):  
        fd = kcs
        kc_name_sec_dict = {kc.soma.name(): kc.soma for kc in model_dict['kc']}
        threshold = -20.0
        if params is not None:
            threshold = cfg.Q_(params['threshold']).to('mV').m
        units = fd[kc_ig_syn_path].attrs['unit'][-1]
        syn_data = fd[kc_ig_syn_path]
        if len(syn_data.shape) == 2:
            syn_data = syn_data[:, 0]
        
        # for row_arr in fd[kc_ig_syn_path].value:
        #     row = row_arr[0]
        for row in syn_data:
            cfg.logger.debug('{}'.format(row))
            presec = kc_name_sec_dict[row['pre']]
            # TODO not using the unit attribute from dataset - assuming NEURON compatible units
            # ignoring 'prepos' field as KCs are single compartmental
            # - important data are weight and delay
            model_dict['ig'].connect_ampa_sec(presec, 0.5, threshold=threshold,
                                              delay=row['delay'],
                                              weight=row['gmax'])
    else: # no explicit synapse info, use constant values from config
        weight = float(params['weight'])
        if weight <= 0:
            cfg.logger.info('No connection from KCs to IG')
            return
        threshold = cfg.Q_(params['threshold']).to('mV').m
        delay = cfg.Q_(params['delay']).to('ms').m
        for kc in kcs:
            model_dict['ig'].connect_ampa_sec(kc.soma, 0.5, threshold=threshold,
                                              delay=delay, weight=weight)
    cfg.logger.info('Completed creation of synapse from KCs to IG')
    
    
def setup_recording(kcs, ggn, n_kc_vm=100, n_ggn_vm=100, t=None):
    """Create vectors for recording data from KCs and GGN sections"""
    global model_dict
    ret = {}
    kc_st = {}
    for kc in kcs:
        nc, stvec = nu.record_spiketimes(kc.soma)
        kc_st[kc.soma.name()] = (nc, stvec)
    ret['kc_spikes'] = kc_st
    if t is None:
        tvec = h.Vector()
        tvec.record(h._ref_t)
    elif isinstance(t, float):
        tvec = h.Vector()
        tvec.record(h._ref_t, t)
    elif isinstance(t, np.ndarray):  # Assume h.Vector
        tvec = h.Vector(t)
    else:
        tvec = t
    ret['time'] = tvec
    if n_kc_vm > kcs:
        n_kc_vm = len(kcs)
    vm_kcs = np.random.choice(kcs, size=n_kc_vm, replace=False)
    if 'test_kc' in model_dict:
        vm_kcs = np.append(vm_kcs, [model_dict['test_kc']])
    ret['kc_vm'] = {kc.soma.name(): nu.setup_recording(kc.soma, 0.5, 'v', t) \
                     for kc in vm_kcs}    
    g = nu.nrngraph(ggn)
    ggn_out_nodes = nu.select_good_nodes_by_sid(g, [ng.name_sid['LCA'], ng.name_sid['MCA']],
                                          [n_ggn_vm, n_ggn_vm])
    ggn_out = [g.node[n]['orig'] for n in ggn_out_nodes]
    ret['ggn_output_vm'] = {sec.name(): nu.setup_recording(sec, 0.5, 'v', t) \
                             for sec in ggn_out}
    ggn_alphaL_input_nodes = nu.select_good_nodes_by_sid(g, [ng.name_sid['alphaL']], [n_ggn_vm])
    ggn_alphaL_input = [g.node[n]['orig'] for n in ggn_alphaL_input_nodes]
    ret['ggn_alphaL_input_vm'] = {sec.name(): nu.setup_recording(sec, 0.5, 'v', t) \
                                   for sec in ggn_alphaL_input}
    ggn_basal_nodes = nu.select_good_nodes_by_sid(g, [ng.name_sid['dend_b']], [n_ggn_vm])
    ggn_basal = [g.node[n]['orig'] for n in ggn_basal_nodes]
    ret['ggn_basal_vm'] = {sec.name(): nu.setup_recording(sec, 0.5, 'v', t) \
                            for sec in ggn_basal}
    if 'ig' in model_dict:
        ret['ig_vm'] = h.Vector()
        ret['ig_vm'].record(model_dict['ig'].izhi._ref_V, t)
    model_dict.update(ret)
    return ret


def save_synapses(fd, model):
    start = timer()
    ###############################################
    # PN->KC synapses
    ###############################################
    conn_dtype = np.dtype([('pre', nsdf.VLENSTR),
                           ('post', nsdf.VLENSTR),
                           ('prepos', np.float64),
                           ('postpos', np.float64),
                           ('gmax', np.float64),
                           ('e', np.float64),
                           ('tau1', np.float64),
                           ('tau2', np.float64),
                           ('delay', np.float64)])
    pn_kc_netcons, pn_kc_syn_dict = model['pn_kc_conn']
    pn_kc_syn_data = np.empty((len(pn_kc_netcons),),
                              dtype=conn_dtype)
    for ii, (pre, post, nc) in enumerate(pn_kc_netcons):
        syn = nc.syn()
        pn_kc_syn_data['pre'][ii] = pre
        pn_kc_syn_data['post'][ii] = post
        pn_kc_syn_data['prepos'][ii] = 0
        pn_kc_syn_data['postpos'][ii] = 0.5  # hard coded for 1-compartment model
        pn_kc_syn_data['gmax'][ii] = nc.weight[0]
        pn_kc_syn_data['e'][ii] = syn.e
        pn_kc_syn_data['tau1'][ii] = syn.tau1
        pn_kc_syn_data['tau2'][ii] = syn.tau2
        pn_kc_syn_data['delay'][ii] = nc.delay
    conn_unit = ['', '', '', '', 'uS', 'mV', 'ms', 'ms', 'ms']
    pn_kc_syn_grp = fd.create_group('/data/static/pn_kc_synapse')
    pn_kc_syn_ds = fd.create_dataset('/data/static/pn_kc_synapse/pn_kc_synapse', data=pn_kc_syn_data)
    pn_kc_syn_ds.attrs['unit'] = conn_unit
    ###############################################
    # Create synapse dataset for KC->GGN
    ###############################################
    ## alphaL - not handling calyx in this case
    kc_ggn_syn_data = np.empty((len(model['kc_ggn_conn']),), dtype=conn_dtype)
    for ii, conn in enumerate(model['kc_ggn_conn']):
        kc_ggn_syn_data['pre'][ii] = conn['pre']
        kc_ggn_syn_data['post'][ii] = conn['post']
        kc_ggn_syn_data['prepos'][ii] = conn['prepos']
        kc_ggn_syn_data['postpos'][ii] = conn['postpos']
        kc_ggn_syn_data['gmax'][ii] = conn['nc'].weight[0]
        kc_ggn_syn_data['e'][ii] = conn['syn'].e
        kc_ggn_syn_data['tau1'][ii] = conn['syn'].tau1
        kc_ggn_syn_data['tau2'][ii] = conn['syn'].tau2
        kc_ggn_syn_data['delay'][ii] = conn['nc'].delay
    kc_ggn_syn_grp = fd.create_group('/data/static/kc_ggn_alphaL_synapse')
    kc_ggn_syn_ds = kc_ggn_syn_grp.create_dataset('kc_ggn_alphaL_synapse',
                                                  data=kc_ggn_syn_data)
    kc_ggn_syn_ds.attrs['unit'] = conn_unit
    ###############################################
    # Create synapse dataset for KC->IG - not handling
    ###############################################
    # if 'ig' in model:
    #     conn_unit = ['', '', '', '', '', 'ms']
    #     conn_dtype = np.dtype([('pre', nsdf.VLENSTR),
    #                            ('post', nsdf.VLENSTR),
    #                            ('prepos', np.float64),
    #                            ('postpos', np.float64),
    #                            ('gmax', np.float64),
    #                            ('delay', np.float64)])
    #     kc_ig_syn_data = np.empty((len(
    #     for row in model.kc_ig_syn_info:
    #         sources.append('{}__{}'.format(row['pre'], row['post']))
    #         kc_ig_syn_data.put_data(sources[-1],
    #                                 (row['pre'], row['post'], row['prepos'], row['postpos'],
    #                                  row['gmax'], row['delay']))
    #     kc_ig_syn_ds = writer.add_static_ds('kc_ig_synapse',
    #                                         sources)
    #     writer.add_static_data(kc_ig_syn_ds, kc_ig_syn_data)

    ###############################################
    # Save GGN->KC synapses
    ###############################################
    conn_dtype = np.dtype([('pre', nsdf.VLENSTR),
                           ('post', nsdf.VLENSTR),
                           ('prepos', np.float64),
                           ('postpos', np.float64),
                           ('vmid', np.float64),
                           ('vslope', np.float64),
                           ('e', np.float64),
                           ('gbar', np.float64),
                           ('tau', np.float64)])
    conn_unit = ['', '', '', '', 'mV', 'mV', 'mV', 'uS', 'ms']
    ggn_kc_conn_data = np.empty((len(model['ggn_kc_conn']),), dtype=conn_dtype)
    for ii, row in enumerate(model['ggn_kc_conn']):
        ggn_kc_conn_data['pre'][ii] = row['pre']
        ggn_kc_conn_data['post'][ii] = row['post']
        ggn_kc_conn_data['prepos'][ii] = row['prepos']
        ggn_kc_conn_data['postpos'][ii] = row['postpos']
        ggn_kc_conn_data['vmid'][ii] = row['vmid']
        ggn_kc_conn_data['vslope'][ii] = row['vslope']
        ggn_kc_conn_data['e'][ii] = row['e']
        ggn_kc_conn_data['gbar'][ii] = row['gbar']
        ggn_kc_conn_data['tau'][ii] = row['tau']
    ggn_kc_conn_grp = fd.create_group('/data/static/ggn_kc_synapse')
    ggn_kc_conn_ds = ggn_kc_conn_grp.create_dataset('ggn_kc_synapse',
                                                    data=ggn_kc_conn_data)
    ggn_kc_conn_ds.attrs['unit'] = conn_unit
    end = timer()
    cfg.logger.info('Saved synapse data in {} s'.format(end - start))
    

# def save(infile, outfile, kc_st, pn_st, kc_vm, ggn_output_vm, ggn_alphaL_input_vm, ggn_basal_vm, ig_vm, time):
def save(infile, outfile, data, model_dict, savesyn=False, config=None):
    """Save the model info and data with infile as reference for synapses
    and original PN spike trains.

    infile: input template file name
    outfile: output data filename
    kc_st: dict kc secname - spike trains
    pn_st: dict pn secname - spike trains
    kc_vm: dict kc secname to vm vec
    ggn_output_vm: dict ggn secname - Vm vec
    ggn_alphaL_input_vm: dict ggn secname - Vm
    ggn_basal_vm: dict ggn secname - Vm
    config: yaml document representing configuration information in template file.
    """
    cfg.logger.info('Start saving data in {}'.format(outfile))
    kc_vm = data['kc_vm']
    ggn_alphaL_input_vm = data['ggn_alphaL_input_vm']
    ggn_basal_vm = data['ggn_basal_vm']
    time = data['time']
    with h5.File(outfile, 'w') as fd:
        fd.attrs['original'] = infile
        fd.attrs['command'] = ' '.join(sys.argv)
        fd.attrs['dt'] = h.dt  # integration time step
        if config is not None:
            fd.attrs['config'] = yaml.dump(config)
        else:
            cfg.logger.info('No config passed')
        # Save synapse info only if requested
        if savesyn:
            save_synapses(fd, model_dict)
        # Save dynamic data
        time_grp = fd.create_group('/map/time')
        time_ds = time_grp.create_dataset('time', data=time)
        kc_st_grp = fd.create_group(kc_st_path)
        for kc, (nc, st) in data['kc_spikes'].items():
            ii = kc.partition('[')[-1].partition(']')[0]
            kc_ds = kc_st_grp.create_dataset(ii, data=np.array(st))
            kc_ds.attrs['source'] = kc
        pn_st_grp = fd.create_group(pn_st_path)
        for pn, st in data['pn_spikes'].items():
            pn_st_grp.create_dataset(pn, data=st)
        alphaL_grp = fd.create_group('/data/uniform/ggn_alphaL_input')
        alphaL_secnames = list(ggn_alphaL_input_vm.keys())
        alphaL_data = list(ggn_alphaL_input_vm.values())
        alphaL_ds = alphaL_grp.create_dataset('GGN_alphaL_input_Vm', data=np.array(alphaL_data),
                                              compression='gzip')
        alphaL_src = fd.create_dataset('/map/uniform/ggn_alphaL_input', data=alphaL_secnames,
                                       compression='gzip')
        alphaL_ds.attrs['dt'] = time[1] - time[0]            
        alphaL_ds.attrs['tunit'] = 'ms'
        alphaL_ds.attrs['unit'] = 'mV'
        alphaL_ds.attrs['field'] = 'v'
        alphaL_ds.attrs['tstart'] = time[0]
        alphaL_ds.dims[0].label = 'source'
        alphaL_ds.dims.create_scale(alphaL_src, 'source')
        alphaL_ds.dims[0].attach_scale(alphaL_src)
        alphaL_ds.dims[1].label = 'time'
        alphaL_ds.dims.create_scale(time_ds, 'time')
        alphaL_ds.dims[1].attach_scale(time_ds)
        basal_grp = fd.create_group('/data/uniform/ggn_basal')
        basal_secnames = list(ggn_basal_vm.keys())
        basal_data = list(ggn_basal_vm.values())
        basal_ds = basal_grp.create_dataset('GGN_basal_Vm', data=np.array(basal_data),
                                              compression='gzip')
        basal_src = fd.create_dataset('/map/uniform/ggn_basal', data=basal_secnames,
                                      compression='gzip')
        basal_ds.attrs['dt'] = time[1] - time[0]            
        basal_ds.attrs['tunit'] = 'ms'
        basal_ds.attrs['unit'] = 'mV'
        basal_ds.attrs['field'] = 'v'
        basal_ds.attrs['tstart'] = time[0]
        basal_ds.dims[0].label = 'source'
        basal_ds.dims.create_scale(basal_src, 'source')
        basal_ds.dims[0].attach_scale(basal_src)
        basal_ds.dims[1].label = 'time'
        basal_ds.dims.create_scale(time_ds, 'time')
        basal_ds.dims[1].attach_scale(time_ds)
        output_grp = fd.create_group('/data/uniform/ggn_output')
        output_secnames = list(data['ggn_output_vm'].keys())
        output_data = list(data['ggn_output_vm'].values())
        output_ds = output_grp.create_dataset('GGN_output_Vm', data=np.array(output_data),
                                              compression='gzip')
        output_src = fd.create_dataset('/map/uniform/ggn_output', data=output_secnames,
                                       compression='gzip')
        output_ds.attrs['dt'] = time[1] - time[0]            
        output_ds.attrs['tunit'] = 'ms'
        output_ds.attrs['unit'] = 'mV'
        output_ds.attrs['field'] = 'v'
        output_ds.attrs['tstart'] = time[0]
        output_ds.dims[0].label = 'source'
        output_ds.dims.create_scale(output_src, 'source')
        output_ds.dims[0].attach_scale(output_src)
        output_ds.dims[1].label = 'time'
        output_ds.dims.create_scale(time_ds, 'time')
        output_ds.dims[1].attach_scale(time_ds)
        kc_vm_grp = fd.create_group('/data/uniform/kc')
        kc_vm_ds = kc_vm_grp.create_dataset('KC_Vm', data=np.array(list(kc_vm.values())))
        kc_vm_src = fd.create_dataset('/map/unifrom/kc', data=list(kc_vm.keys()),
                                      compression='gzip')
        kc_vm_ds.attrs['dt'] = time[1] - time[0]            
        kc_vm_ds.attrs['tunit'] = 'ms'
        kc_vm_ds.attrs['unit'] = 'mV'
        kc_vm_ds.attrs['field'] = 'v'
        kc_vm_ds.attrs['tstart'] = time[0]
        kc_vm_ds.dims[0].label = 'source'
        kc_vm_ds.dims.create_scale(kc_vm_src, 'source')
        kc_vm_ds.dims[0].attach_scale(kc_vm_src)
        kc_vm_ds.dims[1].label = 'time'
        kc_vm_ds.dims.create_scale(time_ds, 'time')
        kc_vm_ds.dims[1].attach_scale(time_ds)
        if data['ig_vm'] is not None:
            ig_grp = fd.create_group('/data/uniform/ig')
            ig_ds = ig_grp.create_dataset('IG_Vm', data=np.array([data['ig_vm']]),
                                       compression='gzip')
            ig_src = fd.create_dataset('/map/uniform/ig', data=['IG'],
                                          compression='gzip')
            ig_ds.attrs['dt'] = time[1] - time[0]            
            ig_ds.attrs['tunit'] = 'ms'
            ig_ds.attrs['unit'] = 'mV'
            ig_ds.attrs['field'] = 'v'
            ig_ds.attrs['tstart'] = time[0]
            ig_ds.dims[0].label = 'source'
            ig_ds.dims.create_scale(ig_src, 'source')
            ig_ds.dims[0].attach_scale(ig_src)
            ig_ds.dims[1].label = 'time'
            ig_ds.dims.create_scale(time_ds, 'time')
            ig_ds.dims[1].attach_scale(time_ds)               
    cfg.logger.info('Finished saving data')


def test(infile='/data/rays3/ggn/olfactory_network/mb_net_UTC2018_09_19__00_49_12-PID4273-JID9673664.h5'):
    with h5.File(infile, 'r') as fd:
        pn_spikes = load_pn_spikes(fd)
        pns, spikes = zip(*pn_spikes)
        spikes[0] = np.array([1.0, 10.0, 20.0])
        stimvecs, vecstims = create_pn_output(spikes)
        kcs = create_kcs(fd)
        kc_name_sec_dict = {kc.soma.name(): kc.soma for kc in kcs}
        ggn = create_ggn()
        ggn_name_sec_dict = {sec.name(): sec for sec in ggn.all}
        nc_pn_kc, syn_pn_kc = create_pn_kc_conn(fd, {pn: vecstim for pn, vecstim in zip(pns, vecstims)},
                                      kc_name_sec_dict)
        for nc in nc_pn_kc:
            assert nc.valid()
                         
        syn_ggn_kc = create_ggn_kc_conn(fd, kc_name_sec_dict, ggn_name_sec_dict)
        syn_kc_ggn, nc_kc_ggn = create_kc_ggn_conn(fd, kc_name_sec_dict, ggn_name_sec_dict, path=kc_ggn_alphaL_syn_path)
        # kc_st = {kc: np.random.random_sample(5) for kc in np.random.choice(kcs, size=5, replace=False)}
        # pn_st = {pn: np.random.random_sample(5) for pn in np.random.choice(pns, size=5, replace=False)}
        # ggn_vm = {sec: np.arange(10) for sec in np.random.choice(list(ggn_name_sec_dict.keys()),
        #                                                          size=5, replace=False)}
        # kc_vm = {kc: np.random.random_sample(5) for kc in np.random.choice(kcs, size=5, replace=False)}
        data = setup_recording(kcs, ggn, n_kc_vm=10, n_ggn_vm=10, t=0.25)
        h.tstop = 100
        h.init()
        h.run()
        kc_st = {kc: st for kc, (nc, st) in data['kc_spikes'].items()}
        save(infile, 'test.h5', kc_st, pn_spikes, data['kc_vm'], data['ggn_output_vm'],
             data['ggn_alphaL_input_vm'], data['ggn_basal_vm'], None, data['time'])
        cfg.logger.info('finished')


def run_model(args):
    """setup and run a model with templates and other parameters specified
    in args (parsed arguments). List of arguments:

    template_filename: HDF5 file containing the network template and
    PN spike trains. These data should be at paths in constants
    specified at the top of this file.

    output_directory: directory to dump simulated data into.

    pn_shift: shift the assignment of spike trains to PNs by this
    amount, i.e., if pn_shift is 2, then the spike train of pn_0 is
    assigned to pn_2, and pn_{i+2} gets the spike train of pn_{i},
    with wrapping around the edge.

    pn_dither: The maximum magnitude of time shift when dithering the
    PN spike times.

    n_kc_vm: number of KCs to record Vm from

    n_ggn_vm: number of GGN sections to record Vm from for each of
    basal, calyceal and alpha lobe regions.

    recstep: number of integration steps between each recording point.

    simtime: total simulation time

    """
    global model_dict
    output_filename = os.path.join(args.output_directory,
                                   'fixed_net_UTC{}-PID{}-JID{}.h5'.format(
                                       cfg.timestamp.strftime('%Y_%m_%d__%H_%M_%S'),
                                       cfg.mypid, cfg.myjobid))
    ggn = create_ggn()
    ggn_name_sec_dict = {sec.name(): sec for sec in ggn.all}
    with h5.File(args.template_filename, 'r') as template_fd:
        config = yaml.load(template_fd.attrs['config'].decode())
        pn_spikes = load_pn_spikes(template_fd)
        pns, spike_trains = zip(*pn_spikes)
        spike_trains = dither_spiketrain(spike_trains, cell_shift=args.pn_shift,
                                         dither=args.pn_dither)
        pn_spike_vecs, vecstims = create_pn_output(spike_trains)
        kcs = create_kcs(template_fd, config=config['kc'])
        if args.test_kc >= 0:
            delay = cfg.Q_(config['stimulus']['onset'])
            if 'delay' in config['pn_kc_syn']:
                delay += cfg.Q_(config['pn_kc_syn']['delay'])
            duration = cfg.Q_(config['stimulus']['duration'])
            amplitude = cfg.Q_(args.kc_current)
            iclamp = ephys.setup_current_clamp(kcs[args.test_kc].soma,
                                               delay=delay, duration=duration, amplitude=amplitude)
            # test_kc_vvec = ephys.setup_sec_rec(kcs[args.test_kc].soma, 'v')[0] # this is added in setup recording
            model_dict['kc_iclamp'] = iclamp
            model_dict['test_kc'] = kcs[args.test_kc]
        kc_name_sec_dict = {kc.soma.name(): kc.soma for kc in kcs}
        nc_pn_kc, syn_pn_kc = create_pn_kc_conn(template_fd,
                                                dict(zip(pns, vecstims)), kc_name_sec_dict,
                                                config=config['pn_kc_syn'])
        syn_ggn_kc = create_ggn_kc_conn(template_fd, kc_name_sec_dict,
                                        ggn_name_sec_dict, config=config['ggn_kc_syn'])
        syn_kc_ggn, nc_kc_ggn = create_kc_ggn_conn(template_fd,
                                                   kc_name_sec_dict,
                                                   ggn_name_sec_dict,
                                                   config=config['kc_ggn_alphaL_syn'])
        setup_ig(template_fd, ggn_name_sec_dict, config=config)
        data = setup_recording(kcs, ggn, n_kc_vm=args.n_kc_vm,
                               n_ggn_vm=args.n_ggn_vm,
                               t=h.dt * args.recstep)
        h.tstop = args.simtime
        start = datetime.utcnow()
        h.init()
        for v in pn_spike_vecs:
            if np.any(np.array(v.x) <= 0):
                print('negative pn spike time', v)
        cfg.logger.info('Finished init. Starting simulation of {} ms'.format(
            h.tstop))
        if h.tstop > 0:
            nu.block_run(logger=cfg.logger)
        else:
            cfg.logger.info('Dry run. Skipping simulation')
        end = datetime.utcnow()
        delta = end - start
        data['tstart'] = start
        data['tend'] = end
        data['pn_spikes'] = dict(zip(pns, spike_trains))
        cfg.logger.info('Finished running simulation in {} s'.format(
            delta.days * 86400 +
            delta.seconds +
            delta.microseconds * 1e-6))
        cfg.logger.info('Starting data save in {}'.format(output_filename))
        ig_vm = model_dict.get('ig_vm',  None)
        data['ig_vm'] = ig_vm
        save(args.template_filename, output_filename, data, model_dict, args.savesyn, config)
        cfg.logger.info('Finished')


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--template_filename', type=str,
                        help='template filename for pn spike trains and network config')
    parser.add_argument('-o', '--output_directory', type=str,
                        help='directory to save simulation data in')
    parser.add_argument('-s', '--pn_shift', type=int, default=0,
                        help='shift the PN spikes over the PN population by this many PNs')
    parser.add_argument('-d', '--pn_dither', type=float, default=0.0,
                        help='dither the PN spikes by maximum of this much time (in ms)')
    parser.add_argument('--n_kc_vm', type=int, default=500,
                        help='number of KCs to record Vm from')
    parser.add_argument('--n_ggn_vm', type=int, default=100,
                        help='number of section of GGN in each region to record from')
    parser.add_argument('--recstep', default=5, type=int,
                        help='number of integration steps per data recording step')
    parser.add_argument('--simtime', type=float, help='Simulation time (ms)')
    parser.add_argument('--savesyn', action='store_true', help='Save the synapse information (these are humongous data)')
    parser.add_argument('--debug', action='store_true',
                        help='Set log level to debug. Default is info')
    parser.add_argument('--test_kc', type=int, default=-1, help='Inject current into KC with this index')
    parser.add_argument('--kc_current', type=str, default='0pA', help='Amount of current (with unit) to inject in test kc')
    return parser


def main():
    cfg.logger.info('COMMAND LINE: {}'.format(' '.join(sys.argv)))
    print('COMMAND LINE: {}'.format(' '.join(sys.argv)))
    sys.stdout.flush()
    parser = make_parser()
    args = parser.parse_args()
    if args.debug:
        cfg.logger.setLevel(10)
    else:
        cfg.logger.setLevel(20)
    run_model(args)
    

if __name__ == '__main__':
    # test(infile='/data/rays3/ggn/fixed_net/mb_net_UTC2018_09_25__22_25_34-PID16017-JID10014019.h5')
    main()



#
# fixed_network_changing_stim.py ends here
