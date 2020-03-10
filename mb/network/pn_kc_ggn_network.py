# pn_kc_ggn_network.py ---
#
# Filename: pn_kc_ggn_network.py
# Description:
# Author: Subhasis Ray
# Maintainer:
# Created: Wed Feb 21 13:17:05 2018 (-0500)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Mar 10 15:42:53 2020 (-0400)
#           By: Subhasis Ray
#     Update #: 2837

# Code:
"""This script if for creating the PN->KC<->GGN network with KC<->GGN
reciprocal synapses are spatially clustered.

The input to KCs may also be spatially correlated (or not).

"""
from __future__ import print_function
from config import Q_, h, logger, timestamp, mypid, myjobid, nrn_version
import sys
import os
import random
import numpy as np
# from matplotlib import pyplot as plt
from timeit import default_timer as timer
from collections import defaultdict
from datetime import datetime
import yaml
import warnings
from sklearn import cluster
import argparse
import itertools as it
import nrnutils as nu
import ephys
import pn_output as pnout
import ig
import neurograph as ng
import nsdf

print('Imported all modules')
sys.stdout.flush()

warnings.filterwarnings('ignore')
with open('network/config.yaml', 'rt') as fd:
    netconf = yaml.load(fd)
print('Loaded config.yaml')
sys.stdout.flush()

# Faking sklearn.KMeans for fake clusters87
class FakeKMeans(object):
    """A dummy class to substitute for sklearn.KMeans when we are creating
    fake clustes"""
    def __init__(self, labels):
        self.labels_ = labels


#################################
# Code for model setup
#################################

class pn_population(object):
    def __init__(self, params):
        """Container for population of PN spike outputs.

        pns: list of strings representing IDs of all the PNs: pn_{ii}
        for ii in range of pn numbers.

        outputs: list of numpy arrays representing spike times of each PN.

        vecstims: list of VecStim objects corresponding to each PN.

        stimves: list of Vector objects storing the spike times for
        each PN to be played to the corresponding VecStim object.

        """
        start = timer()
        # Create the PN spike times
        onset = Q_(params['stimulus']['onset'])
        duration = Q_(params['stimulus']['duration'])
        tail = Q_(params['stimulus']['tail'])
        shifting = params['pn']['shifting']
        newparams = dict(onset=onset,
                         duration=duration,
                         tail=tail,
                         delta=Q_(params['pn']['delta']),
                         npn=params['pn']['number'],
                         stim_rate=Q_(params['pn']['stim_rate']),
                         spont_rate=Q_(params['pn']['spont_rate']),
                         osc_freq=Q_(params['pn']['osc_freq']),
                         osc_amp_scale=params['pn']['osc_amp_scale'])

        if shifting:
            newparams.update(dict(offdur=Q_(params['pn']['offdur']),
                                  start_frac=params['pn']['start_frac'],
                                  shifting_frac=params['pn']['shifting_frac']))

            _output_list, clusters, res_types = pnout.create_shifting_response(newparams)
            self.clusters = clusters
            self.response_types = res_types

        else:
            newparams.update(dict(odor_exc_frac=params['pn']['odor_exc_frac'],
                                  odor_inh_frac=params['pn']['odor_inh_frac'],
                                  spont_exc_frac=params['pn']['spont_exc_frac'],
                                  osc_amp_scale=params['pn']['osc_amp_scale']))
            _output_list = pnout.create_pn_output(newparams)
            self.clusters = {}
            self.response_types = {}
        stabilization_time = Q_(params['stimulus']['stabilization_time'])
        self.pns = []
        self.outputs = []
        self.vecstims = []
        self.stimvecs = []
        # Remove spikes appearing at or before stabilization_time time
        for ii, spiketimes in enumerate(_output_list):
            st = spiketimes[spiketimes > stabilization_time].copy()
            self.pns.append('pn_{}'.format(ii))
            stimvec = h.Vector(st.to('ms').m)
            vecstim = h.VecStim()
            vecstim.play(stimvec)
            self.stimvecs.append(stimvec)
            self.vecstims.append(vecstim)
            self.outputs.append(st)
        end = timer()
        logger.info('Created {} PNs in {} s'.format(
            len(self.outputs), end - start))


class kc_population(object):
    """Container for KC population"""
    def __init__(self, params):
        """Attribute:
        kcs - a list of KC models.

        positions_3d - x, y, z position of the KCs - initially all 0s.
        updated to 3D location of the presynaptic GGN segment
        corresponding to this KC.

        positions_x - 1D section position of the KCs on GGN sections.

        """
        start = timer()
        if not hasattr(h, params['kc']['name']):
            h.xopen(params['kc']['filename'])
        kc_class = eval('h.{}'.format(params['kc']['name']))
        self.kcs = [kc_class() for ii in range(params['kc']['number'])]
        self.positions_3d = np.zeros((params['kc']['number'], 3))
        self.positions_x = np.zeros(params['kc']['number'])
        self.cluster_ids = np.zeros(params['kc']['number'])
        lca_input_count = int(params['kc']['number'] * params['kc']['lca_frac'])
        self.region = ['lca' for ii in range(lca_input_count)] \
                      + ['mca' for ii in range(lca_input_count,
                                               params['kc']['number'])]
        self.lca_kcs = self.kcs[: lca_input_count]
        self.mca_kcs = self.kcs[lca_input_count:]
        end = timer()
        logger.info('Created {} KCs in {} s'.format(
            len(self.kcs), end - start))


class ggn_population(object):
    def __init__(self, params):
        """Container for the GGN. Should have a singleton GGN.

        Based on the configuration we can also add an IG model, either
        as a pre-recorded spike train in an npy file or as an
        Izhikevich model.

        """
        start = timer()
        if not hasattr(h, params['ggn']['name']):
            h.xopen(params['ggn']['filename'])
        self.ggn = eval('h.{}()'.format(params['ggn']['name']))
        for sec in self.ggn.all:
            sec.Ra = Q_(params['ggn']['RA']).to('ohm * cm').m
            sec.g_pas = 1.0 / Q_(params['ggn']['RM']).to('ohm * cm**2').m
        self.ggn.geom_nseg()
        self.ggn_graph = nu.nrngraph(self.ggn)
        self.create_calyceal_clusters(params)
        if params['ggn']['dclamp']:
            self.add_dynamic_clamp(params['ggn']['dclamp_file'],
                                   params['ggn']['dclamp_sec'],
                                   params['ggn']['dclamp_pos'])
        self.ig = None
        if 'ig' in params:
            ig_name = params['ig']['filename'].strip()
            if ig_name == 'IzhiGS':
                logger.info('Creating Izhikevich model of IG')
                inject = params['ig'].get('inject', 70.0)
                gbarGABA = Q_(params['ggn_ig_syn'].get('gmax', '1.0uS')).to('uS').m
                self.ig = ig.IGIzhi(inject=inject, gbarGABA=gbarGABA)
                self.ig.vmidGraded = Q_(params['ggn_ig_syn'].get('vmid','-45.0mV')).to('mV').m
                self.ig.vslopeGraded = Q_(params['ggn_ig_syn'].get('vslope','5.0mV')).to('mV').m
                presec = eval('self.ggn.{}'.format(params['ggn_ig_syn'].get('source', 'soma')))
                h.setpointer(presec(0.5)._ref_v, 'vpreGraded', self.ig.izhi)
                logger.info('Created IG with injection {}'.format(inject))
                self.make_ig_ggn_syn(params['ig_ggn_syn'])
            elif len(ig_name) > 0:  # Spike trains predefined
                logger.info('Creating VecStim model of IG spiketrain')
                stime = Q_(params['stimulus']['stabilization_time'])
                onset = Q_(params['stimulus']['onset']).to('s').m
                duration = Q_(params['stimulus']['duration']).to('s').m
                tail = Q_(params['stimulus']['tail']).to('s').m
                delta = Q_(params['pn']['delta']).to('s').m
                simtime = onset + duration + tail
                self.ig = ig.IGSpikes(params['ig']['filename'], simtime,
                                      delta,
                                      stabilization_time=stime.to('s').m)
                self.ig.connect(self.ggn.soma, 0.5, params['ig_ggn_syn'])
        end = timer()
        logger.info('Created GGN and clusterd its '
                    'KC input segments in {} s'.format(end - start))

    def make_ig_ggn_syn(self, params):
        # This is ugly but only once
        target = params.get('target', 'dend[1]')
        sec = eval('self.ggn.{}'.format(target))
        self.ig_syn = h.Exp2Syn(sec(0.5))
        self.ig_syn.tau1 = Q_(params['tau1']).to('ms').m
        self.ig_syn.tau2 = Q_(params['tau2']).to('ms').m
        self.ig_syn.e = Q_(params['e']).to('mV').m
        self.ig_nc = h.NetCon(self.ig.izhi._ref_V, self.ig_syn,
                              Q_(params['threshold']).to('mV').m,
                              Q_(params['delay']).to('ms').m,
                              Q_(params['gmax']).to('uS').m,
                              sec=self.ig.sec)

    def get_clusters(self, sid, numpos, size, random_state=None, fake=False):
        sec_pos = nu.select_random_segments_by_sid(self.ggn_graph,
                                                   sid, numpos,
                                                   by_length=True,
                                                   replace=True)
        # print('numpos: {}, num secpos: {}'.format(numpos, len(sec_pos)))
        n_clust = int(round(numpos * 1.0 / size))
        if n_clust <= 1:
            logger.info('Number of GGN position clusters <= 1. Returning Fakse clusters')
            fake = True
            n_clust = 1  # avoid 0 for degenrate case
        kmeans, sec_list, x_list, pos_list = cluster_positions(
            sec_pos, n_clust, random_state=random_state, fake=fake)
        return kmeans, sec_list, x_list, pos_list

    def create_calyceal_clusters(self, params):
        """Create the LCA and MCA positional clusters."""
        cluster_size = params['kc']['cluster_size']
        random_state = params['kc'].get('cluster_rs', None)
        lca_input_count = int(params['kc']['number'] * params['kc']['lca_frac'])
        lca_clus = self.get_clusters(ng.custom_types['LCA'], lca_input_count,
                                     cluster_size, random_state=random_state,
                                     fake=params['kc']['fake_clusters'])
        self.lca_kmeans, self.lca_sec, self.lca_x, self.lca_pos3d = lca_clus
        mca_input_count = params['kc']['number'] - lca_input_count
        mca_clus = self.get_clusters(ng.custom_types['MCA'], mca_input_count,
                                     cluster_size, random_state=random_state,
                                     fake=params['kc']['fake_clusters'])
        print('lca_input_count: {}'.format(lca_input_count))
        print('mca_input_count: {}'.format(mca_input_count))
        print('Total: {}'.format(lca_input_count + mca_input_count))
        self.mca_kmeans, self.mca_sec, self.mca_x, self.mca_pos3d = mca_clus

    def add_dynamic_clamp(self, fname, secname='dend[1]', secpos=0.5):
        """Create a dynamic clamp to apply Vm from file fname at secname(secpos).

        file fname should be a recarray (Vm, t) where Vm is in mV and
        t contains time points in ms

        """
        data = np.load(fname)
        sec = eval('self.ggn.{}'.format(secname))
        self.clamp, self.vplay = nu.setup_dynamic_clamp(sec, secpos, data['Vm'], data['t'])
        logger.info('Created dynamic clamp on GGN.{} at {}.'.format(secname, secpos))

    def connect_kcs_to_ig(self, kcs, params):
        """Connect single compartmental KCs to IG"""
        if not isinstance(self.ig, ig.IGIzhi):
            raise Exception('Cannot connect KCs without a real IG model')
        weight = float(params.get('weight'))
        if weight <= 0:
            logger.info('No connection from KCs to IG')
            return
        std = float(params.get('std', 0.0))
        if std > 0:
            logger.info('KC->IG weight distributed lognormally')
            # Use lognormal distribution: specified gbar is the
            # mean and specified std is as a fraction of the
            # mean. Thus the mu and sigma parameters become:
            # mu = ln(mean/sqrt(1 + variance/mean^2)
            # and sigma^2 = ln(1 + variance/mean^2)
            # thus, mu = ln(mean) - 0.5 ln(1 + variance/mean^2) = ln(mean) - 0.5 sigma^2
            sigma2 = np.log(1 + std**2)  # std is specified as a fraction of mean
            mu = np.log(weight) - sigma2 / 2.0
            gmax = np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2),
                                             size=len(kcs))
        else:
            gmax = [weight] * len(kcs)
        threshold = Q_(params.get('threshold', '-20mV')).to('mV').m
        delay = Q_(params.get('delay', '5ms')).to('ms').m
        kc_ig_syn_info = []
        for kc, w in zip(kcs, gmax):
            self.ig.connect_ampa_sec(kc.soma, 0.5, threshold=threshold,
                                     delay=delay, weight=weight)
            kc_ig_syn_info.append({'pre': kc.soma.name(),
                                   'post': 'IG',
                                   'prepos': 0.5,
                                   'postpos': 0.5,
                                   'gmax': weight,
                                   'delay': delay})
        return kc_ig_syn_info

    def connect_pns_to_ig(self, pns, params):
        """Connect VecStim representation of PN spike trains to IG"""
        if not isinstance(self.ig, ig.IGIzhi):
            raise Exception('Cannot connect PNs without a real IG model')
        weight = float(params.get('weight', 1.0))
        if weight <= 0:
            logger.info('No connection from PNs to IG')
            return
        threshold = Q_(params.get('threshold', '-20mV')).to('mV').m
        delay = Q_(params.get('delay', '5ms')).to('ms').m
        if 'tau' in params:
            self.ig.tauAMPA = Q_(params['tau']).to('ms').m
        for pn in pns:
            self.ig.connect_ampa_vecstim(pn, threshold=threshold,
                                         delay=delay, weight=weight)
        logger.info('Connected {} PN->IG synapse'.format(len(self.ig.netcons)))

class pn_kc_ggn_network(object):
    """container for the whole network."""
    def __init__(self, params):
        """Create the network."""
        start = timer()
        self.pns = pn_population(params)
        self.kcs = kc_population(params)
        self.ggn = ggn_population(params)
        self.ig = self.ggn.ig
        self.kc_ig_syn_info = []
        self.connect_ggn_kc(params['ggn_kc_syn'])
        if Q_(params['kc_ggn_CA_syn']['gmax']).to('pS').m > 0.0:
            self.connect_kc_ggn_ca(params['kc_ggn_CA_syn'])
        self.connect_kc_ggn_alpha(params['kc_ggn_alphaL_syn'])
        self.connect_pn_kc(params['pn_kc_syn'],
                           params['kc']['shared_pn_frac'])
        if 'pn_ig_syn' in params:
            logger.info('Connecting PNs to IG')
            self.ggn.connect_pns_to_ig(self.pns.vecstims, params['pn_ig_syn'])
        if 'kc_ig_syn' in params:
            self.kc_ig_syn_info = self.ggn.connect_kcs_to_ig(self.kcs.kcs, params['kc_ig_syn'])
        # self.setup_data_recording(params)
        end = timer()
        logger.info('Created model and data recorders in {} s'.format(end - start))

    def connect_ggn_kc(self, ggn_kc_syn_params):
        """Connect GGN segments to KCs. This also assigns spatial cluster ids
        and ggn segment positions to the post synaptic KC. Note that
        LCA and MCA cluster ids overlap, but they can be separated by
        the fact that the LCA KCs are the first N entries in the kcs
        list where N = len(kcs.lca_kcs), the total number of KCs post
        synaptic to LCA and the entries after that are MCA KCs

        """
        assert len(self.kcs.lca_kcs) == len(self.ggn.lca_sec) == len(self.ggn.lca_x) \
            == len(self.ggn.lca_pos3d)
        self.ggn_kc_adjacency = []
        self.ggn_kc_synapses = []
        if len(self.kcs.kcs) == 0:
            return
        # if Q_(ggn_kc_syn_params['gmax']).to('pS').m <= 0.0:
        #     return
        positions_x = np.concatenate((self.ggn.lca_x, self.ggn.mca_x))
        positions_3d = np.concatenate((self.ggn.lca_pos3d, self.ggn.mca_pos3d))
        sections = np.concatenate((self.ggn.lca_sec, self.ggn.mca_sec))
        cluster_ids = np.concatenate((self.ggn.lca_kmeans.labels_,
                                      self.ggn.mca_kmeans.labels_))
        print('Number of cluster labels: {}'.format(len(cluster_ids)))
        ggn_kc_attrs = {
            'vmid': Q_(ggn_kc_syn_params.get('vmid', '-40mV')).to('mV').m,
            'vslope': Q_(ggn_kc_syn_params.get('vslope', '5.0mV')).to('mV').m,
            'e': Q_(ggn_kc_syn_params.get('e', '-80mV')).to('mV').m,
            'tau': Q_(ggn_kc_syn_params.get('tau', '4.0ms')).to('ms').m}
        frac_weakly_inhibited = ggn_kc_syn_params.get('frac_weakly_inhibited', 0.0)
        wkc_count = int(len(self.kcs.kcs) * frac_weakly_inhibited)  # Weakly inhibited KCs
        nkc_count = len(self.kcs.kcs) - wkc_count                   # Normally inhibited KCs
        mean_gbar = Q_(ggn_kc_syn_params.get('gmax', '1e-3uS')).to('uS').m
        n_gbar = np.ones(nkc_count) * mean_gbar
        mean_w_gbar = Q_(ggn_kc_syn_params.get('gmax_weakly_inhibited', '0.0nS')).to('uS').m
        w_gbar = np.ones(wkc_count) * mean_w_gbar
        logger.info('{} KCs Weak receive inhibitory synapses of '
                    'mean gbar={} from GGN'.format(wkc_count, mean_w_gbar))
        std = float(ggn_kc_syn_params.get('std', 0))
        if std > 0:
            logger.info('GGN->KC gbar distributed lognormally')
            # Use lognormal distribution: specified gbar is the
            # mean and specified std is as a fraction of the
            # mean. Thus the mu and sigma parameters become:
            # mu = ln(mean/sqrt(1 + variance/mean^2)
            # and sigma^2 = ln(1 + variance/mean^2)
            # thus, mu = ln(mean) - 0.5 ln(1 + variance/mean^2) = ln(mean) - 0.5 sigma^2
            sigma2 = np.log(1 + std**2)  # std is specified as a fraction of mean
            mu = np.log(mean_gbar) - sigma2 / 2.0
            n_gbar = np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2),
                                         size=nkc_count)
            w_mu = np.log(mean_w_gbar) - sigma2 / 2.0
            w_gbar = np.random.lognormal(mean=w_mu, sigma=np.sqrt(sigma2),
                                         size=wkc_count)
        gbar = np.concatenate((w_gbar, n_gbar))
        # Randomize the weak connections
        np.random.shuffle(gbar)
        self.kcs.positions_x[:] = positions_x[:]
        self.kcs.positions_3d[:, :] = positions_3d[:, :]
        self.kcs.cluster_ids[:] = cluster_ids[:]
        for ii, kc in enumerate(self.kcs.kcs):
            ggn_kc_syn = h.GradedSyn(kc.soma(0.5))
            for name, value in ggn_kc_attrs.items():
                setattr(ggn_kc_syn, name, value)
            ggn_kc_syn.gbar = gbar[ii]
            h.setpointer(sections[ii](positions_x[ii])._ref_v, 'vpre', ggn_kc_syn)
            self.ggn_kc_adjacency.append((sections[ii].name(),
                                          ggn_kc_syn.get_segment().sec.name(),
                                          positions_x[ii]))
            self.ggn_kc_synapses.append(ggn_kc_syn)

    def connect_kc_ggn_ca(self, kc_ggn_ca_syn_params):
        """three options for connectinng KCs to GGN in calyx:

        1. by spatial clusters - the KCs are assigned cluster ID
        based on their presynaptic GGN segment location. In this
        case, lateral and medial calyx (LCA and MCA respectively)
        are clustered separately. Each KC will synapse onto a
        random segment of the GGN in the same cluster.

        2. By region: the KCs are restricted to either LCA or MCA
        branch of the GGN synapse onto a random segment in the
        corresponding branch of the GGN.

        3. Completely random: the KCs connect to completely random
        segments from LCA + MCA.

        """
        if Q_(kc_ggn_ca_syn_params['gmax']).to('pS').m <= 0.0:
            self.kc_ggn_LCA_synapses = []
            self.kc_ggn_LCA_netcons = []
            self.kc_ggn_LCA_adjacency = []
            self.kc_ggn_MCA_synapses = []
            self.kc_ggn_MCA_netcons = []
            self.kc_ggn_MCA_adjacency = []
            self.kc_ggn_CA_synapses = []
            self.kc_ggn_CA_adjacency = []
            self.kc_ggn_CA_netcons = []
            return
        if len(self.kcs.kcs) == 0:
            return
        if kc_ggn_ca_syn_params['clustered']:
            lca_ret = connect_in_cluster(self.ggn.lca_kmeans.labels_,
                                         self.kcs.lca_kcs,
                                         self.ggn.lca_sec,
                                         self.ggn.lca_x, kc_ggn_ca_syn_params)
            self.lca_kc_ggn_syn = lca_ret['synapses']
            self.lca_kc_ggn_netcon = lca_ret['netcons']
            self.lca_kc_ggn_adjacency = lca_ret['adjacency']
            mca_ret = connect_in_cluster(self.ggn.mca_kmeans.labels_,
                                         self.kcs.mca_kcs,
                                         self.ggn.mca_sec,
                                         self.ggn.mca_x, kc_ggn_ca_syn_params)
            self.mca_kc_ggn_syn = mca_ret['synapses']
            self.mca_kc_ggn_netcon = mca_ret['netcons']
            self.mca_kc_ggn_adjacency = mca_ret['adjacency']
        elif kc_ggn_ca_syn_params['regional']:
            lca_kcs = np.random.permutation(self.kcs.lca_kcs)
            lca_ret = make_exp2syn_one_one(lca_kcs, self.ggn.lca_sec, self.ggn.lca_x,
                                           kc_ggn_ca_syn_params)
            mca_kcs = np.random.permutation(self.kcs.mca_kcs)
            mca_ret = make_exp2syn_one_one(mca_kcs, self.ggn.mca_sec, self.ggn.mca_x,
                                           kc_ggn_ca_syn_params)
        else:
            kcs = np.random.permutation(self.kcs.kcs)
            ret = make_exp2syn_one_one(kcs,
                                       np.concatenate((self.ggn.lca_sec,
                                                       self.ggn.mca_sec)),
                                       np.concatenate((self.ggn.lca_x,
                                                       self.ggn.mca_x)),
                                       kc_ggn_ca_syn_params)
            self.kc_ggn_CA_synapses = ret['synapses']
            self.kc_ggn_CA_netcons = ret['netcons']
            self.kc_ggn_CA_adjacency = ret['adjacency']
            return
        self.kc_ggn_LCA_synapses = lca_ret['synapses']
        self.kc_ggn_LCA_netcons = lca_ret['netcons']
        self.kc_ggn_LCA_adjacency = lca_ret['adjacency']
        self.kc_ggn_MCA_synapses = mca_ret['synapses']
        self.kc_ggn_MCA_netcons = mca_ret['netcons']
        self.kc_ggn_MCA_adjacency = mca_ret['adjacency']
        self.kc_ggn_CA_synapses = lca_ret['synapses'] + mca_ret['synapses']
        self.kc_ggn_CA_adjacency = lca_ret['adjacency'] + mca_ret['adjacency']
        self.kc_ggn_CA_netcons = lca_ret['netcons'] + mca_ret['netcons']

    def connect_kc_ggn_alpha(self, syn_params):
        """Connect all KCs in network to segments of GGN in alpha lobe"""
        if (Q_(syn_params['gmax']).to('pS').m <= 0.0) or (len(self.kcs.kcs) == 0):
            self.kc_ggn_alphaL_synapses = []
            self.kc_ggn_alphaL_netcons = []
            self.kc_ggn_alphaL_adjacency = []
            self.ggn_alphaL_input_sections = []
            return
        alpha_sec_pos = nu.select_random_segments_by_sid(
            self.ggn.ggn_graph,
            ng.custom_types['alphaL'],
            len(self.kcs.kcs),
            by_length=True,
            replace=True)
        sections, xpos = list(zip(*alpha_sec_pos))
        ret = make_exp2syn_one_one(self.kcs.kcs, sections, xpos, syn_params)
        self.kc_ggn_alphaL_synapses = ret['synapses']
        self.kc_ggn_alphaL_netcons = ret['netcons']
        self.kc_ggn_alphaL_adjacency = ret['adjacency']
        self.ggn_alphaL_input_sections = list(set(sections))

    def connect_pn_kc(self, syn_params, shared_pn_frac):
        """Creates one synapse on each KC and connects a bunch of presynaptic
        PN VectStims to this synapse.

        There are several configuration options in syn_params to
        control the connectivity rule.

        'clustered': KCs are clustered by spatial location. KCs in
        same cluster share `presyn_frac` of their presynaptic PNs.

        'clustered_pre': not only are the KCs clustered, but so are
        the presynaptic PNs. The PNs are clustered based on their
        starting time after stimulus. KCs in same clustered share all
        the PNs in the corresponding PN cluster, then unresponsive
        PNs, then the next PN cluster in sequence of priority until
        the required number of presynaptic PNs are found.

        Until 2018-06-30 only the pn/shifting option was used and if
        the PNs were shifting, this resulted in there being shifting
        PNs connecting to clustered KCs.

        """
        self.pn_kc_synapses = []
        self.pn_kc_netcons = []
        self.pn_kc_adjacency = []
        if (len(self.kcs.kcs) == 0) or (len(self.pns.pns) == 0):
            return
        start = timer()
        tau1 = Q_(syn_params.get('tau1', '13.333ms')).to('ms').m
        tau2 = Q_(syn_params.get('tau2', '13.333ms')).to('ms').m
        erev = Q_(syn_params.get('e', '0.0mV')).to('mV').m
        thresh = Q_(syn_params.get('threshold', '-20.0mV')).to('mV').m
        delay = Q_(syn_params.get('delay', '0.0ms')).to('ms').m
        delay_std = Q_(syn_params.get('delay_std', '0.0ms')).to('ms').m
        gmax = Q_(syn_params.get('gmax', '4.5pS')).to('uS').m
        if gmax <= 0:
            print('No synapses from PNs to KCs')
            return
        nsyn = int(len(self.pns.pns) * syn_params.get('presyn_frac', 0.5))
        if 'std' in syn_params:
            std = float(syn_params['std'])
            logger.info('PN->KC gmax lognormally distributed with std={} x mean'.format(std))
        else:
            std = 0.0
        if syn_params['clustered']:  # post synaptic KCs in same cluster share more PNs
            label_kcs = defaultdict(list)
            for label, kc in zip(self.ggn.lca_kmeans.labels_, self.kcs.lca_kcs):
                label_kcs['lca_{}'.format(label)].append(kc)
            for label, kc in zip(self.ggn.mca_kmeans.labels_, self.kcs.mca_kcs):
                label_kcs['mca_{}'.format(label)].append(kc)
            if syn_params['clustered_pre']:  # Select presynaptic PNs from clusters
                self._make_shifting_clustered_pn_kc_conn(label_kcs,
                                                         nsyn, shared_pn_frac,
                                                         tau1, tau2, erev, thresh,
                                                         delay, gmax, std=std, delay_std=delay_std)
            else:  # select presynaptic PNs randomly for each KC cluster
                self._make_random_clustered_pn_kc_conn(label_kcs,
                                                       nsyn,
                                                       shared_pn_frac,
                                                       tau1, tau2,
                                                       erev, thresh,
                                                       delay, gmax, std=std, delay_std=delay_std)
        else:  # PN->KC synapses not clustered
            pn_indices = list(range(len(self.pns.pns)))
            start_a = timer()
            for kc in self.kcs.kcs:
                input_pn_indices = np.random.choice(pn_indices,
                                                    size=nsyn, replace=False)
                self._make_pn_kc_synapses(input_pn_indices, kc,
                                          tau1, tau2, erev, thresh,
                                          delay, gmax, std=std, delay_std=delay_std)
            end_a = timer()
            logger.info('Connected {} kcs in {} s'.format(len(self.kcs.kcs), end_a - start_a))
        end = timer()
        pn_kc_syn_count = sum([len(nc_list) for nc_list in self.pn_kc_netcons])
        logger.info('Created {} PN->KC synapses in {} s'.format(pn_kc_syn_count, end - start))


    def _make_shifting_clustered_pn_kc_conn(self, label_kcs, nsyn, shared_pn_frac,
                                            tau1, tau2, erev, thresh,
                                            delay, gmax, std=0, delay_std=0.0):
        """Connect PN->KC based on PN start time.

        In this case, PNs with the same start time are mapped to the
        same KC cluster. If that does not fulfil the 50% connectivity
        requirement, take the rest from the unresponsive set.

        NOTE: 2018-05-16 - this version ignores the unresponsive PNs altogether.

        NOTE: 2018-06-28 - updated to allow including unresponsive PNs after the
        cluster with same index, then next cluster.

        """
        labels, kcs_list = zip(*list(label_kcs.items()))
        pns_list_by_start = sorted(list(self.pns.clusters.items()))
        nshared = int(nsyn * shared_pn_frac + 0.5)
        for ii, kcs in enumerate(kcs_list):
            pns = []
            jj = ii
            start_clus = jj % len(pns_list_by_start)
            start, pnclus = pns_list_by_start[start_clus]
            pns += pnclus
            required = nshared - len(pns)
            if len(self.pns.response_types['U']) > required:
                unresponsive = np.random.choice(self.pns.response_types['U'],
                                                size=required, replace=False)
                pns += list(unresponsive)
            else:
                pns += self.pns.response_types['U']
            logger.info('KC cluster: {}: after adding unresponsive PNs {} more are required in shared pre'.format(labels[ii], nshared - len(pns)))
            while len(pns) < nshared:
                jj += 1
                start, pnclus = pns_list_by_start[jj % len(pns_list_by_start)]
                # logger.info('start: {} pns: {}, required: {}'.format(
                    # start, len(pnclus), nshared - len(pnclus)))
                pns += pnclus
                # logger.info('Connected PNs starting at {} to {}'.format(start, labels[ii]))
            pns = pns[:nshared]
            for pn in pns:
                logger.debug('kcs: {}, pn: {}'.format(ii, pn))
            rest = list(set(range(len(self.pns.pns))) - set(pns))
            for kc in kcs:
                pre = pns + random.sample(rest, nsyn - nshared)
                self._make_pn_kc_synapses(pre, kc,
                                          tau1, tau2, erev, thresh,
                                          delay, gmax, std=std, delay_std=delay_std)

    def _make_random_clustered_pn_kc_conn(self, label_kcs, nsyn, shared_pn_frac,
                                          tau1, tau2, erev, thresh,
                                          delay, gmax, std=0, delay_std=0.0):
        """Connect PN->KC based on KC clusters where KCs in same clustered
        share a certain fraction of PNs."""
        nshared = int(nsyn * shared_pn_frac + 0.5)
        pn_indices = list(range(len(self.pns.pns)))
        for label, kcs in label_kcs.items():
            shared_pn_indices = np.random.choice(pn_indices, size=nshared, replace=False)
            unshared_pn_indices = list(set(pn_indices) - set(shared_pn_indices))
            for kc in kcs:
                input_pn_indices = np.concatenate(
                    (shared_pn_indices,
                     np.random.choice(unshared_pn_indices,
                                      size=nsyn -
                                      nshared,
                                      replace=False)))
                self._make_pn_kc_synapses(input_pn_indices, kc,
                                          tau1, tau2, erev, thresh,
                                          delay, gmax, std=std, delay_std=delay_std)

    def _make_unclustered_pn_kc_conn(self):
        raise NotImplementedError('')

    def _make_pn_kc_synapses(self, pn_indices, kc, tau1, tau2, erev,
                             thresh, delay, gmax, std=0, delay_std=0.0):
        """Connect the presynaptic PN VecStim objects to target KC with synapses."""
        synapse = h.Exp2Syn(kc.soma(0.5))
        synapse.e, synapse.tau1, synapse.tau2 = erev, tau1, tau2
        self.pn_kc_synapses.append(synapse)
        if std > 0:
            # Use lognormal distribution: specified gmax is the
            # mean and specified std is as a fraction of the
            # mean. Thus the mu and sigma parameters become:
            # mu = ln(mean/sqrt(1 + variance/mean^2)
            # and sigma^2 = ln(1 + variance/mean^2)
            # thus, mu = ln(mean) - 0.5 ln(1 + variance/mean^2) = ln(mean) - 0.5 sigma^2
            sigma2 = np.log(1 + std**2)  # std is specified as a fraction of mean
            mu = np.log(gmax) - sigma2 / 2.0
            _gmax = np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2),
                                        size=len(pn_indices))
        else:
            _gmax = [gmax] * len(pn_indices)
        assert delay >= 0.0
        if delay_std > 0:
            _delay = np.random.normal(delay, delay_std, size=len(pn_indices))
        else:
            _delay = [delay] * len(pn_indices)
        _delay[_delay < 0] = delay    # replace negative delays with mean delay
        netcons = [h.NetCon(self.pns.vecstims[ii],
                            synapse, thresh, _delay[jj], _gmax[jj], sec=kc.soma)
                   for jj, ii in enumerate(pn_indices)]
        self.pn_kc_netcons.append(netcons)
        self.pn_kc_adjacency += [(self.pns.pns[ii], kc.soma.name())
                                 for ii in pn_indices]

    def setup_data_recording(self, params, t=None):
        if hasattr(self, 'data_dict'):
            raise Exception('Data recording already setup')
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
        data = {'time': tvec}
        # PNs and KCs
        if self.ig is not None:
            if isinstance(self.ig, ig.IGSpikes):
                data['ig_spikes'] = Q_(np.array(self.ig.stimvec.x), 'ms')
            elif isinstance(self.ig, ig.IGIzhi):
                data['ig_vm'] = h.Vector()
                data['ig_vm'].record(self.ig.izhi._ref_V, t)
        data['pn_output'] = self.pns.outputs
        data['kc_spikes'] = []
        self.kc_spike_recorders = []
        for kc in self.kcs.kcs:
            nc, stvec = nu.record_spiketimes(kc.soma)
            data['kc_spikes'].append(stvec)
            self.kc_spike_recorders.append(nc)
        if params['kc']['n_vm'] >= len(self.kcs.kcs):
            vm_kcs = self.kcs.kcs
        else:
            vm_kcs = np.random.choice(self.kcs.kcs, size=params['kc']['n_vm'], replace=False)
        data['kc_vm'] = [nu.setup_recording(kc.soma, 0.5, 'v', t)
                         for kc in vm_kcs]
        data['kc_vm_sections'] = [kc.soma for kc in vm_kcs]
        # GGN input in alphaL
        # data['ggn_alphaL_input_vm'] = nu.setup_recording(self.ggn_alphaL_input_sections)
        data['ggn_alphaL_input_vm'] = [nu.setup_recording(sec, 0.5, 'v', t)
                                       for sec in self.ggn_alphaL_input_sections]
        data['ggn_alphaL_input_sections'] = self.ggn_alphaL_input_sections
        # GGN Vm in CA
        sections = self.ggn.lca_sec + self.ggn.mca_sec
        sec_list = list(set(sections))
        data['ggn_output_vm'] = [nu.setup_recording(sec, 0.5, 'v', t) for sec in sec_list]
        data['ggn_output_sections'] = sec_list
        self.data_dict = data
        # GGN Vm in basal dendrite
        data['ggn_basal_vm'] = [nu.setup_recording(sec, 0.5, 'v', t) for sec in self.ggn.ggn.basal]
        data['ggn_basal_sections'] = [sec for sec in self.ggn.ggn.basal]
        return self.data_dict
# class pn_kc_ggn_network ends here


def connect_in_cluster(labels, kcs, sections, xpos, syn_params):
    """Create 1-1 connections from kcs to GGN sections at xpos within same
    cluster as identified by label at corresponding index

    KCs with same cluster label are shuffled before making synapses.

    """
    kc_label_map = defaultdict(list)
    sec_label_map = defaultdict(list)
    pos_label_map = defaultdict(list)
    kc_ggn_adjacency = []
    kc_ggn_synapses = []
    kc_ggn_netcons = []
    for label, sec, x, kc in zip(labels, sections, xpos, kcs):
        sec_label_map[label].append(sec)
        pos_label_map[label].append(x)
        kc_label_map[label].append(kc)
    for label in kc_label_map:
        sec = sec_label_map[label]
        pos = pos_label_map[label]
        kcs = kc_label_map[label]
        # Note: shuffle modifies the sequence in place, OK here
        # but just to be sure use permutation which copies data
        kcs = np.random.permutation(kcs)
        ret = make_exp2syn_one_one(kcs, sec, pos, syn_params)
        kc_ggn_synapses += ret['synapses']
        kc_ggn_netcons += ret['netcons']
        kc_ggn_adjacency += ret['adjacency']
    return {'synapses': kc_ggn_synapses, 'netcons':
            kc_ggn_netcons, 'adjacency': kc_ggn_adjacency}


def make_exp2syn_one_one(kcs, sections, xpos, syn_params):
    adjacency = []
    synapses = []
    netcons = []
    kg_tau1 = Q_(syn_params.get('tau1', '1.0ms')).to('ms').m
    kg_tau2 = Q_(syn_params.get('tau2', '1.0ms')).to('ms').m
    kg_e = Q_(syn_params.get('e', '0.0mV')).to('mV').m
    kg_thresh = Q_(syn_params.get('threshold', '-20.0mV')).to('mV').m
    kg_delay = Q_(syn_params.get('delay', '1.0ms')).to('ms').m
    kg_gmax = Q_(syn_params.get('gmax', '1e-3uS')).to('uS').m
    std = float(syn_params.get('std', 0.0))
    if std > 0:
        logger.info('KC->GGN gmax distributed lognormally')
        # Use lognormal distribution: specified gbar is the
        # mean and specified std is as a fraction of the
        # mean. Thus the mu and sigma parameters become:
        # mu = ln(mean/sqrt(1 + variance/mean^2)
        # and sigma^2 = ln(1 + variance/mean^2)
        # thus, mu = ln(mean) - 0.5 ln(1 + variance/mean^2) = ln(mean) - 0.5 sigma^2
        sigma2 = np.log(1 + std**2)  # std is specified as a fraction of mean
        mu = np.log(kg_gmax) - sigma2 / 2.0
        gmax = np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2),
                                         size=len(kcs))
    else:
        gmax = [kg_gmax] * len(kcs)

    if kg_gmax > 0.0:  # Don't create synapses unless necessary
        for kc, sec, pos, g in zip(kcs, sections, xpos, gmax):
            syn = h.Exp2Syn(sec(pos))
            syn.tau1, syn.tau2, syn.e = kg_tau1, kg_tau2, kg_e
            nc = h.NetCon(kc.soma(0.5)._ref_v, syn, kg_thresh, kg_delay,
                          g, sec=kc.soma)
            synapses.append(syn)
            netcons.append(nc)
            adjacency.append((kc.soma.name(), sec.name(), pos))
    return {'synapses': synapses, 'netcons': netcons,
            'adjacency': adjacency}


def cluster_positions(sec_pos_list, nclust, random_state=None, fake=False):
    """For a list of sections and corresponding segment positions [(sec,
    x),...] `sec_pos_list`, obtain a spatial clustering of the
    segments in 3D.

    """
    print('secpos list len : {}'.format(len(sec_pos_list)))
    print('nclust: {}'.format(nclust))
    sec_pos_dict = defaultdict(list)
    for sec, pos in sec_pos_list:
        sec_pos_dict[sec].append(pos)
    # These are the expanded lists - the original ordering will be
    # destroyed in the dict
    _sec_list = []
    _x_list = []
    _pos_list = []
    for sec, xlist in sec_pos_dict.items():
        p3d_list = nu.get_seg_3dpos(sec, xlist)

        for x, pos in zip(xlist, p3d_list):
            _sec_list.append(sec)
            _x_list.append(x)
            _pos_list.append(pos)
    if (nclust > len(sec_pos_list)) or (nclust == 0):
        nclust = len(sec_pos_list)
    if nclust == 0:
        nclust = 1
    if fake:
        cluster_size = int(round(len(sec_pos_list) * 1.0 / nclust))
        labels = np.concatenate([[ii] * cluster_size for ii in range(nclust)])
        diff = len(labels) - len(sec_pos_list)
        if diff < 0:
            labels = np.concatenate((labels, [nclust] * (-diff)))
        elif diff > 0:
            labels = labels[:-diff].copy()
        np.random.shuffle(labels)
        kmeans = FakeKMeans(labels)
    else:
        kmeans = cluster.KMeans(n_clusters=nclust,
                                random_state=random_state).fit(np.array(_pos_list))
    return kmeans, _sec_list, _x_list, _pos_list


#################################
# Code for saving data
#################################

def get_model_files(model_dir, subdirs):
    """Return a list of model files under model_dir, including only
    directories listed in subdirs."""
    model_files = []
    for subdir in subdirs:
        dirpath = os.path.join(model_dir, subdir)
        for path, dirs, files in os.walk(dirpath):
            if '.git' in path:
                continue
            for filename in files:
                if filename.endswith('.py') or  \
                   filename.endswith('.hoc') or  \
                   filename.endswith('.mod') or  \
                   filename.endswith('.yaml'):
                    model_files.append(os.path.join(path, filename))
    return model_files


def prerun_save(params):
    """Save model files before starting simulation"""
    outfile = os.path.join(
        params['output']['directory'],
        'mb_net_UTC{}-PID{}-JID{}.h5'.format(
            timestamp.strftime('%Y_%m_%d__%H_%M_%S'),
            mypid, myjobid))
    logger.info('Saving model data in {} before starting simulation'.format(outfile))
    writer = nsdf.NSDFWriter(outfile, dialect=nsdf.dialect.ONED, mode='w', compression='gzip')
    model_dir = os.path.join(os.path.dirname(__file__), '..', '..')
    subdirs = ['common', 'mb', 'morphutils', 'nrn']
    model_files = get_model_files(model_dir, subdirs)
    writer.add_model_filecontents(model_files, model_dir)
    writer._fd.attrs['config'] = yaml.dump(params)
    writer.title = 'PN->KC<->GGN network model'
    writer.creator = 'Subhasis Ray (ray dot subhasis at gmail dot com)'
    writer.software = nrn_version
    writer.description = 'PN->KC<->GGN model with KC->GGN' \
                         'synapses in alpha lobe {}.'.format(
                             'and calyx'
                             if Q_(params['kc_ggn_CA_syn']['gmax']).to('pS').m > 0
                             else 'only')
    writer._fd.attrs['command'] = ' '.join(sys.argv)
    logger.info('Model data saved')
    return writer


def save_cluster_labels(writer, region, labels, sections, xpos, pos3d, kcs):
    """Region `lca` or `mca`, everything else corresponds to the region"""
    cluster_info_dtype = np.dtype([('label', np.int32),
                                   ('sec', nsdf.VLENSTR),
                                   ('pos', np.float64),
                                   ('x', np.float64),
                                   ('y', np.float64),
                                   ('z', np.float64)])
    cluster_labels = nsdf.StaticData('{}_cluster_labels'.format(region),
                                     unit='',
                                     dtype=cluster_info_dtype)
    for ii, label in enumerate(labels):
        sec = sections[ii]
        pos = xpos[ii]
        x, y, z = pos3d[ii]
        cluster_labels.put_data('{}_pos_{}'.format(region, ii),
                                (label, sec.name(), pos, x, y, z))
    cluster_labels_ds = writer.add_static_ds(
        '{}_cluster_labels'.format(region),
        cluster_labels.get_sources())
    ds = writer.add_static_data(cluster_labels_ds, cluster_labels)
    ds.attrs['note'] = 'Each entry corresponds to a segment ' \
        '(specified by section and pos) and the ' \
        'nearest preceding 3D point position. These have ' \
        '1-1 correspondence to GGN->KC synapses ' \
        'and the KCs themselves'
    kc_cluster_labels = nsdf.StaticData('{}_kc_cluster_labels'.format(region),
                                        unit='',
                                        dtype=cluster_info_dtype)
    sources = []
    for ii, label in enumerate(labels):
        sec = kcs[ii].soma
        pos = 0.5
        x, y, z = pos3d[ii]
        sources.append('{}_kc_pos_{}'.format(region, ii))
        kc_cluster_labels.put_data(sources[-1],
                                    (label, sec.name(), pos, x, y, z))
    kc_cluster_labels_ds = writer.add_static_ds(
        '{}_kc_cluster_labels'.format(region),
        sources)
    ds = writer.add_static_data(kc_cluster_labels_ds, kc_cluster_labels)
    ds.attrs['note'] = 'Each entry corresponds to a KC ' \
                       '(specified by section and pos) and the ' \
                       'nearest preceding 3D point position of GGN. These have ' \
                       '1-1 correspondence to GGN->KC synapses '


def save_pn_info(writer, pn_population):
    """Save information about the PNs as clustered by start time and
    response class"""
    if len(pn_population.clusters) == 0:
        return
    pn_start_time_dtype = np.dtype([('start_time', np.float64),
                                    ('PN', np.int32)])
    pn_st = nsdf.StaticData('pn_start_time', unit='', dtype=pn_start_time_dtype)
    ii = 0
    for start_time, pns in pn_population.clusters.items():
        for pn in pns:
            pn_st.put_data(str(ii), (start_time, pn))
            ii += 1
    pn_st_ds = writer.add_static_ds('pn_start_time', pn_st.get_sources())
    writer.add_static_data(pn_st_ds, pn_st)
    pn_types_dtype = np.dtype([('PN', np.int32), ('class', nsdf.VLENSTR)])
    pn_type = nsdf.StaticData('pn_class', unit='', dtype=pn_types_dtype)
    for response, pns in pn_population.response_types.items():
        for pn in pns:
            pn_type.put_data(str(pn), (pn, response))
    pn_type_ds = writer.add_static_ds('pn_class', pn_type.get_sources())
    writer.add_static_data(pn_type_ds, pn_type)


def save_synapses(writer, model):
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
    if len(model.pn_kc_synapses) > 0:
        pn_kc_syn_data = np.empty((len(model.pn_kc_adjacency),),
                                  dtype=conn_dtype)
        ii = 0
        for syn, netcons in zip(model.pn_kc_synapses, model.pn_kc_netcons):
            for nc in netcons:
                pn_kc_syn_data['gmax'][ii] = nc.weight[0]
                pn_kc_syn_data['e'][ii] = syn.e
                pn_kc_syn_data['tau1'][ii] = syn.tau1
                pn_kc_syn_data['tau2'][ii] = syn.tau2
                pn_kc_syn_data['pre'][ii] = model.pn_kc_adjacency[ii][0]
                pn_kc_syn_data['post'][ii] = model.pn_kc_adjacency[ii][-1]
                pn_kc_syn_data['prepos'][ii] = 0
                pn_kc_syn_data['postpos'][ii] = 0.5
                pn_kc_syn_data['delay'][ii] = nc.delay
                ii += 1
        pn_kc_syn_src = ['{}__{}'.format(*prepost)
                         for prepost in model.pn_kc_adjacency]
        conn_unit = ['', '', '', '', 'uS', 'mV', 'ms', 'ms', 'ms']
        pn_kc_syn = nsdf.StaticData('pn_kc_synapse',
                                    unit=conn_unit,
                                    dtype=conn_dtype)
        for src, data in zip(pn_kc_syn_src, pn_kc_syn_data):
            pn_kc_syn.put_data(src, data)
            # print(src)
        # print(len(pn_kc_syn_src), len(set(pn_kc_syn_src)))
        pn_kc_syn_ds = writer.add_static_ds('pn_kc_synapse', pn_kc_syn_src)
        writer.add_static_data(pn_kc_syn_ds, pn_kc_syn)
    ###############################################
    # Create synapse dataset for KC->GGN
    ###############################################
    ## alphaL
    if len(model.kc_ggn_alphaL_synapses) > 0:
        kc_ggn_syn_data = nsdf.StaticData('kc_ggn_alphaL_synapse',
                                          unit=conn_unit,
                                          dtype=conn_dtype)
        sources = []
        for prepost, syn, nc in zip(model.kc_ggn_alphaL_adjacency,
                                    model.kc_ggn_alphaL_synapses,
                                    model.kc_ggn_alphaL_netcons):
            sources.append('{}__{}'.format(*prepost))
            kc_ggn_syn_data.put_data(sources[-1],
                                     (prepost[0], prepost[1], 0.5, prepost[2],
                                      nc.weight[0],
                                      syn.e, syn.tau1,
                                      syn.tau2, nc.delay))
        kc_ggn_syn_ds = writer.add_static_ds('kc_ggn_alphaL_synapse',
                                             sources)
        writer.add_static_data(kc_ggn_syn_ds, kc_ggn_syn_data)
    ## CA
    sources = []
    if hasattr(model, 'kc_ggn_CA_synapses') and (len(model.kc_ggn_CA_synapses) > 0):
        kc_ggn_syn_data = nsdf.StaticData('kc_ggn_CA_synapse',
                                          unit=conn_unit,
                                          dtype=conn_dtype)
        for prepost, syn, nc in zip(model.kc_ggn_CA_adjacency,
                                    model.kc_ggn_CA_synapses,
                                    model.kc_ggn_CA_netcons):
            sources.append('{}__{}'.format(*prepost))
            kc_ggn_syn_data.put_data(sources[-1],
                                    (prepost[0], prepost[1], 0.5, prepost[2],
                                     nc.weight[0],
                                     syn.e, syn.tau1,
                                     syn.tau2))
        kc_ggn_syn_ds = writer.add_static_ds('kc_ggn_CA_synapse',
                                             sources)
        writer.add_static_data(kc_ggn_syn_ds, kc_ggn_syn_data)
    ###############################################
    # Create synapse dataset for KC->IG
    ###############################################
    if (model.ig is not None) and (len(model.kc_ig_syn_info) > 0):
        conn_unit = ['', '', '', '', '', 'ms']
        conn_dtype = np.dtype([('pre', nsdf.VLENSTR),
                               ('post', nsdf.VLENSTR),
                               ('prepos', np.float64),
                               ('postpos', np.float64),
                               ('gmax', np.float64),
                               ('delay', np.float64)])
        kc_ig_syn_data = nsdf.StaticData('kc_ig_synapse',
                                     unit=conn_unit,
                                     dtype=conn_dtype)
        sources = []
        for row in model.kc_ig_syn_info:
            sources.append('{}__{}'.format(row['pre'], row['post']))
            kc_ig_syn_data.put_data(sources[-1],
                                    (row['pre'], row['post'], row['prepos'], row['postpos'],
                                     row['gmax'], row['delay']))
        kc_ig_syn_ds = writer.add_static_ds('kc_ig_synapse',
                                            sources)
        writer.add_static_data(kc_ig_syn_ds, kc_ig_syn_data)

    ###############################################
    # Save GGN->KC synapses
    ###############################################
    conn_dtype = np.dtype([('pre', nsdf.VLENSTR),
                           ('post', nsdf.VLENSTR),
                           ('prepos', np.float64),
                           ('postpos', np.float64),
                           ('vmid', np.float64), ('vslope', np.float64),
                           ('e', np.float64), ('gbar', np.float64),
                           ('tau', np.float64)])
    conn_unit = ['', '', '', '', 'mV', 'mV', 'mV', 'uS', 'ms']
    if len(model.ggn_kc_synapses) > 0:
        ggn_kc_syn_data = nsdf.StaticData('ggn_kc_synapse', unit=conn_unit,
                                          dtype=conn_dtype)
        sources = []
        for prepost, syn in zip(model.ggn_kc_adjacency,
                                model.ggn_kc_synapses):
            sources.append('{}__{}'.format(prepost[0], prepost[1]))
            ggn_kc_syn_data.put_data(sources[-1],
                                     (prepost[0], prepost[1], prepost[2], 0.5,
                                      syn.vmid, syn.vslope, syn.e, syn.gbar,
                                      syn.tau))
        ggn_kc_syn_ds = writer.add_static_ds('ggn_kc_synapse',
                                             ggn_kc_syn_data.get_sources())
        writer.add_static_data(ggn_kc_syn_ds, ggn_kc_syn_data)
    end = timer()
    logger.info('Saved synapse data in {} s'.format(end - start))


def save_model(writer, net_model):
    """Save the model information from net_model"""
    logger.info('Saving model')
    start = timer()
    if len(net_model.kcs.lca_kcs) > 0:
        save_cluster_labels(writer, 'lca', net_model.ggn.lca_kmeans.labels_,
                            net_model.ggn.lca_sec, net_model.ggn.lca_x, net_model.ggn.lca_pos3d,
                            net_model.kcs.lca_kcs)
    if len(net_model.kcs.mca_kcs) > 0:
        save_cluster_labels(writer, 'mca', net_model.ggn.mca_kmeans.labels_,
                            net_model.ggn.mca_sec, net_model.ggn.mca_x, net_model.ggn.mca_pos3d,
                            net_model.kcs.mca_kcs)
    model = nsdf.ModelComponent('olf', uid='olf')
    pn_grp = nsdf.ModelComponent('pn', uid='pn', parent=model)
    pns = [nsdf.ModelComponent(pn, uid=pn, parent=pn_grp)
           for pn in net_model.pns.pns]
    kc_grp = nsdf.ModelComponent('kc', uid='kc', parent=model)
    kcs = [nsdf.ModelComponent(kc.soma.name(), uid=kc.soma.name(),
                               parent=kc_grp)
           for kc in net_model.kcs.kcs]
    ggn_name = net_model.ggn.ggn.hname()
    ggn = nsdf.ModelComponent(ggn_name,
                              uid=ggn_name, parent=model)
    for sec in net_model.ggn.ggn.all:
        node = nsdf.ModelComponent(sec.name(), uid=sec.name(), parent=ggn,
                                   attrs={'RA': sec.Ra,
                                          'g_pas': sec.g_pas,
                                          'cm': sec.cm,
                                          'e_pas': sec.e_pas})
    writer.add_modeltree(model)
    save_synapses(writer, net_model)
    save_pn_info(writer, net_model.pns)
    end = timer()
    logger.info('Saved model including synapse data in {} s'.format(end - start))


def save_recorded_data(writer, model, data_dict, params):
    logger.info('Saving data in {}'.format(writer.filename))
    start = timer()
    dt = data_dict['time'].x[1] - data_dict['time'].x[0]
    # PN output spike times
    pn_spikes = data_dict['pn_output']
    if len(pn_spikes) > 0:
        pn_events = nsdf.EventData('pn_spiketime', unit='ms', dtype=np.float64)
        for src, spiketimes in zip(model.pns.pns, pn_spikes):
            pn_events.put_data(src, spiketimes.to('ms').m)
        pn_ds = writer.add_event_ds_1d('pn', 'pn_spiketime',
                                       pn_events.get_sources())
        writer.add_event_1d(pn_ds, pn_events)
    # KC spike times
    kc_spikes = data_dict['kc_spikes']
    if len(kc_spikes) > 0:
        kc_events = nsdf.EventData('kc_spiketime', unit='ms', dtype=np.float64)
        sources = []
        for kc, spiketimes in zip(model.kcs.kcs, kc_spikes):
            sources.append(kc.soma.name())
            kc_events.put_data(sources[-1], np.array(spiketimes.x))
        kc_spikes_ds = writer.add_event_ds_1d('kc', 'kc_spiketime',
                                              sources)
        writer.add_event_1d(kc_spikes_ds, kc_events)
    # KC Vm
    if len(data_dict['kc_vm']) > 0:
        kc_vm = nsdf.UniformData('KC_Vm', unit='mV', field='v')
        kc_vm.set_dt(dt, 'ms')
        sources = []
        for kc_soma, vm_tab in zip(data_dict['kc_vm_sections'], data_dict['kc_vm']):
            sources.append(kc_soma.name())
            kc_vm.put_data(sources[-1], np.array(vm_tab.x))
        kc_vm_ds = writer.add_uniform_ds('kc', sources)
        writer.add_uniform_data(kc_vm_ds, kc_vm)
    # Handle IG data
    if 'ig_spikes' in data_dict:
        ig_spikes = data_dict['ig_spikes']
        ig_events = nsdf.EventData('ig_spiketime', unit='ms', dtype=np.float64)
        ig_events.put_data('IG', ig_spikes.to('ms').m)
        ig_ds = writer.add_event_ds_1d('ig', 'ig_spiketime',
                                       ig_events.get_sources())
        writer.add_event_1d(ig_ds, ig_events)
    elif 'ig_vm' in data_dict:
        ig_vm = data_dict['ig_vm']
        ig_vm_data = nsdf.UniformData('IG_Vm', unit='mV', field='V')
        ig_vm_data.set_dt(dt, 'ms')
        ig_vm_data.put_data('IG', np.array(ig_vm.x))
        ig_vm_ds = writer.add_uniform_ds('ig', ['IG'])
        writer.add_uniform_data(ig_vm_ds, ig_vm_data)
    # GGN alphaLobe Vm
    ggn_input_vm = nsdf.UniformData('GGN_alphaL_input_Vm', unit='mV', field='v')
    ggn_input_vm.set_dt(dt, 'ms')
    sources = []
    for sec, vm_tab in zip(data_dict['ggn_alphaL_input_sections'],
                           data_dict['ggn_alphaL_input_vm']):
        sources.append(sec.name())
        ggn_input_vm.put_data(sources[-1], np.array(vm_tab.x))
    if len(sources) > 0:
        ggn_input_vm_ds = writer.add_uniform_ds('ggn_alphaL_input',
                                                sources)
        writer.add_uniform_data(ggn_input_vm_ds, ggn_input_vm)
    else:
        logger.info('No data recorded from ggn input sections')
    ggn_output_vm = nsdf.UniformData('GGN_output_Vm', unit='mV', field='v')
    ggn_output_vm.set_dt(dt, 'ms')
    sources = []
    for sec, vm_tab in zip(data_dict['ggn_output_sections'], data_dict['ggn_output_vm']):
        sources.append(sec.name())
        ggn_output_vm.put_data(sources[-1], np.array(vm_tab.x))
    if len(sources) > 0:
        ggn_output_vm_ds = writer.add_uniform_ds('ggn_output', ggn_output_vm.get_sources())
        writer.add_uniform_data(ggn_output_vm_ds, ggn_output_vm)
    else:
        logger.info('No data recorded from ggn output sections')
    ggn_basal_vm = nsdf.UniformData('GGN_basal_Vm', unit='mV', field='v')
    ggn_basal_vm.set_dt(dt, 'ms')
    sources = []
    for sec, vm_tab in zip(data_dict['ggn_basal_sections'], data_dict['ggn_basal_vm']):
        sources.append(sec.name())
        ggn_basal_vm.put_data(sources[-1], np.array(vm_tab.x))
    if len(sources) > 0:
        ggn_basal_vm_ds = writer.add_uniform_ds('ggn_basal',
                                                sources)
        writer.add_uniform_data(ggn_basal_vm_ds, ggn_basal_vm)
    else:
        logger.info('No data recorded from ggn basal sections')
    end = timer()
    logger.info('Saved recorded data in {} s'.format(end - start))

#################################
# Main simulation code
#################################


def run_simulation(params, writer, recstep=5, dry=False):
    """Record data every `recstep` to reduce storage space."""
    model = pn_kc_ggn_network(params)
    print('pn_kc_ggn_network built')
    sys.stdout.flush()
    onset = Q_(params['stimulus']['onset'])
    duration = Q_(params['stimulus']['duration'])
    tail = Q_(params['stimulus']['tail'])
    h.tstop = (onset + duration + tail).to('ms').m
    save_model(writer, model)
    print('model saved')
    sys.stdout.flush()
    start = datetime.utcnow()
    data = {}
    if not dry:
        t = np.arange(0, h.tstop, h.dt * recstep)
        tvec = h.Vector(t)
        data = model.setup_data_recording(params, tvec)
        print('data recording set up')
        sys.stdout.flush()
        h.init()
        logger.info('Finished init. Starting simulation of {} ms'.format(
            h.tstop))
        nu.block_run(logger=logger)
    else:
        logger.info('Dry run. Skipping simulation')
    end = datetime.utcnow()
    delta = end - start
    data['tstart'] = start
    data['tend'] = end
    logger.info('Finished running simulation in {} s'.format(
        delta.days * 86400 +
        delta.seconds +
        delta.microseconds * 1e-6))
    return model, data


def run_and_save(params, dry=False):
    logger.info('About to save')
    writer = prerun_save(params)
    logger.info('Done prerun save')
    model, data = run_simulation(params, writer, dry=dry)
    print('finished simulation')
    sys.stdout.flush()
    if not dry:
        save_recorded_data(writer, model, data, params)
    print('saved data')
    sys.stdout.flush()
    writer.tstart = data['tstart']
    writer.tend = data['tend']


def make_parser():
    """

    NOTE: Sat Jun 30 19:50:49 EDT 2018 updates

    Removed --ca option as it is redundant with kc_ggn_CA_syn/gmax is set
              to 0 or positive.
    Renamed --clustered_pn to `--pn_kc_clustered` for setting
              pn_kc_syn/clustered.
    Added   --kc_ggn_ca_gmax and --kc_ggn_alpha_gmax to override config options
              `kc_ggn_CA_syn/gmax` and `kc_ggn_alphaL_syn/gmax`.
    """

    parser = argparse.ArgumentParser(
        description='Run a simulation of with spatially localized' +
        ' connection between KCs and GGN')
    # parser.add_argument('--clustered_pn', action='store_true')
    parser.add_argument('--RA', default=0.0,
                        help='Specific axial resistance in ohm-cm')
    parser.add_argument('--RM', default=0.0,
                        help='Specific membrane resistance in kohm-cm^2')
    # parser.add_argument('--ca', action='store_true',
    #                     help='Create KC->GGN synapses in Calyx')
    parser.add_argument('--regional', action='store_true',
                        help='Keep KC<->GGN connections in calyx '
                        'local in LCA and MCA')
    parser.add_argument('--shifting_pn', dest='shifting_pn', action='store_true',
                        help='If activity start should shift through PN population'
                             ' (overrides config)')
    parser.add_argument('--pn_kc_clustered_pre', dest='pn_kc_clustered_pre',
                        action='store_true', help='whether to respect PN '
                        ' clusters when connecting pn->kc. Has no effect'
                        ' if pn_kc_clustered is False')
    parser.add_argument('--kc_ggn_clustered', dest='kc_ggn_clustered', action='store_true',
                        help='If KCs connect to the same region of GGN they receive inhibition from'
                             ' (overrides config). '
                        'Results in the KC<->GGN connections being regional anyway')
    parser.add_argument('--pn_kc_clustered', dest='pn_kc_clustered', action='store_true',
                        help='If KCs in same cluster should share `kc/shared_pn_frac` of PNs'
                        '. Override config')
    parser.add_argument('--kc_ggn_ca_gmax', dest='kc_ggn_ca_gmax', default='-1pS',
                        help='KC->GGN gmax in CA. Overrides config')
    parser.add_argument('--kc_ggn_alpha_gmax', dest='kc_ggn_alpha_gmax', default='-1pS',
                        help='KC->GGN gmax in alphaL. Overrides config')
    parser.add_argument('--pn_per_kc', dest='pn_per_kc', type=int, default=0,
                        help='number of PN per KC (overrides config)')
    parser.add_argument('--ggn_kc_gmax', type=str, default='-1nS',
                        help='ggn->kc gmax, if negative use config.yaml')
    parser.add_argument('--ggn_kc_gmax_weak', type=str, default='-1nS',
                        help='ggn->kc weak gmax, if negative use config.yaml')
    parser.add_argument('--kc_frac_weak_inh', type=float, default=-1,
                        help='fraction of KCs receiving weak GGN inhibition, '
                        'if negative use config.yaml')
    parser.add_argument('--pn_kc_gmax', type=str, default='-1pS',
                        help='pn->kc gmax, if negative use config.yaml')
    parser.add_argument('--dryrun', dest='dryrun', action='store_true',
                        help='number of PN per KC (overrides config)')
    parser.add_argument('--npseed', type=int, dest='npseed', help='numpy rngseed')
    parser.add_argument('--kmseed', type=int, dest='kmseed', help='kmeans rngseed')
    parser.add_argument('--fake', dest='fake', action='store_true',
                        help='if specified kc clustering is fake')
    parser.add_argument('--ig_ggn_gmax', type=str, default='-1pS',
                        help='ig->ggn gmax, if negative use config.yaml')
    parser.add_argument('--debug', action='store_true',
                        help='Set log level to debug. Default is info')
    parser.add_argument('--kc_ig_weight', dest='kc_ig_weight', type=float, default=-1.0,
                        help='KC->IG synaptic weight. Overrides config')
    parser.add_argument('--kc_clus_size', dest='kc_clus_size', type=int, default=1, help='cluster size for connecting PN->KC in inclusters')
    parser.add_argument('--kc_clus_shared_pn', dest='kc_clus_shared_pn', type=float, default=-1, help='percentage of PNs shared by KCs in the same cluster')
    parser.add_argument('--pn_kc_delay_std', dest='pn_kc_delay_std', type=str, default='-1ms',
                        help='if specified, and of positive magnitude, PN->KC delays will be'
                        ' sampled from a Gaussian distribution with this std')
    parser.add_argument('--pn_kc_delay_mean', dest='pn_kc_delay_mean', type=str, default='-1ms',
                        help='if specified, and of positive magnitude, PN->KC delays will be'
                        ' with this mean value')

    return parser


if __name__ == '__main__':
    args = make_parser().parse_args()
    z_start = timer()
    print('COMMAND LINE: {}'.format(' '.join(sys.argv)))
    print('Logging in', logger.parent.handlers[0].baseFilename)
    print('Parsed arguments', args)
    sys.stdout.flush()
    if args.debug:
        logger.setLevel(10)
    else:
        logger.setLevel(20)
    logger.info('COMMAND LINE: {}'.format(' '.join(sys.argv)))
    # netconf['kc_ggn_CA_syn']['present'] = args.ca
    # netconf['pn_kc_syn']['clustered'] = args.clustered_pn
    if args.RA:
        netconf['ggn']['RA'] = str(Q_(args.RA, 'ohm*cm'))
    if args.RM:
        netconf['ggn']['RM'] = str(Q_(args.RM, 'kohm*cm**2'))
    if Q_(args.kc_ggn_ca_gmax).to('pS').m >= 0.0:
        netconf['kc_ggn_CA_syn']['gmax'] = args.kc_ggn_ca_gmax
    if Q_(args.kc_ggn_alpha_gmax).to('pS').m >= 0.0:
        netconf['kc_ggn_alphaL_syn']['gmax'] = args.kc_ggn_alpha_gmax
    if Q_(args.ggn_kc_gmax).to('nS').m >= 0:
        netconf['ggn_kc_syn']['gmax'] = args.ggn_kc_gmax
    if Q_(args.pn_kc_gmax).to('nS').m >= 0:
        netconf['pn_kc_syn']['gmax'] = args.pn_kc_gmax
    if args.pn_per_kc:
        netconf['pn_kc_syn']['presyn_frac'] = args.pn_per_kc * 1.0 / netconf['pn']['number']
    if args.shifting_pn:
        netconf['pn']['shifting'] = True
    if args.regional:
        netconf['kc_ggn_CA_syn']['regional'] = True
    if args.kc_ggn_clustered:
        netconf['kc_ggn_CA_syn']['clustered'] = True
    if args.pn_kc_clustered:
        netconf['pn_kc_syn']['clustered'] = True
    if args.pn_kc_clustered_pre:
        netconf['pn_kc_syn']['clustered_pre'] = True
    if args.npseed:
        netconf['rngseed'] = args.npseed
    if args.kmseed:
        netconf['kc']['cluster_rs'] = args.kmseed
    if args.kc_frac_weak_inh > 0:
        netconf['ggn_kc_syn']['frac_weakly_inhibited'] = args.kc_frac_weak_inh
    if args.ggn_kc_gmax_weak > 0:
        netconf['ggn_kc_syn']['gmax_weakly_inhibited'] = args.ggn_kc_gmax_weak
    if Q_(args.ig_ggn_gmax).to('pS').m >= 0:
        netconf['ig_ggn_syn']['gmax'] = args.ig_ggn_gmax
    if args.kc_ig_weight >= 0:
        netconf['kc_ig_syn']['weight'] = args.kc_ig_weight
    netconf['kc']['fake_clusters'] = args.fake
    print('Fake? {}'.format(args.fake))
    if args.kc_clus_size > 1:
        netconf['kc']['cluster_size'] = args.kc_clus_size
    else:
        netconf['pn_kc_syn']['clustered'] = False
    if args.kc_clus_shared_pn > 0:
        assert args.kc_clus_shared_pn <= 1.0
        netconf['kc']['shared_pn_frac'] = args.kc_clus_shared_pn
    if Q_(args.pn_kc_delay_mean).to('ms').m >= 0:
        netconf['pn_kc_syn']['delay'] = args.pn_kc_delay_mean
    if Q_(args.pn_kc_delay_std).to('ms').m >= 0:
        netconf['pn_kc_syn']['delay_std'] = args.pn_kc_delay_std

    sys.stdout.flush()
    if 'rngseed' in netconf:
        np.random.seed(netconf['rngseed'])
        logger.info('NUMPY RNG set to {}'.format(netconf['rngseed']))
        logger.info('First rng sample: {}'.format(np.random.sample()))
    run_and_save(netconf, dry=args.dryrun)
    z_end = timer()
    logger.info('Start to finish {} s'.format(z_end - z_start))

#
# pn_kc_ggn_network.py ends here
