# fig_kc_removal.py --- 
# Author: Subhasis  Ray
# Created: Fri Aug 16 11:27:37 2019 (-0400)
# Last-Updated: Tue Aug 20 13:23:15 2019 (-0400)
#           By: Subhasis Ray
# Version: $Id$

# Code:
"""This summarizes the results from series of KC removals.


In these simulations I first removed KCs that generated 6 or more spikes from the model and simulated repeatedly.

Once there were no more KCs to remove, I started removing KCs that
spiked 4 or more spikes.

Finally, I removed KCs that spiked at all.

"""
import os
import numpy as np
import h5py as h5
import pandas as pd
import network_data_analysis as nda
from matplotlib import pyplot as plt


DATA_DIR = 'D:/biowulf_stage/fixed_net/'
TEMPLATE_DIR = 'D:/biowulf_stage/fixed_net_templates/'

SIX_SPIKES = [
    ('22072442', '31370038', '31390354', '31397146', '31519150', '31689528', '31735952', '31766700', '31817932', '+31864090+'),
    ('22087964', '31370044', '31390198', '31401062', '31414919', '31519151', '31689529', '31730432', '31759138', '31793240', '+31832556+'),
    ('22087965', '31370045', '31392169', '31404140', '31419791', '31519152', '31689530', '31741849', '31791436', '31825166', '+31901939+'),
    ('22087966', '31370046', '31396924', '31412598', '31519153', '31689531', '31739668', '31763425', '+31818159+'),
    ('22087967', '31370047', '31397971', '31412573', '31427275', '31519154', '31689533', '31743014', '31763577', '31835183', '+31871621+'),
    ('22087969', '31370049', '31393197', '31404071', '31519155', '31689536', '31741660', '+31767463+'),
    ('22087970', '31370072', '31392820', '31401873', '31519157', '31689537', '31738490', '+31766626+'),
    ('22087971', '31370183', '31396623', '31412233', '31519159', '31689540', '31737800', '+31764868+'),
    ('22087972', '31370186', '31396423', '31412439', '31519161', '31689543', '31729998', '31761039', '+31805984+'),
    ('22087973', '31370242', '31389652', '31519162', '31689544', '31739760', '+31760348+'),]

THREE_SPIKES = [
    ('31817932', '33088890', '33130828', '33150821', '33169363', '33196902', '33230558', '+33248002+'),
    ('31793240', '33088894', '33138100', '33152203', '33201249', '33238084', '+33250972+'),
    ('31825166', '33088901', '33124461', '33149315', '33166758', '33183865', '+33212362+'),
    ('31763425', '33088907', '33131130', '33150972', '33177255', '33205609', '33234799', '33249051', '+33252123+' ),
    ('31835183', '33088915', '33127211', '33157592', '33172116', '33192104', '+33235021+'),
    ('31741660', '33088919', '33147305', '33164595', '33205556', '+33251787+'),
    ('31738490', '33088926', '33138263', '33155039', '33177259', '+33204717+' ),
    ('31737800', '33088930', '33139122', '33153551', '33169010', '33187885', '+33223605+' ),
    ('31761039', '33088933', '33124370', '33148095', '33163998', '33179470', '33217120', '+33239144+' ),
    ('31739760', '33088939', '33131194', '33151255', '+33167279+'),]

ALL_SPIKES = [
    ('33230558','33357781','33413275','33422096','33434464','33447914','33465825','+33477877+'),
    ('33238084','33357783','33393194','33415310','33426780','33437911','33447538','33458282','+33465056+'),
    ('33183865','33357784','33397857','33420873','33431566','33445268','33455574','+33462623+'),
    ('33249051','33357785','33399813','33425820','33436914','33445885','33455651','33462572','+33466633+'),
    ('33192104','33357786','33385823','33412221','33421554','33443595','33452792','+33462456+'),
    ('33205556','33357787','33383346','33407989','33421229','33436076','33444470','+33454014+'),
    ('33177259','33357788','33387409','33409757','33422380','33433787','33447922','+33462482+'),
    ('33187885','33357789','33378257','33400048','33414034','33422392','+33443173+'),
    ('33217120','33357790','33395217','33418472','33432306','33443443','33448552','33461063','33465508','+33485730+'),
    ('33151255','33357791','33380109','33414869','33431578','33441682','33448950','+33460146+'),
]


def make_psth_and_vm(ax_psth, ax_vm, ax_kc_hist):
    binwidth=100
    datalist = (SIX_SPIKES[0][1:-1], THREE_SPIKES[0][1:-1], ALL_SPIKES[0][1:-1])
    colors = ['#e66101', '#5e3c99', '#009292']
    ls = ['-', ':']
    for ii, group in enumerate(datalist):
        print(group)
        for jj, jid in enumerate((group[0], group[-1])):
            print(jid)
            try:
                fname = nda.find_h5_file(jid, DATA_DIR)
            except:
                # First entry is old data and moved to back up
                # Template should still have all the data
                print(f'JID {jid} not in datadir. Looking in template dir')
                fname = nda.find_h5_file(jid, TEMPLATE_DIR)
            with h5.File(fname, 'r') as fd:
                kc_st, kc_id = nda.get_event_times(fd[nda.kc_st_path])
                kc_sc = np.array([len(st) for st in kc_st])
                try:
                    ax_kc_hist.hist(kc_sc, bins=np.arange(1, max(kc_sc)+0.5, 1),
                                    color=colors[ii], ls=ls[jj],                       
                                    label=f'{ii}: {jj}: {jid}',
                                    histtype='step', linewidth=1)
                except IndexError:
                    print(jid, ':', kc_sc, '|')
                
                pop_st = np.concatenate(kc_st)
                try:                    
                    ax_psth.hist(pop_st,
                                 bins=np.arange(500, 2100, binwidth),
                                 color=colors[ii], ls=ls[jj], histtype='step', label=jid)
                except IndexError:
                    print(jid, pop_st)
                ggn_vm, t = nda.get_ggn_vm(fd, 'basal')
                ax_vm.plot(t, ggn_vm[0, :], label=jid, color=colors[ii], ls=ls[jj])
    ax_psth.legend()
    ax_vm.legend()


def plot_spike_counts(ax, fname=None):        
    if fname is None:
        for ii in range(len(SIX_SPIKES)):
            spike_counts = []
            for jj, sim_set in enumerate([SIX_SPIKES, THREE_SPIKES, ALL_SPIKES]):
                sim_list = sim_set[ii][1:-1]
                for kk, jid in enumerate(sim_list):            
                    try:
                        fname = nda.find_h5_file(jid, DATA_DIR)
                    except:
                        # First entry is old data and moved to back up
                        # Template should still have all the data
                        print(f'JID {jid} not in datadir. Looking in template dir')
                        fname = nda.find_h5_file(jid, TEMPLATE_DIR)
                    with h5.File(fname, 'r') as fd:
                        kc_st, kc_id = nda.get_event_times(fd[nda.kc_st_path])                    
                        kc_sc = sum([len(st) for st in kc_st])
                        print('JID:', jid, 'total spikes:', kc_sc)
                        # For all cases, plot the result from the last successful KC removal, and first
                        spike_counts.append(kc_sc)
                ax.plot(len(spike_counts)-1, spike_counts[-1], 'k|')
            ax.plot(spike_counts, 'o-', fillstyle='none')
    else:
        total_spike_count = pd.read_csv(fname, sep=',')
        for ii, (series_id, simgrp) in enumerate(total_spike_count.groupby('series')):
            ax.plot(simgrp['total_spikes'].values, 'o-', fillstyle='none')
            series_df = simgrp.reset_index()
            for removal, remgrp in series_df.groupby('removal'):
                print('#', removal)
                print(remgrp)
#                jj += len(remgrp)
                ax.plot(remgrp.index.values[-1], remgrp['total_spikes'].values[-1], 'k.')
               
def main():
    plt.close('all')
    fig = plt.figure()
    ax_kc_psth = fig.add_subplot(221)
    ax_total_spike_count = fig.add_subplot(222)
    ax_ggn_vm = fig.add_subplot(223, sharex=ax_kc_psth)
    ax_spike_count_hist = fig.add_subplot(224)
    make_psth_and_vm(ax_kc_psth, ax_ggn_vm, ax_spike_count_hist)
    if os.path.exists('total_spike_counts.csv'):
        plot_spike_counts(ax_total_spike_count, 'total_spike_counts.csv')
    ax_total_spike_count.set_xticks([0, 5, 10, 15, 20])
    ax_ggn_vm.set_xlim(500, 2200)
    ax_ggn_vm.set_ylim(-51, -40)
    for ax in fig.axes:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
    fig.set_size_inches(175/25.0, 140/25.0)
    fig.savefig('Figure_kc_removal.svg')
    plt.show()
#
#ax = plt.subplot(111)
#plot_spike_counts(ax, fname='D:/subhasis_ggn/model/analysis/total_spike_counts.csv')
#

if __name__ == '__main__':
    main()

# 
# fig_kc_removal.py ends here
