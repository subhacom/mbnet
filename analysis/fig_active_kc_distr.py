# plot_active_kc_distr.py --- 
# Author: Subhasis  Ray
# Created: Mon Nov 12 15:32:18 2018 (-0500)
# Last-Updated: Fri Dec  7 13:10:59 2018 (-0500)
#           By: Subhasis  Ray
# Version: $Id$

# Code:

import sys
import os
import collections
import itertools as it
import numpy as np
import h5py as h5
import pandas as pd
import matplotlib.pyplot as plt
import yaml 
import network_data_analysis as nda

plt.rc('font', size=11)

os.chdir('/home/rays3/projects/ggn/analysis')

datadir = '/data/rays3/ggn/fixed_net'
trials_file = 'odor_trials_data.csv'


odor_trials_data = pd.read_csv(trials_file)

def find_jid_spiking_kcs(jid_list, datadir):
    ret = []
    for jid in jid_list:
        fname = nda.find_h5_file(str(jid), datadir)
        with h5.File(fname, 'r') as fd:
            ret.append(nda.get_spiking_kcs(fd))
    return ret

jid_spiking_kcs_fname = 'jid_spiking_kcs.yaml'
if not os.path.exists(jid_spiking_kcs_fname):
    spiking_kcs = find_jid_spiking_kcs(odor_trials_data.jid.values, datadir)
    jid_spiking_kcs = dict(zip(odor_trials_data.jid.values, spiking_kcs))
    with open('jid_spiking_kcs.yaml', 'w') as fd:
        yaml.dump(jid_spiking_kcs, fd)                
    print('Dumped spiking KCs for all trials in jid_spiking_kcs.yaml')
else:
    with open(jid_spiking_kcs_fname, 'r') as fd:
        jid_spiking_kcs = yaml.load(fd)
        
jid_kc_count = odor_trials_data.copy()
jid_kc_count['spiking_kcs'] = 0

common_kcs = odor_trials_data[['connection', 'template_jid', 'odor']].copy()
common_kcs.drop_duplicates(inplace=True)
common_kcs.loc[:, 'common_kcs'] = 0
common_kcs.loc[:, 'avg_kcs'] = 0

for odor, ogrp in odor_trials_data.groupby('odor'):
    for conn, cgrp in ogrp.groupby('connection'):
        print('Connection:', conn)
        for template, tgrp in cgrp.groupby('template_jid'):
            print('Template:', template)
            common = None
            counts = []
            for idx, jid in tgrp.jid.iteritems():
                counts.append(len(jid_spiking_kcs[jid]))
                jid_kc_count.loc[jid_kc_count.jid == jid, 'spiking_kcs']   \
                    = counts[-1]
                if common is None:
                    common = set(jid_spiking_kcs[jid])
                else:
                    common = common.intersection(jid_spiking_kcs[jid])
            pos = (common_kcs.odor == odor) & (common_kcs.connection == conn) & (common_kcs.template_jid == template)
            common_kcs.loc[pos, 'common_kcs'] = len(common)
            common_kcs.loc[pos, 'avg_kcs'] = np.mean(counts)


templates = odor_trials_data.template_jid.unique()
# for t in templates:
#     fname = nda.find_h5_file(str(t), datadir)
#     print(fname)
#     with h5.File(fname, 'r') as fd:
#         print(yaml.dump(nda.load_config(fd), default_flow_style=''))

tmp_color_map = {
    10829002: '#e66101',       # iid
    10829014: '#b2abd2',      # clus
    9932209: '#fdb863',        # iid
    9932198: '#5e3c99'       # clus
}

tmp_label_map = {
    10829002: 'Diffuse 2',       # iid
    10829014: 'Clustered 2',      # clus
    9932209: 'Diffuse 1',        # iid
    9932198: 'Clustered 1'       # clus
}


########## Plot PN activity for each odor

ax = None
ii = 0
for odor, ogrp in odor_trials_data.groupby('odor'):
    fig = plt.figure()
    ax = fig.add_subplot(111, sharex=ax)
    print(odor)
    for odata in ogrp.itertuples():
        fname = nda.find_h5_file(str(odata.jid), datadir)
        with h5.File(fname, 'r') as fd:
            pns = list(fd[nda.pn_st_path].keys())
            pns = sorted(pns, key=lambda x: int(x.rpartition('_')[-1]))
            stlist = [fd[nda.pn_st_path][pn].value for pn in pns]
            ylist = [np.ones_like(st) * ii for ii, st in enumerate(stlist)]
            ax.plot(np.concatenate(stlist), np.concatenate(ylist), ',')
            # for sp in ['left', 'right', 'top', 'bottom']:
            #     ax.spines[sp].set_visible(False)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.set_ylim(0, len(pns))
            ax.set_xlim(200, 1500)
        ii += 1
        break
    fig.set_size_inches(0.8, 1)
    fig.set_frameon(False)
    fig.suptitle(odor.upper())
    fig.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.8)
    fig.savefig('fixed_net_odor_{}.png'.format(odor), transparent=True)
# plt.show()    


########## Plot the distribution of spiking KC counts
## Figure 4d
odor_idx = 0
fig, ax = plt.subplots()
odor_ranges = {}
for odor, ogrp in jid_kc_count.groupby('odor'):
    jj = 0    
    for conn, cgrp in ogrp.groupby('connection'):
        print('Connection:', conn)
        for template, tgrp in cgrp.groupby('template_jid'):
            # print('Template:', template)
            counts = tgrp.spiking_kcs.values
            # print(len(counts))
            c = tmp_color_map[template]
            ax.boxplot(counts, positions=[odor_idx + jj], notch=False,
                       boxprops=dict(color=c, facecolor=c),
                       capprops=dict(color=c),
                       whiskerprops=dict(color=c),
                       flierprops=dict(marker='o', color=c,
                                       markerfacecolor=c, markeredgecolor=c,
                                       linestyle='none', markersize=3),
                       medianprops=dict(color='white'),
                       patch_artist=True)
            # print(odor_idx + jj)
            jj += 1
    odor_ranges[odor] = (odor_idx, odor_idx + 4)
    odor_idx += 5

conn_template = jid_kc_count[['connection', 'template_jid']].copy().drop_duplicates()
lines = []
for conn, tgrp in conn_template.groupby('connection'):
    for template in sorted(tgrp.template_jid.values):
        print(template)
        lines.append(plt.Line2D([],[], linewidth=2, color=tmp_color_map[template], label=tmp_label_map[template]))
ax.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), handles=lines, loc=3, ncol=2, mode='expand')        
ax.set_xlim(-1, odor_idx)
texty = 1900.0
texth = 100
for odor, (oleft, oright) in odor_ranges.items():
    patch = plt.Rectangle((oleft-0.5, texty - texth/2.0), oright - oleft + 0.5, texth, edgecolor='black', facecolor='white')
    ax.add_patch(patch)
    bb = ax.text(0.5 * (oleft + oright), texty, odor.upper(), ha="center", va="center")
ax.set_ylabel('# of activated KCs')    
ax.spines['bottom'].set_visible(False)    
ax.spines['top'].set_visible(False)    
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.xaxis.set_visible(False)
fig.set_frameon(False)
fig.set_size_inches(10/2.54, 7/2.54)
fig.tight_layout()
fig.savefig('active_kc_distr.svg')
# plt.show()

############# Plot the distribution of common KCs vs average number of spiking KCs
## Figure 4e
odor_marker = {'a': '*',
               'a1': 's',
               'a2': 'd',
               'b': 'o'}
fig, ax = plt.subplots()
for conn, cgrp in common_kcs.groupby('connection'):
    for template, tgrp in cgrp.groupby('template_jid'):
        for odor, ogrp in tgrp.groupby('odor'):
            assert(len(ogrp) == 1)
            print(ogrp)
            ax.plot(ogrp.avg_kcs.values, ogrp.common_kcs.values, color=tmp_color_map[template], marker=odor_marker[odor], mfc='none')
ax.spines['bottom'].set_visible(False)    
ax.spines['top'].set_visible(False)    
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.set_frameon(False)
lines = []
for odor, marker in odor_marker.items():
    lines.append(plt.Line2D([],[], linestyle='', marker=marker, mec='black', mfc='none', label=odor.upper()))
ax.legend(bbox_to_anchor=(0.6, 0.02, 1.0, 0.5), handles=lines, loc=4, ncol=2, mode='expand')        
# ax.legend(handles=lines, loc=4, ncol=2, mode='expand')        
ax.set_xlabel('Average # of KCs')
ax.set_ylabel('# of common KCs')
fig.set_size_inches(15/2.54, 7/2.54)
fig.tight_layout()
fig.savefig('active_kc_common_vs_average.svg')
# plt.show()

############ Plot distribution of common KCs between different odors
## Figure 4f

spiking_kcs_stats = pd.read_csv('shared_spiking_kc_stats.csv')
spiking_kcs_stats.loc[:, 'odor_distance'] = -1

def get_odist(oleft, oright):
    odist = -1
    if oleft == oright:
        odist = 0
    elif ((oleft == 'a') and (oright == 'a1')) or ((oleft == 'a1') and (oright == 'a2')):
        odist = 15
    elif (oleft == 'a') and (oright == 'a2'):
        odist = 30
    elif (oleft == 'a') and (oright == 'b'):
        odist = 100
    return odist

for (oleft, oright), ogrp in spiking_kcs_stats.groupby(['odor_left', 'odor_right']):
    spiking_kcs_stats.loc[(spiking_kcs_stats.odor_left == oleft) &
                          (spiking_kcs_stats.odor_right == oright),
                          'odor_distance'] = get_odist(oleft, oright)

good_stats = spiking_kcs_stats[spiking_kcs_stats.odor_distance >= 0].copy()


################## Final plot Figure 4 f
## Violin plot of distribution of shared spiking KCs by connection and template
fig, ax = plt.subplots()
pos = 0
grpw = 6
for conn, cgrp in good_stats.groupby('connection'):
    for template, tgrp in cgrp.groupby('template'):
        common_kcs = []
        print('====')
        for odist, distgrp in tgrp.groupby('odor_distance'):
            print(odist)
            common_kcs.append(distgrp.common.values)
        print('---')
        c = tmp_color_map[template]
        vp = ax.violinplot(common_kcs, np.arange(4) * grpw + pos, showmedians=True, points=20)
        plt.setp(vp['bodies'], color=c)
        plt.setp(vp['cmins'], color=c)  #tmp_color_map[template])
        plt.setp(vp['cmedians'], color=c)  #tmp_color_map[template])
        plt.setp(vp['cmaxes'], color=c)  #tmp_color_map[template])
        plt.setp(vp['cbars'], color=c)  #tmp_color_map[template])
        pos += 1
ax.spines['bottom'].set_visible(False)    
ax.spines['top'].set_visible(False)    
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.set_frameon(False)
ax.set_ylabel('# of common KCs')
odist = sorted(good_stats.odor_distance.unique())
texth = 30
texty = -10
for ii, od in enumerate(odist):
    patch = plt.Rectangle((-1 + ii * grpw, texty - texth/2.0), grpw-1, texth, edgecolor='black', facecolor='white')
    ax.add_patch(patch)
    bb = ax.text(ii * grpw + 2.0, texty, str(od), ha="center", va="center")
ax.set_xlabel('Shift in responsive PN population')    
# ax.xaxis.set_visible(False)
ax.set_xticks([])
fig.set_size_inches(10/2.54, 7/2.54)
fig.tight_layout()
fig.savefig('active_common_kc_distr.svg')
plt.show()
######## END Fig 4d

# ##################
# ## Recompute common KCs from scratch and plot for each type of connection
# fig, ax = plt.subplots()
# pos = 0
# for conn, cgrp in odor_trials_data.groupby('connection'):
#     shared_kc_count = collections.defaultdict(list)
#     for template, tgrp in cgrp.groupby('template_jid'):
#         print(template)
#         for (oleft, oright) in it.combinations_with_replacement(tgrp.odor.unique(), 2):
#             odist = get_odist(oleft, oright)
#             print(oleft, oright, odist)
#             if odist < 0:
#                 continue
#             left = np.random.permutation(tgrp.loc[tgrp.odor == oleft, 'jid'].values)
#             right = np.random.permutation(tgrp.loc[tgrp.odor == oright, 'jid'].values)
#             if oleft == oright:            
#                 right = left[len(left) // 2:]
#                 left = left[:len(left) // 2]
#             for jid_left, jid_right in zip(left, right):
#                 left_kcs = set(jid_spiking_kcs[jid_left])
#                 right_kcs = set(jid_spiking_kcs[jid_right])
#                 shared_kc_count[odist].append(len(left_kcs.intersection(right_kcs)))        
            
#     sorted_dist_counts = sorted(list(shared_kc_count.items()), key=lambda x: x[0])
#     print('==========')
#     for k, v in sorted_dist_counts:
#         print(k, len(v))
#     print('------')
#     sorted_counts = [x[1] for x in sorted_dist_counts]
#     vp = ax.violinplot(sorted_counts, np.arange(4)* 4 + pos, showmedians=True, points=20)
#     pos += 1
# fig.suptitle('Confirmation with independent computation of common KC count')
# plt.show()

# ##################
# ## Recompute common KCs from scratch and plot for each type of connection
# ## and for each template
# fig, ax = plt.subplots()
# pos = 0
# for conn, cgrp in odor_trials_data.groupby('connection'):
#     for template, tgrp in cgrp.groupby('template_jid'):
#         print(template)
#         shared_kc_count = collections.defaultdict(list)
#         for (oleft, oright) in it.combinations_with_replacement(tgrp.odor.unique(), 2):
#             odist = get_odist(oleft, oright)
#             print(oleft, oright, odist)
#             if odist < 0:
#                 continue
#             left = np.random.permutation(tgrp.loc[tgrp.odor == oleft, 'jid'].values)
#             right = np.random.permutation(tgrp.loc[tgrp.odor == oright, 'jid'].values)
#             if oleft == oright:            
#                 right = left[len(left) // 2:]
#                 left = left[:len(left) // 2]
#             for jid_left, jid_right in zip(left, right):
#                 left_kcs = set(jid_spiking_kcs[jid_left])
#                 right_kcs = set(jid_spiking_kcs[jid_right])
#                 shared_kc_count[odist].append(len(left_kcs.intersection(right_kcs)))        
            
#         sorted_dist_counts = sorted(list(shared_kc_count.items()), key=lambda x: x[0])
#         print('==========')
#         for k, v in sorted_dist_counts:
#             print(k, len(v))
#         print('------')
#         sorted_counts = [x[1] for x in sorted_dist_counts]
#         vp = ax.violinplot(sorted_counts, np.arange(4)* 4 + pos, showmedians=True, points=20)
#         pos += 1
# plt.show()


# 
# plot_active_kc_distr.py ends here
