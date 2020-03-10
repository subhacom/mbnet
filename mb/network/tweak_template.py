# tweak_template.py --- 
# Author: Subhasis Ray
# Created: Fri Feb 22 16:11:33 2019 (-0500)
# Last-Updated: Tue Mar  5 15:53:53 2019 (-0500)
#           By: Subhasis  Ray
# Version: $Id$

# Code:
"""update synaptic conductances between cell types"""

import sys
import argparse
import shutil
import h5py as h5
import numpy as np
import yaml
from pint import UnitRegistry
from matplotlib import pyplot as plt
ur = UnitRegistry()
Q_ = ur.Quantity


VLENSTR = h5.special_dtype(vlen=bytes)

def make_parser():
    parser = argparse.ArgumentParser('This program updates the specified synaptic conductances in the input file.')
    parser.add_argument('-i', '--input_filename', type=str,
                        help='path to file to be updated.')
    parser.add_argument('-o', '--output_filename', type=str,
                        help='output file name')
    parser.add_argument('--pn_kc_gbar_mean', default='-1.0pS',
                        help='Update PN->KC conductance with this mean')
    parser.add_argument('--pn_kc_gbar_std', type=float, default=-1,
                        help='Update PN->KC conductance with lognorm distr with std.')
    parser.add_argument('--kc_ig_gbar_mean', type=float, default=-1,
                        help='Update KC->IG conductance with this mean')
    parser.add_argument('--kc_ig_gbar_std', type=float, default=-1,
                        help='Update KC->IG weight with lognorm distr with std.')
    parser.add_argument('--kc_ggn_gbar_mean', default='-1.0pS',
                        help='Update KC->GGN conductance with lognorm distr with this mean.')
    parser.add_argument('--kc_ggn_gbar_std', type=float, default=-1,
                        help='Update KC->GGN conductance with lognorm distr with this std.')
    parser.add_argument('--ggn_kc_gbar_mean', default='-1.0pS',
                        help='Update KC<-GGN conductance with lognorm distr with this mean.')
    parser.add_argument('--ggn_kc_gbar_std', type=float, default=-1,
                        help='Update KC<-GGN conductance with lognorm distr with this std.')
    parser.add_argument('--ig_ggn_gbar', default='-1.0pS',
                        help='Update IG->GGN conductance with this value.')
    return parser


def show_distribution(fname, pn_kc=True, kc_ig=True):
    with h5.File(fname, 'r') as fd:
        if pn_kc:
            pn_kc_syn = fd['/data/static/pn_kc_synapse/pn_kc_synapse']
            pn_kc_gmax = pn_kc_syn.value[:, 0]['gmax']
            fig, ax = plt.subplots()
            fig.suptitle('PN->KC gmax distribution')
            ax.hist(pn_kc_gmax)
            ax.set_xlabel('gmax {}'.format(pn_kc_syn.attrs['unit'][4]))
        if kc_ig:
            try:
                kc_ig_syn = fd['/data/static/kc_ig_synapse/kc_ig_synapse']
                kc_ig_gmax = kc_ig_syn.value[:, 0]['gmax']
                fig, ax = plt.subplots()
                fig.suptitle('KC->IG gmax distribution')
                ax.hist(kc_ig_gmax)
                ax.set_xlabel('gmax')
            except KeyError:
                print('NO data for KC->IG synapse')
    plt.show()
    
            
def main():
    parser = make_parser()
    args = parser.parse_args()
    shutil.copyfile(args.input_filename, args.output_filename)
    with h5.File(args.output_filename, 'r+') as fd:
        config = yaml.load(fd.attrs['config'].decode())
        # PN->KC
        std = args.pn_kc_gbar_std
        gmax = Q_(args.pn_kc_gbar_mean).to('uS').m
        if (std >= 0) and (gmax >= 0):
            pn_kc_syn = fd['/data/static/pn_kc_synapse/pn_kc_synapse']
            values = pn_kc_syn.value
            nsyn = values.shape[0]
            sigma2 = np.log(1 + std**2)  # std is specified as a fraction of mean
            mu = np.log(gmax) - sigma2 / 2.0
            _gmax = np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2),
                                        size=nsyn)
            values[:, 0]['gmax'] = _gmax
            pn_kc_syn[...] = values
            config['pn_kc_syn']['gmax'] = args.pn_kc_gbar_mean
            config['pn_kc_syn']['std'] = args.pn_kc_gbar_std
        # KC->IG
        std = args.kc_ig_gbar_std
        gmax = args.kc_ig_gbar_mean
        if (std >= 0) and (gmax >= 0):
            try:
                kc_ig_syn = fd['/data/static/kc_ig_synapse/kc_ig_synapse']
                values = kc_ig_syn.value
            except KeyError:
                kc_ig_syn = None
                dtype = np.dtype([('pre', VLENSTR),
                                  ('post', VLENSTR),
                                  ('prepos', np.float64),
                                  ('postpos', np.float64),
                                  ('gmax', np.float64),
                                  ('delay', np.float64)])
                conn_unit = ['', '', '', '', '', 'ms']
                kcs = list(fd['/model/modeltree/olf/kc'].keys())
                values = np.empty((len(kcs), 1), dtype=dtype)
                values['pre'][:, 0] = kcs
                values['post'][:, 0] = 'IG'
                values['delay'][:, 0] = Q_(config['kc_ig_syn']['delay']).to('ms').m
            sigma2 = np.log(1 + std**2)  # std is specified as a fraction of mean
            mu = np.log(gmax) - sigma2 / 2.0
            _gmax = np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2),
                                        size=len(kcs))
            values['gmax'][:, 0] = _gmax
            if '/data/static/kc_ig_synapse/kc_ig_synapse' in fd:
                fd['/data/static/kc_ig_synapse/kc_ig_synapse'][...] = values
            else:
                grp = fd.create_group('/data/static/kc_ig_synapse')
                ds = grp.create_dataset('kc_ig_synapse', data=values)
                ds.attrs['unit'] = np.array(conn_unit, dtype='S')
            config['kc_ig_syn']['weight'] = args.kc_ig_gbar_mean
            config['kc_ig_syn']['std'] = args.kc_ig_gbar_std
        # KC->GGN
        std = args.kc_ggn_gbar_std
        gmax = Q_(args.kc_ggn_gbar_mean).to('uS').m
        if (std >= 0) and (gmax >= 0):
            try:
                kc_ggn_syn = fd['/data/static/kc_ggn_alphaL_synapse/kc_ggn_alphaL_synapse']
                values = kc_ggn_syn.value
                sigma2 = np.log(1 + std**2)  # std is specified as a fraction of mean
                mu = np.log(gmax) - sigma2 / 2.0
                _gmax = np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2),
                                        size=len(kcs))
                values['gmax'][:, 0] = _gmax
                fd['/data/static/kc_ggn_alphaL_synapse/kc_ggn_alphaL_synapse'][...] = values
                config['kc_ggn_alphaL_syn']['gmax'] = args.kc_ggn_gbar_mean
                config['kc_ggn_alphaL_syn']['std'] = args.kc_ggn_gbar_std
            except KeyError:
                print('No data for KC->GGN alphaL synapses')
        # GGN->KC
        std = args.ggn_kc_gbar_std
        gmax = Q_(args.ggn_kc_gbar_mean).to('uS').m
        if (std >= 0) and (gmax >= 0):
            try:
                ggn_kc_syn = fd['/data/static/ggn_kc_synapse/ggn_kc_synapse']
                values = ggn_kc_syn.value
                sigma2 = np.log(1 + std**2)  # std is specified as a fraction of mean
                mu = np.log(gmax) - sigma2 / 2.0
                _gmax = np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2),
                                        size=len(kcs))
                values['gmax'][:, 0] = _gmax
                fd['/data/static/ggn_kc_synapse/ggn_kc_synapse'][...] = values
                config['ggn_kc_syn']['gmax'] = gmax
                config['ggn_kc_syn']['std'] = std
            except KeyError:
                print('No data for KC<-GGN synapses')
        if Q_(args.ig_ggn_gbar).to('uS').m > 0:
            config['ig_ggn_syn']['gmax'] = args.ig_ggn_gbar
        
                
    show_distribution(args.output_filename)

            
if __name__ == '__main__':
    # test(infile='/data/rays3/ggn/fixed_net/mb_net_UTC2018_09_25__22_25_34-PID16017-JID10014019.h5')
    main()


# 
# tweak_template.py ends here
