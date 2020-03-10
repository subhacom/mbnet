# run_param_sweep_staggered_input.py --- 
# Author: Subhasis Ray
# Created: Mon Sep 10 11:50:27 2018 (-0400)
# Last-Updated: Mon Sep 10 16:28:32 2018 (-0400)
#           By: Subhasis Ray
# Version: $Id$

# Code:
"""Run a parameter sweep for RM and RA values on GGN. RA - specific axial resistance (ohm-cm), RM - specific membrane resistance (kohm-cm2)"""
import os
import argparse
import numpy as np
import subprocess
from datetime import datetime


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest='filename', help='celltemplate filename')
    parser.add_argument('-v', '--vclamp', default=-40.0, help='clamp voltage')
    parser.add_argument('--RAmin', dest='RAmin', type=float, help='minimum value of RA')
    parser.add_argument('--RAmax', dest='RAmax', type=float, help='maximum value of RA')
    parser.add_argument('--RAdivs', dest='RAdivs', type=float, help='divisions in RA')
    parser.add_argument('--RMmin', dest='RMmin', type=float, help='minimum value of RM')
    parser.add_argument('--RMmax', dest='RMmax', type=float, help='maximum value of RM')
    parser.add_argument('--RMdivs', dest='RMdivs', type=float, help='divisions in RM')
    parser.add_argument('--simtime', type=float, default=300.0)
    # parser.add_argument('--rec', type=int, default=0)
    parser.add_argument('--maxrate', type=float)
    parser.add_argument('-n', dest='syncount', type=int, default=10)
    # parser.add_argument('-s', type=float, default=50)
    parser.add_argument('-i', dest='input')
    parser.add_argument('--outdir', default='data')
    return parser

if __name__ == '__main__':
    args = make_parser().parse_args()
    cellname = None
    ts = datetime.now()
    datadir = os.path.join(args.outdir, ts.strftime('%Y_%m_%d__%H_%M_%S'))
    with open(args.filename, 'r') as cellfile:
        for line in cellfile.readlines():
            if 'begintemplate' in line:
                cellname = line.partition('begintemplate')[-1].split()[0]
                break
    if cellname is None:
        raise Exception('Could not find cellname in template file')
    RA = np.linspace(args.RAmin, args.RAmax, args.RAdivs)
    RM = np.linspace(args.RMmin, args.RMmax, args.RMdivs)
    for jj, _rm in enumerate(RM):
        for ii, _ra in enumerate(RA):
            outfilename = 'vm_vclamp_RA{}_RM{}.h5'.format(ii, jj)
            subprocess.Popen(['python', 'ggn_voltage_attenuation_vclamp.py',
                              '-f', args.filename,
                              '-c', cellname,
                              '-v', args.vclamp,
                              '-s', 'dend_8[1]',
                              '--RA', str(_ra),
                              '--RM', str(_rm*1e3),
                              '--datadir', datadir,
                              '--deflection',
                              '--simtime', str(args.simtime),
                              '--outfile', outfilename])
                             
    


# 
# run_param_sweep_staggered_input.py ends here
