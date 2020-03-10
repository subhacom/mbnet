# run_fixed_net_changing_stim.py --- 
# Author: Subhasis  Ray
# Created: Thu Jan  3 18:17:36 2019 (-0500)
# Last-Updated: Thu Jan  3 18:43:20 2019 (-0500)
#           By: Subhasis  Ray
# Version: $Id$

# Code:
"""Simulations of the same template networks with different stimuli, multiple trials."""
import sys
import os
import subprocess

jobids = [
    # clustered PN-KC connection
    '17024513',
    '17024515',
    '17024517',
    '17024519',
    '17024521',
    '17024523',
    '17024526',
    '17024528',
    # diffuse PN-KC connection
    '17024514',
    '17024516',
    '17024518',
    '17024520',
    '17024522',
    '17024525',
    '17024527',
    '17024529']


datadir =  '/data/rays3/ggn/olfactory_network'
outdir = '/data/rays3/ggn/fixed_net/'
trial_count =  5

pn_dither = 5
pn_shifts =  [0,  15,  30,  100]

def find_h5_file(jid, datadir):
    """Copied from network_data_analysis.py"""
    flist = os.listdir(datadir)
    match = [fname for fname in flist if fname.endswith('JID{}.h5'.format(jid))]
    if len(match) > 1:
        raise Exception('Two HDF5 files with same JID: {}'.format(jid))
    elif len(match) == 0:
        raise Exception('No match for jid {} in {}'.format(jid, datadir))
    return os.path.join(datadir, match[0])


def run_trials(templatefile, pn_shift, trials = trial_count, recstep = 5, simtime = 1500):
    for ii in range(trials):
        output =  subprocess.check_output(['sbatch', 'slurm/run_fixed_network_changing_stim.sh', '-f', templatefile, '-o', outdir, '-d', '5', '-s', str(pn_shift), '--n_kc_vm=100', '--n_ggn_vm=100', '--simtime={}'.format(simtime)])
        print('Template {},  PN shift: {},  Run: {}'.format(templatefile,  pn_shift, ii))
        print(output)
        sys.stdout.flush()


if __name__ ==  '__main__':
    for jid in jobids:
        fpath =  find_h5_file(jid, datadir)
        for shift in pn_shifts:
            run_trials(fpath, shift)

# 
# run_fixed_net_changing_stim.py ends here
