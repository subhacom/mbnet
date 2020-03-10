#!/usr/bin/env python

"""Script for running a series of simulations removing high spiking
KCs each time.

Requires Python3.5 or later"""

import os
import subprocess
import re
import time
import argparse
import logging

logging.basicConfig(format='%(asctime)s %(funcName)s: %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


TEMPLATE_DIR = '/data/rays3/ggn/fixed_net_templates'
DATA_DIR = '/data/rays3/ggn/fixed_net'

LIMIT = 5

def find_template(jid, template_dir, limit):
    pattern = re.compile(f'.+-JID{jid}_kc{limit}.h5')
    for fname in os.listdir(template_dir):
        if pattern.match(fname) is not None:
            return fname
    return None


def find_data(jid, data_dir, retries=3, sleep_time=600):
    pattern = re.compile(f'.+-JID{jid}.h5')
    trial = 0
    for fname in os.listdir(data_dir):
        if pattern.match(fname) is not None:
            return fname
        if (trial < retries) and (sleep_time > 0):
            time.sleep(sleep_time)
            trial += 1
            logger.info('Retry %s', trial)
            continue
    return None
    

def check_template(jid, template_dir, data_dir, limit):
    template = find_template(jid, template_dir, limit)
    if template is None:
        data_file = find_data(jid, data_dir)
        if data_file is None:
            logger.info(f'Could not find data file for jid {jid} in {data_dir}')
            return
        proc = subprocess.run(['python',
                               '/home/rays3/projects/ggn/analysis/remove_high_firing_kcs.py',
                               '--limit', str(limit),
                               '--sdir', data_dir,
                               '--tdir', template_dir,
                               '--jid', jid],
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        logger.info(proc.stdout)
        if proc.returncode != 0:
            logger.info(proc.stderr)
            return
        template = find_template(jid, template_dir, limit)        
    return template


def run_template(jid, template_dir, data_dir, limit, sleep_time):
    logger.info(f'===== Run template {jid} {template_dir} {data_dir} {limit}')
    template = check_template(jid, template_dir, data_dir, limit)
    if template is None:
        return
    template_path = os.path.join(template_dir, template)
    proc = subprocess.run(['sbatch',
                           '/home/rays3/projects/ggn/mb/slurm/run_fixed_network_changing_stim.sh',
                           '-f', template_path,
                           '-o', data_dir,
                           '-d', '0',
                           '--n_kc_vm', '100',
                           '--n_ggn_vm', '100',
                           '--simtime', '2500', '--savesyn'],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if proc.returncode != 0:
        logger.info(proc.stderr)
        return
    logger.info(f'Started run with template for {jid}')
    logger.info(proc.stdout)
    new_jid = proc.stdout.strip()
    while True:
        # logger.info('Sleeping for', sleep_time)
        time.sleep(sleep_time)
        proc = subprocess.run(['squeue', '-h', '-u', 'rays3', '-j', new_jid], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        # logger.info('squeue shows:')
        # logger.info(proc.stdout)
        # logger.info('-'*10)
        if len(proc.stdout.strip()) == 0:
            logger.info(f'job {new_jid} finished.')
            break
    logger.info('-'*80)
    run_template(new_jid, template_dir, data_dir, limit, sleep_time)
        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser('remove high spiking kcs, create template and run the template until there are no more highspiking kcs to remove')
    parser.add_argument('--data', type=str, help='data directory')
    parser.add_argument('--template', type=str, help='template directory, ')
    parser.add_argument('--jid', type=str, help='starting jid')
    parser.add_argument('--limit', type=int, help='kcs firing more than this many spikes will be removed')
    parser.add_argument('--sleep', type=int, help='time between checks for end of simulation')    
    args = parser.parse_args()
    run_template(args.jid, args.template, args.data, args.limit, sleep_time=args.sleep)
