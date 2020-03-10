#!/bin/bash
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
set -e
# LIMIT=5
# JID_FILE=start_jid.txt

DATA_DIR='/data/rays3/ggn/fixed_net'
TEMPLATE_DIR='/data/rays3/ggn/fixed_net_templates'
GGNDIR=$HOME/projects/ggn
SLURMDIR=$HOME/projects/ggn/mb/slurm/circular_run

LIMIT=0
JID_FILE=start_jid_0_spikes.txt
SIMTIME=2500

for jid in `cat ${MYDIR}/${JID_FILE}`; do
    sbatch ${SLURMDIR}/run_remove_kcs_dep_fixed_net.sh --jid=${jid} --limit=${LIMIT} --sdir=${DATA_DIR} --tdir=${TEMPLATE_DIR}
done
