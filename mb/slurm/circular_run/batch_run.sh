#!/bin/bash
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
set -e
# LIMIT=5
# JID_FILE=start_jid.txt

LIMIT=0
JID_FILE=start_jid_0_spikes.txt
SIMTIME=2500

for jid in `cat ${MYDIR}/${JID_FILE}`; do    
    sbatch ${MYDIR}/run_fixed_net_dep_remove_kcs.sh -j ${jid} -l ${LIMIT} -t ${SIMTIME}
done
