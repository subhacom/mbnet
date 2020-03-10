#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --partition=norm

# This script runs a simulatin with a template of a given JID
# And schedules removal of KCs from the resulting dataset via run_remove_kcs_dep_fixed_net.sh
# And scedules itself conditional on completion of the latter
# Usage: run_fixed_net_dep_remove_kcs.sh -j jid -l limit -t simtime

set -e
source $HOME/miniconda2/etc/profile.d/conda.sh
module load neuron/7.5
conda deactivate
conda activate rhel7
echo `which python`

DATA_DIR='/data/rays3/ggn/fixed_net'
TEMPLATE_DIR='/data/rays3/ggn/fixed_net_templates'
GGNDIR=$HOME/projects/ggn
SLURMDIR=$HOME/projects/ggn/mb/slurm/circular_run

OLD_PYTHONPATH=$PYTHONPATH

export PYTHONPATH=$PYTHONPATH:$GGNDIR/mb:$GGNDIR/nrn:$GGNDIR/common:$GGNDIR/morphutils:/data/rays3/python/envs/rhel7/lib/python2.7/site-packages/
echo "PYTHONPATH=$PYTHONPATH"
pushd $GGNDIR/mb
unset DISPLAY
echo "Starting fixed_net run. Arguments: $*"

while getopts "j:l:t:" opt; do
    case $opt in
	j)
	    echo "Template jid $OPTARG"
	    template_jid=$OPTARG
	    ;;
	l)
	    echo "Limit $OPTARG"
	    limit=$OPTARG
	    ;;
	t)
	    echo "Simulation time $OPTARG"
	    simtime=$OPTARG
	    ;;
	:)
	    echo "Option -$OPTARG requires positional arg"
	    exit 1
	    ;;
    esac
done

template_file=`find $TEMPLATE_DIR -name "*${template_jid}*kc${limit}.h5"`
echo "$SLURM_JOB_ID: running simulation with template=$template_file"
time python network/fixed_network_changing_stim.py -f $template_file -o $DATA_DIR -s 0 -d 0  --n_kc_vm=100 --n_ggn_vm=100 --simtime=$simtime --savesyn
popd
export PYTHONPATH=$OLD_PYTHONPATH

mydir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
myfname="$(basename $0)"
jobid=`sbatch --depend=afterany:$SLURM_JOB_ID $SLURMDIR/run_remove_kcs_dep_fixed_net.sh --jid=$SLURM_JOB_ID --limit=$limit --sdir=$DATA_DIR --tdir=$TEMPLATE_DIR`
nextrun=`sbatch --depend=afterany:$jobid $SLURMDIR/run_fixed_net_dep_remove_kcs.sh -j $SLURM_JOB_ID -l $limit -t $simtime`
echo "Finished $myfname running from $mydir"
