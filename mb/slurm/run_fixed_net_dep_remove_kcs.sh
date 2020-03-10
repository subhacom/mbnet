#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --partition=norm
set -e
module load neuron/7.5
source deactivate rhel7
source activate rhel7
echo `which python`
GGNDIR=$HOME/projects/ggn
OLD_PYTHONPATH=$PYTHONPATH
export PYTHONPATH=$PYTHONPATH:$GGNDIR/mb:$GGNDIR/nrn:$GGNDIR/common:$GGNDIR/morphutils:/data/rays3/python/envs/rhel7/lib/python2.7/site-packages/
echo "PYTHONPATH=$PYTHONPATH"
pushd $GGNDIR/mb
unset DISPLAY
echo "Arguments: $*" 
time python network/fixed_network_changing_stim.py  $*
popd
export PYTHONPATH=$OLD_PYTHONPATH
sbatch --depend=afterany:$SLURM_JOB_ID $GGNDIR/mb/slurm/run_remove_kcs_dep_simulation.sh -j $SLURM_JOB_ID
