#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --partition=norm

module load neuron/7.5
source deactivate rhel7
source activate rhel7
echo `which python`
GGNDIR=$HOME/projects/ggn
export PYTHONPATH=$PYTHONPATH:$GGNDIR/mb:$GGNDIR/nrn:$GGNDIR/common:$GGNDIR/morphutils:/data/rays3/python/envs/rhel7/lib/python2.7/site-packages/
echo "PYTHONPATH=$PYTHONPATH"
pushd $GGNDIR/mb
unset DISPLAY
echo "Arguments: $*" 
time python network/pn_kc_ggn_network.py  $*
popd

