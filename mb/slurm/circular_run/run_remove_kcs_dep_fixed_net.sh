#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --partition=norm

# This script is for removing KCs from  dataset with specified JID and creating a template

set -e
source $HOME/miniconda2/etc/profile.d/conda.sh

# DATA_DIR='/data/rays3/ggn/fixed_net'
# TEMPLATE_DIR='/data/rays3/ggn/fixed_net_templates'
GGNDIR=$HOME/projects/ggn
SLURMDIR=$HOME/projects/ggn/mb/slurm/circular_run


# set PATH="$PATH":$HOME/miniconda2/condabin:$HOME/miniconda2/bin
conda deactivate
conda activate py3
echo `which python`
pushd $GGNDIR/mb
OLD_PYTHONPATH=$PYTHONPATH

export PYTHONPATH=$PYTHONPATH:$GGNDIR/analysis
echo "$SLURM_JOB_ID: removing KCs for args $*"
echo "PYTHONPATH=$PYTHONPATH"
time python $SLURMDIR/remove_high_firing_kcs.py $*
conda deactivate
popd
export PYTHONPATH=OLD_PYTHONPATH
