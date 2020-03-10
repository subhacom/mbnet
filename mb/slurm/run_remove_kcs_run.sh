#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --partition=norm

set -e

source $HOME/miniconda2/etc/profile.d/conda.sh
GGNDIR=$HOME/projects/ggn
# set PATH="$PATH":$HOME/miniconda2/condabin:$HOME/miniconda2/bin
conda deactivate
conda activate py3
echo `which python`
pushd $GGNDIR/mb
python slurm/run_remove_kcs_run.py $*
conda deactivate
popd
