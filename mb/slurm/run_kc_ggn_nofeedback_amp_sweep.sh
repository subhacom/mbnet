#!/bin/bash                                                                                 
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --partition=norm

module load neuron/7.5
conda deactivate rhel7
conda activate rhel7
echo `which python`
GGNDIR=$HOME/projects/ggn
export PYTHONPATH=$PYTHONPATH:$GGNDIR/mb:$GGNDIR/nrn:$GGNDIR/common:$GGNDIR/morphutils:/data/rays3/python/envs/rhel7/lib/python2.7/site-packages/
echo "PYTHONPATH=$PYTHONPATH"
pushd $GGNDIR/mb
unset DISPLAY
echo "Arguments: $*" 
python network/kc_ggn_nofeedback.py --kc-file=cell_templates/kc_1_comp.hoc --kc=KC --ggn-file=cell_templates/GGN_20170309_sc.hoc --ggn=GGN_20170309_sc --amp 0.01 0.035 0.001 --out /data/rays3/ggn/kc_ggn_feedback/kc_ggn_amp_sweep_nofeedback
popd
