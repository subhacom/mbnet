#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --partition=norm

# Arguments:
# '-f', '--template_filename', type=str,
#                         help='template filename for pn spike trains and network config'
# '-o', '--output_directory', type=str,
#                     help='directory to save simulation data in'
# '-s', '--pn_shift', type=int,
#                     help='shift the PN spikes over the PN population by this many PNs'
# '-d', '--pn_dither', type=float, default=5.0,
#                     help='dither the PN spikes by maximum of this much time (in ms)'
# '--n_kc_vm', type=int, default=500,
#                     help='number of KCs to record Vm from'
# '--n_ggn_vm', type=int, default=100,
#                     help='number of section of GGN in each region to record from'
# '--recstep', default=5, type=int,
#                     help='number of integration steps per data recording step'
# '--simtime', type=float, help='Simulation time (ms)'

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

