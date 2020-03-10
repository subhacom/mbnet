#!/bin/bash                                                                                 
module load neuron/7.4                                                                      
module load python                                                                          
source activate subha_py2.7                                                                 
GGNDIR=$HOME/projects/ggn                                                                   
export PYTHONPATH=$PYTHONPATH:$GGNDIR/mb:$GGNDIR/nrn:$GGNDIR/common:$GGNDIR/morphutils      
echo "PYTHONPATH=$PYTHONPATH"                                                               
pushd $GGNDIR/mb                                                                            
unset DISPLAY                                                                               
python network/kc_ggn_feedback_frequency_sweep.py --kc-file=cell_templates/kc_1_comp.hoc --kc=KC --ggn-file=cell_templates/GGN_20170309_sc.hoc --ggn=GGN_20170309_sc --amp 0.01 0.035 0.001 --out /data/rays3/ggn/kc_ggn_feedback/kc_ggn_amp_sweep --jitter 0.0 1 --nkc 1000  # nkc=1 avoid scaling gsyn by lumping kcs of synapses (each synapse has standard conductance)
popd

