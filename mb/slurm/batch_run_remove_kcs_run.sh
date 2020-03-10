#!/bin/bash


## First round from original simulations
# jids=('22072442'
#       '22087964'
#       '22087965'
#       '22087966'
#       '22087967'
#       '22087969'
#       '22087970'
#       '22087971'
#       '22087972'
# for jid in "${jids[@]}"
# do
#     echo "Run starting with $jid"
#     sbatch  slurm/run_remove_kcs_run.sh --jid=$jid --template=/data/rays3/ggn/fixed_net_templates/ --data=/data/rays3/ggn/fixed_net/ --limit=5 --sleep=1800
# done
#       '22087973')

## second round - retry the failed ones
LIMIT=5
jids=('30977902'
      '30987321'
      '30978680'
      '30984569')

for jid in "${jids[@]}"
do
    echo "Run starting with $jid, LIMIT=$LIMIT"
    sbatch  slurm/run_remove_kcs_run.sh --jid=$jid --template=/data/rays3/ggn/fixed_net_templates/ --data=/data/rays3/ggn/fixed_net/ --limit=$LIMIT --sleep=1800
done

## Third round - for runs where no more > 5 spiking kcs were left, lower the spike limit
LIMIT=3
jids=('30983028'
      '31004736'
      '30985103'
      '30988926'
      '30977839')

for jid in "${jids[@]}"
do
    echo "Run starting with JID=$jid, LIMIT=$LIMIT"
    sbatch  slurm/run_remove_kcs_run.sh --jid=$jid --template=/data/rays3/ggn/fixed_net_templates/ --data=/data/rays3/ggn/fixed_net/ --limit=$LIMIT --sleep=1800
done
