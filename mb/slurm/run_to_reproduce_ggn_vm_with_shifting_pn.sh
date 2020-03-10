#!/bin/bash

echo "Starting to queue jobs"
for igg in $(seq 50.0 100.0 255.0); do
    echo "A. IGG: ${igg}"
    for gkg in $(seq 1 0.2 1.6); do
	echo "B. GKG: ${gkg}"
	for pkg in $(seq 2.5 0.2 3.1); do
	    echo "C. PKG: ${pkg}"
	    for kgg in $(seq 15.0 5.0 25.0); do
		echo "D. KG: ${kgg}"
		echo "ig_ggn_gmax=${igg}, ggn_kc_gmax=${gkg}, kc_ggn_gmax=${kgg} pn_kc_gmax=${pkg}"
		sbatch slurm/run_mb_net.sh --shifting_pn \
		       --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=${kgg}pS \
		       --ggn_kc_gmax=${gkg}nS --pn_kc_gmax=${pkg}pS \
		       --ig_ggn_gmax=${igg}nS \
		       --kc_frac_weak_inh=0.0
		sleep 1
	    done
	done
    done
done
