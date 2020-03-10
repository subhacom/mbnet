#!/bin/bash
for i in "${!pi1[@]}"; do
    echo "Repeat $i, NP seed ${pi1[$i]}, KMeans seed: ${pi2[$i]}"
    echo "1. Shifting PNs (in temporal clusters)"
    echo "1.2. No KC->GGN synapse in CA"
    echo "1.2._.1. KC clusters share more PNs"
    echo "1.2._.1.1. KC clusters share more PNs from same PN cluster"
    echo "KCs in same cluster share more PNs, from simultaneously active PN groups"
    sbatch slurm/run_mb_net.sh --shifting_pn  --pn_kc_clustered --pn_kc_clustered_pre \
	   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=10.0pS \
	   --ggn_kc_gmax=1.0nS --pn_kc_gmax=3.0pS \
	   --kc_frac_weak_inh=0.0
    sleep 3

    # sbatch slurm/run_mb_net.sh --shifting_pn  --pn_kc_clustered --pn_kc_clustered_pre \
    # 	   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=10.0pS \
    # 	   --ggn_kc_gmax=0.7nS --pn_kc_gmax=3.0pS \
    # 	   --kc_frac_weak_inh=0.0  \
    # 	   --npseed=${pi1[$i]} --kmseed=${pi2[$i]} --fake
    # sleep 3
    echo "1.2._.1.2. KC clusters share more PNs but random ones"
    echo "KCs in same cluster share more PNs, but random set of PNs"
    sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered  \
	   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=10.0pS \
	   --ggn_kc_gmax=1.0nS --pn_kc_gmax=3pS \
	   --kc_frac_weak_inh=0.0 \
	   --npseed=${pi1[$i]} --kmseed=${pi2[$i]}
    sleep 3
    # sbatch slurm/run_mb_net.sh --shifting_pn --pn_kc_clustered  \
    # 	   --kc_ggn_ca_gmax=0.0pS  --kc_ggn_alpha_gmax=10.0pS \
    # 	   --ggn_kc_gmax=0.7nS --pn_kc_gmax=3.0pS \
    # 	   --kc_frac_weak_inh=0.0  \
    # 	   --npseed=${pi1[$i]} --kmseed=${pi2[$i]} --fake
    # sleep 3
    echo "1.2._.2. KCs in same cluster do not share more PNs than with other clusters"
    echo "KCs in same cluster do not share more PNs"
    sbatch slurm/run_mb_net.sh --shifting_pn \
	   --kc_ggn_ca_gmax=0.0pS --kc_ggn_alpha_gmax=10.0pS \
	   --ggn_kc_gmax=1.0nS --pn_kc_gmax=3.0pS \
	   --kc_frac_weak_inh=0.0 \
	   --npseed=${pi1[$i]} --kmseed=${pi2[$i]}
    sleep 3
    echo "2. No structure in PN activity (not shifting)"
    echo "2.2. No KC->GGN connection in CA"
    echo "2.2._.1. KCs in same cluster share more PNs"
    sbatch slurm/run_mb_net.sh --pn_kc_clustered  \
	   --kc_ggn_ca_gmax=0.0pS   --kc_ggn_alpha_gmax=10.0pS \
	   --ggn_kc_gmax=0.9nS --pn_kc_gmax=3.5pS  \
	   --kc_frac_weak_inh=0.0 \
	   --npseed=${pi1[$i]} --kmseed=${pi2[$i]}
    sleep 3
    # sbatch slurm/run_mb_net.sh --pn_kc_clustered   \
    # 	   --kc_ggn_ca_gmax=0.0pS   --kc_ggn_alpha_gmax=10.0pS \
    # 	   --ggn_kc_gmax=0.9nS --pn_kc_gmax=3.7pS \
    # 	   --kc_frac_weak_inh=0.0 \
    # 	   --npseed=${pi1[$i]} --kmseed=${pi2[$i]} --fake
    # sleep 3
    echo "2.2._.2. KCs in same cluster do not share more PNs"
    sbatch slurm/run_mb_net.sh  --kc_ggn_ca_gmax=0.0pS  --kc_ggn_alpha_gmax=10.0pS \
	   --ggn_kc_gmax=0.9nS --pn_kc_gmax=3.5pS \
	   --kc_frac_weak_inh=0.0 \
	   --npseed=${pi1[$i]} --kmseed=${pi2[$i]}
    sleep 3
done
