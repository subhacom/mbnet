#!/bin/bash

# templates=(/data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42990-JID13081216_20190225_2.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_20190225_2.h5)

# templates=(/data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_20190226.h5)
# templates=(/data/rays3/ggn/olfactory_network/mb_net_UTC2019_02_27__22_54_44-PID52151-JID21337799_edited_1.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2019_02_27__22_54_44-PID52151-JID21337799_edited_2.h5)

# templates=(/data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_new_input.h5)
# templates=(/data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_new_input_20190304_1.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_new_input_20190304_2.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_new_input_20190304_3.h5)


# templates=(/data/rays3/ggn/olfactory_network/mb_net_UTC2019_02_28__16_26_41-PID18990-JID21392633_fixed.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2019_02_28__16_26_41-PID18990-JID21392633_fixed_A2.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2019_02_28__16_26_41-PID18990-JID21392633_fixed_B1.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2019_02_28__16_26_41-PID18990-JID21392633_fixed_B2.h5)
# templates=(/data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_new_input_20190304_1B.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_new_input_20190304_2B.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_new_input_20190304_3B.h5)
# templates=(/data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_4B.h5 /data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_4.h5)
templates=(/data/rays3/ggn/olfactory_network/mb_net_UTC2018_11_08__17_07_09-PID42992-JID13081218_4E.h5)
shifts=(0 50 100 200 300)
for template in "${templates[@]}"; do
    for s in "${shifts[@]}"; do
	echo "Template=${template}, Shift=${s}"
	sbatch slurm/run_fixed_network_changing_stim.sh -f ${template} -o /data/rays3/ggn/fixed_net/ -s ${s} --simtime=2500.0
	echo
	sleep 1
    done
done
