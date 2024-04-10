# test_kc.py --- 
# 
# Filename: test_kc.py
# Description: 
# Author: Subhasis Ray
# Created: Wed Apr 10 16:56:23 2024 (+0530)
# Last-Updated: Wed Apr 10 17:48:41 2024 (+0530)
#           By: Subhasis Ray
# 

# Code:

"""Create model from cell template file `filename`. The cell name in
    the template should be specified in `cellname`.

"""
from config import Q_, ur, h
import numpy as np
from matplotlib import pyplot as plt
import ephys

# filename and cellname specify cell template file and the name of the
# cell in the template respectively
filename='mb/cell_templates/kc_1_comp.hoc'
cellname='KC'


Em = Q_('-70mV')        # Wustenberg, et al. 2004, TABLE 2 legend

h.xopen(filename)
kc = eval(f'h.{cellname}()')
delay = Q_(100, 'ms')
duration=Q_(500.0, 'ms')
inj_current = Q_(16, 'pA')
clamp = ephys.setup_current_clamp(kc.soma, delay=delay, duration=duration, amplitude=inj_current)

# Recording data
vvec = ephys.setup_sec_rec(kc.soma, 'v')[0]
ik_vec = ephys.setup_sec_rec(kc.soma, 'ik')[0]
ina_vec = ephys.setup_sec_rec(kc.soma, 'ina')[0]
im_vec = h.Vector()
im_vec.record(clamp._ref_i)
tvec = h.Vector()
tvec.record(h._ref_t)

h.tstop = delay.to('ms').m + duration.to('ms').m
h.v_init = Em.to('mV').m
h.init()
h.run()
# Now plot the data
t = Q_(np.asarray(tvec.x), 'ms')
fig, axes = plt.subplots(nrows=2, sharex='all')
axes[0].plot(t, np.array(vvec.x))
axes[1].plot(t, np.array(im_vec.x))
axes[0].set_ylabel('Vm (mV)')
axes[1].set_ylabel('Im (nA)')

data = np.array([t.to('s').m, np.array(vvec.x)], dtype=[('t', float), ('v', float)])
data_file = 'kc_vm.npy'
np.save(data_file, data)

print(f'Simulated {cellname} from template file {filename}.')
print(f'Settling time {delay} followed by {inj_current} current injection for {duration}')
print(f'Saved Vm in {data_file}')
plt.show()


# 
# test_kc.py ends here
