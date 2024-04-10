# config.py --- 
# 
# Filename: config.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Jun 30 10:27:20 2016 (-0400)
# Version: 
# Package-Requires: ()
# Last-Updated: Wed Apr 10 16:35:58 2024 (+0530)
#           By: Subhasis Ray
#     Update #: 103
# URL: 
# Doc URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# 
# 
# 

# Change Log:
# 
# 
# 
# 

# Code:


import sys
import os
from datetime import datetime
import logging
import subprocess

from neuron import h

nrn_version = subprocess.check_output(['nrniv', '--version']).strip()
timestamp = datetime.utcnow()
mypid = os.getpid()
myjobid = os.environ.get('SLURM_JOBID', '0')

# Create logging directory
if not os.path.exists('log'):
    os.mkdir('log')
    print(f'Created log directory "{os.getcwd()}/log"')

logfilename = 'log/mb_model_UTC{}-PID{}-JID{}.log'.format(
    timestamp.strftime('%Y_%m_%d__%H_%M_%S'), mypid, myjobid)
logging.basicConfig(filename=logfilename, level=logging.DEBUG,
                    format='%(asctime)s %(funcName)s: %(message)s',
                    filemode='w')

logger = logging.getLogger('mb_model')

try:
    nrnhome = os.environ['NEURONHOME']
    sys.path.append(os.path.join(nrnhome, 'lib', 'python'))
except KeyError as err:
    raise KeyError('Set envrionment variable NEURONHOME to the path of NEURON installation') from err

if sys.platform == 'win32':
    # The mod files are in the directory `mod` and nrnmech.dll is created there
    dll = os.path.join(os.path.dirname(__file__), 'mod', 'nrnmech.dll')
    print('nrnmech', dll)
    h.nrn_load_dll(dll)

logger.info('sys.path={}'.format(str(sys.path)))
h.load_file('stdrun.hoc')

from pint import UnitRegistry
ur = UnitRegistry()
Q_ = ur.Quantity


# 
# config.py ends here
