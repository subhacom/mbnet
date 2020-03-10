# config.py --- 
# 
# Filename: config.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Jun 30 10:27:20 2016 (-0400)
# Version: 
# Package-Requires: ()
# Last-Updated: Fri Feb  9 11:49:50 2018 (-0500)
#           By: Subhasis Ray
#     Update #: 84
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

from __future__ import print_function
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

logfilename = 'log/mb_model_UTC{}-PID{}-JID{}.log'.format(
    timestamp.strftime('%Y_%m_%d__%H_%M_%S'), mypid, myjobid)
logging.basicConfig(filename=logfilename, level=logging.DEBUG,
                    format='%(asctime)s %(funcName)s: %(message)s',
                    filemode='w')

logger = logging.getLogger('mb_model')

if sys.platform == 'win32':
    os.environ['NEURONHOME'] = 'c:\\nrn'
    sys.path += ['c:\\nrn\\lib\\python', 'd:\\subhasis_ggn\\model\\common', 'd:\\subhasis_ggn\\model\\morphutils', 'd:\\subhasis_ggn\\model\\mb\\network']
    # The mod files are in the directory `mod` and nrnmech.dll is created there
    dll = os.path.join(os.path.dirname(__file__), 'mod', 'nrnmech.dll')
    print('nrnmech', dll)
    h.nrn_load_dll(dll)
else:
    # os.environ['NEURONHOME'] = '/usr/local/apps/neuron/nrn-7.4/'
    sys.path += ['/home/rays3/projects/ggn/mb', '/home/rays3/projects/ggn/common', '/home/rays3/projects/ggn/morphutils', '/home/rays3/projects/ggn/nrn', '/home/rays3/projects/ggn/mb/network']

logger.info('sys.path={}'.format(str(sys.path)))
h.load_file('stdrun.hoc')

from pint import UnitRegistry
ur = UnitRegistry()
Q_ = ur.Quantity


# 
# config.py ends here
