This directory contains scripts for simulating the GGN model and the
mushroom body olfactory circuit around it in NEURON 7.4 with Python
2.7 related to the article:

"Feedback inhibition and its control in an insect olfactory circuit".
Subhasis Ray, Zane Aldworth, and Mark Stopfer; 2020. eLife.


<a id="org2f03dda"></a>

# analysis

-   various example scripts for analyzing simulated data.


<a id="org66870b1"></a>

# common

-   nhpp.py : Generate times for nonhomogeneous Poisson process


<a id="org436e58b"></a>

# mb


<a id="org18dfcff"></a>

## cell<sub>templates</sub>

-   GGN<sub>20170309</sub><sub>sc.hoc</sub> : cell template for GGN model
-   kc<sub>1</sub><sub>comp.hoc</sub> : single compartmental KC model


<a id="org2d6beb8"></a>

## mod: contains mechanism files


<a id="org97a33c4"></a>

## network

-   change<sub>pn</sub><sub>spikes.py</sub> : take a simulated data file and create
    another one after modifying the PN spike times, so that same model
    is simulated with modified spike input.
-   config.yaml : configuration file for setting network model
    parameters
-   fixed<sub>network</sub><sub>changing</sub><sub>stim.py</sub> : simulate a given network model
    with different PN spike inputs.
-   kc<sub>ggn</sub><sub>feedback</sub><sub>dclamp.py</sub> : test of single KC with GGN inhibition
    when the GGN is driven by a dynamic clamp.
-   kc<sub>ggn</sub><sub>feedback</sub><sub>frequency</sub><sub>sweep.py</sub> : amplitude and frequency
    sweeps for testing effect of GGN feedback on a single KC.
-   kc<sub>ggn</sub><sub>nofeedback</sub> : script to simulate KC with no feedback
    alongside KC with feedback inhibition from GGN.
-   pn<sub>kc</sub><sub>ggn</sub><sub>network.py</sub> : script to setup and simulate the mushroom
    body network model (uses config.yaml for parameters).
-   pn<sub>output.py</sub> : script to setup PN spike trains
-   tweak<sub>template.py</sub> : script to modify an existing network template
    (in a data file dumped by an earlier simulation).


<a id="org87bf80e"></a>

## slurm : utility scripts for running simulations in batch mode under slurm (on NIH biowulf).

-   batch<sub>run</sub><sub>remove</sub><sub>kcs</sub><sub>run.sh</sub> : example script for running successive
    simulations after removing high spiking KCs.
-   circular<sub>run</sub> : scripts for running running successive simulations
    after removing high spiking KCs. These scripts read last job ids
    from a specified file to identify the corresponding data dumps,
    remove high spiking kcs from those model templates, and run the
    simulation again until no more highspiking KC is left.
-   run<sub>fixed</sub><sub>net</sub><sub>changing</sub><sub>stim.py</sub>: example utility script to run
    simulation of a given network model with changed PN inputs as a
    subprocess.
-   run<sub>fixed</sub><sub>net</sub><sub>with</sub><sub>ig.sh</sub> : sample script to run a given network
    template including IG with different PN inputs.
-   run<sub>fixed</sub><sub>network</sub><sub>changing</sub><sub>stim.sh</sub> : script to run a fixed network
    template while changing the PN input pattern.
-   run<sub>kc</sub><sub>ggn</sub><sub>feedback</sub><sub>amp</sub><sub>sweep.sh</sub>,
    run<sub>kc</sub><sub>ggn</sub><sub>nofeedback</sub><sub>amp</sub><sub>sweep.sh</sub> : scripts to test single KC with
    and without GGN feedback.
-   run<sub>mb</sub><sub>net.sh</sub> : script to run mushroom body network model in batch mode.
-   run<sub>to</sub><sub>reproduce</sub><sub>ggn</sub><sub>vm</sub><sub>no</sub><sub>seed.sh</sub>, run<sub>to</sub><sub>reproduce</sub><sub>ggn</sub><sub>vm.sh</sub>,
    run<sub>to</sub><sub>reproduce</sub><sub>ggn</sub><sub>vm</sub><sub>with</sub><sub>shifting</sub><sub>pn.sh</sub>: scripts that run
    multiple simulations while changing parameters to reproduce
    realistic GGN voltage trace.
-   run<sub>with</sub><sub>ig.sh</sub> : script to simulate network model with IG
    included.


<a id="org2e77656"></a>

## test<sub>cell</sub>

-   ggn<sub>voltage</sub><sub>attenuation</sub><sub>vclamp.py</sub> : check voltage attenuation
    along GGN arbor when one branch is voltage clamped.
-   run<sub>param</sub><sub>sweep</sub><sub>vclamp.py</sub> : test voltage attenuation along GGN
    arbor with voltage clamp while changing passive properties.


<a id="org0fc9385"></a>

# morphutils

-   cellmovie.py : dump a video of rotating neuron in 3D
-   displaycell.py : display neuronal morphology from SWC file
-   morph3d<sub>\*</sub>.py : display neuronal morphology in 3D using
    corresponding module.  Use vtk/vispy/matplotlib to display 3D
    morphology.
-   neurograph.py : handle morphology as a graph.  This file contains
    data type definitions and functions to read an SWC file and turn
    it into a graph (using networkx).  There is also GGN specific
    mapping where I assign custom types to specific branches based on
    anatomical location.


<a id="orge51345c"></a>

# nrn

-   nrnutils.py : Utilities for handling NEURON model
    -   convert a NEURON cell model into a networkx graph
    -   insert alpha synapses
    -   insert ion channel mechanisms
    -   set up recording of Vm
-   localized<sub>input</sub><sub>output.py</sub> : apply synaptic inputs at specified
    branches.  This scripts runs simulation with synchronous synaptic
    inputs at multiple compartments on specific branches.
-   localized<sub>input</sub><sub>output</sub><sub>passive</sub><sub>sweep.py</sub> : apply synaptic inputs at
    specified branches of GGN models with different passive
    properties.
-   nrninit.bat, nrninit.sh : batch file (Windows) and shell script
    (Linux).  Initialize PYTHONPATH to include various model related
    Python scripts
-   staggered<sub>input.py</sub> : simulate synaptic input at random
    compartments in specified branch
