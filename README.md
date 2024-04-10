This directory contains scripts for simulating the GGN model and the
mushroom body olfactory circuit around it (originally in NEURON 7.4 with Python
2.7) related to the article:

"Feedback inhibition and its control in an insect olfactory circuit".
Subhasis Ray, Zane Aldworth, and Mark Stopfer; 2020. eLife.

- April 2024: Updated to run with NEURON 8.2.4 on Python 3.12 during COMBINE/Harmony-2024 Hackathon.

# Running the simulations
  Required packages: 
  
  * networkx (graph data structures)
  * neuron (simulation)
  * sklearn (spatial clustering based connections)
  * pint (units)
  * pyyaml (configuration)
  * nsdf (saving data): `pip install git+https://github.com/nsdf/nsdf.git`
  * vtk (for 3D visualization)
  * numpy (numerics)
  * matplotlib (visualization)
  
  For testing, run simulation scripts from the root directory of this
  repo (where this `README.md` file resides). Edit the `nrn/nrninit.sh`
  (or `nrninit.bat` on Windows) with the path for this directory.
  
  If you are running from another directory, edit `nrn/nrninit.sh` and
  set `MODEL_DIR` to the path of that directory.
  * To initialize the environment variables run the `nrn/nrninit.sh`
    file. For Linux/OSX with bash shell: `source nrn/nrninit.sh`
  * You also have to build the NEURON mechanisms (.mod) files:

    `nrnivmodl mb/mod`
    
  * To simulate GGN in isolation to study signal attenuation with
    distance, run:
  
      `python mb/test_cell/ggn_voltage_attenuation_vclamp.py`
      
  * To simulate the PN->KC<->GGN network model, create a directory
    called `data` for dumping simulated data and run:

    `python mb/network/pn_kc_ggn_network.py`
    
    You may want to edit the configuration file
    `mb/network/config.yaml` to reduce the number of cells in the
    network and the duration of the simulation.


# analysis

-   various example scripts for analyzing simulated data.



# common

-   `nhpp.py` : Generate times for nonhomogeneous Poisson process



# mb

## `cell_templates`
-   `GGN_20170309_sc.swc` : morphology trace of GGN shrinkage corrected for
	methylsalicylate fixing
-   `GGN_20170309_sc.hoc` : cell template for GGN model
-   `kc_1_comp.hoc` : single compartmental KC model



## mod: contains mechanism files

## network

-   `change_pn_spikes.py` : take a simulated data file and create
    another one after modifying the PN spike times, so that same model
    is simulated with modified spike input.
-   `config.yaml` : configuration file for setting network model
    parameters
-   `fixed_network_changing_stim.py` : simulate a given network model
    with different PN spike inputs.
-   `kc_ggn_feedback_dclamp.py` : test of single KC with GGN
    inhibition when the GGN is driven by a dynamic clamp.
-   `kc_ggn_feedback_frequency_sweep.py` : amplitude and frequency
    sweeps for testing effect of GGN feedback on a single KC.
-   `kc_ggn_nofeedback.py` : script to simulate KC with no feedback
    alongside KC with feedback inhibition from GGN.
-   `pn_kc_ggn_network.py` : script to setup and simulate the mushroom
    body network model (uses config.yaml for parameters).
-   `pn_output.py` : script to setup PN spike trains
-   `tweak_template.py` : script to modify an existing network
    template (in a data file dumped by an earlier simulation).



## slurm : utility scripts for running simulations in batch mode under slurm (on NIH biowulf).
    *NOTE: Many of these scripts have hard coded absolute paths which
    need to be updated for running in a new environment.*
    
-   `batch_run_remove_kcs_run.sh` : example script for running successive
    simulations after removing high spiking KCs.
-   `circular_run` : scripts for running running successive simulations
    after removing high spiking KCs. These scripts read last job ids
    from a specified file to identify the corresponding data dumps,
    remove high spiking kcs from those model templates, and run the
    simulation again until no more highspiking KC is left.
-   `run_fixed_net_changing_stim.py` : example utility script to run
    simulation of a given network model with changed PN inputs as a
    subprocess.
-   `run_fixed_net_with_ig.sh` : sample script to run a given network
    template including IG with different PN inputs.
-   `run_fixed_network_changing_stim.sh` : script to run a fixed network
    template while changing the PN input pattern.
-   `run_kc_ggn_feedback_amp_sweep.sh`,
    `run_kc_ggn_nofeedback_amp_sweep.sh` : scripts to test single KC with
    and without GGN feedback.
-   `run_mb_net.sh` : script to run mushroom body network model in batch mode.
-   `run_to_reproduce_ggn_vm_no_seed.sh`, `run_to_reproduce_ggn_vm.sh`,
    `run_to_reproduce_ggn_vm_with_shifting_pn.sh`: scripts that run
    multiple simulations while changing parameters to reproduce
    realistic GGN voltage trace.
-   `run_with_ig.sh` : script to simulate network model with IG
    included.



## `test_cell`

-   `ggn_voltage_attenuation_vclamp.py` : check voltage attenuation
    along GGN arbor when one branch is voltage clamped.
-   `run_param_sweep_vclamp.py` : test voltage attenuation along GGN
    arbor with voltage clamp while changing passive properties.



# morphutils

-   `cellmovie.py` : dump a video of rotating neuron in 3D
-   `displaycell.py` : display neuronal morphology from SWC file
-   `morph3d_*.py` : display neuronal morphology in 3D using
    corresponding module.  Use vtk/vispy/matplotlib to display 3D
    morphology.
-   `neurograph.py` : handle morphology as a graph.  This file contains
    data type definitions and functions to read an SWC file and turn
    it into a graph (using networkx).  There is also GGN specific
    mapping where I assign custom types to specific branches based on
    anatomical location.



# nrn

-   `nrnutils.py` : utilities for handling NEURON model
    -   convert a NEURON cell model into a networkx graph
    -   insert alpha synapses
    -   insert ion channel mechanisms
    -   set up recording of Vm
- `display_celltemplate.py`: 3D visualization of the morphology of a
  neuron specified in NEURON cell template
-   `localized_input_output.py` : apply synaptic inputs at specified
    branches.  This scripts runs simulation with synchronous synaptic
    inputs at multiple compartments on specific branches.
-   `localized_input_output_passive_sweep.py` : apply synaptic inputs at
    specified branches of GGN models with different passive
    properties.
-   `nrninit.bat`, `nrninit.sh` : batch file (Windows) and shell script
    (Linux).  Initialize PYTHONPATH to include various model related
    Python scripts
-   `staggered_input.py` : simulate synaptic input at random
    compartments in specified branch
