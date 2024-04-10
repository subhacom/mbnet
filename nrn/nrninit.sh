#!/bin/bash

# To initialize PYTHONPATH variable before running a simulation, run this script from a bash shell like this:
# `source nrn/nrninit.sh`

if [[ -z "${MODEL_DIR}" ]]; then
    export MODEL_DIR=$(pwd)
    echo "Setting model root directory to $MODEL_DIR"
    echo "To use a custom path, set the MODEL_DIR environment variable to the desired path."
fi

export PYTHONPATH="${MODEL_DIR}/mb:${MODEL_DIR}/nrn:${MODEL_DIR}/morphutils:${MODEL_DIR}/common"
echo "PYTHONPATH=$PYTHONPATH"

