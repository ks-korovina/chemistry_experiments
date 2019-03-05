# DOE framework for joint optimization and synthesis

## Upcoming TODOs:
* Implement the simplest Explorer, with random mutations, not necessarily anything synthesizeable
* Deal with data

## Requirements:

* rdkit for working with the ChEMBL dataset
* networkx for working with graphs

## Structure of the repo:

* **data** contains datasets, data structures and helpers for handling molecules.
* **synth** contains models for performing synthesis (such as forward reaction prediction) in separation of the remaining framework. Data for training the model is stored in 'data/datasets/synth_data'. A better option would be training the predictor jointly with DOE, but that may be too difficult.

## Getting started:

Set up environment for RDKit:

```
conda create -c rdkit -n my-rdkit-env rdkit
conda activate my-rdkit-env
```

Set PYTHONPATH for imports:

```
export PYTHONPATH=":$PWD/rdkit_contrib:$PWD/synth/:$PWD/synth/rexgen_direct"
```
