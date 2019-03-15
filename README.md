# DOE framework for joint optimization and synthesis

This project aims at performing Bayesian optimization on organic molecules 

## Requirements:

* See `requirements.txt`

## Datasets:

ChEMBL data can be found here: https://github.com/kevinid/molecule_generator/releases/download/1.0/datasets.tar.gz

## Structure of the repo:

* `datasets` contains datasets
* `mols` isolates all domain-related things: data structures for molecules, molecular kernels, molecular GP implementation, and various helpers for handling molecules.
* `synth` contains models for performing synthesis (such as forward reaction prediction) in separation of the remaining framework. Data for training the model is stored in 'datasets/synth_data'. A better option would be training the predictor jointly with DOE, but that may be too difficult.
* `gp` contains basic GP for any domain. It is completely borrowed from [Nasbot](https://github.com/kirthevasank/nasbot)
* `explore` implements the search
* `chemist_opt` performs BO with custom kernels and exploration (acquisition optimization) strategies. It also follows the [Nasbot](https://github.com/kirthevasank/nasbot "Nasbot")
* `utils` contains misc utils, mostly from [Nasbot](https://github.com/kirthevasank/nasbot) at the moment being

## Getting started:

### Installing dependencies

Set up environment for RDKit:

```bash
conda create -c rdkit -n my-rdkit-env rdkit
conda activate my-rdkit-env
```

Install requirements with conda:

```bash
while read requirement; do conda install --yes $requirement; done < requirements.txt
```

or pip:

```bash
pip install -r requirements.txt
```

In addition to these requirements, a `graphkernels` package should be installed. It automatically installs `igraph` and other dependencies. However, it does not install `eigen3`, `pkg-config`, therefore those are included into requirements. Install `graphkernels` via pip (on Mac):

```bash
pip install graphkernels
```

If the above fails on MacOS (see [stackoverflow](https://stackoverflow.com/questions/16229297/why-is-the-c-standard-library-not-working)), the simplest solution is

```bash
MACOSX_DEPLOYMENT_TARGET=10.9 pip install graphkernels
```

### Environment

Next, set PYTHONPATH for imports:

```bash
export PYTHONPATH=":$PWD/rdkit_contrib:$PWD/synth/:$PWD/synth/rexgen_direct"
```

Or run setup script:

```bash
source setup.sh 
```

## Running tests

TBA

## Credits

TBA

