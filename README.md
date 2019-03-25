# DOE framework for joint optimization and synthesis

This project aims at performing Bayesian optimization on organic molecules 

## Requirements:

* See `requirements.txt`

## Datasets:

ChEMBL data can be found here: https://github.com/kevinid/molecule_generator/releases/download/1.0/datasets.tar.gz

## Structure of the repo:

* `datasets` contains datasets and the `loaders.py` module, which provides dataset loading, initial pool for GA/Explorer, etc.
* `mols` isolates all domain-related things: data structures for molecules, molecular kernels, molecular GP implementation, and various helpers for handling molecules.
* `synth` contains models for performing synthesis (such as forward reaction prediction) in separation of the remaining framework. Data for training the model is stored in 'datasets/synth_data'. A better option would be training the predictor jointly with DOE, but that may be too difficult.
* `gp` contains basic GP for any domain. It is completely borrowed from [Nasbot](https://github.com/kirthevasank/nasbot)
* `explore` implements the search
* `chemist_opt` performs BO with custom kernels and exploration (acquisition optimization) strategies. It also follows the [Nasbot](https://github.com/kirthevasank/nasbot "Nasbot")
* `utils` contains misc utils, mostly from [Nasbot](https://github.com/kirthevasank/nasbot) at the moment being

## Getting started:

### Installing dependencies

**Fortran compiler.**

Computations are faster with `direct fortran` library installed (moreover, some of the computations seem to depend on having direct fortran). For this `cd` into `utils/direct_fortran` and run `bash make_direct.sh`. You will need a fortran compiler such as gnu95, which in the case of MacOS can be installed from [here](https://github.com/fxcoudert/gfortran-for-macOS). Once this is done, you can run `simple_direct_test.py` to make sure that it was installed correctly.

**Python packages.** 

Next, set up environment for RDKit:

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

* [Nasbot](https://github.com/kirthevasank/nasbot) - BO
* Rexgen - forward synthesis
* Datasets?


