"""
Class that performs molecule space traversal
Currently the interface is only targeting EA algorithms.
Author: kkorovin@cs.cmu.edu


TODO:
* Make a superclass for all exporatory algorithms
* Think about what a smarter explorer would do

"""

import numpy as np
import pandas as pd
from tqdm import tqdm

from collections import defaultdict
from time import time

from chemist_opt.blackbox_optimiser import blackbox_opt_args
from synth.forward_synth import RexgenForwardSynthesizer
from rdkit import Chem
from rdkit_contrib.sascorer import calculateScore as calculateSAScore
from mols.data_struct import Molecule
from datasets.loaders import get_chembl_prop, get_initial_pool

ga_opt_args = [
  # get_option_specs('num_mutations_per_epoch', False, 50,
  #   'Number of mutations per epoch.'),
  # get_option_specs('num_candidates_to_mutate_from', False, -1,
  #   'The number of candidates to choose the mutations from.'),
  # get_option_specs('fitness_sampler_scaling_const', False, 2,
  #   'The scaling constant for sampling according to exp_probs.'),
]


class Explorer:
    """
    TODO: common interface to-be
    """
    pass


class RandomExplorer(Explorer):
    """
    Implements a random evolutionary algorithm
    for exploring molecule space.
    """
    def __init__(self, fitness_func, initial_pool=None, max_pool_size=None):
        """
        Params:
        :fitness_func: function to optimize over evolution
        :initial_pool: just what it says
        :max_pool_size: int or None

        TODO:
        :mutation_op: mutates a given Molecule
        :crossover_op: takes two Molecules
                    and returns one new Molecule
        """
        self.fitness_func = fitness_func
        self.synth = RexgenForwardSynthesizer()
        if initial_pool is None:
            initial_pool = get_initial_pool()
        self.pool = initial_pool
        self.max_pool_size = max_pool_size

        # TODO: think whether to add additional *synthesized* pool

    def evolve_step(self):
        """
        TODO docs
        """

        # choose molecules to cross-over
        r_size = np.random.randint(2,3)
        mols = np.random.choice(self.pool, size=r_size)
        # print(mols)

        # evolve
        outcomes = self.synth.predict_outcome(mols)
        top_outcome = sorted(outcomes, key=lambda mol: self.fitness_func(mol))[-1]
        print("Newly generated mol value:", self.fitness_func(top_outcome))
        self.pool.append(top_outcome)

        # filter
        if self.max_pool_size is not None:
            self.pool = sorted(self.pool, key=lambda mol: self.fitness_func(mol))[-self.max_pool_size:]

    def evolve(self, capital):
        """
        Params:
        :data: start dataset (list of Molecules)
        :capital: number of steps or other cost of exploration
        """
        for _ in range(capital):
            self.evolve_step()

    def get_best(self, k):
        top = sorted(self.pool, key=lambda mol: self.fitness_func(mol))[-k:]
        return top


# APIs
# ======================================================================================
def ga_optimise_from_args(func, max_capital):
    # the func may accept iterable or a single Molecule
    mol = Molecule(smiles="c1cc(OCCCN2CCCCC2)ccc1CN1CCC2(CC1)OCCO2")
    try:
        func(mol)
        func_ = func
    except Exception as e:
        # print("Failed,", e)
        func_ = lambda m: func([m])

    explorer = RandomExplorer(func_)
    explorer.evolve(max_capital)
    top = explorer.get_best(k=1)
    val = func(top)
    return top, val

