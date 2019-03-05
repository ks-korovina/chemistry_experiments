"""

Class that performs molecule space traversal
Currently the interface is only targeting EA algorithms.

TODO:
* Maybe the Explorer should behave as follows:
  starting up, construct a database to which the operations
  would be able to consult (like a chemistry knowledge base)

* Make a superclass for all exporatory algorithms

"""

import numpy as np
from time import time
from synth.forward_synth import RexgenForwardSynthesizer
from rdkit import Chem
from rdkit_contrib.sascorer import calculateScore as calculateSAScore
from data.data_struct import Molecule

class Explorer:
	"""
	What should the common interface be?
	"""
	pass


class RandomExplorer(Explorer):
	"""
	Implements a random evolutionary algorithm
	for exploring molecule space.
	"""
	def __init__(self, fitness_func, initial_pool, max_pool_size=None):
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
		self.pool = initial_pool
		self.max_pool_size = max_pool_size

		# TODO: how to keep synthesis paths?


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

	# not necessarily a good idea
	# def __next__(self):
	# 	pass



############## Mutation ops ####################
# Maybe this should be factored out of here    #
# Also perhaps these should be functor classes #
################################################

def mutate_random(mol, database):
	"""
	Mutate a molecule by randomly adding a valid component.
	Database knows which components are ok to add.
	"""
	
	pass

def mutate_forward_synth_random(mol, database):
	"""
	Takes a molecule and randomly "bump" into each other.
	Database has available molecules inside.
	"""
	pass

# mutation ops that are guided:
def mutate_random_guided(mol, eval_func, database):
	pass

def mutate_forward_synth_guided(mol, eval_func, database):
	pass


if __name__=="__main__":
	# Temporarily funcs operate on SMILES strings
	dummy_func = lambda mol: len(mol.smiles)
	sas_func = lambda mol: calculateSAScore(Chem.MolFromSmiles(mol.smiles))
	print(sas_func(Molecule("CC")))

	test_pool = ["CC", "O=C=O", "C#N", "CCN(CC)CC", "CC(=O)O", "C1CCCCC1", "c1ccccc1"]
	test_pool = [Molecule(smiles) for smiles in test_pool]

	###### Stupid example with len-of-smiles fitness #
	exp = RandomExplorer(dummy_func, initial_pool=test_pool)
	print("Starting len of SMILES optimization")
	exp.evolve(2)

	#check
	print(exp.pool)

	###### SA optimization example #
	exp = RandomExplorer(sas_func, initial_pool=test_pool)
	print("Starting SA score optimization")
	t0 = time()
	exp.evolve(10)

	#check
	print("Completed SA score optimization, time elapsed: %.3fs" % (time()-t0))
	print(exp.pool)
	top = exp.get_best(1)[0]
	print(top.get_synthesis_path())






