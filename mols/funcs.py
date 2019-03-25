"""

Functions that operate on Molecules
TODO: write a better wrapper

"""

# import functools
from rdkit import Chem
from rdkit_contrib.sascorer import calculateScore as calculateSAScore


def SAScore(mol):
    """ Synthetic accessibility score """
    try:
    	rdkit_mol = Chem.MolFromSmiles(mol.smiles)
    except:
    	rdkit_mol = Chem.MolFromSmiles(mol[0].smiles)
    return calculateSAScore(rdkit_mol)

def SMILES_len(mol):
	return len(mol.smiles)