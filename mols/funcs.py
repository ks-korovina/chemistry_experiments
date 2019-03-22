"""

Functions that operate on Molecules

"""

from rdkit import Chem
from rdkit_contrib.sascorer import calculateScore as calculateSAScore

def SAScore(mol):
    """ Synthetic accessibility score """
    rdkit_mol = Chem.MolFromSmiles(mol.smiles)
    return calculateSAScore(rdkit_mol)

def SMILES_len(mol):
	return len(mol.smiles)