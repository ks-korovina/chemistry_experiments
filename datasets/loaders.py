"""
Loading and handling chemical data
Author: kkorovin@cs.cmu.edu

"""

import pandas as pd
from mols.data_struct import Molecule


def get_chembl_prop(n_mols=None, as_mols=False):
    """ Returns (pool, smile->prop mappings) """
    path = "./datasets/ChEMBL_prop.txt"
    df = pd.read_csv(path, sep="\t", header=None)
    # smile: v for the first of two properties
    smile_to_prop = {s: v for (s, v) in zip(df[0], df[1])}
    smile_to_prop = defaultdict(int, smile_to_prop)
    return df[0].values, smile_to_prop


def get_chembl(n_mols=None, as_mols=True):
    """ Return list of SMILES """
    path = "./datasets/ChEMBL.txt"
    with open(path, "r") as f:
        if n_mols is None:
            res = [line.strip() for line in f]
        else:
            res = [f.readline().strip() for _ in range(n_mols)]
    return [Molecule(smile) for smile in res]


if __name__ == "__main__":
    print(get_chembl()[:4])
