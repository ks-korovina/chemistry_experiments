"""

Implements forward synthesis

TODO:
* convert from string representations to general Molecule class

Notes:
* Using pretrained models - are there any?
* There is code for MoleculeTransformers:
  see https://github.com/pschwllr/MolecularTransformer
  so train it?
* Another option:
  https://github.com/connorcoley/rexgen_direct

"""

import sys
# sys.path.append("/Users/ksk/Desktop/CMU/DOE/chemistry")
sys.path.append("/Users/ksk/Desktop/CMU/DOE/chemistry/synth")

from rexgen_direct.core_wln_global.directcorefinder import DirectCoreFinder 
from rexgen_direct.scripts.eval_by_smiles import edit_mol
from rexgen_direct.rank_diff_wln.directcandranker import DirectCandRanker

TEMP = ["[CH3:26][c:27]1[cH:28][cH:29][cH:30][cH:31][cH:32]1", 
        "[Cl:18][C:19](=[O:20])[O:21][C:22]([Cl:23])([Cl:24])[Cl:25]",
        "[NH2:1][c:2]1[cH:3][cH:4][c:5]([Br:17])[c:6]2[c:10]1[O:9][C:8]([CH3:11])([C:12](=[O:13])[O:14][CH2:15][CH3:16])[CH2:7]2"
        ]


class ForwardSynthesizer:
    """
    Class for answering forward prediction queries.

    """
    def __init__(self):
        # load trained model
        pass

    def predict_outcome(self, list_of_mols):
        """
        Using a predictor, produce the most likely reaction
        
        Params:
        :list_of_mols: list of reactants and reagents
                       (former contribute atoms, latter don't)
        TODO: what else?
        """
        pass


class RexgenForwardSynthesizer:
    def __init__(self):
        # load trained model
        self.directcorefinder = DirectCoreFinder()
        self.directcorefinder.load_model()
        self.directcandranker = DirectCandRanker()
        self.directcandranker.load_model()

    def predict_outcome(self, list_of_mols=TEMP, k=1):
        """
        Using a predictor, produce top-k most likely reactions
        
        Params:
        :list_of_mols: list of reactants and reagents (smiles format)
                       (former contribute atoms, latter don't)
        """
        react = ".".join(list_of_mols)
        (react, bond_preds, bond_scores, cur_att_score) = self.directcorefinder.predict(react)
        outcomes = self.directcandranker.predict(react, bond_preds, bond_scores)
        # for outcome in outcomes:
        #     print(outcome)
        res = [out["smiles"][0] for out in outcomes[:k]]
        return res


if __name__=="__main__":
    t = RexgenForwardSynthesizer()
    t.predict_outcome()


