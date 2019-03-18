"""

Molecular kernels unit tests.
A:  kkorovin@cs.cmu.edu

"""

from mols.data_struct import Molecule
from mols import mol_kernels
from utils.base_test_class import BaseTestClass, execute_tests

S1, S2, S3 = "Cc1ccccc1", "C1OC1", "CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H](N)C1"

class MolKernelsTestCase(BaseTestClass):
    def __init__(self, *args, **kwargs):
        super(MolKernelsTestCase, self).__init__(*args, **kwargs)

    def test_visualize(self):
        ## Plotting this test graph
        # plt.subplot(121)
        # nx.draw(graph)
        # plt.subplot(122)
        # nx.draw(graph, pos=nx.circular_layout(graph), nodecolor='r', edge_color='b')
        # plt.show()

        # TODO: how to visualize igraph mol?
        pass

    def test_conversions(self):
        mol = Molecule(S1)
        graph = mol_kernels.mol2graph_igraph(mol)
        print(graph)

    def test_edgehist_kernel(self):
        mols = [Molecule(S1), Molecule(S2)]
        print(mol_kernels.compute_edgehist_kernel(mols))

    def test_wl_kernel(self):
        mols = [Molecule(S1), Molecule(S2), Molecule(S3)]
        print(mol_kernels.compute_wl_kernel(mols))


if __name__=="__main__":
    execute_tests()