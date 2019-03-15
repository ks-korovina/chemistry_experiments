"""

Kernels on molecules:

1) Graph-based kernels: operate on graph representations of molecules
  Some graph kernel implementations are taken from https://github.com/jajupmochi/py-graph

2) String-based kernels: operate on SMILES strings
   TBA

A:  kkorovin@cs.cmu.edu

TODO:
    - Issue: for some reason graphkernels fails when graphs are attributed,
      this seems to be a serious drawback (that good are these kernels then?)

"""

GRAPH_LIB = "igraph"  # depending on package for graph kernels

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops
if GRAPH_LIB == "igraph":
    import igraph
else:
    import networkx
import graphkernels.kernels as gk


# Graph-based kernels ---------------------------------------------------------
# TODO: already in the Random.ipynb notebook, transfer here


def mol2graph_igraph(mol):
    """
    Convert molecule to nx.Graph

    Adapted from
    https://iwatobipen.wordpress.com/2016/12/30/convert-rdkit-molecule-object-to-igraph-graph-object/
    """
    mol = mol.to_rdkit()
    admatrix = rdmolops.GetAdjacencyMatrix(mol)
    bondidxs = [(b.GetBeginAtomIdx(),b.GetEndAtomIdx() ) for b in mol.GetBonds()]
    adlist = np.ndarray.tolist(admatrix)
    graph = igraph.Graph()
    g = graph.Adjacency(adlist).as_undirected()

    ## set properties
    # for idx in g.vs.indices:
    #     g.vs[idx][ "AtomicNum" ] = mol.GetAtomWithIdx(idx).GetAtomicNum()
    #     g.vs[idx][ "AtomicSymbole" ] = mol.GetAtomWithIdx(idx).GetSymbol()
    
    # for bd in bondidxs:
    #     btype = mol.GetBondBetweenAtoms(bd[0], bd[1]).GetBondTypeAsDouble()
    #     g.es[g.get_eid(bd[0], bd[1])]["BondType"] = btype
    #     print( bd, mol.GetBondBetweenAtoms(bd[0], bd[1]).GetBondTypeAsDouble() )
    return g


def mol2graph_networkx(mol):
    """
    Convert molecule to nx.Graph

    Adapted from
    https://iwatobipen.wordpress.com/2016/12/30/convert-rdkit-molecule-object-to-igraph-graph-object/
    """
    mol = mol.to_rdkit()
    admatrix = Chem.rdmolops.GetAdjacencyMatrix(mol)
    bondidxs = [(b.GetBeginAtomIdx(),b.GetEndAtomIdx() ) for b in mol.GetBonds()]
    graph = nx.Graph(admatrix)

    for idx in graph.nodes:
        graph.nodes[idx]["AtomicNum"] = mol.GetAtomWithIdx(idx).GetAtomicNum()
        graph.nodes[idx]["AtomicSymbol"] = mol.GetAtomWithIdx(idx).GetSymbol()

    for bd in bondidxs:
        btype = mol.GetBondBetweenAtoms(bd[0], bd[1]).GetBondTypeAsDouble()
        graph.edges[bd[0], bd[1]]["BondType"] = str(int(btype))
        # print(bd, m1.GetBondBetweenAtoms(bd[0], bd[1]).GetBondTypeAsDouble())
    return graph


"""
Kernels available in graphkernels: TODO into functions

    K1 = gk.CalculateEdgeHistKernel(graph_list)
    K2 = gk.CalculateVertexHistKernel(graph_list) 
    K3 = gk.CalculateVertexEdgeHistKernel(graph_list)
    K4 = gk.CalculateVertexVertexEdgeHistKernel(graph_list)
    K5 = gk.CalculateEdgeHistGaussKernel(graph_list)
    K6 = gk.CalculateVertexHistGaussKernel(graph_list)
    K7 = gk.CalculateVertexEdgeHistGaussKernel(graph_list)
    K8 = gk.CalculateGeometricRandomWalkKernel(graph_list)
    K9 = gk.CalculateExponentialRandomWalkKernel(graph_list)
    K10 = gk.CalculateKStepRandomWalkKernel(graph_list)
    K11 = gk.CalculateWLKernel(graph_list)
    K12 = gk.CalculateConnectedGraphletKernel(graph_list, 4)
    K13 = gk.CalculateGraphletKernel(graph_list, 4)
    K14 = gk.CalculateShortestPathKernel(graph_list)

"""


def compute_edgehist_kernel(mols):
    """
    Compute edge hist kernel
    Arguments:
        mols {list[Molecule]} -- [description]
    """
    mol_graphs_list = [mol2graph_igraph(m) for m in mols]
    return gk.CalculateEdgeHistKernel(mol_graphs_list)


def compute_wl_kernel(mols):
    """
    Compute edge hist kernel
    Arguments:
        mols {list[Molecule]} -- [description]
    """
    mol_graphs_list = [mol2graph_igraph(m) for m in mols]
    return gk.CalculateWLKernel(mol_graphs_list)



# String-based kernels ---------------------------------------------------------






