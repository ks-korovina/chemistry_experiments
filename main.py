"""

Main experiment script
Author: kkorovin@cs.cmu.edu

"""

from chemist_opt.chemist import optimize_chemist
from mols.funcs import SAScore as mol_func
from chemist_opt.domains import MolDomain
from chemist_opt.function_caller import FunctionCaller
from chemist_opt.worker_manager import SyntheticWorkerManager
from chemist_opt.utils import prep_optimiser_args
from chemist_opt.chemist import all_chemist_args


if __name__ == "__main__":
    worker_manager = SyntheticWorkerManager(1, time_distro='const')
    func_caller = FunctionCaller(mol_func, MolDomain)
    options, options_clone, reporter = prep_optimiser_args(mol_func, all_chemist_args)
    opt_val, opt_pt, history = optimize_chemist(func_caller, worker_manager, 10,
                                                options=options, reporter=reporter)
    print(opt_val, opt_pt)
    print("History")
    print(history)
