"""

Chemist unit tests

TODO: adapt all of these from NASBOT to Chemist

"""

from argparse import Namespace
from copy import deepcopy

from chemist_opt.function_caller import FunctionCaller
from chemist_opt.worker_manager import SyntheticWorkerManager
from utils.reporters import get_reporter
from utils.option_handler import load_options
from datasets.loaders import get_initial_pool

from chemist_opt import chemist
from mols.funcs import SAScore as mol_func
from chemist_opt.domains import MolDomain
from utils.base_test_class import BaseTestClass, execute_tests
from chemist_opt.worker_manager import SyntheticWorkerManager


def get_mol_opt_arguments():
    """ Returns arguments for NN Optimisation. """
    ret = Namespace()
    ret.ga_init_pool = get_initial_pool()
    ret.ga_init_vals = [mol_func(mol) for mol in ret.ga_init_pool]
    ret.mol_domain = MolDomain()
    return ret


class ChemistTestCase(BaseTestClass):
    def setUp(self):
        """ Set up. """
        ret = get_mol_opt_arguments()
        for key, val in ret.__dict__.items():
            setattr(self, key, val)

    def test_instantiation(self):
        """ Test Creation of Chemist object. """
        self.report('Testing Random Optimiser instantiation.')
        func_caller = FunctionCaller(mol_func, MolDomain)
        worker_manager = SyntheticWorkerManager(1, time_distro='const')
        optimiser = chemist.Chemist(func_caller, worker_manager, reporter='silent')
        for attr in dir(optimiser):
            if not attr.startswith('_'):
                self.report('optimiser.%s = %s'%(attr, str(getattr(optimiser, attr))),
                                                'test_result')

    def _get_optimiser_args(self, optimiser_args=chemist.all_chemist_args):
        """ Returns the options and reporter. """
        reporter = get_reporter('default')
        options = load_options(optimiser_args, reporter=reporter)
        options.pre_eval_points = self.ga_init_pool
        options.pre_eval_vals = self.ga_init_vals
        options.pre_eval_true_vals = self.ga_init_vals
        options_clone = deepcopy(options)
        return options, options_clone, reporter

    @classmethod
    def _test_optimiser_results(cls, opt_val, _, history, options, options_clone):
        """ Tests optimiser results. """
        pass

        # assert opt_val == history.curr_opt_vals[-1]
        # test_if_pre_eval_networks_have_changed(options, options_clone)

    def test_nasbot_optimisation_single(self):
        """ Tests optimisation with a single worker. """
        pass

        self.report('Testing Chemist with a single worker.')
        worker_manager = SyntheticWorkerManager(1, time_distro='const')
        func_caller = FunctionCaller(mol_func, MolDomain)
        options, options_clone, reporter = self._get_optimiser_args()
        opt_val, opt_pt, history = chemist.chemist(func_caller, worker_manager, 10,
                                                    options=options, reporter=reporter)
        # self._test_optimiser_results(opt_val, opt_pt, history, options, options_clone)
        # self.report('')

if __name__ == "__main__":
    execute_tests()

