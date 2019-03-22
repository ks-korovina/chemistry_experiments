"""

Chemist unit tests

TODO: adapt all of these from NASBOT to Chemist

"""

from chemist_opt.function_caller import FunctionCaller
from chemist_opt.worker_manager import SyntheticWorkerManager

from chemist_opt import chemist
from utils.base_test_class import BaseTestClass, execute_tests
from chemist_opt.worker_manager import SyntheticWorkerManager


class ChemistTestCase(BaseTestClass):
    def setUp(self):
        """ Set up. """
        pass

    #     ret = get_mol_opt_arguments()
    #     for key, val in ret.__dict__.iteritems():
    #       setattr(self, key, val)

    def test_instantiation(self):
        """ Test Creation of NASBOT object. """
        pass

        # self.report('Testing Random Optimiser instantiation.')
        # tp_comp = get_tp_comp('cnn')
        # func_caller = FunctionCaller(cnn_syn_func1, self.cnn_domain)
        # worker_manager = SyntheticWorkerManager(1, time_distro='const')
        # optimiser = nasbot.NASBOT(func_caller, worker_manager, tp_comp, reporter='silent')
        # for attr in dir(optimiser):
        #     if not attr.startswith('_'):
        #         self.report('optimiser.%s = %s'%(attr, str(getattr(optimiser, attr))),
        #                                         'test_result')

    def _get_optimiser_args(self):
        """ Returns the options and reporter. """
        pass

        # return get_optimiser_args(self, nn_type, nasbot.all_nasbot_args)

    @classmethod
    def _test_optimiser_results(cls, opt_val, _, history, options, options_clone):
        """ Tests optimiser results. """
        pass

        # assert opt_val == history.curr_opt_vals[-1]
        # test_if_pre_eval_networks_have_changed(options, options_clone)

    def test_nasbot_optimisation_single(self):
        """ Tests optimisation with a single worker. """
        pass

        # self.report('Testing NASBOT with a single worker.')
        # worker_manager = SyntheticWorkerManager(1, time_distro='const')
        # func_caller = FunctionCaller(cnn_syn_func1, self.cnn_domain)
        # tp_comp = get_tp_comp('cnn')
        # options, options_clone, reporter, _, _ = self._get_optimiser_args('cnn')
        # opt_val, opt_pt, history = nasbot.nasbot(
        #                       func_caller, worker_manager, 10, tp_comp,
        #                       options=options, reporter=reporter)
        # self._test_optimiser_results(opt_val, opt_pt, history, options, options_clone)
        # self.report('')

if __name__ == "__main__":
    execute_tests()

