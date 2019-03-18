"""
Module for BO with graph kernel and synthesizeable exploration.

Core class is Chemist

TODO:
* complete all places with ImplementMe

"""

import numpy as np

# Local imports
from opt.blackbox_optimiser import blackbox_opt_args
from mols.mol_gp import mol_gp_args, MolGPFitter
from mols.mol_kernels import *  # kernel names
from opt.nn_opt_utils import get_initial_pool
from opt.gp_bandit import GPBandit, gp_bandit_args
from utils.general_utils import block_augment_array
from utils.reporters import get_reporter
from utils.option_handler import get_option_specs, load_options


chemist_specific_args = []  # TODO

all_chemist_args = chemist_specific_args + gp_bandit_args + \
                     blackbox_opt_args + nn_gp_args


class Chemist(GPBandit):
    """
    Analog of NASBOT class.
    To not have it inherit from any GPBandit,
    must merge and simplify all functionality.
    """
    def __init__(self, func_caller, worker_manager, tp_comp, options=None, reporter=None):
        # self.options.acq_opt_method = self.options.chemist_acq_opt_method
        # super(Chemist, self)._child_set_up()
        # self.list_of_dists = None
        # self.already_evaluated_dists_for = None
        # # Create a GP fitter with no data and use its tp_comp as the bandit's tp_comp
        # init_gp_fitter = MolGPFitter([], [], self.domain.get_type(), tp_comp=self.tp_comp,
        #                             list_of_dists=None, options=self.options,
        #                             reporter=self.reporter)
        # self.tp_comp = init_gp_fitter.tp_comp
        # self.mislabel_coeffs = init_gp_fitter.mislabel_coeffs
        # self.struct_coeffs = init_gp_fitter.struct_coeffs
        
        raise NotImplementedError("ImplementMe")

    def _child_set_up(self):
        """ Child up. """
        raise NotImplementedError("ImplementMe")

    def _set_up_acq_opt_ga(self):
        raise NotImplementedError("ImplementMe")

    def _compute_list_of_dists(self, X1, X2):
        raise NotImplementedError("ImplementMe")

    def _get_gp_fitter(self, reg_X, reg_Y):
        raise NotImplementedError("ImplementMe")

    def _add_data_to_gp(self, new_points, new_vals):
        raise NotImplementedError("ImplementMe")

    def _child_set_gp_data(self, reg_X, reg_Y):
        raise NotImplementedError("ImplementMe")


# APIs ---------------------------------------------------------

# TODO

