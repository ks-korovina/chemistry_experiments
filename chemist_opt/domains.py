"""
    Harness to manage optimisation domains.
    -- kandasamy@cs.cmu.edu,
       kkorovin@cs.cmu.edu
"""

# pylint: disable=no-member
# pylint: disable=invalid-name
# pylint: disable=arguments-differ
# pylint: disable=abstract-class-not-used

import numpy as np
# Local
from explore.explorer import ga_opt_args, ga_optimise_from_args
from gp.kernel import SEKernel
from utils.oper_utils import random_maximise, direct_ft_maximise
from utils.option_handler import load_options
from utils.reporters import get_reporter
from chemist_opt.worker_manager import SyntheticWorkerManager


_EUCLIDEAN_DFLT_OPT_METHOD = 'rand'
_NN_DFLT_OPT_METHOD = 'ga'


class Domain(object):
    """ Domain class. An abstract class which implements domains. """

    def __init__(self, dflt_domain_opt_method):
        """ Constructor. """
        super(Domain, self).__init__()
        self.dflt_domain_opt_method = dflt_domain_opt_method

    def maximise_obj(self, opt_method, obj, num_evals, *args, **kwargs):
        """ Optimises the objective and returns it. """
        if opt_method == 'dflt_domain_opt_method':
            opt_method = self.dflt_domain_opt_method
        print("MYSELF:", self)
        return self._child_maximise_obj(opt_method, obj, num_evals, *args, **kwargs)

    def _child_maximise_obj(self, opt_method, obj, num_evals, *args, **kwargs):
        """ Child class implementation for optimising an objective. """
        raise NotImplementedError('Implement in a child class.')

    def get_default_kernel(self, *args, **kwargs):
        """ Get the default kernel for this domain. """
        raise NotImplementedError('Implement in a child class.')

    def get_type(self):
        """ Returns the type of the domain. """
        raise NotImplementedError('Implement in a child class.')

    def get_dim(self):
        """ Returns the dimension of the space. """
        raise NotImplementedError('Implement in a child class.')


# For euclidean spaces ---------------------------------------------------------------
class EuclideanDomain(Domain):
    """ Domain for Euclidean spaces. """

    def __init__(self, bounds):
        """ Constructor. """
        self.bounds = np.array(bounds)
        self._dim = len(bounds)
        super(EuclideanDomain, self).__init__('rand')

    def _child_maximise_obj(self, opt_method, obj, num_evals):
        """ Child class implementation for optimising an objective. """
        if opt_method == 'rand':
            return self._rand_maximise_obj(obj, num_evals)
        elif opt_method == 'direct':
            return self._direct_maximise_obj(obj, num_evals)
        else:
            raise ValueError('Unknown opt_method=%s for EuclideanDomain'%(opt_method))

    def _rand_maximise_obj(self, obj, num_evals):
        """ Maximise with random evaluations. """
        if num_evals is None:
            lead_const = 10 * min(5, self.dim)**2
            num_evals = lambda t: np.clip(lead_const * np.sqrt(min(t, 1000)), 2000, 3e4)
        opt_val, opt_pt = random_maximise(obj, self.bounds, num_evals)
        return opt_val, opt_pt

    def _direct_maximise_obj(self, obj, num_evals):
        """ Maximise with direct. """
        if num_evals is None:
            lead_const = 10 * min(5, self.dim)**2
            num_evals = lambda t: np.clip(lead_const * np.sqrt(min(t, 1000)), 2000, 3e4)
        lb = self.bounds[:, 0]
        ub = self.bounds[:, 1]
        opt_val, opt_pt, _ = direct_ft_maximise(obj, lb, ub, num_evals)
        return opt_val, opt_pt

    def get_default_kernel(self, range_Y):
        """ Returns the default (SE) kernel. """
        return SEKernel(self.dim, range_Y/4.0, dim_bandwidths=0.05*np.sqrt(self.dim))

    def get_type(self):
        """ Returns the type of the domain. """
        return 'euclidean'

    def get_dim(self):
        """ Return the dimensions. """
        return self._dim


# For molecules. ---------------------------------------------------------------
class MolDomain(Domain):
    """ Domain for Molecules. """
    def __init__(self, constraint_checker=None):
        # this may check validity of molecules:
        self.constraint_checker = constraint_checker
        super(MolDomain, self).__init__('ga')

    @staticmethod
    def maximise_obj(opt_method, obj, num_evals, *args, **kwargs):
        """ Optimises the objective and returns it. """
        opt_pt, opt_val = ga_optimise_from_args(obj, num_evals)
        return opt_val, opt_pt

    def _child_maximise_obj(self, opt_method, obj, num_evals, *args, **kwargs):
        """ Child class implementation for optimising an objective. """
        if opt_method == 'ga':
            return self._ga_maximise(obj, num_evals, *args, **kwargs)
        elif opt_method == 'rand_ga':
            return self._rand_ga_maximise(obj, num_evals)
        else:
            raise ValueError('Unknown method=%s for MolDomain'%(opt_method))

    def _ga_maximise(self, obj, num_evals, mutation_op, 
                     init_pool, init_pool_vals=None,
                     expects_inputs_to_be_iterable=False):
        """ Maximise with genetic algorithms.
            if expects_inputs_as_list is True it means the function expects the inputs to
            be iterable by default.
            Arguments:
                obj - target function to optimize (*objective*)
        """
        # Optimization happens here:
        opt_pt, opt_val = ga_optimise_from_args(obj, num_evals)
        return opt_val, opt_pt

    def _rand_ga_maximise(self, obj, num_evals):
        """ Maximise over the space of neural networks via rand_ga. """
        raise NotImplementedError('Not implemented rand_ga for MolDomain yet.')

    def get_default_kernel(self, tp_comp, mislabel_coeffs, struct_coefs, powers,
                                                 dist_type, range_Y):
        """ Returns the default (SE) kernel. """
        raise NotImplementedError("TODO.")

    def get_dim(self):
        """ Return the dimensions. """
        raise NotImplementedError("TODO.")

