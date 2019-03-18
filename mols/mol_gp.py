"""

GP implementation working on molecular domains.
Kernels are in mol_kernels.py

Author:  kkorovin@cs.cmu.edu

*******************************************************************************

TODO:
* complete all places with ImplementMe

Notes:
self.tp_comp is an object that has .evaluate_single(x1, x2) method
it is a component of kernels in nasbot, so I shouldn't need it

"""

from mol_kernels import MolKernel
from utils.ancillary_utils import get_list_of_floats_as_str
from utils.reporters import get_reporter
from utils.option_handler import get_option_specs, load_options


mol_gp_specific_args = []  # check what these should be

mol_gp_args = gp_core.mandatory_gp_args + basic_gp_args + mol_gp_specific_args


# GP implementation for molecules ---------------------------------------------

class MolGP(gp_core.GP):
""" A Gaussian process for Neural Networks. """

    def __init__(self, X, Y, kernel_type, kernel_hyperparams, mean_func, noise_var, *args, **kwargs):
        """ Constructor.
            kernel_type: Should be one of [TODO: kernel names]
            kernel_hyperparams: is a python dictionary specifying the hyper-parameters for the
                                kernel. Will have parameters [TODO: which hyperparams of kernels]
            list_of_dists: Is a list of distances for the training set.
        """
        kernel = self._get_kernel_from_type(kernel_type, kernel_hyperparams)
        self.list_of_dists = list_of_dists
        # Call super constructor
        super(MolGP, self).__init__(X, Y, kernel, mean_func, noise_var,
                                     handle_non_psd_kernels='project_first', *args, **kwargs)

    # TODO: list_of_dists?
    def build_posterior(self):
        """ Checks if the sizes of list of distances and the data are consistent before
            calling build_posterior from the super class.
        """
        if self.list_of_dists is not None:
            assert self.list_of_dists[0].shape == (len(self.X), len(self.Y))
        super(MolGP, self).build_posterior()

    # TODO: list_of_dists?
    def set_list_of_dists(self, list_of_dists):
        """ Updates the list of distances. """
        self.list_of_dists = list_of_dists

    # TODO: list_of_dists/eval_...?
    def _get_training_kernel_matrix(self):
        """ Compute the kernel matrix from distances if they are provided. """
        if self.list_of_dists is not None:
            # TODO: see if this method is needed and makes sense
            return self.kernel.evaluate_from_dists(self.list_of_dists)
        else:
            return self.kernel(self.X, self.X)

    @classmethod
    def _get_kernel_from_type(cls, kernel_type, kernel_hyperparams):
        """ Returns the kernel. """
        return MolKernel(kernel_type)


# GP fitter for molecules -----------------------------------------------------

class MolGPFitter(gp_core.GPFitter):
    """
    Fits a GP by tuning the kernel hyper-params.
    This is the interface for MolGP.
    """
    def __init__(self, X, Y, options=None, reporter=None, *args, **kwargs):
        """ Constructor. """
        self.X = X
        self.Y = Y
        self.reporter = get_reporter(reporter)
        # self.list_of_dists = list_of_dists  # TODO: is this needed?
        self.num_data = len(X)
        if options is None:
            options = load_options(mol_gp_args, 'GPFitter', reporter=reporter)
        super(MolGPFitter, self).__init__(options, *args, **kwargs)

    def _child_set_up(self):
        """ technical method that sets up attributes, mostly """
        raise NotImplementedError("ImplementMe")

    def _preprocess_struct_mislabel_coeffs(self):
        """ also technical """
        raise NotImplementedError("ImplementMe")

    def _child_set_up_ml_tune(self):
        """
        Sets up tuning for Maximum likelihood.
        TODO: choose one kernel and add its params here.
        """
        raise NotImplementedError("ImplementMe")

    def _child_set_up_post_sampling(self):
        raise NotImplementedError("Not implemented post sampling yet.")

    def _child_build_gp(self, gp_hyperparams, build_posterior=True):
        """ Builds the GP from the hyper-parameters. """
        gp_hyperparams = gp_hyperparams[:]  # create a copy of the list
        kernel_hyperparams = {}
        # mean functions
        mean_func = gp_core.get_mean_func_from_options(self.options, self.Y)
        # extract GP hyper-parameters
        # 1. Noise variance ------------------------------------------------------
        noise_var = gp_core.get_noise_var_from_options_and_hyperparams(self.options,
                                                          gp_hyperparams, self.Y, 0)
        if self.options.noise_var_type == 'tune':
            gp_hyperparams = gp_hyperparams[1:]

        # SETTING kernel_hyperparams BASED ON KERNEL TYPE--------------------------

        # 2. kernel_hyperparams['scale'] and kernel_hyperparams['alphas']
        # 3. Setting kernel_hyperparams['betas']
        # 4 & 5. Set kernel_hyperparams['mislabel_coeffs']
        #        and kernel_hyperparams['struct_coeffs']
        # 6. Set kernel_hyperparams['non_assignment_penalty']
        #    = self.options.non_assignment_penalty

        return MolGP(self.X, self.Y, self.kernel_type, 
                     kernel_hyperparams, mean_func, noise_var,
                     build_posterior=build_posterior)


