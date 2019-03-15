"""

GP implementation working on molecular domains.
Kernels are in mol_kernels.py

TODO:
* complete all places with ImplementMe

"""

from mol_kernels import *
from utils.ancillary_utils import get_list_of_floats_as_str
from utils.reporters import get_reporter
from utils.option_handler import get_option_specs, load_options


mol_gp_specific_args = []  # check what these should be

mol_gp_args = gp_core.mandatory_gp_args + basic_gp_args + mol_gp_specific_args


# GP implementation for molecules ---------------------------------------------

class MolGP(gp_core.GP):
""" A Gaussian process for Neural Networks. """

  def __init__(self, X, Y, kernel_type, kernel_hyperparams, mean_func, noise_var,
        list_of_dists=None, *args, **kwargs):
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

  def build_posterior(self):
    """ Checks if the sizes of list of distances and the data are consistent before
      calling build_posterior from the super class.
    """
    if self.list_of_dists is not None:
      assert self.list_of_dists[0].shape == (len(self.X), len(self.Y))
    super(MolGP, self).build_posterior()

  def set_list_of_dists(self, list_of_dists):
    """ Updates the list of distances. """
    self.list_of_dists = list_of_dists

  def _get_training_kernel_matrix(self):
    """ Compute the kernel matrix from distances if they are provided. """
    if self.list_of_dists is not None:
      return self.kernel.evaluate_from_dists(self.list_of_dists)
    else:
      return self.kernel(self.X, self.X)

  @classmethod
  def _get_kernel_from_type(cls, kernel_type, kernel_hyperparams):
    """ Returns the kernel. """
    if kernel_type == "SOME_NAME":
      return mol_kernels.SomeKernel()
    else:
      raise ValueError('Unknown kernel_type %s.'%(kernel_type))


# GP fitter for molecules -----------------------------------------------------

class MolGPFitter(gp_core.GPFitter):
  """ Fits a GP by tuning the kernel hyper-params. """
  def __init__(self, X, Y, nn_type, tp_comp=None, list_of_dists=None,
         options=None, reporter=None, *args, **kwargs):
    """ Constructor. """
    self.X = X
    self.Y = Y
    self.reporter = get_reporter(reporter)
    self.list_of_dists = list_of_dists
    self.num_data = len(X)
    if options is None:
      options = load_options(mol_gp_args, 'GPFitter', reporter=reporter)
    super(MolGPFitter, self).__init__(options, *args, **kwargs)

  def _child_set_up(self):
    raise NotImplementedError("ImplementMe")

  def _preprocess_struct_mislabel_coeffs(self):
    raise NotImplementedError("ImplementMe")

  def _child_set_up_ml_tune(self):
    """ Sets up tuning for Maximum likelihood. """
    raise NotImplementedError("ImplementMe")

  def _child_set_up_post_sampling(self):
    raise NotImplementedError("Not implemented in NASBOT")

  def _child_build_gp(self, gp_hyperparams, build_posterior=True):
    """ Builds the GP from the hyper-parameters. """
    raise NotImplementedError("ImplementMe")

    return MolGP(self.X, self.Y, self.kernel_type, 
           kernel_hyperparams,mean_func, noise_var,
           list_of_dists=self.list_of_dists,
           build_posterior=build_posterior)


