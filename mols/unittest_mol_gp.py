"""
Unit tests for MolGP and MolGPFitter.
Author: kkorovin@cs.cmu.edu

NOTES:
*  basic test shows worse performance
    on WL kernel than on edgehist kernel

"""

import numpy as np

# Local imports
from mols import mol_gp
from datasets.loaders import get_chembl
from utils.base_test_class import BaseTestClass, execute_tests
from mols.funcs import SAScore


def gen_gp_test_data():
    """ Xs are molecules, Ys are some numeric value """
    n_all = 2000

    mols = get_chembl(n_all)
    n_train = int(n_all * 0.8)
    ys = np.array([SAScore(m) for m in mols])
    X1_tr, X1_te = mols[:n_train], mols[n_train:]
    Y1_tr, Y1_te = ys[:n_train], ys[n_train:]
    return [(X1_tr, Y1_tr, X1_te, Y1_te)]


# Some utilities we will need for testing ------------------------------------------------
def build_molgp_with_dataset(dataset, kernel_type):
  """ Builds a GP using the training set in dataset. """
  X, Y = dataset[:2]
  mean_func = lambda x: np.array([np.median(Y)] * len(x))
  noise_var = (Y.std() ** 2)/20
  kernel_hyperparams = get_kernel_hyperparams(kernel_type)
  return mol_gp.MolGP(X, Y, kernel_type,
                      kernel_hyperparams,
                      mean_func, noise_var)


def get_kernel_hyperparams(kernel_type):
    """ Returns the kernel hyperparams for the unit-tests below. """
    kernel_hyperparams = {}
    kernel_hyperparams["cont_par"] = 2.0
    kernel_hyperparams["int_par"] = 3
    return kernel_hyperparams


def fit_molgp_with_dataset(dataset, kernel_type):
    """ Fits an MolGP to this dataset. """
    options = load_options(mol_gp.mol_gp_args, '')
    options.kernel_type = kernel_type
    gp_fitter = mol_gp.MolGPFitter(dataset[0], dataset[1],
                                   options=options, reporter=None)
    _, fitted_gp, _ = gp_fitter.fit_gp()
    return fitted_gp


# TODO
def compute_average_prediction_error(dataset, preds, true_labels_idx=None):
    """ Computes the prediction error. """
    true_labels_idx = 3 if true_labels_idx is None else true_labels_idx
    return (np.linalg.norm(dataset[true_labels_idx] - preds)**2)/ \
            len(dataset[true_labels_idx])



class MolGPTestCase(BaseTestClass):
    def __init__(self, *args, **kwargs):
        """ Constructor. """
        super(MolGPTestCase, self).__init__(*args, **kwargs)
        self.datasets = gen_gp_test_data()
        self.kernel_types = ['edgehist_kernel', 'wl_kernel']

    def test_basics(self):
        """ Tests for adding data, evaluation and marginal likelihood. """
        self.report('Testing adding data, evaluation and marginal likelihood.' +
                    ' Probabilistic test, might fail.')
        num_tests = 0
        num_successes = 0
        for dataset in self.datasets:
            for kernel_type in self.kernel_types:
                curr_gp = build_molgp_with_dataset(dataset, kernel_type)
                # Predictions & Marginal likelihood
                curr_preds, _ = curr_gp.eval(dataset[2], 'std')
                curr_gp_err = compute_average_prediction_error(dataset, curr_preds)
                const_err = compute_average_prediction_error(dataset, dataset[1].mean())
                lml = curr_gp.compute_log_marginal_likelihood()
                is_success = curr_gp_err < const_err
                num_tests += 1
                num_successes += is_success
                self.report(('(%s, ntr=%d, nte=%d):: GP-lml=%0.4f, GP-err=%0.4f, ' +
                           'Const-err=%0.4f.  succ=%d')%(kernel_type, len(dataset[0]),
                           len(dataset[2]), lml, curr_gp_err, const_err, is_success),
                           'test_result')
        
        succ_frac = num_successes / float(num_tests)
        self.report('Summary: num_successes / num_floats = %d/%d = %0.4f'%(num_successes,
                    num_tests, succ_frac), 'test_result')
        assert succ_frac > 0.5


class MolGPFitterTestCase(BaseTestClass):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def test1(self):
        pass


if __name__=="__main__":
    execute_tests()
