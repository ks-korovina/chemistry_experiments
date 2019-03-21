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
from utils.option_handler import load_options
from mols.funcs import SAScore

_TOL = 1e-5

func = SAScore
func = lambda mol: len(mol.smiles)


def gen_gp_test_data():
    """ Xs are molecules, Ys are some numeric value """
    n_all = 100

    mols = get_chembl(n_all * 3)
    ys = np.array([func(m) for m in mols])

    mols1, mols2, mols3 = mols[:n_all], mols[n_all:2*n_all], mols[2*n_all:3*n_all]
    ys1, ys2, ys3 = ys[:n_all], ys[n_all:2*n_all], ys[2*n_all:3*n_all]

    n_train = int(n_all * 0.8)
    ys = np.array([SAScore(m) for m in mols])

    X1_tr, X1_te = mols1[:n_train], mols1[n_train:]
    Y1_tr, Y1_te = ys1[:n_train], ys1[n_train:]
    X2_tr, X2_te = mols2[:n_train], mols2[n_train:]
    Y2_tr, Y2_te = ys2[:n_train], ys2[n_train:]
    X3_tr, X3_te = mols3[:n_train], mols3[n_train:]
    Y3_tr, Y3_te = ys3[:n_train], ys3[n_train:]

    return [(X1_tr, Y1_tr, X1_te, Y1_te),
            (X2_tr, Y2_tr, X2_te, Y2_te),
            (X3_tr, Y3_tr, X3_te, Y3_te)]


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


    def test_hallucinated_predictions(self):
        """ Testing hallucinated predictions for NNGP. """
        self.report('Testing hallucinated predictions for MolGP.')
        for dataset in self.datasets:
            for kernel_type in self.kernel_types:
                curr_gp = build_molgp_with_dataset(dataset, kernel_type)
                curr_preds, curr_stds = curr_gp.eval(dataset[2][1:], 'std')
                ha_preds, ha_stds = curr_gp.eval_with_hallucinated_observations(
                                    dataset[2][1:], dataset[0][4:] + dataset[2][:1], 'std')
                assert np.linalg.norm(curr_preds - ha_preds) < _TOL
                assert np.all(curr_stds > ha_stds)


class MolGPFitterTestCase(BaseTestClass):
    """ Contains unit tests for the TransportNNDistanceComputer class. """
    def __init__(self, *args, **kwargs):
        """ Constructor. """
        super(MolGPFitterTestCase, self).__init__(*args, **kwargs)
        self.datasets = gen_gp_test_data()
        self.kernel_types = ['wl_kernel', "edgehist_kernel"]

    @classmethod
    def _read_kernel_hyperparams(cls, kernel_hyperparams):
        """ Returns the kernel hyper-params as a string. """
        return ",".join([k + "_" + str(v)
            for k, v in kernel_hyperparams.items()])

    def test_marg_likelihood_and_prediction(self):
        """ Tests marginal likelihood and prediction. """
        # pylint: disable=too-many-locals
        self.report('Testing evaluation and marginal likelihood on Fitted GP.' +
                    ' Probabilistic test, might fail.')
        num_tests = 0
        num_lml_successes = 0
        num_const_err_successes = 0
        num_naive_err_successes = 0
        for dataset_idx, dataset in enumerate(self.datasets):
            for kernel_type in self.kernel_types:
                # Build the naive GP
                naive_gp = build_molgp_with_dataset(dataset, kernel_type)
                # Obtain a fitted GP
                fitted_gp = fit_molgp_with_dataset(dataset, kernel_type)
                # Obtain log marginal likelihoods
                naive_lml  = naive_gp.compute_log_marginal_likelihood()
                fitted_lml = fitted_gp.compute_log_marginal_likelihood()
                lml_succ = naive_lml < fitted_lml
                num_lml_successes += lml_succ
                self.report(('(%s, ntr=%d, nte=%d):: naive-gp-lml=%0.4f, ' +
                           'fitted-gp-lml=%0.4f, lml-succ=%d.')%(
                           kernel_type, len(dataset[0]), len(dataset[2]), naive_lml,
                           fitted_lml, lml_succ), 'test_result')
                # Predictions & Marginal likelihood
                naive_preds, _ = naive_gp.eval(dataset[2], 'std')
                naive_gp_err = compute_average_prediction_error(dataset, naive_preds)
                fitted_preds, _ = fitted_gp.eval(dataset[2], 'std')
                fitted_gp_err = compute_average_prediction_error(dataset, fitted_preds)
                const_err = compute_average_prediction_error(dataset, dataset[1].mean())
                const_err_succ = const_err > fitted_gp_err
                naive_err_succ = naive_gp_err > fitted_gp_err
                num_const_err_successes += const_err_succ
                num_naive_err_successes += naive_err_succ
                self.report(('  dataset #%d: const-err: %0.4f (%d), naive-gp-err=%0.4f (%d), ' +
                           'fitted-gp-err=%0.4f.')%(
                           dataset_idx, const_err, const_err_succ,
                           naive_gp_err, naive_err_succ, fitted_gp_err),
                          'test_result')
                # Print out betas
                self.report('  fitted kernel %s hyper-params: %s'%(kernel_type,
                          self._read_kernel_hyperparams(fitted_gp.kernel.hyperparams)),
                          'test_result')
                num_tests += 1
        # Print out some statistics
        lml_frac = float(num_lml_successes) / float(num_tests)
        const_err_frac = float(num_const_err_successes) / float(num_tests)
        naive_err_frac = float(num_naive_err_successes) / float(num_tests)
        self.report('LML tests success fraction = %d/%d = %0.4f'%(
                    num_lml_successes, num_tests, lml_frac), 'test_result')
        self.report('Const-err success fraction = %d/%d = %0.4f'%(
                    num_const_err_successes, num_tests, const_err_frac), 'test_result')
        self.report('Naive-GP-err success fraction = %d/%d = %0.4f'%(
                    num_naive_err_successes, num_tests, naive_err_frac), 'test_result')
        assert num_lml_successes == num_tests
        assert const_err_frac >= 0.5
        assert naive_err_frac >= 0.3


if __name__=="__main__":
    execute_tests()
