"""

Chemist utils
Author: kkorovin@cs.cmu.edu

"""

from copy import deepcopy
from utils.reporters import get_reporter
from utils.option_handler import load_options
from datasets.loaders import get_initial_pool


def prep_optimiser_args(obj_func, optimiser_args):
    """ Returns the options and reporter. """
    reporter = get_reporter('default')
    options = load_options(optimiser_args, reporter=reporter)
    options.pre_eval_points = get_initial_pool()
    options.pre_eval_vals = [obj_func(mol) for mol in options.pre_eval_points]
    options.pre_eval_true_vals = options.pre_eval_vals
    options_clone = deepcopy(options)
    return options, options_clone, reporter
