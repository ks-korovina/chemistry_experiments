"""
  Unit tests for explorer.py
"""

import numpy as np

# Local imports
from explore import explorer
from utils.base_test_class import BaseTestClass, execute_tests


class ExplorerTestCase(BaseTestClass):
  def __init__(self, *args, **kwargs):
    super(ExplorerTestCase, self).__init__(*args, **kwargs)

  def test_testcase_pass(self):
    print("Nothing to do, pass")


if __name__ == '__main__':
  execute_tests()