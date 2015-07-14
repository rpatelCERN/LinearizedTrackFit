from unittest import TestCase

__author__ = 'demattia'

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__).rstrip("unittests")))
import combinations


class TestCombinationsGenerator(TestCase):
    def test_combinations_generator(self):
        c_gen = combinations.CombinationsGenerator()
        self.assertEqual(c_gen.combinations(6), {6: [[0, 1, 2, 3, 4], [0, 1, 2, 3, 5], [0, 1, 2, 4, 5], [0, 1, 3, 4, 5], [0, 2, 3, 4, 5], [1, 2, 3, 4, 5]]})
