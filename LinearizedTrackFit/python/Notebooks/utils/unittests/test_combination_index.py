import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__).rstrip("unittests")))
import combinations

from unittest import TestCase

__author__ = 'demattia'


class TestCombinationIndex(TestCase):
    def test_combination_index(self):
        # Test a simple combination in the barrel
        layers = [5, 6, 7, 8, 9, 10]
        radius = [0]*6
        self.assertEquals(combinations.combination_index(layers, radius), 2016)
        # Test one with a disk (PS)
        layers = [5, 6, 7, 8, 9, 11]
        self.assertEquals(combinations.combination_index(layers, radius), 3040)
        # Test one with a disk (2S)
        radius[5] = 70
        self.assertEquals(combinations.combination_index(layers, radius), 2100192)
