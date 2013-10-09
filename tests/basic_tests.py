# -*- coding: utf-8 -*-

# Tests for Mosaic utility functions

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import unittest

import numpy as N

from mosaic.api import validate_units

def valid_units(units):
    try:
        validate_units(units, "")
        return True
    except (ValueError, TypeError):
        return False

class UnitStringTest(unittest.TestCase):

    def test_units(self):
        self.assertTrue(valid_units("nm"))
        self.assertTrue(valid_units("nm ps-1"))
        self.assertTrue(valid_units("g mol-1"))
        self.assertTrue(valid_units("60 s"))
        self.assertTrue(valid_units("1000 m"))
        self.assertTrue(valid_units("0.001 m"))
        # Not a string
        self.assertFalse(valid_units(42))
        # Invalid symbols
        self.assertFalse(valid_units("foo"))
        self.assertFalse(valid_units("m foo-2"))
        # Two numerical factors
        self.assertFalse(valid_units("10 10 m"))
        # Numerical factor not in initial position
        self.assertFalse(valid_units("m 10 s-1"))
        # Non-integer exponent
        self.assertFalse(valid_units("m2.5"))
        # Numerical factor neither integer nor decimal fraction
        self.assertFalse(valid_units("1.2e5 s"))

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(UnitStringTest))
    return s

if __name__ == '__main__':
    unittest.main()


