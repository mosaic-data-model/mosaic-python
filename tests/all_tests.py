# -*- coding: utf-8 -*-

# Run the complete test suite

#-----------------------------------------------------------------------------
#       Copyright (C) 2013-2014 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import unittest

import basic_tests
import mutable_model_tests
import immutable_model_tests
import array_model_tests
import equivalence_tests
import hdf5_tests
import h5md_tests
import xml_tests
import molecular_structure_tests
import pdb_tests

def suite():
    test_suite = unittest.TestSuite()
    test_suite.addTests(basic_tests.suite())
    test_suite.addTests(mutable_model_tests.suite())
    test_suite.addTests(immutable_model_tests.suite())
    test_suite.addTests(array_model_tests.suite())
    test_suite.addTests(equivalence_tests.suite())
    test_suite.addTests(hdf5_tests.suite())
    test_suite.addTests(h5md_tests.suite())
    test_suite.addTests(xml_tests.suite())
    test_suite.addTests(molecular_structure_tests.suite())
    test_suite.addTests(pdb_tests.suite())
    return test_suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())

