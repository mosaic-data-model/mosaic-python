# -*- coding: utf-8 -*-

# Tests for module mosaic.h5md

#-----------------------------------------------------------------------------
#       Copyright (C) 2014 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import copy
import os
import tempfile
import unittest

import numpy as N
import numpy.random as NR
import immutable.np as IN

import mosaic.immutable_model as M

import h5py

from mosaic.h5md import Trajectory

class TestHDF5File(object):

    @classmethod
    def initialize(cls):
        cls.tempdir = tempfile.mkdtemp()
        cls.filename = os.path.join(cls.tempdir, 'test.h5')

    @classmethod
    def cleanup(cls):
        os.remove(cls.filename)
        os.removedirs(cls.tempdir)

    def __init__(self, mode='r'):
        self.file = h5py.File(self.filename, mode)

    def close(self):
        self.file.close()

    def __enter__(self):
        return self.file

    def __exit__(self, type, value, traceback):
        self.close()
        return False

water = M.fragment("water", (),
                   (("H1", M.atom(M.element("H"), 1)),
                    ("H2", M.atom(M.element("H"), 1)),
                    ("O",  M.atom(M.element("O"), 1))),
                   (("H1", "O", "single"), ("H2", "O", "single")))

class WriteReadTest(unittest.TestCase):

    def setUp(self):
        TestHDF5File.initialize()

    def tearDown(self):
        TestHDF5File.cleanup()

    def test_infinite(self):
        universe = M.universe('infinite', [(water, 'water', 5)])
        nsites = universe.number_of_sites
        conf1 = M.Configuration(universe,
                                IN.array(NR.normal(size=(nsites, 3))),
                                None)
        conf2 = M.Configuration(universe,
                                IN.array(NR.normal(size=(nsites, 3))),
                                None)
        self._write_and_read_back(universe, conf1, conf2)

    def test_cubic(self):
        universe = M.universe('cube', [(water, 'water', 10)])
        nsites = universe.number_of_sites
        conf1 = M.Configuration(universe,
                                IN.array(NR.normal(size=(nsites, 3))),
                                IN.array(10., N.float64))
        conf2 = M.Configuration(universe,
                                IN.array(NR.normal(size=(nsites, 3))),
                                IN.array(11., N.float64))
        self._write_and_read_back(universe, conf1, conf2)

    def test_cuboid(self):
        universe = M.universe('cuboid', [(water, 'water', 1)])
        nsites = universe.number_of_sites
        conf1 = M.Configuration(universe,
                                IN.array(NR.normal(size=(nsites, 3))),
                                IN.array(NR.normal(size=(3, ))))
        conf2 = M.Configuration(universe,
                                IN.array(NR.normal(size=(nsites, 3))),
                                IN.array(NR.normal(size=(3, ))))
        self._write_and_read_back(universe, conf1, conf2)

    def test_parallelepiped(self):
        universe = M.universe('parallelepiped', [(water, 'water', 100)])
        nsites = universe.number_of_sites
        conf1 = M.Configuration(universe,
                                IN.array(NR.normal(size=(nsites, 3))),
                                IN.array(NR.normal(size=(3, 3))))
        conf2 = M.Configuration(universe,
                                IN.array(NR.normal(size=(nsites, 3))),
                                IN.array(NR.normal(size=(3, 3))))
        self._write_and_read_back(universe, conf1, conf2)

    def _write_and_read_back(self, universe, conf1, conf2):
        with TestHDF5File('w') as h5file:
            tr = Trajectory(h5file, universe=universe, float_type=N.float64)
            tr.write_step(0, 1., conf1)
            tr.write_step(1, 2., conf2)
            tr.close()
        with TestHDF5File('r') as h5file:
            tr = Trajectory(h5file)
            r_step1, r_time1, r_conf1 = tr.read_step(0)
            r_step2, r_time2, r_conf2 = tr.read_step(1)
            tr.close()
        self.assertEqual(r_step1, 0)
        self.assertEqual(r_time1, 1.)
        self.assertTrue(conf1.is_equivalent(r_conf1))
        self.assertEqual(r_step2, 1)
        self.assertEqual(r_time2, 2.)
        self.assertTrue(conf2.is_equivalent(r_conf2))


def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(WriteReadTest))
    return s

if __name__ == '__main__':
    unittest.main()
