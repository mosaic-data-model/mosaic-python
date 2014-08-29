# -*- coding: utf-8 -*-

# Tests for module mosaic.array_model

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import unittest

import numpy as N

import mosaic.mutable_model as M
import mosaic.array_model as AM
from mosaic.api import is_valid

class UniverseTest(unittest.TestCase):

    def setUp(self):
        self.factory = AM.Factory()
        mol = M.Fragment("water", "water", (),
                         (M.Atom("H1", M.Element("H")),
                          M.Atom("H2", M.Element("H")),
                          M.Atom("O",  M.Element("O"), 2)),
                         (("H1", "O", "single"), ("H2", "O", "single")))
        self.universe = M.Universe('infinite', [(mol, 10)],
                                   convention='my_own')
        self.am_universe = self.factory(self.universe)

    def test_basics(self):
        self.assertTrue(self.am_universe.is_equivalent(self.universe))
        self.assertEqual(self.am_universe.number_of_molecules, 10)
        self.assertEqual(self.am_universe.number_of_atoms, 30)
        self.assertEqual(self.am_universe.number_of_sites, 40)
        self.assertEqual(self.am_universe.number_of_bonds, 20)
        self.assertEqual(self.am_universe.cell_shape, "infinite")
        self.assertEqual(self.am_universe.convention, "my_own")

    def test_return_types(self):
        self.assertTrue(type(self.am_universe.number_of_molecules) is int)
        self.assertTrue(type(self.am_universe.number_of_atoms) is int)
        self.assertTrue(type(self.am_universe.number_of_sites) is int)
        self.assertTrue(type(self.am_universe.number_of_bonds) is int)
        fragment = self.am_universe.molecules[0][0]
        self.assertTrue(type(fragment.number_of_atoms) is int)
        self.assertTrue(type(fragment.number_of_sites) is int)
        self.assertTrue(type(fragment.number_of_bonds) is int)
        
    def test_properties(self):
        mm_masses = M.TemplateAtomProperty(self.universe, "mass", "amu",
                                           N.array([1., 1., 16.], N.float32))
        masses = self.factory(mm_masses)
        self.assertTrue(masses.is_equivalent(mm_masses))
        self.assertEqual(masses.type, 'template_atom')
        self.assertTrue(masses.universe is self.am_universe)
        self.assertEqual(masses.element_shape, ())
        self.assertEqual(masses.data.shape, (3,))
        mm_masses = M.TemplateSiteProperty(self.universe, "mass", "amu",
                                           N.array([1., 1., 8., 8.], N.float32))
        bead_masses = self.factory(mm_masses)
        self.assertTrue(bead_masses.is_equivalent(mm_masses))
        self.assertEqual(bead_masses.type, 'template_site')
        self.assertTrue(bead_masses.universe is self.am_universe)
        self.assertEqual(bead_masses.element_shape, ())
        self.assertEqual(bead_masses.data.shape, (4,))
        mm_velocities = M.SiteProperty(self.universe, "velocity", "nm ps-1",
                                       dtype=N.float64, element_shape=(3,))
        velocities = self.factory(mm_velocities)
        self.assertTrue(velocities.is_equivalent(mm_velocities))
        self.assertEqual(velocities.type, 'site')
        self.assertTrue(velocities.universe is self.am_universe)
        self.assertEqual(velocities.data.shape, (40, 3))
        self.assertEqual(velocities.element_shape, (3,))
        mm_ap = M.AtomProperty(self.universe, "foo", "",
                               dtype=N.int16, element_shape=(2, 2))
        ap = self.factory(mm_ap)
        self.assertTrue(ap.is_equivalent(mm_ap))
        self.assertEqual(ap.type, 'atom')
        self.assertTrue(ap.universe is self.am_universe)
        self.assertEqual(ap.data.shape, (30, 2, 2))
        self.assertEqual(ap.element_shape, (2, 2))

    def test_labels(self):
        labels = tuple(a.name
                       for f, n in self.universe.molecules
                       for a in f.recursive_atom_iterator())
        mm_el = M.TemplateAtomLabel(self.universe, "element", labels)
        el = self.factory(mm_el)
        self.assertTrue(is_valid(el))
        self.assertTrue(el.is_equivalent(mm_el))
        self.assertEqual(el.name, "element")
        self.assertTrue(el.universe.is_equivalent(self.universe))
        self.assertTrue(len(el.strings)
                        == self.universe.number_of_template_atoms)
        for s1, s2 in zip(labels, el.strings):
            self.assertEqual(s1, s2)

        labels = tuple(a.name
                       for f, n in self.universe.molecules
                       for _ in range(n)
                       for a in f.recursive_atom_iterator())
        mm_el = M.AtomLabel(self.universe, "element", labels)
        el = self.factory(mm_el)
        self.assertTrue(is_valid(el))
        self.assertTrue(el.is_equivalent(mm_el))
        self.assertEqual(el.name, "element")
        self.assertTrue(el.universe.is_equivalent(self.universe))
        self.assertTrue(len(el.strings)
                        == self.universe.number_of_atoms)
        for s1, s2 in zip(labels, el.strings):
            self.assertEqual(s1, s2)

        labels = tuple(a.name
                       for f, n in self.universe.molecules
                       for a in f.recursive_atom_iterator()
                       for _ in range(a.number_of_sites))
        mm_el = M.TemplateSiteLabel(self.universe, "element", labels)
        el = self.factory(mm_el)
        self.assertTrue(is_valid(el))
        self.assertTrue(el.is_equivalent(mm_el))
        self.assertEqual(el.name, "element")
        self.assertTrue(el.universe.is_equivalent(self.universe))
        self.assertTrue(len(el.strings)
                        == self.universe.number_of_template_sites)
        for s1, s2 in zip(labels, el.strings):
            self.assertEqual(s1, s2)

        labels = tuple(a.name
                       for f, n in self.universe.molecules
                       for _ in range(n)
                       for a in f.recursive_atom_iterator()
                       for __ in range(a.number_of_sites))
        mm_el = M.SiteLabel(self.universe, "element", labels)
        el = self.factory(mm_el)
        self.assertTrue(is_valid(el))
        self.assertTrue(el.is_equivalent(mm_el))
        self.assertEqual(el.name, "element")
        self.assertTrue(el.universe.is_equivalent(self.universe))
        self.assertTrue(len(el.strings)
                        == self.universe.number_of_sites)
        for s1, s2 in zip(labels, el.strings):
            self.assertEqual(s1, s2)
        
    def test_bonds(self):
        bonds = self.universe.bond_index_array()
        self.assertEqual(len(bonds), self.universe.number_of_bonds)
        self.assertTrue((bonds >= 0).all())
        self.assertTrue((bonds < self.universe.number_of_atoms).all())
        for i in range(10):
            self.assertEqual(bonds[2*i, 0], 3*i)
            self.assertEqual(bonds[2*i, 1], 3*i+2)
            self.assertEqual(bonds[2*i+1, 0], 3*i+1)
            self.assertEqual(bonds[2*i+1, 1], 3*i+2)

    def test_polymers(self):
        peptide_chain = M.Polymer('chain', 'ALA,LYS',
                                  (M.Fragment('ALA1', 'ALA', (),
                                              (M.Atom('CA',
                                                      M.CGParticle('CA'), 1),),
                                              ()),
                                   M.Fragment('LYS2', 'LYS', (),
                                              (M.Atom('CA',
                                                      M.CGParticle('CA'), 1),),
                                              ())),
                                  (('ALA1.CA', 'LYS2.CA', ''),),
                                  'polypeptide')
        protein = M.Fragment('protein', 'protein', (peptide_chain,), (), ())
        universe = M.Universe('infinite', [(protein, 1)])
        array_universe = self.factory(universe)
        array_peptide_chain = array_universe.molecules[0][0].fragments[0]
        self.assert_(array_peptide_chain.is_polymer)
        self.assertEqual(array_peptide_chain.polymer_type, 'polypeptide')

    def test_selections(self):
        universe = self.am_universe
        s = AM.Selection(universe, N.array([0, 2], N.uint8), 'atom')
        self.assertEqual(s.number_of_atoms, 2)
        self.assertEqual(s.number_of_sites, 3)
        s = AM.Selection(universe, N.array([1], N.uint8), 'template_atom')
        self.assertEqual(s.number_of_atoms, 10)
        self.assertEqual(s.number_of_sites, 10)
        s = AM.Selection(universe, N.array([2], N.uint8), 'template_atom')
        self.assertEqual(s.number_of_atoms, 10)
        self.assertEqual(s.number_of_sites, 20)
        s = AM.Selection(universe, N.array([2], N.uint8), 'site')
        self.assertEqual(s.number_of_sites, 1)
        s = AM.Selection(universe, N.array([0], N.uint8), 'template_site')
        self.assertEqual(s.number_of_sites, 10)

    def test_optimizations(self):
        import mosaic.api
        universe = self.am_universe
        for a1, a2 in zip(universe.recursive_atom_iterator(),
                          mosaic.api.MosaicUniverse.recursive_atom_iterator(universe)):
            self.assertTrue(a1.is_equivalent(a2))
        bonds1 = universe.bond_index_array()
        bonds2 = mosaic.api.MosaicUniverse.bond_index_array(universe)
        self.assertTrue((bonds1 == bonds2).all())

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(UniverseTest))
    return s

if __name__ == '__main__':
    unittest.main()
