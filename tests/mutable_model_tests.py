# -*- coding: utf-8 -*-

# Tests for module mosaic.mutable_model

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import copy
import unittest

import numpy as N

import mosaic.mutable_model as M
from mosaic.api import is_valid

def make_water_fragment(nsites=1):
    return M.Fragment("solvent", "water", (),
                      (M.Atom("H1", M.Element("H"), nsites),
                       M.Atom("H2", M.Element("H"), nsites),
                       M.Atom("O",  M.Element("O"), nsites)),
                      (("H1", "O", "single"), ("H2", "O", "single")))

class AtomDescriptorTest(unittest.TestCase):

    def test_singleton(self):
        self.assertTrue(M.Dummy() is M.Dummy())
        self.assertTrue(M.Dummy('a') is M.Dummy('a'))
        self.assertTrue(M.Dummy('a') is not M.Dummy('b'))
        self.assertTrue(M.Dummy('a') is not M.Unknown('a'))
        self.assertTrue(M.Dummy('C') is not M.Element('C'))
        self.assertTrue(M.Dummy('D') is not M.CGParticle('D'))
        self.assertTrue(M.Element('C') is M.Element('C'))

    def test_name(self):
        for name in ['a', 'b', 'c']:
            self.assertEqual(M.Unknown(name).name, name)

    def test_type(self):
        self.assertEqual(M.Dummy().type, "dummy")
        self.assertEqual(M.Unknown().type, "")
        self.assertEqual(M.Element('O').type, "element")
        self.assertEqual(M.CGParticle('ala').type, "cgparticle")

class WaterTest(unittest.TestCase):

    def setUp(self):
        self.mol = make_water_fragment()

    def test_basics(self):
        self.assertEqual(self.mol.number_of_atoms, 3)
        self.assertEqual(self.mol.number_of_sites, 3)
        self.assertEqual(self.mol.number_of_bonds, 2)
        atoms = [a for a in self.mol.recursive_atom_iterator()]
        self.assertEqual(self.mol.number_of_atoms, len(atoms))
        self.assertTrue(self.mol['H1'] is atoms[0])
        self.assertTrue(self.mol['H2'] is atoms[1])
        self.assertTrue(self.mol['O'] is atoms[2])
        for a in atoms:
            self.assertTrue(a.parent is self.mol)
        self.assertEqual(self.mol.label, "solvent")
        self.assertEqual(self.mol.species, "water")
        
    def test_copy(self):
        mol_copy = copy.copy(self.mol)
        self.assertEqual(mol_copy.number_of_atoms, 3)
        atoms = [a for a in self.mol.recursive_atom_iterator()]
        copy_atoms = [a for a in mol_copy.recursive_atom_iterator()]
        self.assertEqual(len(copy_atoms), len(atoms))
        self.assertTrue(mol_copy['H1'] is copy_atoms[0])
        self.assertTrue(mol_copy['H2'] is copy_atoms[1])
        self.assertTrue(mol_copy['O'] is copy_atoms[2])
        for name in ['H1', 'H2', 'O']:
            self.assertFalse(self.mol[name] is mol_copy[name])
        for a in copy_atoms:
            self.assertTrue(a.parent is mol_copy)

    def test_equality(self):
        mol_copy = copy.copy(self.mol)
        equal_mol = make_water_fragment()
        changed_bond_order = M.Fragment("solvent", "water", (),
                                        (M.Atom("H1", M.Element("H")),
                                         M.Atom("H2", M.Element("H")),
                                         M.Atom("O",  M.Element("O"))),
                                        (("O", "H2", "single"),
                                         ("O", "H1", "single")))
        changed_atom_order = M.Fragment("solvent", "water", (),
                                        (M.Atom("O",  M.Element("O")),
                                         M.Atom("H1", M.Element("H")),
                                         M.Atom("H2", M.Element("H"))),
                                        (("O", "H1", "single"),
                                         ("O", "H2", "single")))
        self.assertEqual(self.mol, self.mol)
        self.assertTrue(self.mol.is_equivalent(self.mol))
        self.assertEqual(self.mol, mol_copy)
        self.assertTrue(self.mol.is_equivalent(mol_copy))
        self.assertEqual(self.mol, equal_mol)
        self.assertEqual(self.mol, equal_mol)
        self.assertTrue(self.mol.is_equivalent(equal_mol))
        self.assertEqual(self.mol, changed_bond_order)
        self.assertTrue(self.mol.is_equivalent(changed_bond_order))
        self.assertNotEqual(self.mol, changed_atom_order)
        self.assertFalse(self.mol.is_equivalent(changed_atom_order))

class PeptideTest(unittest.TestCase):

    def _make_molecule(self, label):
        C = M.Element('C')
        H = M.Element('H')
        N = M.Element('N')
        O = M.Element('O')
        peptide_group = M.Fragment('peptide', 'peptide',
                                   (),
                                   (M.Atom('CA', C),
                                    M.Atom('HA', H),
                                    M.Atom('H', H),
                                    M.Atom('N', N),
                                    M.Atom('C', C),
                                    M.Atom('O', O)),
                                   (('N', 'H', "single"),
                                    ('N', 'CA', "single"),
                                    ('CA', 'HA', "single"),
                                    ('CA', 'C', "single"),
                                    ('C', 'O', "double")))
        ala_sidechain = M.Fragment('sidechain', 'ala_sidechain',
                                   (),
                                   (M.Atom('CB', C),
                                    M.Atom('HB1', H),
                                    M.Atom('HB2', H),
                                    M.Atom('HB3', H)),
                                   (('CB', 'HB1', "single"),
                                    ('CB', 'HB2', "single"),
                                    ('CB', 'HB3', "single"),))
        ala = lambda label: M.Fragment(label, 'ala',
                                       (copy.copy(peptide_group),
                                        copy.copy(ala_sidechain)),
                                       (),
                                       (('peptide.CA', 'sidechain.CB',
                                         "single"),))
        return M.Polymer(label, 'alanine_dipeptide',
                         (ala('ALA1'), ala('ALA2')),
                         (('ALA1.peptide.C', 'ALA2.peptide.N', "single"),),
                         'polypeptide')

    def test_basic(self):
        mol = self._make_molecule('di_ala')
        self.assertTrue(is_valid(mol))
        self.assertEqual(mol.number_of_atoms, 20)
        self.assertEqual(mol.number_of_sites, 20)
        self.assertEqual(mol.number_of_bonds, 19)
        atoms = [a for a in mol.recursive_atom_iterator()]
        self.assertEqual(mol.number_of_atoms, len(atoms))
        self.assertEqual(mol.polymer_type, "polypeptide")

    def test_equality(self):
        self.assertEqual(self._make_molecule('di_ala'),
                         self._make_molecule('di_ala'))
        self.assertTrue(self._make_molecule('di_ala')
                        .is_equivalent(self._make_molecule('di_ala')))
        self.assertNotEqual(self._make_molecule('di_ala_1'),
                            self._make_molecule('di_ala_2'))
        self.assertFalse(self._make_molecule('di_ala_1')
                        .is_equivalent(self._make_molecule('di_ala_2')))

    def test_iterators(self):
        mol = self._make_molecule('di_ala')
        atoms = tuple(mol.recursive_atom_iterator())
        self.assertEqual(len(atoms), mol.number_of_atoms)
        bonds = tuple(mol.recursive_bond_iterator())
        self.assertEqual(len(bonds), mol.number_of_bonds)
        for a1, a2, order in bonds:
            for a in a1, a2:
                node = mol
                for p in a.split('.'):
                    node = node[p]
                self.assertTrue(isinstance(node, M.Atom))
        paths = tuple(mol.recursive_atom_path_iterator())
        self.assertEqual(len(paths), mol.number_of_atoms)
        for ap in paths:
            node = mol
            for p in ap.split('.'):
                node = node[p]
            self.assertTrue(isinstance(node, M.Atom))

class ErrorCheckingTest(unittest.TestCase):

    def test_atom_descriptor(self):
        self.assertRaises(TypeError, lambda: M.Dummy(42))
        self.assertRaises(TypeError, lambda: M.Element(42))
        self.assertRaises(ValueError, lambda: M.Element("X"))
        
    def test_atom(self):
        carbon = M.Element("C")
        self.assertRaises(TypeError, lambda: M.Atom(42, carbon, 1))
        self.assertRaises(ValueError, lambda: M.Atom('2.', carbon, 1))
        self.assertRaises(TypeError, lambda: M.Atom('some_label', 'C', 1))
        self.assertRaises(ValueError, lambda: M.Atom('some_label', carbon, 0))
        
    def test_fragment(self):
        carbon = M.Element("C")
        # Illegal labels
        atom = M.Atom("C", carbon)
        self.assertRaises(TypeError,
                          lambda: M.Fragment(42, '', (), (atom,), ()))
        atom = M.Atom("C", carbon)
        self.assertRaises(ValueError,
                          lambda: M.Fragment('2.', '', (), (atom,), ()))
        # Illegal fragments
        atom = M.Atom("C", carbon)
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', None, (atom,), ()))
        atom = M.Atom("C", carbon)
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', [1, 2], (atom,), ()))
        atom = M.Atom("C", carbon)
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', (atom,), (), ()))
        # Illegal atoms
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', (), None, ()))
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', (), [1, 2], ()))
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', (), (carbon,), ()))
        self.assertRaises(ValueError,
                          lambda: M.Fragment('m', 'm', (),
                                             (M.Atom("C", carbon),
                                              M.Atom("C", carbon)),
                                             ()))
        # Illegal bond lists
        atom = M.Atom("C", carbon)
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', (), (atom,), None))
        atom = M.Atom("C", carbon)
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', (), (atom,), [1, 2, 3]))
        atom = M.Atom("C", carbon)
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', (), (atom,),
                                             (('X', 'X'))))
        atom = M.Atom("C", carbon)
        self.assertRaises(TypeError,
                          lambda: M.Fragment('m', 'm', (), (atom,),
                                             (['X', 'X', 'single'])))
        
    def test_bonds(self):
        carbon = lambda label: M.Atom(label, M.Element("C"))
        # Bond specified by only one atom
        self.assertRaises(ValueError,
                          lambda: M.Fragment('m', 'm', (),
                                             (carbon('C1'), carbon('C2')),
                                             (('C1', ),)))
        # Bond specified by two atoms but no bond order
        self.assertRaises(ValueError,
                          lambda: M.Fragment('m', 'm', (),
                                             (carbon('C1'), carbon('C2')),
                                             (('C1', 'C2'),)))
        # Bond specified by two identical atoms
        self.assertRaises(ValueError,
                          lambda: M.Fragment('m', 'm', (),
                                             (carbon('C1'), carbon('C2')),
                                             (('C1', 'C1', ''),)))
        # Bond specified by an atom outside the fragment
        self.assertRaises(ValueError,
                          lambda: M.Fragment('m', 'm', (),
                                             (carbon('C1'), carbon('C2')),
                                             (('C1', 'C3', ''),)))
        # Bond specified by an atom name that is undefined
        self.assertRaises(ValueError,
                          lambda: M.Fragment('m', 'm', (),
                                             (carbon('C1'), carbon('C2')),
                                             (('C1', 'C3', ''),)))
        # Bond specified at the wrong fragment level
        f = M.Fragment('x', 'x', (), (carbon('C1'), carbon('C2')), ())
        self.assertRaises(ValueError,
                          lambda: M.Fragment('m', 'm', (f,),
                                             (carbon('C3'),),
                                             (('x.C1', 'x.C2', ''),)))

    def test_universe(self):
        mol = make_water_fragment()
        self.assertRaises(TypeError,
                          lambda: M.Universe(0, [(mol, 10)]))
        self.assertRaises(ValueError,
                          lambda: M.Universe('strange', [(mol, 10)]))
        self.assertRaises(TypeError,
                          lambda: M.Universe('infinite', mol))
        self.assertRaises(ValueError,
                          lambda: M.Universe('infinite', [("water", 10)]))
        self.assertRaises(TypeError,
                          lambda: M.Universe('infinite', [mol]))
        self.assertRaises(ValueError,
                          lambda: M.Universe('infinite', [(10, mol)]))
        self.assertRaises(ValueError,
                          lambda: M.Universe('infinite', [(mol, 10)],
                                             [(N.zeros((3,3), N.float64),
                                               N.zeros((3,), N.float64))]))

    def test_configuration(self):
        mol = make_water_fragment()
        universe = M.Universe('cube', [(mol, 10)])
        # Wrong dtype
        self.assertRaises(ValueError,
                          lambda: M.Configuration(universe, dtype=N.int))
        # Positions but no cell parameters
        self.assertRaises(ValueError,
                          lambda: M.Configuration(universe,
                                                  N.zeros((30, 3), N.float32)))
        # Positions and cell parameters of different dtype
        self.assertRaises(ValueError,
                          lambda: M.Configuration(universe,
                                                  N.zeros((30, 3), N.float32),
                                                  N.float64(10.)))
        # Positions not an array
        self.assertRaises(TypeError,
                          lambda: M.Configuration(universe,
                                                  list(N.zeros((30, 3),
                                                               N.float32)),
                                                  N.float32(10.)))
        # Positions of wrong shape
        self.assertRaises(ValueError,
                          lambda: M.Configuration(universe,
                                                  N.zeros((25, 3), N.float32),
                                                  N.float32(10.)))
        # Cell parameters of wrong shape
        self.assertRaises(ValueError,
                          lambda: M.Configuration(universe,
                                                  N.zeros((30, 3), N.float32),
                                                  N.zeros((3,), N.float32)))
        # Different explicit dtype
        self.assertRaises(ValueError,
                          lambda: M.Configuration(universe,
                                                  N.zeros((30, 3), N.float32),
                                                  N.float32(10.),
                                                  N.float64))
        
    def test_selection(self):
        mol = make_water_fragment(nsites=2)
        universe = M.Universe('cube', [(mol, 5)])
        # Index occurs twice
        self.assertRaises(ValueError,
                          lambda: M.AtomSelection(universe,
                                                  N.zeros((2,), N.uint16)))
        # Atom index too large
        self.assertRaises(ValueError,
                          lambda: M.AtomSelection(universe,
                                                  N.array([20], N.uint16)))
        # Template atom index too large
        self.assertRaises(ValueError,
                          lambda: M.TemplateAtomSelection(universe,
                                                          N.array([3], N.uint8)))
        # Site index too large
        self.assertRaises(ValueError,
                          lambda: M.SiteSelection(universe,
                                                  N.array([40], N.uint16)))
        # Template site index too large
        self.assertRaises(ValueError,
                          lambda: M.TemplateSiteSelection(universe,
                                                          N.array([8], N.uint8)))

class UniverseTest(unittest.TestCase):

    def setUp(self):
        mol = make_water_fragment(2)
        self.universe = M.Universe('infinite', [(mol, 10)],
                                   convention='my_own')

    def test_basics(self):
        self.assertTrue(is_valid(self.universe))
        self.assertEqual(self.universe.number_of_molecules, 10)
        self.assertEqual(self.universe.number_of_atoms, 30)
        self.assertEqual(self.universe.number_of_sites, 60)
        self.assertEqual(self.universe.number_of_bonds, 20)
        self.assertEqual(self.universe.cell_shape, "infinite")
        self.assertEqual(self.universe.convention, "my_own")

    def test_properties(self):
        masses = M.TemplateAtomProperty(self.universe, "mass", "amu",
                                        N.array([1., 1., 16.], N.float32))
        self.assertTrue(is_valid(masses))
        self.assertEqual(masses.type, 'template_atom')
        self.assertTrue(masses.universe is self.universe)
        self.assertEqual(masses.element_shape, ())
        self.assertEqual(masses.data.shape, (3,))
        bead_masses = M.TemplateSiteProperty(self.universe, "mass", "amu",
                                             N.array([1., 1.,
                                                      1., 1.,
                                                      8., 8.], N.float32))
        self.assertTrue(is_valid(bead_masses))
        self.assertEqual(bead_masses.type, 'template_site')
        self.assertTrue(bead_masses.universe is self.universe)
        self.assertEqual(bead_masses.element_shape, ())
        self.assertEqual(bead_masses.data.shape, (6,))
        velocities = M.SiteProperty(self.universe, "velocity", "nm ps-1",
                                    dtype=N.float64, element_shape=(3,))
        self.assertTrue(is_valid(velocities))
        self.assertEqual(velocities.type, 'site')
        self.assertTrue(velocities.universe is self.universe)
        self.assertEqual(velocities.data.shape, (60, 3))
        self.assertEqual(velocities.element_shape, (3,))
        foo = M.AtomProperty(self.universe, "foo", "",
                             dtype=N.int16, element_shape=(2, 2))
        self.assertTrue(is_valid(foo))
        self.assertEqual(foo.type, 'atom')
        self.assertTrue(foo.universe is self.universe)
        self.assertEqual(foo.data.shape, (30, 2, 2))
        self.assertEqual(foo.element_shape, (2, 2))

    def test_labels(self):
        labels = tuple(a.name
                       for f, n in self.universe.molecules
                       for a in f.recursive_atom_iterator())
        el = M.TemplateAtomLabel(self.universe, "element", labels)
        self.assertTrue(is_valid(el))
        self.assertEqual(el.name, "element")
        self.assertTrue(el.universe == self.universe)
        self.assertTrue(len(el.strings)
                        == self.universe.number_of_template_atoms)
        for s1, s2 in zip(labels, el.strings):
            self.assertEqual(s1, s2)

        labels = tuple(a.name
                       for f, n in self.universe.molecules
                       for _ in range(n)
                       for a in f.recursive_atom_iterator())
        el = M.AtomLabel(self.universe, "element", labels)
        self.assertTrue(is_valid(el))
        self.assertEqual(el.name, "element")
        self.assertTrue(el.universe == self.universe)
        self.assertTrue(len(el.strings)
                        == self.universe.number_of_atoms)
        for s1, s2 in zip(labels, el.strings):
            self.assertEqual(s1, s2)

        labels = tuple(a.name
                       for f, n in self.universe.molecules
                       for a in f.recursive_atom_iterator()
                       for _ in range(a.number_of_sites))
        el = M.TemplateSiteLabel(self.universe, "element", labels)
        self.assertTrue(is_valid(el))
        self.assertEqual(el.name, "element")
        self.assertTrue(el.universe == self.universe)
        self.assertTrue(len(el.strings)
                        == self.universe.number_of_template_sites)
        for s1, s2 in zip(labels, el.strings):
            self.assertEqual(s1, s2)

        labels = tuple(a.name
                       for f, n in self.universe.molecules
                       for _ in range(n)
                       for a in f.recursive_atom_iterator()
                       for __ in range(a.number_of_sites))
        el = M.SiteLabel(self.universe, "element", labels)
        self.assertTrue(is_valid(el))
        self.assertEqual(el.name, "element")
        self.assertTrue(el.universe == self.universe)
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

    def test_index_mappings(self):
        mol = self.universe.molecules[0][0]
        s2a = mol.site_to_atom_index_mapping()
        self.assertTrue((s2a == N.array([0, 0, 1, 1, 2, 2])).all())

        s2a = self.universe.site_to_atom_index_mapping()
        s2a_ref = N.repeat(N.arange(30), 2)
        self.assertTrue((s2a == s2a_ref).all())

        st2at = self.universe.template_site_to_template_atom_index_mapping()
        st2at_ref = N.array([0, 0, 1, 1, 2, 2])
        self.assertTrue((st2at == st2at_ref).all())

        s2t = self.universe.site_to_template_index_mapping()
        s2t_ref = N.resize(N.arange(mol.number_of_sites),
                           (self.universe.number_of_sites,))
        self.assertTrue((s2t == s2t_ref).all())

        a2t = self.universe.atom_to_template_index_mapping()
        a2t_ref = N.resize(N.arange(mol.number_of_atoms),
                           (self.universe.number_of_atoms,))
        self.assertTrue((a2t == a2t_ref).all())

        a2s = self.universe.atom_to_site_index_mapping()
        a2s_ref = 2*N.arange(self.universe.number_of_atoms)
        self.assertTrue((a2s == a2s_ref).all())

    def test_selections(self):
        s = M.AtomSelection(self.universe, N.array([0, 2], N.uint8))
        self.assertEqual(s.number_of_atoms, 2)
        self.assertEqual(s.number_of_sites, 4)
        s = M.TemplateAtomSelection(self.universe, N.array([1], N.uint8))
        self.assertEqual(s.number_of_atoms, 10)
        self.assertEqual(s.number_of_sites, 20)
        s = M.SiteSelection(self.universe, [2])
        self.assertEqual(s.number_of_sites, 1)
        s = M.TemplateSiteSelection(self.universe, [0])
        self.assertEqual(s.number_of_sites, 10)

class PBCTest(unittest.TestCase):

    def setUp(self):
        self.infinite = M.Universe('infinite', [])
        self.cube = M.Universe('cube', [])
        self.cuboid = M.Universe('cuboid', [])
        self.parallelepiped = M.Universe('parallelepiped', [])
    
    def test_lattice_vectors(self):
        conf = M.Configuration(self.infinite,
                               N.zeros((0, 3), N.float32))
        self.assertEqual(conf.lattice_vectors(), ())
        self.assertEqual(conf.cell_volume(), None)
        conf = M.Configuration(self.cube,
                               N.zeros((0, 3), N.float32),
                               N.array(1., N.float32))
        lv = conf.lattice_vectors()
        self.assertTrue((lv[0] == N.array([1., 0., 0.])).all())
        self.assertTrue((lv[1] == N.array([0., 1., 0.])).all())
        self.assertTrue((lv[2] == N.array([0., 0., 1.])).all())
        self.assertEqual(conf.cell_volume(), 1.)
        conf = M.Configuration(self.cuboid,
                               N.zeros((0, 3), N.float32),
                               N.array([1., 2., 4.], N.float32))
        lv = conf.lattice_vectors()
        self.assertTrue((lv[0] == N.array([1., 0., 0.])).all())
        self.assertTrue((lv[1] == N.array([0., 2., 0.])).all())
        self.assertTrue((lv[2] == N.array([0., 0., 4.])).all())
        self.assertEqual(conf.cell_volume(), 8.)
        conf = M.Configuration(self.parallelepiped,
                               N.zeros((0, 3), N.float32),
                               N.array([[1., 2., 4.],
                                        [8., 4., 2.],
                                        [16., 4., 8.]], N.float32))
        lv = conf.lattice_vectors()
        self.assertTrue((lv[0] == N.array([1., 2., 4.])).all())
        self.assertTrue((lv[1] == N.array([8., 4., 2.])).all())
        self.assertTrue((lv[2] == N.array([16., 4., 8.])).all())
        self.assertAlmostEqual(conf.cell_volume(), 168.)
        
def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(AtomDescriptorTest))
    s.addTest(loader.loadTestsFromTestCase(WaterTest))
    s.addTest(loader.loadTestsFromTestCase(PeptideTest))
    s.addTest(loader.loadTestsFromTestCase(ErrorCheckingTest))
    s.addTest(loader.loadTestsFromTestCase(UniverseTest))
    s.addTest(loader.loadTestsFromTestCase(PBCTest))
    return s

if __name__ == '__main__':
    unittest.main()
