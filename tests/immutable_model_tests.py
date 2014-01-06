# -*- coding: utf-8 -*-

# Tests for module mosaic.immutable_model

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import unittest

import numpy as N
import immutable.np as IN

import mosaic.immutable_model as M
from mosaic.api import is_valid

def make_water_fragment(nsites=1):
    return M.fragment("water", (),
                      (("H1", M.atom(M.element("H"), nsites)),
                       ("H2", M.atom(M.element("H"), nsites)),
                       ("O",  M.atom(M.element("O"), nsites))),
                      (("H1", "O", "single"), ("H2", "O", "single")))


class AtomDescriptorTest(unittest.TestCase):

    def test_singleton(self):
        self.assertTrue(M.dummy() is M.dummy())
        self.assertTrue(M.dummy('a') is M.dummy('a'))
        self.assertTrue(M.dummy('a') is not M.dummy('b'))
        self.assertTrue(M.dummy('a') is not M.unknown('a'))
        self.assertTrue(M.dummy('C') is not M.element('C'))
        self.assertTrue(M.element('C') is M.element('C'))

    def test_name(self):
        for name in ['a', 'b', 'c']:
            self.assertEqual(M.unknown(name).name, name)

    def test_type(self):
        self.assertEqual(M.dummy().type, "dummy")
        self.assertEqual(M.unknown().type, "")
        self.assertEqual(M.element('O').type, "element")
        self.assertEqual(M.cgparticle('ala').type, "cgparticle")

class WaterTest(unittest.TestCase):

    def setUp(self):
        self.mol = make_water_fragment()

    def test_basics(self):
        self.assertEqual(self.mol.number_of_atoms, 3)
        self.assertEqual(self.mol.number_of_sites, 3)
        self.assertEqual(self.mol.number_of_bonds, 2)
        self.assertEqual(self.mol.species, "water")

    def test_equality(self):
        same_mol = make_water_fragment()
        changed_bond_order = M.fragment("water", (),
                                        (("H1", M.atom(M.element("H"))),
                                         ("H2", M.atom(M.element("H"))),
                                         ("O",  M.atom(M.element("O")))),
                                        (("O", "H2", "single"),
                                         ("O", "H1", "single")))
        changed_atom_order = M.fragment("water", (),
                                        (("O",  M.atom(M.element("O"))),
                                         ("H1", M.atom(M.element("H"))),
                                         ("H2", M.atom(M.element("H")))),
                                        (("O", "H1", "single"),
                                         ("O", "H2", "single")))
        self.assertEqual(self.mol, self.mol)
        self.assertEqual(self.mol, same_mol)
        self.assertEqual(self.mol, changed_bond_order)
        self.assertNotEqual(self.mol, changed_atom_order)
        
class PeptideTest(unittest.TestCase):

    def _make_molecule(self):
        C = M.element('C')
        H = M.element('H')
        N = M.element('N')
        O = M.element('O')
        peptide_group = M.fragment('peptide',
                                   (),
                                   (('CA', M.atom(C)),
                                    ('HA', M.atom(H)),
                                    ('H', M.atom(H)),
                                    ('N', M.atom(N)),
                                    ('C', M.atom(C)),
                                    ('O', M.atom(O))),
                                   (('N', 'H', "single"),
                                    ('N', 'CA', "single"),
                                    ('CA', 'HA', "single"),
                                    ('CA', 'C', "single"),
                                    ('C', 'O', "double")))
        ala_sidechain = M.fragment('ala_sidechain',
                                   (),
                                   (('CB', M.atom(C)),
                                    ('HB1', M.atom(H)),
                                    ('HB2', M.atom(H)),
                                    ('HB3', M.atom(H))),
                                   (('CB', 'HB1', "single"),
                                    ('CB', 'HB2', "single"),
                                    ('CB', 'HB3', "single"),))
        ala = M.fragment('alanine',
                         (('peptide', peptide_group),
                          ('sidechain', ala_sidechain)),
                         (),
                         (('peptide.CA', 'sidechain.CB', "single"),))
        return M.polymer('alanine_dipeptide',
                         (('ALA1', ala),
                          ('ALA2', ala)),
                         (('ALA1.peptide.C', 'ALA2.peptide.N', "single"),),
                         'polypeptide')

    def test_basic(self):
        mol = self._make_molecule()
        self.assertEqual(mol.number_of_atoms, 20)
        self.assertEqual(mol.number_of_sites, 20)
        self.assertEqual(mol.number_of_bonds, 19)
        self.assertEqual(mol.polymer_type, "polypeptide")

    def test_equality(self):
        self.assertEqual(self._make_molecule(),
                         self._make_molecule())

    def test_iterators(self):
        mol = self._make_molecule()
        mol_ref = M.FragmentRef('x', mol)
        atoms = tuple(mol_ref.recursive_atom_iterator())
        self.assertEqual(len(atoms), mol.number_of_atoms)
        bonds = tuple(mol_ref.recursive_bond_iterator())
        self.assertEqual(len(bonds), mol.number_of_bonds)
        for a1, a2, order in bonds:
            for a in a1, a2:
                node = mol
                for p in a.split('.'):
                    node = node[p]
                self.assertTrue(isinstance(node, M.Atom))
        paths = tuple(mol_ref.recursive_atom_path_iterator())
        self.assertEqual(len(paths), mol.number_of_atoms)
        for ap in paths:
            node = mol
            for p in ap.split('.'):
                node = node[p]
            self.assertTrue(isinstance(node, M.Atom))

class ErrorCheckingTest(unittest.TestCase):

    def test_atom_descriptor(self):
        self.assertRaises(TypeError, lambda: M.dummy(42))
        self.assertRaises(ValueError, lambda: M.element(42))
        self.assertRaises(ValueError, lambda: M.element("X"))
        
    def test_atom(self):
        carbon = M.element("C")
        self.assertRaises(TypeError, lambda: M.atom('C', 1))
        self.assertRaises(ValueError, lambda: M.atom(carbon, 0))
        
    def test_fragment(self):
        carbon = M.atom(M.element("C"))
        # Illegal fragments
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', None, (("C", carbon),), ()))
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', [1, 2], (("C", carbon),), ()))
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', (("C", carbon),), (), ()))
        # Illegal atoms
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', (), None, ()))
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', (), [1, 2], ()))
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', (), (carbon,), ()))
        self.assertRaises(ValueError,
                          lambda: M.fragment('m', (),
                                             (("C", carbon),
                                              ("C", carbon)),
                                             ()))
        # Illegal bond lists
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', (), (("C", carbon),), None))
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', (), (("C", carbon),),
                                             [1, 2, 3]))
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', (), (("C", carbon),),
                                             (('X', 'X'))))
        self.assertRaises(TypeError,
                          lambda: M.fragment('m', (), (("C", carbon),),
                                             (['X', 'X', 'single'])))
        
    def test_bonds(self):
        carbon = M.atom(M.element("C"))
        # Bond specified by only one atom
        self.assertRaises(ValueError,
                          lambda: M.fragment('m', (),
                                             (('C1', carbon), ('C2', carbon)),
                                             (('C1', ),)))
        # Bond specified by two atoms but no bond order
        self.assertRaises(ValueError,
                          lambda: M.fragment('m', (),
                                             (('C1', carbon), ('C2', carbon)),
                                             (('C1', 'C2'),)))
        # Bond specified by two identical atoms
        self.assertRaises(ValueError,
                          lambda: M.fragment('m', (),
                                             (('C1', carbon), ('C2', carbon)),
                                             (('C1', 'C1', ''),)))
        # Bond specified by an atom name that is undefined
        self.assertRaises(ValueError,
                          lambda: M.fragment('m', (),
                                             (('C1', carbon), ('C2', carbon)),
                                             (('C1', 'C3', ''),)))
        # Bond specified at the wrong fragment level
        f = M.fragment('x', (), (('C1', carbon), ('C2', carbon)), ())
        self.assertRaises(ValueError,
                          lambda: M.fragment('m', (('x', f),),
                                             (('C3', carbon),),
                                             (('x.C1', 'x.C2', ''),)))

    def test_universe(self):
        mol = M.fragment("water", (),
                         (("H1", M.atom(M.element("H"), 8)),
                          ("H2", M.atom(M.element("H"), 8)),
                          ("O",  M.atom(M.element("O"), 2))),
                         (("H1", "O", "single"), ("H2", "O", "single")))
        self.assertRaises(TypeError,
                          lambda: M.universe(0, [(mol, 'water', 10)]))
        self.assertRaises(ValueError,
                          lambda: M.universe('strange', [(mol, 'water', 10)]))
        self.assertRaises(ValueError,
                          lambda: M.universe('strange', [(mol, 10)]))
        self.assertRaises(TypeError,
                          lambda: M.universe('infinite', mol))
        self.assertRaises(ValueError,
                          lambda: M.universe('infinite', [("water", 10)]))
        self.assertRaises(TypeError,
                          lambda: M.universe('infinite', [mol]))
        self.assertRaises(ValueError,
                          lambda: M.universe('infinite', [(10, mol)]))
        self.assertRaises(ValueError,
                          lambda: M.universe('infinite', [(mol, 'water', 10)],
                                             [(IN.zeros((3,3), N.float64),
                                               IN.zeros((3,), N.float64))]))

    def test_configuration(self):
        mol = make_water_fragment()
        universe = M.universe('cube', [(mol, 'water', 10)])
        # Missing data
        self.assertRaises(TypeError,
                          lambda: M.Configuration(universe))
        # Positions but no cell parameters
        self.assertRaises(TypeError,
                          lambda: M.Configuration(universe,
                                                  IN.zeros((30, 3), N.float32)))
        # Positions and cell parameters of different dtype
        self.assertRaises(ValueError,
                          lambda: M.Configuration(universe,
                                                  IN.zeros((30, 3), N.float32),
                                                  N.float64(10.)))
        # Positions not an array
        self.assertRaises(TypeError,
                          lambda: M.Configuration(universe,
                                                  list(IN.zeros((30, 3),
                                                               N.float32)),
                                                  N.float32(10.)))
        # Positions of wrong shape
        self.assertRaises(ValueError,
                          lambda: M.Configuration(universe,
                                                  IN.zeros((25, 3), N.float32),
                                                  N.float32(10.)))
        # Cell parameters of wrong shape
        self.assertRaises(ValueError,
                          lambda: M.Configuration(universe,
                                                  IN.zeros((30, 3), N.float32),
                                                  IN.zeros((3,), N.float32)))
        
    def test_selection(self):
        mol = make_water_fragment(nsites=2)
        universe = M.universe('cube', [(mol, 'water', 5)])
        # Index occurs twice
        self.assertRaises(ValueError,
                          lambda: M.AtomSelection(universe,
                                                  IN.zeros((2,), N.uint16)))
        # Atom index too large
        self.assertRaises(ValueError,
                          lambda: M.AtomSelection(universe,
                                                  IN.array([20], N.uint16)))
        # Template atom index too large
        self.assertRaises(ValueError,
                          lambda: M.TemplateAtomSelection(universe,
                                                          IN.array([3], N.uint8)))
        # Site index too large
        self.assertRaises(ValueError,
                          lambda: M.SiteSelection(universe,
                                                  IN.array([40], N.uint16)))
        # Template site index too large
        self.assertRaises(ValueError,
                          lambda: M.TemplateSiteSelection(universe,
                                                          IN.array([8], N.uint8)))

class UniverseTest(unittest.TestCase):

    def setUp(self):
        mol = make_water_fragment(2)
        self.universe = M.universe('infinite', [(mol, 'water', 10)],
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
        masses = M.TemplateAtomProperty(self.universe,
                                        "masses", "amu",
                                        IN.array([1., 1., 16.], N.float32))
        self.assertTrue(is_valid(masses))
        self.assertEqual(masses.type, 'template_atom')
        self.assertTrue(masses.universe == self.universe)
        self.assertEqual(masses.element_shape, ())
        self.assertEqual(masses.data.shape, (3,))
        bead_masses = M.TemplateSiteProperty(self.universe,
                                             "mass", "amu",
                                             IN.array([1., 1.,
                                                       1., 1.,
                                                       8., 8.], N.float32))
        self.assertTrue(is_valid(bead_masses))
        self.assertEqual(bead_masses.type, 'template_site')
        self.assertTrue(bead_masses.universe is self.universe)
        self.assertEqual(bead_masses.element_shape, ())
        self.assertEqual(bead_masses.data.shape, (6,))
        velocities = M.SiteProperty(self.universe,
                                    "velocity", "nm ps-1",
                                    IN.zeros((60, 3), dtype=N.float64))
        self.assertTrue(is_valid(velocities))
        self.assertEqual(velocities.type, 'site')
        self.assertTrue(velocities.universe is self.universe)
        self.assertEqual(velocities.data.shape, (60, 3))
        self.assertEqual(velocities.element_shape, (3,))
        foo = M.AtomProperty(self.universe,
                             "foo", "",
                             IN.zeros((30, 2, 2), dtype=N.int16))
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
        s = M.AtomSelection(self.universe, IN.array([0, 2], N.uint8))
        self.assertEqual(s.number_of_atoms, 2)
        self.assertEqual(s.number_of_sites, 4)
        s = M.TemplateAtomSelection(self.universe, IN.array([1], N.uint8))
        self.assertEqual(s.number_of_atoms, 10)
        self.assertEqual(s.number_of_sites, 20)
        s = M.SiteSelection(self.universe, [2])
        self.assertEqual(s.number_of_sites, 1)
        s = M.TemplateSiteSelection(self.universe, [0])
        self.assertEqual(s.number_of_sites, 10)

class PBCTest(unittest.TestCase):

    def setUp(self):
        self.infinite = M.universe('infinite', ())
        self.cube = M.universe('cube', ())
        self.cuboid = M.universe('cuboid', ())
        self.parallelepiped = M.universe('parallelepiped', ())

    def test_lattice_vectors(self):
        conf = M.Configuration(self.infinite,
                               IN.zeros((0, 3), N.float32),
                               None)
        self.assertEqual(conf.lattice_vectors(), ())
        self.assertEqual(conf.cell_volume(), None)
        conf = M.Configuration(self.cube,
                               IN.zeros((0, 3), N.float32),
                               IN.array(1., N.float32))
        lv = conf.lattice_vectors()
        self.assertTrue((lv[0] == N.array([1., 0., 0.], N.float32)).all())
        self.assertTrue((lv[1] == N.array([0., 1., 0.], N.float32)).all())
        self.assertTrue((lv[2] == N.array([0., 0., 1.], N.float32)).all())
        self.assertEqual(conf.cell_volume(), 1.)
        conf = M.Configuration(self.cuboid,
                               IN.zeros((0, 3), N.float32),
                               IN.array([1., 2., 4.], N.float32))
        lv = conf.lattice_vectors()
        self.assertTrue((lv[0] == N.array([1., 0., 0.], N.float32)).all())
        self.assertTrue((lv[1] == N.array([0., 2., 0.], N.float32)).all())
        self.assertTrue((lv[2] == N.array([0., 0., 4.], N.float32)).all())
        self.assertEqual(conf.cell_volume(), 8.)
        conf = M.Configuration(self.parallelepiped,
                               IN.zeros((0, 3), N.float32),
                               IN.array([[1., 2., 4.],
                                         [8., 4., 2.],
                                         [16., 4., 8.]], N.float32))
        lv = conf.lattice_vectors()
        self.assertTrue((lv[0] == N.array([1., 2., 4.], N.float32)).all())
        self.assertTrue((lv[1] == N.array([8., 4., 2.], N.float32)).all())
        self.assertTrue((lv[2] == N.array([16., 4., 8.], N.float32)).all())
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
