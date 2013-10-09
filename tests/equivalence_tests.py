# -*- coding: utf-8 -*-

# Equivalence tests between different models

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import unittest
import copy
import itertools as IT

import numpy as N
import immutable.np as IN

import mosaic.immutable_model as IM
import mosaic.mutable_model as MM
import mosaic.array_model as AM
from mosaic.api import is_valid

class WaterTest(unittest.TestCase):

    def make_mutable_universe(self):
        water = MM.Fragment("solvent", "water", (),
                            (MM.Atom("H1", MM.Element("H"), 8),
                             MM.Atom("H2", MM.Element("H"), 8),
                             MM.Atom("O",  MM.Element("O"), 2)),
                            (("H1", "O", "single"), ("H2", "O", "single")))
        return MM.Universe('infinite', [(water, 10)],
                           convention='my_own')

    def make_immutable_universe(self):
        water = IM.fragment("water", (),
                            (("H1", IM.atom(IM.element("H"), 8)),
                             ("H2", IM.atom(IM.element("H"), 8)),
                             ("O",  IM.atom(IM.element("O"), 2))),
                            (("H1", "O", "single"), ("H2", "O", "single")))
        return IM.universe('infinite', [(water, "solvent", 10)],
                           convention='my_own')

    def test_equivalence(self):
        models = [self.make_mutable_universe(),
                  self.make_immutable_universe()]
        models.extend([AM.Universe(m) for m in models])
        for m1 in models:
            for m2 in models:
                self.assertTrue(m1.is_equivalent(m2))


class PeptideTest(unittest.TestCase):

    def make_mutable_universe(self):
        C = MM.Element('C')
        H = MM.Element('H')
        N = MM.Element('N')
        O = MM.Element('O')
        peptide_group = MM.Fragment('peptide', 'peptide',
                                    (),
                                    (MM.Atom('CA', C),
                                     MM.Atom('HA', H),
                                     MM.Atom('H', H),
                                     MM.Atom('N', N),
                                     MM.Atom('C', C),
                                     MM.Atom('O', O)),
                                    (('N', 'H', "single"),
                                     ('N', 'CA', "single"),
                                     ('CA', 'HA', "single"),
                                     ('CA', 'C', "single"),
                                     ('C', 'O', "double")))
        ala_sidechain = MM.Fragment('sidechain', 'ala_sidechain',
                                    (),
                                    (MM.Atom('CB', C),
                                     MM.Atom('HB1', H),
                                     MM.Atom('HB2', H),
                                     MM.Atom('HB3', H)),
                                    (('CB', 'HB1', "single"),
                                     ('CB', 'HB2', "single"),
                                     ('CB', 'HB3', "single"),))
        ala = lambda label: MM.Fragment(label, 'alanine',
                                        (copy.copy(peptide_group),
                                         copy.copy(ala_sidechain)),
                                        (),
                                        (('peptide.CA', 'sidechain.CB',
                                          "single"),))
        di_ala = MM.Polymer('di_ala', 'alanine_dipeptide',
                            (ala('ALA1'), ala('ALA2')),
                            (('ALA1.peptide.C', 'ALA2.peptide.N', "single"),),
                            'polypeptide')
        water = MM.Fragment("solvent", "water", (),
                            (MM.Atom("H1", MM.Element("H"), 8),
                             MM.Atom("H2", MM.Element("H"), 8),
                             MM.Atom("O",  MM.Element("O"), 2)),
                            (("H1", "O", "single"), ("H2", "O", "single")))
        return MM.Universe('cube', [(di_ala, 1), (water, 10)],
                           convention='my_own')

    def make_immutable_universe(self):
        C = IM.element('C')
        H = IM.element('H')
        N = IM.element('N')
        O = IM.element('O')
        peptide_group = IM.fragment('peptide',
                                    (),
                                    (('CA', IM.atom(C)),
                                     ('HA', IM.atom(H)),
                                     ('H', IM.atom(H)),
                                     ('N', IM.atom(N)),
                                     ('C', IM.atom(C)),
                                     ('O', IM.atom(O))),
                                    (('N', 'H', "single"),
                                     ('N', 'CA', "single"),
                                     ('CA', 'HA', "single"),
                                     ('CA', 'C', "single"),
                                     ('C', 'O', "double")))
        ala_sidechain = IM.fragment('ala_sidechain',
                                    (),
                                    (('CB', IM.atom(C)),
                                     ('HB1', IM.atom(H)),
                                     ('HB2', IM.atom(H)),
                                     ('HB3', IM.atom(H))),
                                    (('CB', 'HB1', "single"),
                                     ('CB', 'HB2', "single"),
                                     ('CB', 'HB3', "single"),))
        ala = IM.fragment('alanine',
                          (('peptide', peptide_group),
                           ('sidechain', ala_sidechain)),
                          (),
                          (('peptide.CA', 'sidechain.CB', "single"),))
        di_ala = IM.polymer('alanine_dipeptide',
                            (('ALA1', ala),
                             ('ALA2', ala)),
                            (('ALA1.peptide.C', 'ALA2.peptide.N', "single"),),
                            'polypeptide')
        water = IM.fragment("water", (),
                            (("H1", IM.atom(H, 8)),
                             ("H2", IM.atom(H, 8)),
                             ("O",  IM.atom(O, 2))),
                            (("H1", "O", "single"), ("H2", "O", "single")))
        return IM.universe('cube',
                           [(di_ala, "di_ala", 1),
                            (water, "solvent", 10)],
                           convention='my_own')

    def make_configuration(self, M, universe):
        return M.Configuration(universe,
                               IN.zeros((universe.number_of_sites, 3),
                                        N.float32),
                               IN.array(1, N.float32))

    def make_masses(self, M, universe):
        mass = {'C': 12., 'H': 1., 'N': 14., 'O': 16.}
        m = [mass[a.name]
             for a in IT.chain.from_iterable(f.recursive_atom_iterator()
                                             for f, _ in universe.molecules)]
        m = IN.array(m, N.float32)
        try:
            return M.TemplateAtomProperty(universe, "mass", "g mol-1", m)
        except AttributeError:
            return M.Property(universe, "mass", "g mol-1", m, 'template_atom')

    def test_equivalence(self):
        models = [self.make_immutable_universe(), self.make_mutable_universe()]
        models.extend([AM.Universe(m) for m in models])
        for m1 in models:
            for m2 in models:
                self.assertTrue(m1.is_equivalent(m2))
        confs = [self.make_configuration(M, u)
                 for M, u in zip((IM, MM, AM), models)]
        for c1 in confs:
            for c2 in confs:
                self.assertTrue(c1.is_equivalent(c2))
        props = [self.make_masses(M, u)
                 for M, u in zip((IM, MM, AM), models)]
        for p1 in props:
            for p2 in props:
                self.assertTrue(p1.is_equivalent(p2))

    def test_factory(self):
        models = [self.make_immutable_universe(),
                  self.make_mutable_universe()]
        models.extend([AM.Universe(m) for m in models])
        confs = [self.make_configuration(M, u)
                 for M, u in zip((IM, MM, AM), models)]
        factories = [M.Factory() for M in (IM, MM, AM)]
        for f in factories:
            for c in confs:
                c.validate()
                self.assertTrue(is_valid(c))
                fc = f(c)
                self.assertTrue(c.is_equivalent(fc))
                self.assertTrue(is_valid(fc))

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(WaterTest))
    s.addTest(loader.loadTestsFromTestCase(PeptideTest))
    return s

if __name__ == '__main__':
    unittest.main()
