# -*- coding: utf-8 -*-

# Tests for module mosaic.xml_io

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import unittest
import copy
import io
import sys

import numpy as N
import numpy.random as NR
import immutable.np as IN

import mosaic.mutable_model
import mosaic.immutable_model
import mosaic.array_model
from mosaic.api import is_valid
from mosaic.xml_io import XMLWriter, XMLReader

if sys.version_info[0] == 2:
    def StringIO(buffer=None):
        return io.BytesIO(buffer)
else:
    def StringIO(buffer=None):
        return io.StringIO(buffer)

class PeptideTest(object):

    def test_universe(self):
        file = StringIO()
        with XMLWriter(file) as writer:
            writer.store("universe", self.universe)
        xml_string = file.getvalue()
        file = StringIO(xml_string)
        with XMLReader(file) as reader:
            data = {xml_id: obj for xml_id, obj in reader}
        universe = data['universe']
        self.assertTrue(universe.is_equivalent(self.universe))
        self.assertTrue(is_valid(universe))

    def test_configuration(self):
        M = self.M
        for dtype in (N.float32, N.float64):
            conf = M.Configuration(self.universe,
                                   IN.array(NR.normal(size=(40,3)),
                                            dtype=dtype),
                                   IN.array(10., dtype))
            self.assertTrue(is_valid(conf))
            file = StringIO()
            with XMLWriter(file) as writer:
                writer.store("universe", self.universe)
                writer.store("configuration", conf)
            xml_string = file.getvalue()
            file = StringIO(xml_string)
            with XMLReader(file) as reader:
                data = {xml_id: obj for xml_id, obj in reader}
            conf_r = data['configuration']
            self.assertTrue(conf_r.is_equivalent(conf))
            self.assertTrue(is_valid(conf_r))

    def test_property(self):
        M = self.M
        universe = self.universe
        for dtype in (N.float32, N.float64):
            for p in [M.AtomProperty(universe, "mass", "amu",
                         IN.array(NR.normal(size=(universe.number_of_atoms,)))),
                      M.SiteProperty(universe, "foo", "",
                         IN.array(NR.normal(size=(universe.number_of_sites, 3)))),
                      M.TemplateAtomProperty(universe, "mass", "amu",
                         IN.array(NR.normal(size=(universe.number_of_template_atoms,
                                                  1)))),
                      M.TemplateSiteProperty(universe, "foo", "bar",
                         IN.array(NR.normal(size=(universe.number_of_template_sites,
                                                  2, 2)))),
                      ]:
                self.assertTrue(is_valid(p))
                file = StringIO()
                with XMLWriter(file) as writer:
                    writer.store("universe", self.universe)
                    writer.store("property", p)
                xml_string = file.getvalue()
                file = StringIO(xml_string)
                with XMLReader(file) as reader:
                    data = {xml_id: obj for xml_id, obj in reader}
                p_r = data['property']
                self.assertTrue(p_r.is_equivalent(p))
                self.assertTrue(is_valid(p_r))

    def test_label(self):
        M = self.M
        universe = self.universe
        for l in [M.AtomLabel(universe, "element",
                              tuple(a.name
                                    for f, n in universe.molecules
                                    for _ in range(n)
                                    for a in f.recursive_atom_iterator())),
                  M.TemplateAtomLabel(universe, "element",
                                      tuple(a.name
                                            for f, n in universe.molecules
                                            for a in f.recursive_atom_iterator())),
                  M.SiteLabel(universe, "element",
                              tuple(a.name
                                    for f, n in self.universe.molecules
                                    for _ in range(n)
                                    for a in f.recursive_atom_iterator()
                                    for __ in range(a.number_of_sites))),
                  M.TemplateSiteLabel(universe, "element",
                                      tuple(a.name
                                            for f, n in self.universe.molecules
                                            for a in f.recursive_atom_iterator()
                                            for _ in range(a.number_of_sites))),
                  ]:
            self.assertTrue(is_valid(l))
            file = StringIO()
            with XMLWriter(file) as writer:
                writer.store("universe", self.universe)
                writer.store("label", l)
            xml_string = file.getvalue()
            file = StringIO(xml_string)
            with XMLReader(file) as reader:
                data = {xml_id: obj for xml_id, obj in reader}
            l_r = data['label']
            self.assertTrue(l_r.is_equivalent(l))
            self.assertTrue(is_valid(l_r))

    def test_selection(self):
        M = self.M
        universe = self.universe
        for s in [M.AtomSelection(universe, IN.array([0], N.uint8)),
                  M.SiteSelection(universe, IN.array([1, 2], N.uint8)),
                  M.TemplateAtomSelection(universe, IN.array([0, 2], N.uint8)),
                  M.TemplateSiteSelection(universe, IN.array([0, 3], N.uint8))]:
            self.assertTrue(is_valid(s))
            file = StringIO()
            with XMLWriter(file) as writer:
                writer.store("universe", self.universe)
                writer.store("selection", s)
            xml_string = file.getvalue()
            file = StringIO(xml_string)
            with XMLReader(file) as reader:
                data = {xml_id: obj for xml_id, obj in reader}
            s_r = data['selection']
            self.assertTrue(s_r.is_equivalent(s))
            self.assertTrue(is_valid(s_r))


class MMPeptideTest(unittest.TestCase,
                    PeptideTest):
    M = mosaic.mutable_model

    def setUp(self):
        M = self.M
        C = M.Element('C')
        H = M.Element('H')
        N = M.Element('N')
        O = M.Element('O')
        peptide_group = M.Fragment('peptide', 'peptide_group',
                                   (),
                                   (M.Atom('CA', C),
                                    M.Atom('HA', H),
                                    M.Atom('H', H),
                                    M.Atom('N', N),
                                    M.Atom('C', C),
                                    M.Atom('O', O)),
                                    (('N', 'H', 'single'),
                                     ('N', 'CA', 'single'),
                                     ('CA', 'HA', 'single'),
                                     ('CA', 'C', 'single'),
                                     ('C', 'O', 'double')))
        ala_sidechain = M.Fragment('sidechain', 'ala_sidechain',
                                   (),
                                   (M.Atom('CB', C),
                                    M.Atom('HB1', H),
                                    M.Atom('HB2', H),
                                    M.Atom('HB3', H)),
                                   (('CB', 'HB1', 'single'),
                                    ('CB', 'HB2', 'single'),
                                    ('CB', 'HB3', 'single'),))
        ala = lambda label: M.Fragment(label, 'alanine',
                                       (copy.copy(peptide_group),
                                        copy.copy(ala_sidechain)),
                                       (),
                                       (('peptide.CA',
                                         'sidechain.CB',
                                         'single'),))
        ala_dipeptide = M.Polymer('di_ala', 'alanine_dipeptide',
                                  (ala('ALA1'), ala('ALA2')),
                                  (('ALA1.peptide.C',
                                    'ALA2.peptide.N',
                                    'single'),),
                                  'polypeptide')
        self.universe = M.Universe('cube', [(ala_dipeptide, 2)],
                                   convention='my_own')
        self.assertTrue(is_valid(self.universe))


class IMPeptideTest(unittest.TestCase,
                    PeptideTest):
    M = mosaic.immutable_model

    def setUp(self):
        M = self.M
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
        ala_dipeptide = M.polymer('alanine_dipeptide',
                                  (('ALA1', ala),
                                   ('ALA2', ala)),
                                  (('ALA1.peptide.C', 'ALA2.peptide.N', "single"),),
                                  'polypeptide')
        self.universe = M.universe('cube', ((ala_dipeptide, 'di_ala', 2),),
                                   convention='my_own')
        self.assertTrue(is_valid(self.universe))

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(MMPeptideTest))
    s.addTest(loader.loadTestsFromTestCase(IMPeptideTest))
    return s

if __name__ == '__main__':
    unittest.main()
