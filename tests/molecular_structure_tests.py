# -*- coding: utf-8 -*-

# Tests for the molecular structure utilities

#-----------------------------------------------------------------------------
#       Copyright (C) 2014 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import unittest

from mosaic.xml_io import XMLReader
import mosaic.mutable_model as M
from mosaic.molecular_structure import FragmentIterator


def peptide_predicate(fragment):
    return fragment.is_polymer and fragment.polymer_type == 'polypeptide'

def two_atoms(fragment):
    return fragment.number_of_atoms == 2

def three_atoms(fragment):
    return fragment.number_of_atoms == 3

class PDB2ONXTest(unittest.TestCase):

    def setUp(self):
        data = dict(XMLReader("2ONX.xml"))
        self.universe = data['universe']

    def test_peptide_chains(self):
        peptide_chains = list(FragmentIterator(self.universe,
                                               peptide_predicate))
        self.assertEqual(len(peptide_chains), 1)
        chain, atom_index, site_index = peptide_chains[0]
        self.assertEqual(chain.label, 'A')
        self.assertEqual(atom_index, 0)
        self.assertEqual(site_index, 0)


class CalphaModelTest(unittest.TestCase):

    def setUp(self):
        water = M.Fragment("water", "water", (),
                           (M.Atom("H1", M.Element("H")),
                            M.Atom("H2", M.Element("H")),
                            M.Atom("O",  M.Element("O"))),
                           (("H1", "O", "single"), ("H2", "O", "single")))
        peptide_chain = M.Polymer('chain', 'ALA,LYS',
                                  (M.Fragment('ALA1', 'ALA', (),
                                              (M.Atom('CA',
                                                      M.CGParticle('CA'), 1),),
                                              ()),
                                   M.Fragment('LYS2', 'LYS', (),
                                              (M.Atom('CA',
                                                      M.CGParticle('CA'), 2),),
                                              ())),
                                  (('ALA1.CA', 'LYS2.CA', ''),),
                                  'polypeptide')
        protein = M.Fragment('protein', 'protein', (peptide_chain,), (), ())
        self.universe = M.Universe('infinite', [(protein, 2), (water, 5)])

    def test_peptide_chains(self):
        peptide_chains = list(FragmentIterator(self.universe,
                                               peptide_predicate))
        self.assertEqual(len(peptide_chains), 2)
        chain, atom_index, site_index = peptide_chains[0]
        self.assertEqual(chain.label, 'chain')
        self.assertEqual(atom_index, 0)
        self.assertEqual(site_index, 0)
        chain, atom_index, site_index = peptide_chains[1]
        self.assertEqual(chain.label, 'chain')
        self.assertEqual(atom_index, 2)
        self.assertEqual(site_index, 3)

    def test_template_peptide_chains(self):
        peptide_chains = list(FragmentIterator(self.universe, peptide_predicate,
                                               template=True))
        self.assertEqual(len(peptide_chains), 1)
        chain, atom_index, site_index = peptide_chains[0]
        self.assertEqual(chain.label, 'chain')
        self.assertEqual(atom_index, 0)
        self.assertEqual(site_index, 0)

    def test_three_atoms(self):
        matches = list(FragmentIterator(self.universe, three_atoms))
        self.assertEqual(len(matches), 5)
        for i, (fragment, atom_index, site_index) in enumerate(matches):
            self.assertEqual(fragment.label, 'water')
            self.assertEqual(atom_index, 4+3*i)
            self.assertEqual(site_index, 6+3*i)

    def test_template_three_atoms(self):
        matches = list(FragmentIterator(self.universe, three_atoms,
                                    template=True))
        self.assertEqual(len(matches), 1)
        fragment, atom_index, site_index = matches[0]
        self.assertEqual(fragment.label, 'water')
        self.assertEqual(atom_index, 2)
        self.assertEqual(site_index, 3)

    def test_two_atoms(self):
        matches = list(FragmentIterator(self.universe, two_atoms))
        self.assertEqual(len(matches), 2)
        for i, (fragment, atom_index, site_index) in enumerate(matches):
            self.assertEqual(fragment.label, 'protein')
            self.assertEqual(atom_index, 2*i)
            self.assertEqual(site_index, 3*i)

    def test_two_atoms_descend(self):
        matches = list(FragmentIterator(self.universe, two_atoms,
                                        descend_on_match=True))
        self.assertEqual(len(matches), 4)
        for i in range(2):
            fragment, atom_index, site_index = matches[2*i]
            self.assertEqual(fragment.label, 'protein')
            self.assertEqual(atom_index, 2*i)
            self.assertEqual(site_index, 3*i)
            fragment, atom_index, site_index = matches[2*i+1]
            self.assertEqual(fragment.label, 'chain')
            self.assertEqual(atom_index, 2*i)
            self.assertEqual(site_index, 3*i)

    def test_template_two_atoms_descend(self):
        matches = list(FragmentIterator(self.universe, two_atoms,
                                        template=True,
                                        descend_on_match=True))
        self.assertEqual(len(matches), 2)
        fragment, atom_index, site_index = matches[0]
        self.assertEqual(fragment.label, 'protein')
        self.assertEqual(atom_index, 0)
        self.assertEqual(site_index, 0)
        fragment, atom_index, site_index = matches[1]
        self.assertEqual(fragment.label, 'chain')
        self.assertEqual(atom_index, 0)
        self.assertEqual(site_index, 0)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(CalphaModelTest))
    s.addTest(loader.loadTestsFromTestCase(PDB2ONXTest))
    return s

if __name__ == '__main__':
    unittest.main()
