# -*- coding: utf-8 -*-

# Tests for the PDB importer

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import unittest

from mosaic.xml_io import XMLReader
from mosaic_pdb.import_structure import make_models, MMCIFStructure

class PDB2ONXTest(unittest.TestCase):

    def setUp(self):
        self.pdb_model = dict(make_models("2ONX.cif"))
        self.xml_model = dict(XMLReader("2ONX.xml"))

    def test_equivalence(self):
        for name, data in self.pdb_model.items():
            self.assertTrue(name in self.xml_model)
            self.assertTrue(data.is_equivalent(self.xml_model[name]))

class AuthSpecTest(unittest.TestCase):

    def setUp(self):
        self.structure = MMCIFStructure(structure_file="2ONX.cif",
                                        with_auth_spec=True)

    def test_auth_spec(self):
        auth_spec = self.structure.auth_spec
        self.assertFalse(auth_spec is None)
        model = self.structure.models[0]
        for asym_id, sites in model.items():
            for residue in sites:
                for atom in residue:
                    unique_id, atom_id, alt_id, comp_id, seq_id, \
                        ins_code, atom_type, x, y, z, occupancy, \
                        u_iso = atom
                    auth = auth_spec[unique_id]
                    self.assertEqual(comp_id, auth['comp_id'])
                    self.assertEqual(seq_id, auth['seq_id'])
                    self.assertEqual(asym_id, auth['asym_id'])
                    self.assertEqual(atom_id, auth['atom_id'])

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(PDB2ONXTest))
    s.addTest(loader.loadTestsFromTestCase(AuthSpecTest))
    return s

if __name__ == '__main__':
    unittest.main()
