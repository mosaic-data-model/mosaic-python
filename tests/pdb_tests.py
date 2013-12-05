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
from mosaic_pdb.import_structure import make_models

class PDB2ONXTest(unittest.TestCase):

    def setUp(self):
        self.pdb_model = dict(make_models("2ONX.cif"))
        self.xml_model = dict(XMLReader("2ONX.xml"))

    def test_equivalence(self):
        for name, data in self.pdb_model.items():
            self.assertTrue(name in self.xml_model)
            self.assertTrue(data.is_equivalent(self.xml_model[name]))

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(PDB2ONXTest))
    return s

if __name__ == '__main__':
    unittest.main()
