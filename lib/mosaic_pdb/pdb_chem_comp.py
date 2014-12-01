# -*- coding: utf-8 -*-

"""PDB Chemical Component Dictionary

The PDB uses a dictionary defining all the chemical components from
its structure files, each identified by a three-letter code. This
dictionary permits in particular the identification of chemical bonds.
The dictionary consists of a huge mmCIF file plus a smaller one.
Downloading and parsing these files takes a lot of time, therefore
the result is cached in a pickle file in $HOME/.mosaic.

The dictionary is updated as new components are added with new
structures being published. If a conversion fails, it's a good idea to
delete the file
$HOME/.mosaic/pdb_chemical_components_dictionary.pickle which is
then regenerated automatically on the next run.

.. moduleauthor:: Konrad Hinsen <konrad.hinsen@cnrs-orleans.fr>

"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import gzip
import os
import sys

# Python 2/3 compatibility issues
if sys.version_info[0] == 2:
    from StringIO import StringIO
    from urllib2 import urlopen
    import cPickle as pickle
    bytes2text = lambda x: x
else:
    from io import BytesIO as StringIO
    from urllib.request import urlopen
    import pickle
    from io import TextIOWrapper
    def bytes2text(stream):
        return TextIOWrapper(stream, encoding="utf8")

from mosaic_pdb import mmcif

# For now, MOSAIC_DIR is used only for the PDB Chemical Components Dictionary,
# which is why it is handled here.

MOSAIC_DIR = os.path.join(os.environ['HOME'], '.mosaic')
if not os.path.exists(MOSAIC_DIR):
    os.makedirs(MOSAIC_DIR)

# PDB Chemical Components Dictionary
#
# ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz
# ftp://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif.gz
#
# These big files are downloaded once and converted to a pickle
# file stored in PDB_CHEM_COMP_DICT that is loaded much faster.

if sys.version_info[0] == 2:
    PDB_CHEM_COMP_DICT = os.path.join(MOSAIC_DIR,
                                      'pdb_chemical_components_dictionary.pickle')
else:
    PDB_CHEM_COMP_DICT = os.path.join(MOSAIC_DIR,
                                      'pdb_chemical_components_dictionary_py3.pickle')

chem_comp_urls = ['ftp://ftp.wwpdb.org/pub/pdb/data/monomers/'
                  'components.cif.gz',
                  'ftp://ftp.wwpdb.org/pub/pdb/data/monomers/'
                  'aa-variants-v1.cif.gz']


# The in-memory representation of the dictionary

class ChemComp(object):

    def __repr__(self):
        return "ChemComp('%s')" % self.id


class Table(object):

    def __init__(self, name, header, conversions):
        self.name = name
        self.fields = [f for f in header if f in conversions]
        self.select = [(i, conversions[f])
                       for i, f in enumerate(header)
                       if f in conversions]
        self.entries = []

    def addEntry(self, data):
        self.entries.append(tuple(c(data[i]) for i, c in self.select))

    def __repr__(self):
        return "Table('%s')" % self.name

    def __getitem__(self, item):
        if isinstance(item, int):
            data = self.entries[item]
            return dict((k, data[i]) for i, k in enumerate(self.fields))
        elif isinstance(item, str):
            try:
                i = self.fields.index(item)
            except ValueError:
                raise KeyError(item)
            return [e[i] for e in self.entries]


class SingleEntryTable(Table):

    def __init__(self, name, conversions):
        self.name = name
        self.conversions = conversions
        self.fields = []
        self.entries = [[]]

    def __setitem__(self, item, value):
        if item in self.conversions:
            self.fields.append(item)
            self.entries[0].append(self.conversions[item](value))


def string_item(s):
    if s == '?':
        return None
    else:
        return s


def float_item(s):
    if s == '?':
        return None
    else:
        return float(s)


# Try to load the pickle file. If it doesn't exist, load
# the data from the official URLs and create the pickle file.

try:

    with open(PDB_CHEM_COMP_DICT, 'rb') as pickled_dict:
        components, variants = pickle.load(pickled_dict)

except IOError:

    sys.stderr.write("Downloading the PDB Chemical Component Dictionary.\n")
    sys.stderr.write("This make take a few minutes.\n")

    table_data = {'chem_comp_atom': 'atoms',
                  'chem_comp_bond': 'bonds'}

    table_items = {'chem_comp_atom': {'atom_id': string_item,
                                      'type_symbol': string_item,
                                      'charge': float_item},
                   'chem_comp_bond': {'atom_id_1': string_item,
                                      'atom_id_2': string_item,
                                      'value_order': string_item}}

    components = {}
    variants = {}

    for url in chem_comp_urls:
        handle = urlopen(url)
        buffer = StringIO(handle.read())
        handle.close()
        f = bytes2text(gzip.GzipFile(fileobj=buffer))
        parser = mmcif.MMCIFParser(file_object=f)
        for item_type, item in parser.parse():
            if item_type is mmcif.KEYWORD:
                if item[0] == 'data':
                    comp_id = item[1]
                    comp = ChemComp()
                    components[comp_id] = comp
                else:
                    raise ValueError("Keyword %s not yet implemented"
                                     % item[0])
            elif item_type is mmcif.TABLE_HEADER:
                table_name = item[0][0]
                if table_name in table_data:
                    table = Table(table_name, [i[1] for i in item],
                                  table_items[table_name])
                    setattr(comp, table_data[table_name], table)
                else:
                    table = None
            elif item_type is mmcif.TABLE_DATA:
                if table is not None:
                    table.addEntry(item)
            elif item_type is mmcif.DATA:
                if item[0] == 'chem_comp':
                    setattr(comp, item[1], item[2])
                    if item[1] == 'three_letter_code':
                        if item[2] != comp_id:
                            vs = variants.get(item[2], [])
                            vs.append(comp_id)
                            variants[item[2]] = vs
                elif item[0] in ['chem_comp_atom', 'chem_comp_bond']:
                    try:
                        table = getattr(comp, table_data[item[0]])
                    except AttributeError:
                        table = SingleEntryTable(item[0], table_items[item[0]])
                        setattr(comp, table_data[item[0]], table)
                    table[item[1]] = item[2]
                else:
                    if item[0] not in ['pdbx_chem_comp_descriptor',
                                       'pdbx_chem_comp_identifier',
                                       'pdbx_chem_comp_feature',
                                       'pdbx_chem_comp_audit']:
                        raise ValueError("Unexpected data item %s" % item[0])
            else:
                raise ValueError("Unexpected item type %s"
                                 % str(item_type),
                                 parser.line_number)

    f = open(PDB_CHEM_COMP_DICT, 'wb')
    pickle.dump((components, variants), f, protocol=-1)
    f.close()
