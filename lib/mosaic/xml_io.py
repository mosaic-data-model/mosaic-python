# -*- coding: utf-8 -*-
"""XML I/O

.. moduleauthor:: Konrad Hinsen

The XML representation of Mosaic data consists of a single element
``<mosaic version="x.y">...</mosaic>`` that can contain any number of
Mosaic data items, each of which is identified by a unique id.
The classes :class:`XMLWriter` and :class:`XMLReader` handle the
translation between Mosaic data items in memory that can reference
each other and XML elements that reference each other by id.

A common pattern for generating an XML file is:

::

    with XMLWriter('molecule.xml') as writer:
       writer.store("universe", universe)
       writer.store("configuration1", configuration1)
       writer.store("configuration2", configuration2)

A common pattern for reading from an XML file is:

::

    items = {}
    for id, data in XMLReader('molecule.xml'):
       items[id] = data

"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

try:
    import xml.etree.cElementTree as ET
except IOError:
    import xml.etree.ElementTree as ET

import numpy as N
import immutable.np as IN

import mosaic.api as api
import mosaic.immutable_model as im
from mosaic.utility import MethodRegister
from mosaic.utility import uint_for_max_value
from mosaic.utility import isstring
from mosaic.utility import xml_encoding


# Number formatting

def n2s(x, dp=False):
    if isinstance(x, int):
        return str(x)
    else:
        if dp:
            return "{0:.20g}".format(x)
        else:
            return "{0:.10g}".format(x)


# Data type formatting

def t2s(t):
    return {N.dtype(N.int8): "int8",
            N.dtype(N.int16): "int16",
            N.dtype(N.int32): "int32",
            N.dtype(N.int32): "int32",
            N.dtype(N.uint8): "uint8",
            N.dtype(N.uint16): "uint16",
            N.dtype(N.uint32): "uint32",
            N.dtype(N.uint32): "uint32",
            N.dtype(N.float32): "float32",
            N.dtype(N.float64): "float64",
            N.dtype(N.bool): "boolean"
            }[t]


def s2t(s):
    return {"int8": N.int8,
            "int16": N.int16,
            "int32": N.int32,
            "int32": N.int32,
            "uint8": N.uint8,
            "uint16": N.uint16,
            "uint32": N.uint32,
            "uint32": N.uint32,
            "float32": N.float32,
            "float64": N.float64,
            "boolean": N.bool
            }[s]

class XMLStore(object):

    def __init__(self):
        self._id_map = {}
        self._data_map = {}

    def _register_data_item(self, xml_id, data_item):
        self._data_map[xml_id] = data_item
        self._id_map[data_item] = xml_id

    def _get_id(self, data_item):
        return self._id_map.get(data_item, None)

    def _get_data(self, xml_id):
        return self._data_map.get(xml_id, None)

    # Make XMLStores work as context managers
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


class XMLWriter(XMLStore):

    """Output handler for XML files

    This class handles the translation of references between Mosaic
    data items in memory to XML references by id. References are
    allowed only to data items that have been written earlier using
    the same XMLWriter instance.

    An XMLWriter can be used like a file, in which case it must be
    closed after the last ``store`` operation, or as a context manager
    in a with-statement.
    """

    def __init__(self, xml_file):
        """
        :param xml_file: a writeable file object, or a string
                         interpreted as a file name
        :type xml_file: str or file-like
        """
        XMLStore.__init__(self)
        if isstring(xml_file):
            self.file = open(xml_file, 'w')
            self._close_file = True
        else:
            # assume it is a file-like object
            self.file = xml_file
            self._close_file = False
        # Keep a set of already used ids in order to check they remain unique
        self._xml_ids = set()
        # Start with the XML declaration
        self.file.write('<?xml version="1.0" encoding="utf-8"?>')
        # Open the top-level element
        self.file.write('<mosaic version="%d.%d">' % api.MOSAIC_VERSION)

    def close(self):
        # Close the top-level element
        self.file.write("</mosaic>")
        # Close the underlying file if we own it
        if self._close_file:
            self.file.close()
        # No more output after this
        self.file = None
        # Clear cache
        self._id_map = None
        self._data_map = None

    storage_handler = MethodRegister()

    def store(self, xml_id, data):
        """
        :param xml_id: the id of the XML element representing the data item
        :type xml_id:  str
        :param data:   a Mosaic data item
        :type data:    :class:`mosaic.api.MosaicDataItem`
        """
        api.validate_type(xml_id, str, "xml_id")
        api.validate_type(data, api.MosaicDataItem, "data")
        xml_id = xml_id
        if xml_id in self._xml_ids:
            raise ValueError("XML ID %s has already been used" % xml_id)
        handler = self.storage_handler[type(data)]
        if handler is None:
            raise TypeError("Storage of %s not yet implemented"
                            % str(type(data)))
        handler(self, xml_id, data)
        self._register_data_item(xml_id, data)
        self._xml_ids.add(xml_id)

    def _write(self, element):
        ET.ElementTree(element).write(self.file,
                                      encoding=xml_encoding)

    #
    # Storage handlers
    #
    @storage_handler(api.MosaicUniverse)
    def _store_universe(self, xml_id, universe):
        el = ET.Element('universe')
        el.set("id", xml_id)
        el.set("cell_shape", universe.cell_shape)
        el.set("convention", universe.convention)
        st = universe.symmetry_transformations
        if len(st) > 0:
            el_st = ET.SubElement(el, 'symmetry_transformations')
            for rot, trans in st:
                el_st.append(self._symmetry_el(rot, trans))
        el_molecules = ET.SubElement(el, 'molecules')
        for fragment, count in universe.molecules:
            el_m = ET.Element('molecule', count=str(count))
            el_m.append(self._fragment_element(fragment))
            el_molecules.append(el_m)
        self._write(el)

    def _symmetry_el(self, rot, trans):
        el = ET.Element('transformation')
        el_rot = ET.SubElement(el, 'rotation')
        el_rot.text = ' '.join(n2s(x) for x in rot.flat)
        el_trans = ET.SubElement(el, 'translation')
        el_trans.text = ' '.join(n2s(x) for x in trans.flat)
        return el

    def _fragment_element(self, fragment):
        el = ET.Element('fragment',
                        label=fragment.label,
                        species=fragment.species)
        if fragment.is_polymer:
            el.set('polymer_type', fragment.polymer_type)
        if len(fragment.fragments) > 0:
            el_fragments = ET.SubElement(el, 'fragments')
            for f in fragment.fragments:
                el_fragments.append(self._fragment_element(f))
        if len(fragment.atoms) > 0:
            el_atoms = ET.SubElement(el, 'atoms')
            for a in fragment.atoms:
                el_atom = ET.Element('atom',
                                     label=a.label,
                                     type=a.type,
                                     name=a.name)
                if a.number_of_sites != 1:
                    el_atom.set('nsites', str(a.number_of_sites))
                el_atoms.append(el_atom)
        if len(fragment.bonds) > 0:
            el_bonds = ET.SubElement(el, 'bonds')
            for a1, a2, order in fragment.bonds:
                el_bonds.append(ET.Element('bond',
                                           atoms=a1 + ' ' + a2,
                                           order=order))
        return el

    @storage_handler(api.MosaicConfiguration)
    def _store_configuration(self, xml_id, configuration):
        universe = configuration.universe
        universe_id = self._get_id(universe)
        if universe_id is None:
            raise IOError("universe must be stored first")
        el = ET.Element('configuration')
        el.set("id", xml_id)
        el.append(ET.Element('universe', ref=universe_id))
        float_type = t2s(configuration.positions.dtype)
        dp = float_type == "float64"
        if N.product(configuration.cell_parameters.shape) != 0:
            shape_str = ' '.join(str(x)
                                 for x in configuration.cell_parameters.shape)
            el_cp = ET.Element('cell_parameters', shape=shape_str)

            el_cp.text = ' '.join(n2s(p, dp)
                                  for p in configuration.cell_parameters.flat)
            el.append(el_cp)
        el_pos = ET.Element('positions', type=float_type)
        el_pos.text = ' '.join(' '.join(n2s(x, dp) for x in p)
                               for p in configuration.positions)
        el.append(el_pos)
        self._write(el)

    @storage_handler(api.MosaicProperty)
    def _store_property(self, xml_id, property):
        universe = property.universe
        universe_id = self._get_id(universe)
        if universe_id is None:
            raise IOError("universe must be stored first")
        el = ET.Element(property.type + '_property')
        el.set("id", xml_id)
        el.set("name", property.name)
        el.set("units", property.units)
        el.append(ET.Element('universe', ref=universe_id))
        shape_str = ' '.join(str(x) for x in property.element_shape)
        el_data = ET.Element('data',
                             shape=shape_str,
                             type=t2s(property.data.dtype))
        dp = property.data.dtype == N.dtype(N.float64)
        el_data.text = ' '.join(' '.join(n2s(x, dp) for x in v.flat)
                                for v in property.data)
        el.append(el_data)
        self._write(el)

    @storage_handler(api.MosaicLabel)
    def _store_label(self, xml_id, label):
        universe = label.universe
        universe_id = self._get_id(universe)
        if universe_id is None:
            raise IOError("universe must be stored first")
        el = ET.Element(label.type + '_label')
        el.set("id", xml_id)
        el.set("name", label.name)
        el.append(ET.Element('universe', ref=universe_id))
        el_strings = ET.Element("strings")
        el_strings.text = ' '.join(label.strings)
        el.append(el_strings)
        self._write(el)

    @storage_handler(api.MosaicSelection)
    def _store_selection(self, xml_id, selection):
        universe = selection.universe
        universe_id = self._get_id(universe)
        if universe_id is None:
            raise IOError("universe must be stored first")
        el = ET.Element(selection.type + '_selection')
        el.set("id", xml_id)
        el.append(ET.Element('universe', ref=universe_id))
        el_indices = ET.Element("indices")
        el_indices.text = ' '.join(str(i) for i in selection.indices)
        el.append(el_indices)
        self._write(el)


class XMLReader(XMLStore):

    """Input handler for XML files

    This class handles the translation of references by id in an XML
    file to in-memory references between data items.

    An XMLReader is used as an iterator over ``(id, data_item)`` pairs.

    The current implementation does not do any validation, it assumes
    the XML input to be correct.
    """


    def __init__(self, xml_file):
        """
        :param xml_file: a file object, or a string
                         interpreted as a file name
        :type xml_file: str or file-like
        """
        XMLStore.__init__(self)
        if isstring(xml_file):
            # file name given: open file and close it at the end
            self.file = open(xml_file, 'r')
            self._close_file = True
        else:
            # assume it is a file-like object, don't close at the end
            self.file = xml_file
            self._close_file = False

    def close(self):
        if self._close_file:
            self.file.close()
        # No more reading after this
        self.file = None
        self._close_file = False
        # Clear cache
        self._id_map = None
        self._data_map = None

    def __del__(self):
        self.close()

    def __iter__(self):
        for event, el in ET.iterparse(self.file):
            if el.tag == 'mosaic':
                # check version number
                version = tuple(int(s) for s in el.get("version").split('.'))
                if version[0] > api.MOSAIC_VERSION[0] \
                   or version[1] > api.MOSAIC_VERSION[1]:
                    raise ValueError("XML data is for version %d.%d, "
                                     "software is version %d.%d"
                                     % (version + api.MOSAIC_VERSION))
                continue
            xml_id = el.get("id")
            if xml_id is not None:
                handler = self.data_handler[el.tag]
                if handler is None:
                    raise ValueError("Unknown element type %s" % el.tag)
                data = handler(self, el)
                el.clear()
                self._register_data_item(xml_id, data)
                yield xml_id, data

    data_handler = MethodRegister()

    @data_handler("universe")
    def _read_universe(self, el):
        ref = el.get('ref', None)
        if ref is not None:
            return self._get_data(ref)
        molecules = \
          tuple(self._parse_molecule(el_m)
                for el_m in el.iter('molecule'))
        symmetry_transformations = \
          tuple((self._array((3, 3), N.float64, el_t.find('rotation').text),
                 self._array((3,), N.float64, el_t.find('translation').text))
                for el_t in el.iter('transformation'))
        return im.universe(el.get('cell_shape'),
                           molecules,
                           symmetry_transformations,
                           el.get('convention'))

    def _parse_molecule(self, el):
        label, fragment = self._parse_fragment_tree(el.find('fragment'))
        count = int(el.get('count'))
        return fragment, label, count

    def _parse_fragment_tree(self, el):
        fragments = []
        for el_fragments in el.findall('fragments'):
            for el_f in el_fragments:
                fragments.append(self._parse_fragment_tree(el_f))
        atoms = []
        for el_atoms in el.findall('atoms'):
            for el_a in el_atoms:
                atom_classes = {"dummy": im.dummy,
                                "unknown": im.unknown,
                                "element": im.element}
                descr = atom_classes[el_a.get('type')](el_a.get('name'))
                atoms.append((el_a.get('label'),
                              im.atom(descr, int(el_a.get('nsites', '1')))))
        bonds = []
        for el_bonds in el.findall('bonds'):
            for el_b in el_bonds:
                a1, a2 = el_b.get('atoms').split()
                bonds.append((a1, a2, el_b.get('order')))
        label = el.get('label')
        species = el.get('species')
        polymer_type = el.get('polymer_type', None)
        if polymer_type is None:
            return label, im.fragment(species, fragments, atoms, bonds)
        else:
            assert len(atoms) == 0
            return label, im.polymer(species, fragments, bonds, polymer_type)

    @data_handler("configuration")
    def _read_configuration(self, el):
        universe = self._read_universe(el.find('universe'))
        el_pos = el.find('positions')
        dtype = s2t(el_pos.get('type'))
        pos = self._property((3,), dtype, el_pos.text)
        el_cp = el.find('cell_parameters')
        if el_cp is None:
            cp = None
        else:
            cp = self._array(el_cp.get('shape'), dtype, el_cp.text)
        return im.Configuration(universe, pos, cp)

    @data_handler("atom_property")
    @data_handler("site_property")
    @data_handler("template_atom_property")
    @data_handler("template_site_property")
    def _read_property(self, el):
        universe = self._read_universe(el.find('universe'))
        ptype = el.tag[:-9]   # strip off '_property'
        klass = {"atom": im.AtomProperty,
                 "site": im.SiteProperty,
                 "template_atom": im.TemplateAtomProperty,
                 "template_site": im.TemplateSiteProperty}[ptype]
        el_d = el.find('data')
        data = self._property(el_d.get('shape'),
                              s2t(el_d.get('type')),
                              el_d.text)
        return klass(universe, el.get('name'), el.get('units'), data)

    @data_handler("atom_label")
    @data_handler("site_label")
    @data_handler("template_atom_label")
    @data_handler("template_site_label")
    def _read_label(self, el):
        universe = self._read_universe(el.find('universe'))
        ptype = el.tag[:-6]   # strip off '_label'
        klass = {"atom": im.AtomLabel,
                 "site": im.SiteLabel,
                 "template_atom": im.TemplateAtomLabel,
                 "template_site": im.TemplateSiteLabel}[ptype]
        el_strings = el.find('strings')
        strings = tuple(el_strings.text.split())
        return klass(universe, el.get('name'), strings)

    @data_handler("atom_selection")
    @data_handler("site_selection")
    @data_handler("template_atom_selection")
    @data_handler("template_site_selection")
    def _read_selection(self, el):
        universe = self._read_universe(el.find('universe'))
        ptype = el.tag[:-10]   # strip off '_selection'
        klass = {"atom": im.AtomSelection,
                 "site": im.SiteSelection,
                 "template_atom": im.TemplateAtomSelection,
                 "template_site": im.TemplateSiteSelection}[ptype]
        el_d = el.find('indices')
        indices = [int(s) for s in el_d.text.split()]
        indices = IN.array(indices, uint_for_max_value(max(indices)))
        return klass(universe, indices)

    # Parse array data

    def _array(self, shape, dtype, text):
        if isstring(shape):
            shape = tuple(int(s) for s in shape.split())
        return IN.array([dtype(s) for s in text.split()]).reshape(shape)

    def _property(self, el_shape, dtype, text):
        data = IN.array([dtype(s) for s in text.split()])
        if isstring(el_shape):
            el_shape = tuple(int(s) for s in el_shape.split())
        n_els = N.product(el_shape, dtype=N.int)
        assert len(data) % n_els == 0
        n_entries = len(data) // n_els
        return data.reshape((n_entries,)+el_shape)
