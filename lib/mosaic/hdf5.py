# -*- coding: utf-8 -*-
"""HDF5 I/O

.. moduleauthor:: Konrad Hinsen

Any model implementing the Mosaic API can be written to a HDF5 file,
but will be implicitly converted to a :mod:`mosaic.array_model`
representation.  Data items read from a HDF5 file are always
represented by an :mod:`mosaic.array_model`.

The HDF5 storage layout is designed for clarity rather than
simplicity, the files should be as self-documenting as possible.

A common pattern for generating an HDF5 file is:

::

    with HDF5Store('molecule.h5', 'w') as writer:
       writer.store("universe", universe)
       writer.store("configuration1", configuration1)
       writer.store("configuration2", configuration2)

A common pattern for reading from an HDF5 file is:

::

    items = {}
    for id, data in HDF5Store('molecule.h5'):
       items[id] = data

"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import copy
import operator
import weakref

import h5py
import h5py.h5t
import h5py.h5d
import h5py.h5s
import h5py.version
import numpy as N

import mosaic.array_model
import mosaic.api as api
from mosaic.utility import MethodRegister
from mosaic.utility import ascii
from mosaic.utility import isstring
from mosaic.utility import py_str

import sys
del sys


if h5py.version.hdf5_version_tuple < (1, 8, 7):
    raise ValueError("HDF5 version >= 1.8.7 required")


def create_dataset(group, path, data):
    if isstring(data):
        ds = group.require_dataset(path, shape=(),
                                   dtype=h5py.special_dtype(vlen=bytes))
        ds[...] = data
    else:
        ds = group.require_dataset(path, shape=data.shape, dtype=data.dtype)
        if 0 not in data.shape:
            ds[...] = data
    return ds


class HDF5Store(object):

    """I/O handler for HDF5 files

    This class handles the translation of references between Mosaic
    data items in memory and HDF5 object references in a file. When
    writing data, references are allowed only to data items that have
    been written earlier using the same HDF5Store instance.

    An HDF5Store can be used like a file, in which case it must be
    closed after the last ``store`` or ``retrieve`` operation, or as a
    context manager in a with-statement. For reading, it also behaves
    like an iterator yielding ``(name, data_item)`` pairs.
    """

    def __init__(self, hdf5_group_or_filename, file_mode=None):
        """
        :param hdf5_group_or_filename: a group in an open HDF5 file,
                                       or a string interpreted as a file name
        :type hdf5_group_or_filename: str or ``h5py.Group``
        :param file_mode: ``'r'`` or ``'w'``, used only with a filename argument
        :type file_mode:  str
        """
        if isstring(hdf5_group_or_filename):
            if file_mode is None:
                file_mode = 'r'
            self.root = h5py.File(hdf5_group_or_filename, file_mode)
            self._close = True
        else:
            api.validate_value(file_mode, [None], "file_mode")
            self.root = hdf5_group_or_filename
            self._close = False
        self._path_map = weakref.WeakKeyDictionary()
        self._data_map = weakref.WeakValueDictionary()
        self._factory = mosaic.array_model.Factory()

    def close(self):
        if self._close:
            self.root.file.close()
        self.root = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def flush(self):
        self.root.file.flush()

    storage_handler = MethodRegister()

    def store(self, path, data):
        """
        :param path: a HDF5 path, relative to the group used by the HDF5Store
        :type path:  str
        :param data: a Mosaic data item
        :type data:  :class:`mosaic.api.MosaicDataItem`
        """
        api.validate_type(path, str, "path")
        if path[0] == "/":
            raise ValueError("absolute paths not allowed")
        api.validate_type(data, api.MosaicDataItem, "data")
        handler = self.storage_handler[type(data)]
        if handler is None:
            raise TypeError("Storage of %s not yet implemented"
                            % str(type(data)))
        handler(self, path, data)
        self._register_data_item(path, data)

    def retrieve(self, path_or_node):
        """
        :param path_or_node: a HDF5 path, relative to the group used by
                             the HDF5Store, or an HDF5 node
        :type path_or_node:  str or h5py.Node
        :returns: a Mosaic data item
        :rtype:  :class:`mosaic.api.MosaicDataItem`
        """
        if isinstance(path_or_node, h5py.h5r.Reference):
            path_or_node = self.root[path_or_node]
        if isinstance(path_or_node, str):
            path = path_or_node
        else: # Group, Dataset, or classes with the same interface
            path = path_or_node.name
        if path[0] == '/':
            if path.startswith(self.root.name):
                path = path[len(self.root.name):]
                if path[0] == '/':
                    path = path[1:]
            else:
                raise ValueError("path name not inside root group")
        data = self._get_data(path)
        if data is not None:
            return data
        node = self.root[path]
        label = py_str(self._read_stamp(node))
        if label not in self._node_types:
            raise ValueError("Undefined node type %s" % label)
        handler = self.storage_handler[label]
        if handler is None:
            raise ValueError("Retrieval of node type %s not yet implemented"
                             % label)
        data = handler(self, node)
        self._register_data_item(path, data)
        return data

    _node_types = ['universe', 'configuration', 'property', 'label',
                   'selection']

    def __iter__(self):
        # First pass: universes
        for name in self.root:
            node = self.root[name]
            if node.attrs.get('DATA_MODEL', None) == ascii('MOSAIC') \
               and node.attrs['MOSAIC_DATA_TYPE'] == ascii('universe'):
                yield name, self.retrieve(node)
        # Second pass: anything else
        for name in self.root:
            node = self.root[name]
            if node.attrs.get('DATA_MODEL', None) == ascii('MOSAIC') \
               and node.attrs['MOSAIC_DATA_TYPE'] != ascii('universe'):
                yield name, self.retrieve(node)

    def recursive_iterator(self, group=None):
        if group is None:
            group = self.root
        # First pass: universes
        for name in group:
            node = group[name]
            if node.attrs.get('DATA_MODEL', None) == ascii('MOSAIC') \
               and node.attrs['MOSAIC_DATA_TYPE'] == ascii('universe'):
                yield name, self.retrieve(node)
            if isinstance(node, h5py.Group):
                for gname, gnode in self.recursive_iterator(node):
                    yield name+'/'+gname, gnode
        # Second pass: anything else
        for name in group:
            node = group[name]
            if node.attrs.get('DATA_MODEL', None) == ascii('MOSAIC') \
               and node.attrs['MOSAIC_DATA_TYPE'] != ascii('universe'):
                yield name, self.retrieve(node)
            if isinstance(node, h5py.Group):
                for gname, gnode in self.recursive_iterator(node):
                    yield name+'/'+gname, gnode

    def _stamp(self, node, label):
        node.attrs['DATA_MODEL'] = ascii('MOSAIC')
        node.attrs['DATA_MODEL_MAJOR_VERSION'] = api.MOSAIC_VERSION[0]
        node.attrs['DATA_MODEL_MINOR_VERSION'] = api.MOSAIC_VERSION[1]
        node.attrs['MOSAIC_DATA_TYPE'] = ascii(label)

    def _read_stamp(self, node):
        if node.attrs['DATA_MODEL'] != ascii('MOSAIC'):
            raise ValueError("Node %s not an Mosaic data item" % str(node))
        version = (node.attrs['DATA_MODEL_MAJOR_VERSION'],
                   node.attrs['DATA_MODEL_MINOR_VERSION'])
        if version[0] > api.MOSAIC_VERSION[0] \
           or version[1] > api.MOSAIC_VERSION[1]:
            raise ValueError("HDF5 data is for version %d.%d, "
                             "software is version %d.%d"
                             % (version + api.MOSAIC_VERSION))
        return node.attrs['MOSAIC_DATA_TYPE']

    def _register_data_item(self, path, data_item):
        self._data_map[path] = data_item
        self._path_map[data_item] = path

    def _get_path(self, data_item):
        return self._path_map.get(data_item, None)

    def _get_data(self, path):
        return self._data_map.get(path, None)

    #
    # Store data in HDF5 file
    #
    @storage_handler(api.MosaicUniverse)
    def _store_universe(self, path, universe):
        group = self.root.require_group(path)
        try:
            # ActivePapers support
            group.mark_as_data_item()
        except AttributeError:
            pass
        self._stamp(group, 'universe')

        # Convert the universe to an array_model
        universe = self._factory(universe)

        # Create HDF5 datasets
        create_dataset(group, 'cell_shape',
                       data=universe.cell_shape)
        create_dataset(group, 'symmetry_transformations',
                       data=N.array(universe._symmetry_transformations))
        create_dataset(group, 'convention',
                       data=universe.convention)
        create_dataset(group, 'fragments', data=universe._fragment_array)
        create_dataset(group, 'atoms', data=universe._atom_array)
        create_dataset(group, 'bonds', data=universe._bond_array)
        create_dataset(group, 'molecules', data=universe._molecule_array)
        if len(universe._polymer_array) > 0:
            create_dataset(group, 'polymers', data=universe._polymer_array)
        symbols = universe._symbols
        symbol_table = group.require_dataset('symbols',
                                shape=(len(symbols),),
                                dtype=h5py.special_dtype(vlen=bytes))
        for i, s in enumerate(symbols):
            symbol_table[i] = ascii(s)

    @storage_handler(api.MosaicConfiguration)
    def _store_configuration(self, path, configuration):
        universe = configuration.universe
        universe_path = self._get_path(universe)
        if universe_path is None:
            raise IOError("universe must be stored first")
        group = self.root.require_group(path)
        try:
            # ActivePapers support
            group.mark_as_data_item()
        except AttributeError:
            pass
        self._stamp(group, 'configuration')
        group.attrs['universe'] = self.root[universe_path].ref
        ncellparams = N.multiply.reduce(universe.cell_parameter_array_shape)
        if ncellparams > 0:
            create_dataset(group, 'cell_parameters',
                           data=configuration.cell_parameters)
        arr = N.ascontiguousarray(configuration.positions)
        positions = group.require_dataset('positions',
                                          shape=(len(arr),),
                                          dtype='(3,)f%d' % arr.dtype.itemsize)
        # There doesn't seem to be any way to write this array
        # using high-level operations, so we use the low-level access.
        mtype = h5py.h5t.py_create(positions.id.dtype)
        mspace = h5py.h5s.create_simple(positions.shape)
        fspace = positions.id.get_space()
        positions.id.write(mspace, fspace, arr, mtype)

    @storage_handler(api.MosaicProperty)
    def _store_property(self, path, property):
        universe = property.universe
        universe_path = self._get_path(universe)
        if universe_path is None:
            raise IOError("universe must be stored first")
        element_shape = property.data.shape[1:]
        if element_shape:
            dtype = N.dtype((property.data.dtype, element_shape))
        else:
            dtype = property.data.dtype
        arr = N.ascontiguousarray(property.data)
        ds = self.root.require_dataset(path,
                                       shape=(len(arr),),
                                       dtype=dtype)
        self._stamp(ds, 'property')
        ds.attrs['universe'] = self.root[universe_path].ref
        ds.attrs['name'] = ascii(property.name)
        ds.attrs['units'] = ascii(property.units)
        ds.attrs['property_type'] = ascii(property.type)
        # There doesn't seem to be any way to write this array
        # using high-level operations, so we use the low-level access.
        mtype = h5py.h5t.py_create(ds.id.dtype)
        mspace = h5py.h5s.create_simple(ds.shape)
        fspace = ds.id.get_space()
        ds.id.write(mspace, fspace, arr, mtype)

    @storage_handler(api.MosaicLabel)
    def _store_label(self, path, label):
        universe = label.universe
        label = self._factory(label)
        universe_path = self._get_path(universe)
        if universe_path is None:
            raise IOError("universe must be stored first")
        ds = create_dataset(self.root, path, data=label._string_array)
        self._stamp(ds, 'label')
        ds.attrs['universe'] = self.root[universe_path].ref
        ds.attrs['name'] = ascii(label.name)
        ds.attrs['label_type'] = ascii(label.type)

    @storage_handler(api.MosaicSelection)
    def _store_selection(self, path, selection):
        universe = selection.universe
        universe_path = self._get_path(universe)
        if universe_path is None:
            raise IOError("universe must be stored first")
        ds = create_dataset(self.root, path,
                            N.ascontiguousarray(selection.indices))
        ds.attrs['universe'] = self.root[universe_path].ref
        ds.attrs['selection_type'] = ascii(selection.type)
        self._stamp(ds, 'selection')

    #
    # Retrieve data from HDF5 file
    #
    @storage_handler('universe')
    def _retrieve_universe(self, node):
        # Work around a restriction/bug in h5py that prevents reading
        # arrays of length zero.
        if node['symmetry_transformations'].shape[0] == 0:
            st = N.zeros((0,), dtype=node['symmetry_transformations'].dtype)
        else:
            st = node['symmetry_transformations'][...]
        return mosaic.array_model.Universe.from_arrays(
                        node['atoms'][...],
                        node['bonds'][...],
                        node['fragments'][...],
                        node['polymers'][...]
                            if 'polymers' in node else None,
                        node['molecules'][...],
                        tuple(py_str(s) for s in node['symbols'][...]),
                        st,
                        py_str(node['cell_shape'][()]),
                        py_str(node['convention'][()]))

    @storage_handler('configuration')
    def _retrieve_configuration(self, node):
        universe = self.retrieve(node.attrs['universe'])
        positions = node['positions'][...]
        if 'cell_parameters' in node:
            ds = node['cell_parameters']
            cell_parameters = N.empty(ds.shape, ds.dtype)
            cell_parameters[...] = ds[...]
        else:
            assert N.product(universe.cell_parameter_array_shape) == 0
            cell_parameters = N.empty(universe.cell_parameter_array_shape,
                                      positions.dtype)
        return mosaic.array_model.Configuration(universe, positions,
                                                cell_parameters)

    @storage_handler('property')
    def _retrieve_property(self, node):
        universe = self.retrieve(node.attrs['universe'])
        data = node[...]
        name = py_str(node.attrs['name'])
        units = py_str(node.attrs['units'])
        property_type = py_str(node.attrs['property_type'])
        return mosaic.array_model.Property(universe, name, units,
                                           data, property_type)

    @storage_handler('label')
    def _retrieve_label(self, node):
        universe = self.retrieve(node.attrs['universe'])
        data = node[...]
        name = py_str(node.attrs['name'])
        label_type = py_str(node.attrs['label_type'])
        return mosaic.array_model.Label.from_arrays(universe, name,
                                                    data, label_type)

    @storage_handler('selection')
    def _retrieve_selection(self, node):
        universe = self.retrieve(node.attrs['universe'])
        indices = node[...]
        selection_type = py_str(node.attrs['selection_type'])
        return mosaic.array_model.Selection(universe, indices,
                                            selection_type)
