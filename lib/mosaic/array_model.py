# -*- coding: utf-8 -*-
"""Array model

.. moduleauthor:: Konrad Hinsen

An array model stores all data in NumPy arrays, in particular the tree
structure defining the fragments in universes. The array model is an
in-memory version of the HDF5 storage layout. It is used for HDF5 I/O,
but can also be of use for interfacing with C or Fortran code.

"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import operator

import numpy as N

import mosaic.api as api
from mosaic.utility import SymbolDict
from mosaic.utility import AbstractFactory
from mosaic.utility import MethodRegister
from mosaic.utility import isascii
from mosaic.utility import py_str

# dtypes for the various structured arrays
fragment_array_dtype = N.dtype([('parent_index', N.uint32),
                                ('label_symbol_index', N.uint32),
                                ('species_symbol_index', N.uint32),
                                # redundant:
                                ('number_of_fragments', N.uint32)])

atom_array_dtype = N.dtype([('parent_index', N.uint32),
                            ('label_symbol_index', N.uint32),
                            ('type_symbol_index', N.uint32),
                            ('name_symbol_index', N.uint32),
                            ('number_of_sites', N.uint32)])

bond_array_dtype = N.dtype([('atom_index_1', N.uint32),
                            ('atom_index_2', N.uint32),
                            ('bond_order_symbol_index', N.uint32)])

molecule_array_dtype = N.dtype([('fragment_index', N.uint32),
                                ('number_of_copies', N.uint32),
                                # redundant:
                                ('first_atom_index', N.uint32),
                                ('number_of_atoms', N.uint32),
                                ('first_bond_index', N.uint32),
                                ('number_of_bonds', N.uint32),
                                ('first_site_index', N.uint32),
                                ('number_of_sites', N.uint32)])

polymer_array_dtype = N.dtype([('fragment_index', N.uint32),
                               ('polymer_type_symbol_index', N.uint32)])

st_array_dtype = N.dtype([("rotation", N.float64, (3, 3)),
                          ("translation", N.float64, (3,))])


class Universe(api.MosaicUniverse):

    def __init__(self, universe):
        def count_fragments(f):
            return 1 + sum(count_fragments(ff) for ff in f.fragments)

        mol_fragments = tuple(fragment for fragment, _ in universe.molecules)
        natoms = sum(f.number_of_atoms for f in mol_fragments)
        nfragments = sum(count_fragments(f) for f in mol_fragments)
        nbonds = sum(f.number_of_bonds for f in mol_fragments)
        nmols = len(mol_fragments)

        symbol_id = SymbolDict()

        self._atom_array = N.empty((natoms,), atom_array_dtype)
        self._fragment_array = N.empty((nfragments + 1,),
                                       fragment_array_dtype)
        self._fragment_array[0] = (0, symbol_id[''], symbol_id[''], 1)
        self._bond_array = N.empty((nbonds,), bond_array_dtype)
        self._molecule_array = N.empty((nmols,), molecule_array_dtype)

        def store_fragment(f, iparent, ifrag, iatom, ibond, isite,
                           polymers, path, atom_index):
            if f.is_polymer:
                polymers.append((ifrag, symbol_id[f.polymer_type]))
            self._fragment_array[ifrag]['parent_index'] = iparent
            self._fragment_array[ifrag]['label_symbol_index'] = \
                                                         symbol_id[f.label]
            self._fragment_array[ifrag]['species_symbol_index'] = \
                                                         symbol_id[f.species]
            iparent = ifrag
            ifrag += 1
            for ff in f.fragments:
                ifrag, iatom, ibond, isite = \
                       store_fragment(ff, iparent, ifrag, iatom, ibond, isite,
                                      polymers, '.'.join([path, ff.label]),
                                      atom_index)
            for a in f.atoms:
                atom_index['.'.join([path, a.label])] = iatom
                self._atom_array[iatom]['parent_index'] = iparent
                self._atom_array[iatom]['label_symbol_index'] = \
                                                symbol_id[a.label]
                self._atom_array[iatom]['type_symbol_index'] = \
                                                symbol_id[a.type]
                self._atom_array[iatom]['name_symbol_index'] = \
                                                symbol_id[a.name]
                self._atom_array[iatom]['number_of_sites'] = a.number_of_sites
                iatom += 1
                isite += a.number_of_sites
            bonds = []
            for a1, a2, order in f.bonds:
                i1 = atom_index['.'.join([path, a1])]
                i2 = atom_index['.'.join([path, a2])]
                if i1 > i2:
                    i1, i2 = i2, i1
                bonds.append((i1, i2, symbol_id[order]))
            bonds.sort()
            self._bond_array[ibond:ibond + len(bonds)] = bonds
            ibond += len(bonds)

            # ifrag now points to the first fragment that is not
            # a child of the current one. ifrag-iparent is thus
            # the number of fragments that are part of the current one.
            self._fragment_array[iparent]['number_of_fragments'] = \
                                                          ifrag - iparent
            return ifrag, iatom, ibond, isite

        ifrag = 1
        iatom = 0
        ibond = 0
        isite = 0
        imol = 0
        polymers = []
        for fragment, count in universe.molecules:
            self._molecule_array[imol]['fragment_index'] = ifrag
            self._molecule_array[imol]['number_of_copies'] = count
            self._molecule_array[imol]['first_atom_index'] = iatom
            self._molecule_array[imol]['first_bond_index'] = ibond
            self._molecule_array[imol]['first_site_index'] = isite
            ifrag, iatom, ibond, isite = \
                   store_fragment(fragment, 0, ifrag, iatom, ibond, isite,
                                  polymers, '', {})
            self._molecule_array[imol]['number_of_atoms'] = \
                       iatom - self._molecule_array[imol]['first_atom_index']
            self._molecule_array[imol]['number_of_bonds'] = \
                       ibond - self._molecule_array[imol]['first_bond_index']
            self._molecule_array[imol]['number_of_sites'] = \
                       isite - self._molecule_array[imol]['first_site_index']
            imol += 1
        if len(polymers) > 0:
            self._polymer_array = N.array(polymers, dtype=polymer_array_dtype)
        else:
            self._polymer_array = N.zeros((0,), dtype=polymer_array_dtype)
        assert ifrag == nfragments + 1
        assert iatom == natoms
        assert ibond == nbonds

        # Other universe attributes
        self._cell_shape = universe.cell_shape
        self._convention = universe.convention
        # Numpy raises an exception if the transformation list is a
        # Python tuple AND it happens to be empty, so we convert
        # it to a list first.
        self._symmetry_transformations = \
                N.array(list(universe.symmetry_transformations),
                        dtype=st_array_dtype)
        symbols = list(symbol_id.items())
        symbols.sort(key=operator.itemgetter(1))
        self._symbols = tuple(s for s, i in symbols)

    # Create from existing arrays, typically read from HDF5
    @classmethod
    def from_arrays(cls, atom_array, bond_array, fragment_array,
                    polymer_array, molecule_array, symbols,
                    symmetry_transformations, cell_shape, convention):
        assert atom_array.dtype == atom_array_dtype
        assert bond_array.dtype == bond_array_dtype
        assert fragment_array.dtype == fragment_array_dtype
        if polymer_array is None:
            polymer_array = N.zeros((0,), dtype=polymer_array_dtype)
        assert polymer_array.dtype == polymer_array_dtype
        assert molecule_array.dtype == molecule_array_dtype
        assert isinstance(symbols, tuple)
        assert symmetry_transformations.dtype == st_array_dtype
        assert isascii(cell_shape)
        assert isascii(convention)

        universe = cls.__new__(cls)
        universe._atom_array = atom_array
        universe._bond_array = bond_array
        universe._fragment_array = fragment_array
        universe._polymer_array = polymer_array
        universe._molecule_array = molecule_array
        universe._symbols = tuple(symbols)
        universe._symmetry_transformations = symmetry_transformations
        universe._cell_shape = cell_shape
        universe._convention = convention
        return universe

    # API implementation
    @property
    def cell_shape(self):
        return self._cell_shape

    @property
    def symmetry_transformations(self):
        return tuple(tuple(t) for t in self._symmetry_transformations)

    @property
    def convention(self):
        return self._convention

    @property
    def molecules(self):
        return tuple((Fragment(self, m['fragment_index']),
                      int(m['number_of_copies']))
                     for m in self._molecule_array)

    # More efficient implementations of the number properties

    @property
    def number_of_molecules(self):
        return int(N.sum(self._molecule_array['number_of_copies']))

    @property
    def number_of_atoms(self):
        return int(N.sum(self._molecule_array['number_of_copies']
                         * self._molecule_array['number_of_atoms']))

    @property
    def number_of_sites(self):
        return int(N.sum(self._molecule_array['number_of_copies']
                         * self._molecule_array['number_of_sites']))

    @property
    def number_of_bonds(self):
        return int(N.sum(self._molecule_array['number_of_copies']
                         * self._molecule_array['number_of_bonds']))

    @property
    def number_of_template_atoms(self):
        return int(N.sum(self._molecule_array['number_of_atoms']))

    @property
    def number_of_template_sites(self):
        return int(N.sum(self._molecule_array['number_of_sites']))

    # More efficient implementation of various methods

    def recursive_atom_iterator(self):
        for mol in self._molecule_array:
            first = mol['first_atom_index']
            natoms = mol['number_of_atoms']
            for _ in range(mol['number_of_copies']):
                for iatom in range(first, first+natoms):
                    yield Atom(self, iatom)

    def bond_index_array(self):
        atoms = 0
        fatoms = 0
        bonds = []
        for mol in self._molecule_array:
            first = mol['first_bond_index']
            nbonds = mol['number_of_bonds']
            natoms = mol['number_of_atoms']
            a1 = self._bond_array[first:first+nbonds]['atom_index_1']
            a2 = self._bond_array[first:first+nbonds]['atom_index_2']
            fbonds = N.concatenate([a1[:, N.newaxis], a2[:, N.newaxis]],
                                    axis=1)
            for _ in range(mol['number_of_copies']):
                bonds.append(fbonds - fatoms + atoms)
                atoms += natoms
            fatoms += natoms
        return N.concatenate(bonds, axis=0)

class Fragment(api.MosaicFragment):

    def __init__(self, universe, ifrag):
        self._universe = universe
        self._ifrag = ifrag

    def __repr__(self):
        return "Fragment('%s', '%s')" % (self.label, self.species)

    @property
    def parent(self):
        ifrag = self._universe._fragment_array[self._ifrag]['parent_index']
        if ifrag == 0:
            return None
        else:
            return Fragment(self._universe, ifrag)

    # API implementation
    @property
    def label(self):
        isym = \
            self._universe._fragment_array[self._ifrag]['label_symbol_index']
        return self._universe._symbols[isym]

    @property
    def species(self):
        isym = \
            self._universe._fragment_array[self._ifrag]['species_symbol_index']
        return self._universe._symbols[isym]

    @property
    def fragments(self):
        mask = self._universe._fragment_array['parent_index'] == self._ifrag
        return [Fragment(self._universe, ifrag)
                for ifrag in N.repeat(N.arange(len(mask)), mask)]

    @property
    def atoms(self):
        mask = self._universe._atom_array['parent_index'] == self._ifrag
        return [Atom(self._universe, iatom)
                for iatom in N.repeat(N.arange(len(mask)), mask)]

    @property
    def bonds(self):

        universe = self._universe

        def ancestors(ip):
            anc = []
            while True:
                if ip == 0:
                    return None
                if ip == self._ifrag:
                    break
                anc.insert(0, ip)
                ip = universe._fragment_array[ip]['parent_index']
            return anc

        bonds = []
        for ia1, ia2, iorder in universe._bond_array[
                                ['atom_index_1', 'atom_index_2',
                                 'bond_order_symbol_index']]:
            anc1 = ancestors(universe._atom_array[ia1]['parent_index'])
            if anc1 is None:
                continue
            anc2 = ancestors(universe._atom_array[ia2]['parent_index'])
            if anc2 is None:
                continue
            if len(anc1) == 0 or len(anc2) == 0 or anc1[0] != anc2[0]:
                # bond belongs to the current fragment
                p1 = [universe._fragment_array[if1]['label_symbol_index']
                      for if1 in anc1]
                p1.append(universe._atom_array[ia1]['label_symbol_index'])
                p2 = [universe._fragment_array[if2]['label_symbol_index']
                      for if2 in anc2]
                p2.append(universe._atom_array[ia2]['label_symbol_index'])
                bonds.append(('.'.join(universe._symbols[l] for l in p1),
                              '.'.join(universe._symbols[l] for l in p2),
                              universe._symbols[iorder]))

        return tuple(bonds)

    @property
    def is_polymer(self):
        return N.sum(self._universe._polymer_array['fragment_index']
                     == self._ifrag) == 1

    @property
    def polymer_type(self):
        mask = self._universe._polymer_array['fragment_index'] == self._ifrag
        if mask.any():
            index = self._universe._polymer_array['polymer_type_symbol_index']
            isym = N.repeat(index, mask)[0]
            return self._universe._symbols[isym]
        else:
            raise ValueError("polymer_type not defined in a non-polymer")

    # More efficient implementation of the number properties

    @property
    def number_of_atoms(self):
        ifrag = self._ifrag
        nfrag = ifrag + \
                self._universe._fragment_array[ifrag]['number_of_fragments']
        return int(N.sum(N.logical_and(
                           self._universe._atom_array['parent_index'] >= ifrag,
                           self._universe._atom_array['parent_index'] < nfrag),
                         dtype=N.int64))

    @property
    def number_of_sites(self):
        ifrag = self._ifrag
        nfrag = ifrag + \
                self._universe._fragment_array[ifrag]['number_of_fragments']
        mask = N.logical_and(
                    self._universe._atom_array['parent_index'] >= ifrag,
                    self._universe._atom_array['parent_index'] < nfrag)
        return int(N.sum(N.repeat(self._universe._atom_array['number_of_sites'],
                                  mask),
                         dtype=N.int64))

    @property
    def number_of_bonds(self):
        ifrag = self._ifrag
        nfrag = ifrag + \
                self._universe._fragment_array[ifrag]['number_of_fragments']
        a1 = N.take(self._universe._atom_array['parent_index'],
                    self._universe._bond_array['atom_index_1'])
        a2 = N.take(self._universe._atom_array['parent_index'],
                    self._universe._bond_array['atom_index_2'])
        return int(N.sum(N.logical_and(N.logical_and(a1 >= ifrag, a1 < nfrag),
                                       N.logical_and(a2 >= ifrag, a2 < nfrag)),
                         dtype=N.int64))


class Atom(api.MosaicAtom):

    def __init__(self, universe, iatom):
        self._universe = universe
        self._iatom = iatom

    def __repr__(self):
        return "Atom('%s')" % self.label

    @property
    def parent(self):
        ifrag = self._universe._atom_array[self._iatom]['parent_index']
        assert ifrag > 0
        return Fragment(self._universe, ifrag)

    # API implementation
    @property
    def label(self):
        isym = self._universe._atom_array[self._iatom]['label_symbol_index']
        return self._universe._symbols[isym]

    @property
    def type(self):
        isym_type = \
            self._universe._atom_array[self._iatom]['type_symbol_index']
        return self._universe._symbols[isym_type]

    @property
    def name(self):
        isym_name = \
            self._universe._atom_array[self._iatom]['name_symbol_index']
        return self._universe._symbols[isym_name]

    @property
    def number_of_atoms(self):
        return 1

    @property
    def number_of_sites(self):
        return int(self._universe._atom_array[self._iatom]['number_of_sites'])


class Property(api.MosaicProperty):

    def __init__(self, universe, name, units, data, property_type):
        api.validate_type(universe, Universe, "universe")
        api.validate_label(name, "name")
        api.validate_units(units, "units")
        api.validate_array(data, None, self._allowed_dtypes, "data")
        api.validate_value(property_type, self._allowed_types, "property_type")
        self._universe = universe
        self._name = name
        self._units = units
        self._data = data
        self._type = property_type

    # Data Model API implementation
    @property
    def universe(self):
        return self._universe

    @property
    def name(self):
        return self._name

    @property
    def units(self):
        return self._units

    @property
    def data(self):
        return self._data

    @property
    def type(self):
        return self._type


class Label(api.MosaicLabel):

    def __init__(self, universe, name, strings, label_type):
        api.validate_type(universe, Universe, "universe")
        api.validate_label(name, "name")
        api.validate_sequence(strings, str, "strings")
        api.validate_value(label_type, self._allowed_types, "label_type")
        self._universe = universe
        self._name = name
        self._string_array = N.array(list('\0'.join(strings) + '\0'),
                                     dtype='S')
        self._type = label_type

    @classmethod
    def from_arrays(cls, universe, name, array, label_type):
        label = cls.__new__(cls)
        label._universe = universe
        label._name = name
        label._string_array = array
        label._type = label_type
        return label

    # Data Model API implementation
    @property
    def universe(self):
        return self._universe

    @property
    def name(self):
        return self._name

    @property
    def strings(self):
        eos = N.repeat(N.arange(self._string_array.shape[0]),
                       self._string_array.view(dtype=N.uint8) == 0)
        strings = []
        f = 0
        for l in eos:
            strings.append(py_str(self._string_array[f:l].tostring()))
            f = l + 1
        return strings

    @property
    def type(self):
        return self._type


class Configuration(api.MosaicConfiguration):

    def __init__(self, universe, positions, cell_parameters):
        api.validate_type(universe, Universe, "universe")
        api.validate_array(positions, (universe.number_of_sites, 3),
                           self._allowed_dtypes, "positions")
        api.validate_array(cell_parameters,
                           universe._cell_parameter_array_shapes
                                                      [universe.cell_shape],
                           [positions.dtype], "cell_parameters")
        self._universe = universe
        self._positions = positions
        self._cell_parameters = cell_parameters

    # Data Model API implementation
    @property
    def universe(self):
        return self._universe

    @property
    def positions(self):
        return self._positions

    @property
    def cell_parameters(self):
        return self._cell_parameters


class Selection(api.MosaicSelection):

    def __init__(self, universe, indices, selection_type):
        api.validate_type(universe, Universe, "universe")
        api.validate_array(indices, None,
                           [N.uint8, N.uint16, N.uint32, N.uint64], "indices")
        api.validate_value(selection_type, self._allowed_types,
                           "selection_type")
        self._universe = universe
        self._indices = indices
        self._type = selection_type

    # Data Model API implementation
    @property
    def universe(self):
        return self._universe

    @property
    def indices(self):
        return self._indices

    @property
    def type(self):
        return self._type


class Factory(AbstractFactory):

    handler = MethodRegister()

    @handler(api.MosaicUniverse)
    def _make_universe(self, universe):
        return Universe(universe)

    @handler(api.MosaicProperty)
    def _make_property(self, property):
        universe = self(property.universe)
        return Property(universe, property.name, property.units,
                        property.data, property.type)

    @handler(api.MosaicLabel)
    def _make_label(self, label):
        universe = self(label.universe)
        return Label(universe, label.name, label.strings, label.type)

    @handler(api.MosaicConfiguration)
    def _make_configuration(self, configuration):
        universe = self(configuration.universe)
        return Configuration(universe,
                             configuration.positions,
                             configuration.cell_parameters)

    @handler(api.MosaicSelection)
    def _make_selection(self, selection):
        universe = self(selection.universe)
        return Selection(universe, selection.indices, selection.type)
