# -*- coding: utf-8 -*-
"""Immutable model

.. moduleauthor:: Konrad Hinsen

This experimental module implements the Mosaic API in terms of
immutable data structures as implemented in ImmutablePy.  The use of
immutable data structures avoids many nasty bugs, such as modifying a
universe after creating a configuration object that refer to it,
invalidating the configuration. However, the immutable data structures
from ImmutablePy are not very efficient at the moment.

"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import copy
import collections
import itertools as IT
import types

import numpy as N

from immutable import Immutable, ImmutableTuple, ImmutableDict, immutable
import immutable.np as IN

import mosaic.api as api
from mosaic.utility import AbstractFactory
from mosaic.utility import MethodRegister
from mosaic.utility import isascii
from mosaic.utility import uint_for_max_value


# Look up template objects in the class' instance cache.
# If not found, create new object.

def _lookup(klass, key):
    key = immutable(key)
    try:
        # Under Python 3, this raises a TypeError
        # for ImmutableNDArray.
        return klass._instances[key]
    except (KeyError, TypeError):
        obj = klass(key)
        try:
            klass._instances[key] = obj
        except TypeError:
            pass
        return obj


# Atom descriptors

class AtomDescriptor(Immutable):

    def __init__(self, name=""):
        self._validate_name(name)
        self._name = name

    def _validate_name(cls, name):
        api.validate_label(name, "name")

    def __repr__(self):
        return "%s('%s')" % (self.__class__.__name__, self._name)

    @property
    def name(self):
        return self._name


class Element(AtomDescriptor):

    _instances = {}

    def _validate_name(cls, name):
        api.MosaicAtom.validate_element_name(name)

    def __repr__(self):
        return "Element('%s')" % self._name

    # Data Model API implementation
    @property
    def type(self):
        return "element"


class CGParticle(AtomDescriptor):

    _instances = {}

    def __repr__(self):
        if self._name:
            return "CGParticle('%s')" % self._name
        else:
            return "CGParticle()"

    @property
    def type(self):
        return "cgparticle"


class Dummy(AtomDescriptor):

    _instances = {}

    def __repr__(self):
        if self._name:
            return "Dummy('%s')" % self._name
        else:
            return "Dummy()"

    @property
    def type(self):
        return "dummy"


class Unknown(AtomDescriptor):

    _instances = {}

    def __repr__(self):
        if self._name:
            return "Unknown('%s')" % self._name
        else:
            return "Unknown()"

    @property
    def type(self):
        return ""


element = lambda name: _lookup(Element, name)
cgparticle = lambda name: _lookup(CGParticle, name)
dummy = lambda name="": _lookup(Dummy, name)
unknown = lambda name="": _lookup(Unknown, name)


# Atom and fragment templates

class Atom(Immutable):

    _instances = {}

    def __init__(self, args):
        descriptor, nsites = args
        api.validate_type(descriptor, AtomDescriptor, "atom descriptor")
        self.descriptor = descriptor
        api.validate_type(nsites, int, "number of sites")
        if nsites < 1:
            raise ValueError("number of sites must be >= 1")
        self.nsites = nsites

    def __repr__(self):
        return "Atom(%s, %d)" % (repr(self.descriptor), self.nsites)

    # API-compatible properties

    @property
    def number_of_atoms(self):
        return 1

    @property
    def number_of_sites(self):
        return self.nsites

class AtomRef(Immutable, api.MosaicAtom):

    def __init__(self, label, atom):
        self._label = label
        self._atom = atom

    # Data Model API implementation

    @property
    def label(self):
        return self._label

    @property
    def name(self):
        return self._atom.descriptor.name

    @property
    def type(self):
        return self._atom.descriptor.type

    @property
    def number_of_sites(self):
        return self._atom.nsites

atom = lambda descriptor, nsites=1: _lookup(Atom, (descriptor, nsites))


class Fragment(Immutable, collections.Mapping):

    _instances = {}

    def __init__(self, args):
        species, fragments, atoms, bonds = args
        api.validate_label(species, "species")
        self.species = species

        labels = set()

        api.validate_sequence(fragments, ImmutableTuple, "fragments",
                              ((lambda p: len(p) == 2, "must have length 2"),))
        for label, fragment in fragments:
            api.validate_label(label, "fragment label")
            if label in labels:
                raise ValueError("label %s occurs more than once" % label)
            labels.add(label)
            api.validate_type(fragment, Fragment, "fragment template")
        self.fragments = fragments

        api.validate_sequence(atoms, ImmutableTuple, "atoms",
                              ((lambda p: len(p) == 2, "must have length 2"),))
        for label, atom in atoms:
            api.validate_label(label, "atom label")
            if label in labels:
                raise ValueError("label %s occurs more than once" % label)
            labels.add(label)
            api.validate_type(atom, Atom, "atom")
        self.atoms = atoms

        self.attrs = ImmutableDict(IT.chain(fragments, atoms))

        api.validate_sequence(bonds, tuple, "bonds",
                              ((lambda p: len(p) == 3, "must have length 3"),))
        for a1, a2, order in bonds:
            self._validate_bond_atom(a1)
            self._validate_bond_atom(a2)
            if a1.split('.')[0] == a2.split('.')[0]:
                raise ValueError("bond between %s and %s must be defined "
                                 "in fragment %s" % (a1, a2, a1.split('.')[0]))
            api.validate_value(order, api.MosaicFragment._bond_orders,
                               "bond order")
        self.bonds = bonds

    def _validate_bond_atom(self, atom_label):
        api.validate_type(atom_label, str, "atom label")
        obj = self
        for item in atom_label.split('.'):
            if isinstance(obj, Atom):
                raise ValueError("invalid atom reference "
                                 "%s (refers to child object of an atom)"
                                 % atom_label)
            try:
                obj = obj.attrs[item]
            except KeyError:
                raise ValueError("invalid atom reference %s "
                                 "(child %s not found)" % (atom_label, item))

    def __repr__(self):
        return "Fragment(%s, ...)" % self.species

    def __eq__(self, other):
        return isinstance(other, Fragment) \
               and self.species == other.species \
               and all(s == o
                       for s, o in zip(self.fragments, other.fragments)) \
               and all(s == o
                       for s, o in zip(self.atoms, other.atoms)) \
               and frozenset([(frozenset([a1, a2]), o)
                              for a1, a2, o in self.bonds]) == \
                   frozenset([(frozenset([a1, a2]), o)
                              for a1, a2, o in other.bonds])

    def __hash__(self):
        # Don't use bonds in hash because the bond list can be
        # different for fragments that test equal.
        return hash((self.species, self.fragments, self.atoms))

    # Mapping interface

    def __getitem__(self, item):
        return self.attrs[item]

    def __len__(self):
        return len(self.attrs)

    def __iter__(self):
        return iter(self.attrs)

    # API-compatible properties

    @property
    def number_of_atoms(self):
        return sum(f.number_of_atoms for l, f in self.fragments) \
               + len(self.atoms)

    @property
    def number_of_sites(self):
        return sum(f.number_of_sites for l, f in self.fragments) \
               + sum(a.number_of_sites for l, a in self.atoms)

    @property
    def number_of_bonds(self):
        return sum(f.number_of_bonds for l, f in self.fragments) \
               + len(self.bonds)

class Polymer(Fragment, collections.Sequence):

    _instances = {}

    def __init__(self, args):
        species, fragments, bonds, polymer_type = args
        Fragment.__init__(self, (species, fragments, ImmutableTuple(), bonds))
        self.polymer_type = polymer_type

    def __repr__(self):
        return "Polymer(%s, ...)" % self.species

    # Sequence and mapping interface
    def __getitem__(self, item):
        if isinstance(item, int):
            return self.fragments[int]
        else:
            return Fragment.__getitem__(self, item)

class FragmentRef(Immutable, api.MosaicFragment):

    def __init__(self, label, fragment):
        self._label = label
        self._fragment = fragment

    # Data Model API implementation

    @property
    def label(self):
        return self._label

    @property
    def species(self):
        return self._fragment.species

    @property
    def fragments(self):
        return [FragmentRef(l, f) for l, f in self._fragment.fragments]

    @property
    def atoms(self):
        return [AtomRef(l, a) for l, a in self._fragment.atoms]

    @property
    def bonds(self):
        return self._fragment.bonds

    @property
    def is_polymer(self):
        return isinstance(self._fragment, Polymer)

    @property
    def polymer_type(self):
        if isinstance(self._fragment, Polymer):
            return self._fragment.polymer_type
        else:
            raise ValueError("polymer_type not defined in a non-polymer")

    # Sequence and mapping interface

    def __getitem__(self, item):
        return self._fragment[item]

fragment = lambda species, fragments, atoms, bonds: \
            _lookup(Fragment, (species, fragments, atoms, bonds))

polymer = lambda species, fragments, bonds, polymer_type: \
            _lookup(Polymer, (species, fragments, bonds, polymer_type))


# Universe

class Universe(Immutable, api.MosaicUniverse):

    _instances = {}

    def __init__(self, args):
        cell_shape, molecules, symmetry_transformations, convention \
                = args
        api.validate_type(cell_shape, str, "cell_shape")
        api.validate_value(cell_shape,
                           list(self._cell_parameter_array_shapes.keys()),
                           "cell_shape")
        api.validate_sequence(molecules, ImmutableTuple, "molecules",
                              ((lambda p: len(p) == 3, "must have length 3"),
                               (lambda p: isinstance(p[0], Fragment)
                                          and isascii(p[1])
                                          and isinstance(p[2], int),
                                "elements must be (fragment, label, count) "
                                "triples")))
        api.validate_sequence(symmetry_transformations, ImmutableTuple,
                              "symmetry_transformations",
                              ((lambda p: len(p) == 2, "must have length 2"),
                               (lambda p: hasattr(p[0], 'shape')
                                          and p[0].shape == (3, 3)
                                          and p[0].dtype == N.float64,
                                "rotation matrix must be float64 "
                                "and have shape (3,3)"),
                               (lambda p: hasattr(p[1], 'shape')
                                          and p[1].shape == (3,)
                                          and p[0].dtype == N.float64,
                                "translation vector must be float64 "
                                "and have shape (3,)"),))
        if cell_shape == "infinite" \
           and len(symmetry_transformations) > 0:
            raise ValueError("Symmetry transformations are allowed "
                             "only in periodic universes")
        api.validate_label(convention, "Universe.convention")

        self._cell_shape = cell_shape
        self._molecules = ImmutableTuple(
                ImmutableTuple((FragmentRef(label, fragment), count))
                for fragment, label, count in molecules)
        self._symmetry_transformations = symmetry_transformations
        self._convention = convention

        self._fragments = ImmutableTuple(ImmutableTuple((f, l))
                                         for f, l, c in molecules)
        self._molecule_counts = ImmutableTuple(c for f, l, c in molecules)
        self._atom_counts = ImmutableTuple(f.number_of_atoms
                                           for f, l, c in molecules)
        self._site_counts = ImmutableTuple(f.number_of_sites
                                           for f, l, c in molecules)
        self._bond_counts = ImmutableTuple(f.number_of_bonds
                                           for f, l, c in molecules)

    def __eq__(self, other):
        return isinstance(other, Universe) \
               and self.cell_shape == other.cell_shape \
               and self.symmetry_transformations == \
                       other.symmetry_transformations \
               and self.convention == other.convention \
               and all(c1 == c2
                       for c1, c2 in zip(self._molecule_counts,
                                         other._molecule_counts)) \
               and all(f1 == f2
                       for f1, f2 in zip(self._fragments, other._fragments))

    def __hash__(self):
        # a cheap hash that doesn't look at the fragments
        return hash(self._cell_shape) + len(self._fragments)

    def __repr__(self):
        mol_str = ', '.join("(Fragment('%s'), %d)" % (f.label, c)
                            for f, c in self._molecules)
        return "Universe('%s', [%s])" % (self.cell_shape, mol_str)

    # Data Model API implementation

    @property
    def cell_shape(self):
        return self._cell_shape

    @property
    def symmetry_transformations(self):
        return self._symmetry_transformations

    @property
    def convention(self):
        return self._convention

    @property
    def molecules(self):
        return self._molecules

    # More efficient implementations of the number properties

    @property
    def number_of_molecules(self):
        return sum(self._molecule_counts)

    @property
    def number_of_atoms(self):
        return sum(c * na for c, na in zip(self._molecule_counts,
                                           self._atom_counts))

    @property
    def number_of_sites(self):
        return sum(c * ns for c, ns in zip(self._molecule_counts,
                                           self._site_counts))

    @property
    def number_of_bonds(self):
        return sum(c * nb for c, nb in zip(self._molecule_counts,
                                           self._bond_counts))

    @property
    def number_of_template_atoms(self):
        return sum(self._atom_counts)

    @property
    def number_of_template_sites(self):
        return sum(self._site_counts)

universe = lambda cell_shape, molecules, \
                    symmetry_transformations=(), convention='': \
            _lookup(Universe, (cell_shape, molecules,
                               symmetry_transformations, convention))


# Properties

class Property(Immutable, api.MosaicProperty):

    def __init__(self, universe, name, units, data):
        api.validate_type(universe, Universe, "universe")
        self._universe = universe
        api.validate_label(name, "name")
        self._name = name
        api.validate_units(units, "units")
        self._units = units
        api.validate_array(data, None, self._allowed_dtypes, "data")
        if data.shape[0] != self._number_of_values():
            raise ValueError("data array has incorrect shape")
        self._data = data

    def __repr__(self):
        return "%s()" % (self.__class__.__name__)

    def __eq__(self, other):
        return isinstance(other, Property) \
               and self.type == other.type \
               and self.name == other.name \
               and self.units == other.units \
               and self.universe == other.universe \
               and (self.data == other.data).all()

    def __hash__(self):
        # cheap hash function that doesn't look at the data
        return hash((self.type, self.universe, self.data.shape))

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


class AtomProperty(Property):

    def _number_of_values(self):
        return self._universe.number_of_atoms

    # Data Model API implementation
    @property
    def type(self):
        return 'atom'


class SiteProperty(Property):

    def _number_of_values(self):
        return self._universe.number_of_sites

    # Data Model API implementation
    @property
    def type(self):
        return 'site'


class TemplateAtomProperty(Property):

    def _number_of_values(self):
        return self._universe.number_of_template_atoms

    # Data Model API implementation
    @property
    def type(self):
        return 'template_atom'


class TemplateSiteProperty(Property):

    def _number_of_values(self):
        return self._universe.number_of_template_sites

    # Data Model API implementation
    @property
    def type(self):
        return 'template_site'


# Labels

class Label(Immutable, api.MosaicLabel):

    def __init__(self, universe, name, strings):
        api.validate_type(universe, Universe, "universe")
        self._universe = universe
        api.validate_label(name, "name")
        self._name = name
        api.validate_sequence(strings, str, "strings")
        if len(strings) != self._number_of_values():
            raise ValueError("incorrect number of strings")
        for s in strings:
            api.validate_ascii_string(s, "label")
        self._strings = ImmutableTuple(strings)

    def __repr__(self):
        return "%s()" % (self.__class__.__name__)

    def __eq__(self, other):
        return isinstance(other, Label) \
               and self.type == other.type \
               and self.name == other.name \
               and self.universe == other.universe \
               and self.strings == other.strings

    def __hash__(self):
        return hash((self.type, self.name, self.universe, self.strings))

    # Data Model API implementation
    @property
    def universe(self):
        return self._universe

    @property
    def name(self):
        return self._name

    @property
    def strings(self):
        return self._strings


class AtomLabel(Label):

    def _number_of_values(self):
        return self._universe.number_of_atoms

    # Data Model API implementation
    @property
    def type(self):
        return 'atom'


class SiteLabel(Label):

    def _number_of_values(self):
        return self._universe.number_of_sites

    # Data Model API implementation
    @property
    def type(self):
        return 'site'


class TemplateAtomLabel(Label):

    def _number_of_values(self):
        return self._universe.number_of_template_atoms

    # Data Model API implementation
    @property
    def type(self):
        return 'template_atom'


class TemplateSiteLabel(Label):

    def _number_of_values(self):
        return self._universe.number_of_template_sites

    # Data Model API implementation
    @property
    def type(self):
        return 'template_site'


# Configuration

class Configuration(Immutable, api.MosaicConfiguration):

    def __init__(self, universe, positions, cell_parameters):
        api.validate_type(universe, Universe, "universe")
        cell_param_shape = universe.cell_parameter_array_shape
        nsites = universe.number_of_sites
        # Allow cell_parameters=None for universes that don't
        # need cell parameters.
        if cell_parameters is None and cell_param_shape == (0,):
            cell_parameters = IN.zeros(cell_param_shape,
                                       positions.dtype)
        # At this point, positions and cell_parameters must be arrays
        # of the required shapes.
        api.validate_array(positions, (nsites, 3),
                           self._allowed_dtypes, "positions")
        api.validate_array(cell_parameters, cell_param_shape,
                           self._allowed_dtypes, "cell_parameters")
        if positions.dtype != cell_parameters.dtype:
            raise ValueError("positions and cell parameters must have"
                             " the same element type")

        self._universe = universe
        self._positions = positions
        self._cell_parameters = cell_parameters

    def __repr__(self):
        return "Configuration()"

    def __eq__(self, other):
        return isinstance(other, Configuration) \
               and self.universe == other.universe \
               and (self.cell_parameters == other.cell_parameters).all() \
               and (self.positions == other.positions).all()

    def __hash__(self):
        # cheap hash function that doesn't look at the positions
        return hash((self.universe, self.positions.shape,
                     tuple(self.cell_parameters.flat)))

    # Reimplementation that returns immutable data
    # Inactive for now because this leads to strange bugs with various
    # versions of NumPy.
    if False:
        def lattice_vectors(self):
            return immutable(api.MosaicConfiguration.lattice_vectors(self))

    # Data Model API implementation
    @property
    def universe(self):
        return self._universe

    @property
    def cell_parameters(self):
        return self._cell_parameters

    @property
    def positions(self):
        return self._positions


class Selection(Immutable, api.MosaicSelection):

    def __init__(self, universe, indices):
        api.validate_type(universe, Universe, "universe")
        self._universe = universe
        indices = IN.array(indices, uint_for_max_value(max(indices)))
        api.validate_indices(indices, self._max_index(), "indices")
        self._indices = indices

    def __str__(self):
        return "%s()" % (self.__class__.__name__)

    def __eq__(self, other):
        return isinstance(other, Selection) \
               and self.type == other.type \
               and self.universe == other.universe \
               and (self.indices == other.indices).all()

    def __hash__(self):
        return hash(self.type) + hash(self.universe) + hash(len(self.indices))

    # Data Model API implementation
    @property
    def universe(self):
        return self._universe

    @property
    def indices(self):
        return self._indices


class AtomSelection(Selection):

    def _max_index(self):
        return self._universe.number_of_atoms

    # Data Model API implementation
    @property
    def type(self):
        return 'atom'


class SiteSelection(Selection):

    def _max_index(self):
        return self._universe.number_of_sites

    # Data Model API implementation
    @property
    def type(self):
        return 'site'


class TemplateAtomSelection(Selection):

    def _max_index(self):
        return self._universe.number_of_template_atoms

    # Data Model API implementation
    @property
    def type(self):
        return 'template_atom'


class TemplateSiteSelection(Selection):

    def _max_index(self):
        return self._universe.number_of_template_sites

    # Data Model API implementation
    @property
    def type(self):
        return 'template_site'


class Factory(AbstractFactory):

    handler = MethodRegister()

    @handler(api.MosaicUniverse)
    def _make_universe(self, mosaic_universe):

        def make_atom(a):
            cons = {"dummy": dummy,
                    "unknown": unknown,
                    "element": element}[a.type]
            return atom(cons(a.name), a.number_of_sites)

        def make_fragment(f):
            fragments = [(sub.label, make_fragment(sub))
                         for sub in f.fragments]
            atoms = [(a.label, make_atom(a)) for a in f.atoms]
            if f.is_polymer:
                return polymer(f.species, fragments, f.bonds, f.polymer_type)
            else:
                return fragment(f.species, fragments, atoms, f.bonds)

        return universe(mosaic_universe.cell_shape,
                        tuple((make_fragment(f), f.label, n)
                              for f, n in mosaic_universe.molecules),
                        mosaic_universe.symmetry_transformations,
                        mosaic_universe.convention)

    @handler(api.MosaicProperty)
    def _make_property(self, mosaic_property):
        universe = self(mosaic_property.universe)
        klass = {"atom": AtomProperty,
                 "site": SiteProperty,
                 "template_atom": TemplateAtomProperty,
                 "template_site": TemplateSiteProperty}[mosaic_property.type]
        return klass(universe, mosaic_property.name, mosaic_property.units,
                     mosaic_property.data)

    @handler(api.MosaicLabel)
    def _make_label(self, mosaic_label):
        universe = self(mosaic_label.universe)
        klass = {"atom": AtomLabel,
                 "site": SiteLabel,
                 "template_atom": TemplateAtomLabel,
                 "template_site": TemplateSiteLabel}[mosaic_label.type]
        return klass(universe, mosaic_label.name, mosaic_label.strings)

    @handler(api.MosaicConfiguration)
    def _make_configuration(self, mosaic_configuration):
        universe = self(mosaic_configuration.universe)
        return Configuration(universe,
                             mosaic_configuration.positions,
                             mosaic_configuration.cell_parameters)

    @handler(api.MosaicSelection)
    def _make_selection(self, mosaic_selection):
        universe = self(mosaic_selection.universe)
        klass = {"atom": AtomSelection,
                 "site": SiteSelection,
                 "template_atom": TemplateAtomSelection,
                 "template_site": TemplateSiteSelection}[mosaic_selection.type]
        return klass(universe, mosaic_selection.indices)
