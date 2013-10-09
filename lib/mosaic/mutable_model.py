# -*- coding: utf-8 -*-
"""Mutable model

.. moduleauthor:: Konrad Hinsen

This is the main in-memory representation of Mosaic data.  Like most
Python code, it uses mutable data structures, which are efficient but
can become inconsistent if not used carefully.

"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import copy
import itertools as IT
import types

import numpy as N

import mosaic.api as api
from mosaic.utility import AbstractFactory
from mosaic.utility import MethodRegister
from mosaic.utility import uint_for_max_value


class MosaicObject(object):

    def __repr__(self):
        return str(self)

    def __ne__(self, other):
        return not self.__eq__(other)


class AtomDescriptor(MosaicObject):

    def __new__(cls, name=""):
        try:
            return cls._instances[name]
        except KeyError:
            ob = super(AtomDescriptor, cls).__new__(cls)
            api.validate_type(name, str, "name")
            cls._validate_name(name)
            ob._name = name
            cls._instances[name] = ob
            return ob

    @classmethod
    def _validate_name(cls, name):
        pass

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __eq__(self, other):
        return self is other

    def __hash__(self):
        return id(self)

    def __str__(self):
        return "%s('%s')" % (self.__class__.__name__, self._name)

    @property
    def name(self):
        return self._name


class Element(AtomDescriptor):

    _instances = {}

    @classmethod
    def _validate_name(cls, name):
        api.MosaicAtom.validate_element_name(name)

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __str__(self):
        return "Element(%s)" % self._name

    @property
    def type(self):
        return "element"


class CGParticle(AtomDescriptor):

    _instances = {}

    def __str__(self):
        if self._name:
            return "CGParticle('%s')" % self._name
        else:
            return "CGParticle()"

    @property
    def type(self):
        return "cgparticle"


class Dummy(AtomDescriptor):

    _instances = {}

    def __str__(self):
        if self._name:
            return "Dummy('%s')" % self._name
        else:
            return "Dummy()"

    @property
    def type(self):
        return "dummy"


class Unknown(AtomDescriptor):

    _instances = {}

    def __str__(self):
        if self._name:
            return "Unknown('%s')" % self._name
        else:
            return "Unknown()"

    @property
    def type(self):
        return ""


class Path(tuple):

    def __str__(self):
        return '.'.join(self)


class MolecularStructureObject(MosaicObject):

    def __copy__(self):
        return copy.deepcopy(self)

    def has_ancestor(self, obj):
        return self is obj or \
               (self._parent is not None and self._parent.has_ancestor(obj))

    def toplevel_ancestor(self):
        if self._parent is None:
            return self
        else:
            return self._parent.toplevel_ancestor()

    def ancestors(self, target=None):
        if self._parent is target:
            return (self,)
        else:
            if self._parent is None:
                raise ValueError("%s is not an ancestor of %s"
                                 % (str(target), str(self)))
            return self._parent.ancestors(target) + (self,)

    def common_ancestor(self, objs):
        ancestors = [obj.ancestors() for obj in objs]
        common = None
        while all(ancestors):
            top = set(a[0] for a in ancestors)
            if len(top) == 1:
                common = top.pop()
                ancestors = [a[1:] for a in ancestors]
            else:
                return common

    def path_to(self, obj):
        return Path(a._label for a in obj.ancestors(self))

    def full_path(self):
        return Path(obj._label for obj in self.ancestors())


class Atom(MolecularStructureObject, api.MosaicAtom):

    def __init__(self, label, descriptor, nsites=1):
        self.label = label

        api.validate_type(descriptor, AtomDescriptor, "atom descriptor")
        api.validate_type(nsites, int, "number of sites")
        if nsites < 1:
            raise ValueError("number of sites must be >= 1")

        self._descriptor = descriptor
        self.nsites = nsites
        self._parent = None

    def __deepcopy__(self, memo):
        return Atom(self._label, self._descriptor, self.nsites)

    def __getitem__(self, item):
        assert isinstance(item, str)
        raise KeyError(item)

    def __eq__(self, other):
        return isinstance(other, Atom) \
               and self.type == other.type \
               and self.name == other.name \
               and self.label == other.label \
               and self.nsites == other.nsites

    def __hash__(self):
        return hash(self.name + self.label + str(self.nsites))

    def __str__(self):
        return "%s('%s')" % (self.__class__.__name__, self._label)

    @property
    def parent(self):
        return self._parent

    # Data Model API implementation
    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, new_label):
        api.validate_label(new_label, "Atom.label")
        self._label = new_label

    @property
    def name(self):
        return self._descriptor.name

    @property
    def type(self):
        return self._descriptor.type

    @property
    def number_of_sites(self):
        return self.nsites


class Fragment(MolecularStructureObject, api.MosaicFragment):

    def __init__(self, label, species, fragments=[], atoms=[], bonds=[]):
        self.label = label
        self.species = species

        self._parent = None
        self._fragments = []
        self._atoms = []
        self._bonds = []
        self._attrs = {}

        self.add_fragments(fragments)
        self.add_atoms(atoms)
        self.add_bonds(bonds)

    def __str__(self):
        return "%s('%s', '%s')" % (self.__class__.__name__,
                                   self._label, self._species)

    def add_fragments(self, fragments):
        api.validate_sequence(fragments, Fragment, "fragments")
        new_attrs = set()
        for f in fragments:
            self._check_child_obj(f, new_attrs)
        for f in fragments:
            self._add_child_obj(f)
        self._fragments.extend(fragments)

    def add_atoms(self, atoms):
        api.validate_sequence(atoms, Atom, "atoms")
        new_attrs = set()
        for a in atoms:
            self._check_child_obj(a, new_attrs)
        for a in atoms:
            self._add_child_obj(a)
        self._atoms.extend(atoms)

    def add_bonds(self, bonds):
        api.validate_sequence(bonds, tuple, "bonds",
                              ((lambda p: len(p) == 3,
                                "must have length  3"),
                               (lambda p: (isinstance(p[0], str)
                                           or isinstance(p[0], Atom))
                                          and (isinstance(p[1], str)
                                               or isinstance(p[1], Atom)),
                                "elements must be strings or atoms"),
                               (lambda p: p[2] in self._bond_orders,
                                "bond order must be one of "
                                + str(self._bond_orders))))
        for a1, a2, order in bonds:
            atom1 = self._get_atom(a1)
            atom2 = self._get_atom(a2)
            if atom1 is atom2:
                raise ValueError("bond between %s and itself"
                                 % str(atom1.full_path()))
            p = self.common_ancestor((atom1, atom2))
            if p is self:
                self._bonds.append((atom1, atom2, order))
            else:
                raise ValueError("bond %s-%s must be defined in fragment %s"
                                 % (str(atom1.full_path()),
                                    str(atom2.full_path()),
                                    p._label))

    def _check_child_obj(self, obj, new_attrs):
        if obj._parent is not None:
            raise ValueError("%s is already part of another fragment (%s)"
                             % (str(obj), obj._parent._label))
        l = obj._label
        if l in self._attrs or l in new_attrs:
            raise ValueError("Label %s occurs more than once" % l)
        new_attrs.add(l)

    def _add_child_obj(self, obj):
        obj._parent = self
        self._attrs[obj._label] = obj

    def __deepcopy__(self, memo):
        return Fragment(self._label, self._species,
                        copy.deepcopy(self._fragments, memo),
                        copy.deepcopy(self._atoms, memo),
                        tuple((str(self.path_to(a1)),
                               str(self.path_to(a2)),
                               order)
                              for a1, a2, order in self._bonds))

    def eq_structure(self, other):
        def bond_path_set(frag):
            set(tuple(sorted((a1, a2)) + [order])
                for a1, a2, order in frag.bonds)
        return self.species == other.species \
               and len(self.fragments) == len(other.fragments) \
               and len(self.atoms) == len(other.atoms) \
               and len(self.bonds) == len(other.bonds) \
               and all(f1 == f2
                       for f1, f2 in zip(self.fragments, other.fragments)) \
               and all(a1 == a2 for a1, a2 in zip(self.atoms, other.atoms)) \
               and bond_path_set(self) == bond_path_set(other)

    def __eq__(self, other):
        return isinstance(other, Fragment) \
               and self.label == other.label \
               and self.eq_structure(other)

    def __hash__(self):
        return hash(self.label + self.species)

    def _get_atom(self, atom_spec):
        if isinstance(atom_spec, str):
            try:
                return self[atom_spec]
            except KeyError:
                raise ValueError("fragment %s does not contain atom %s"
                                 % (self._label, atom_spec))
        else:
            if not atom_spec.has_ancestor(self):
                raise ValueError("fragment %s does not contain atom %s"
                                 % (self._label, str(atom_spec)))
            return atom_spec

    def __getitem__(self, path):
        assert isinstance(path, str) or isinstance(path, Path)
        obj = self._attrs
        if isinstance(path, str):
            path = path.split('.')
        for item in path:
            obj = obj[item]
        return obj

    @property
    def parent(self):
        return self._parent

    # Data Model API implementation
    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, new_label):
        api.validate_label(new_label, "Fragment.label")
        self._label = new_label

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, new_species):
        api.validate_label(new_species, "Fragment.species")
        self._species = new_species

    @property
    def is_polymer(self):
        return False

    @property
    def polymer_type(self):
        raise ValueError("Non-polymer fragment does not have a polymer type")

    @property
    def fragments(self):
        return self._fragments

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return tuple((str(self.path_to(a1)),
                      str(self.path_to(a2)),
                      order)
                     for a1, a2, order in self._bonds)


class Polymer(Fragment):

    def __init__(self, label, species, fragments, bonds, polymer_type=''):
        Fragment.__init__(self, label, species, fragments, (), bonds)
        self.polymer_type = polymer_type

    def __deepcopy__(self, memo):
        return Polymer(self._label, self._species,
                       copy.deepcopy(self._fragments, memo),
                       tuple((str(self.path_to(a1)),
                              str(self.path_to(a2)),
                              order)
                             for a1, a2, order in self._bonds),
                       self._polymer_type)

    def __getitem__(self, item):
        if isinstance(item, str):
            return Fragment.__getitem__(self, item)
        assert isinstance(item, int)
        return self._fragments[item]

    def __len__(self):
        return len(self._fragments)

    def __eq__(self, other):
        return isinstance(other, Polymer) \
               and Fragment.__eq__(self, other) \
               and self.polymer_type == other.polymer_type

    __hash__ = Fragment.__hash__

    # Data Model API implementation
    @property
    def is_polymer(self):
        return True

    @property
    def polymer_type(self):
        return self._polymer_type

    @polymer_type.setter
    def polymer_type(self, new_polymer_type):
        api.validate_value(new_polymer_type, self._polymer_types,
                           'polymer_type')
        self._polymer_type = new_polymer_type


class Universe(MosaicObject, api.MosaicUniverse):

    def __init__(self, cell_shape, molecules,
                 symmetry_transformations=(),
                 convention=''):
        self.cell_shape = cell_shape
        self.symmetry_transformations = symmetry_transformations
        self.convention = convention
        self._molecules = []
        self._templates = []
        self._molecule_counts = []
        self._atom_counts = []
        self._site_counts = []
        self._bond_counts = []
        self.add_molecules(molecules)

    def add_molecules(self, molecules):
        api.validate_sequence(molecules, tuple, "molecules",
                              ((lambda p: len(p) == 2, "must have length 2"),
                               (lambda p: isinstance(p[0], Fragment)
                                          and isinstance(p[1], int),
                                "elements must be (fragment, count) pairs")))
        self._molecules.extend(molecules)
        self._templates.extend([f for f, c in molecules])
        self._molecule_counts.extend([c for f, c in molecules])
        self._atom_counts.extend([f.number_of_atoms for f, c in molecules])
        self._site_counts.extend([f.number_of_sites for f, c in molecules])
        self._bond_counts.extend([f.number_of_bonds for f, c in molecules])

    def __eq__(self, other):
        return isinstance(other, Universe) \
               and self.cell_shape == other.cell_shape \
               and self.symmetry_transformations == \
                              other.symmetry_transformations \
               and self.convention == other.convention \
               and all(f1 == f2 and c1 == c2
                       for (f1, c1), (f2, c2)
                            in zip(self.molecules, other.molecules))

    def __hash__(self):
        return hash(self.cell_shape) \
               + len(self.molecules)

    def __str__(self):
        mol_str = ', '.join("(%s('%s'), %d)"
                            % (f.__class__.__name__, f._label, c)
                            for f, c in self._molecules)
        return "Universe('%s', [%s])" % (self.cell_shape, mol_str)

    # Data Model API implementation

    @property
    def cell_shape(self):
        return self._cell_shape

    @cell_shape.setter
    def cell_shape(self, new_cell_shape):
        api.validate_type(new_cell_shape, str, "cell_shape")
        api.validate_value(new_cell_shape,
                           self._cell_parameter_array_shapes.keys(),
                           "cell_shape")
        self._cell_shape = new_cell_shape

    @property
    def symmetry_transformations(self):
        return self._symmetry_transformations

    @symmetry_transformations.setter
    def symmetry_transformations(self, symmetry_transformations):
        if self.cell_shape == "infinite" \
           and len(symmetry_transformations) > 0:
            raise ValueError("Symmetry transformations are allowed "
                             "only in periodic universes")
        api.validate_sequence(symmetry_transformations, tuple,
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
        self._symmetry_transformations = symmetry_transformations

    @property
    def convention(self):
        return self._convention

    @convention.setter
    def convention(self, new_convention):
        api.validate_label(new_convention, "convention")
        self._convention = new_convention

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


class Property(MosaicObject, api.MosaicProperty):

    def __init__(self, universe, name, units,
                 data=None, dtype=None, element_shape=None):
        api.validate_type(universe, Universe, "universe")
        self._universe = universe
        api.validate_label(name, "name")
        self._name = name
        api.validate_units(units, "units")
        self._units = units
        if data is None:
            if dtype is None:
                dtype = N.float64
            if dtype not in self._allowed_dtypes:
                raise ValueError("dtype must be " +
                                 " or ".join(str(t)
                                             for t in self._allowed_dtypes))
            if element_shape is None:
                element_shape = ()
            array_shape = (self._number_of_values(),) + element_shape

            self._data = N.empty(array_shape, dtype)
        else:
            api.validate_array(data, None, self._allowed_dtypes, "data")
            if dtype is not None and data.dtype != dtype:
                raise ValueError("elements of data array do not have "
                                 "the requested type")
            if data.shape[0] != self._number_of_values():
                raise ValueError("data array has incorrect shape")
            if element_shape is not None and data.shape[1:] != element_shape:
                raise ValueError("elements of data array do not have "
                                 "the requested shape")
            self._data = data

    def __str__(self):
        return "%s()" % (self.__class__.__name__)

    def __eq__(self, other):
        return isinstance(other, Property) \
               and self.type == other.type \
               and self.universe == other.universe \
               and self.name == other.name \
               and self.units == other.units \
               and (self.data == other.data).all()

    def __hash__(self):
        return hash(self.type) + hash(self.universe) + hash(self.data.shape)

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


class Label(MosaicObject, api.MosaicLabel):

    def __init__(self, universe, name, strings=None):
        api.validate_type(universe, Universe, "universe")
        self._universe = universe
        api.validate_label(name, "name")
        self._name = name
        if strings is None:
            self._strings = self._number_of_values() * ['']
        else:
            api.validate_sequence(strings, str, "strings")
            if len(strings) != self._number_of_values():
                raise ValueError("incorrect number of strings")
            for s in strings:
                api.validate_ascii_string(s, "label")
            self._strings = strings

    def __str__(self):
        return "%s()" % (self.__class__.__name__)

    def __eq__(self, other):
        return isinstance(other, Label) \
               and self.type == other.type \
               and self.name == other.name \
               and self.universe == other.universe \
               and self.strings == other.strings

    def __hash__(self):
        return hash(self.type) + hash(self.name) + \
               hash(self.universe) + hash(self.strings)

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


class Configuration(MosaicObject, api.MosaicConfiguration):

    def __init__(self, universe, positions=None, cell_parameters=None,
                 dtype=None):
        api.validate_type(universe, Universe, "universe")
        cell_param_shape = universe.cell_parameter_array_shape
        nsites = universe.number_of_sites
        if positions is None and cell_parameters is None:
            if dtype is None:
                dtype = N.float64
            if dtype not in self._allowed_dtypes:
                raise ValueError("dtype must be " +
                                 " or ".join(str(t)
                                             for t in self._allowed_dtypes))
            positions = N.empty((nsites, 3), dtype)
            cell_parameters = N.empty(cell_param_shape, dtype)
        else:
            # Allow cell_parameters=None for universes that don't
            # need cell parameters.
            if positions is not None and cell_parameters is None \
               and cell_param_shape == (0,):
                cell_parameters = N.empty(cell_param_shape,
                                          positions.dtype)
            # Require both positions and cell parameters, or neither
            if positions is None or cell_parameters is None:
                raise ValueError("configuration requires both "
                                 "positions and cell parameters")
            # At this point, positions and cell_parameters must be arrays
            # of the required shapes.
            api.validate_array(positions, (nsites, 3),
                               self._allowed_dtypes, "positions")
            api.validate_array(cell_parameters, cell_param_shape,
                               self._allowed_dtypes, "cell_parameters")
            if positions.dtype != cell_parameters.dtype:
                raise ValueError("positions and cell parameters must have"
                                 " the same element type")
            # If an explicit dtype is given, check arrays for conformance
            if dtype is not None:
                if positions.dtype != dtype:
                    raise ValueError("arrays don't have the requested"
                                     " element type " + str(dtype))

        self._universe = universe
        self._positions = positions
        self._cell_parameters = cell_parameters

    def __str__(self):
        return "Configuration()"

    def __eq__(self, other):
        return isinstance(other, Configuration) \
               and self.universe == other.universe \
               and (self.cell_parameters == other.cell_parameters).all() \
               and (self.positions == other.positions).all()

    def __hash__(self):
        return hash(self.universe) \
               + hash(self.positions.shape) + hash(self.cell_parameters.shape)

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


class Selection(MosaicObject, api.MosaicSelection):

    def __init__(self, universe, indices):
        api.validate_type(universe, Universe, "universe")
        self._universe = universe
        indices = N.array(indices, uint_for_max_value(max(indices)))
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
    def _make_universe(self, universe):

        def atom(a):
            klass = {"dummy": Dummy,
                     "unknown": Unknown,
                     "element": Element}[a.type]
            return Atom(a.label, klass(a.name), a.number_of_sites)

        def fragment(f):
            fragments = [fragment(sub) for sub in f.fragments]
            atoms = [atom(a) for a in f.atoms]
            if f.is_polymer:
                return Polymer(f.label, f.species,
                               fragments, f.bonds, f.polymer_type)
            else:
                return Fragment(f.label, f.species,
                                fragments, atoms, f.bonds)

        return Universe(universe.cell_shape,
                        [(fragment(f), n) for f, n in universe.molecules],
                        universe.symmetry_transformations,
                        universe.convention)

    @handler(api.MosaicProperty)
    def _make_property(self, property):
        universe = self(property.universe)
        klass = {"atom": AtomProperty,
                 "site": SiteProperty,
                 "template_atom": TemplateAtomProperty,
                 "template_site": TemplateSiteProperty}[property.type]
        return klass(universe, property.name, property.units,
                     property.data)

    @handler(api.MosaicLabel)
    def _make_label(self, label):
        universe = self(label.universe)
        klass = {"atom": AtomLabel,
                 "site": SiteLabel,
                 "template_atom": TemplateAtomLabel,
                 "template_site": TemplateSiteLabel}[label.type]
        return klass(universe, label.name, label.strings)

    @handler(api.MosaicConfiguration)
    def _make_configuration(self, configuration):
        universe = self(configuration.universe)
        return Configuration(universe,
                             configuration.positions,
                             configuration.cell_parameters)

    @handler(api.MosaicSelection)
    def _make_selection(self, selection):
        universe = self(selection.universe)
        klass = {"atom": AtomSelection,
                 "site": SiteSelection,
                 "template_atom": TemplateAtomSelection,
                 "template_site": TemplateSiteSelection}[selection.type]
        return klass(universe, selection.indices)

