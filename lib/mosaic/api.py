# -*- coding: utf-8 -*-
"""The Mosaic Python API

.. moduleauthor:: Konrad Hinsen

This module provides abstract base classes that define the Mosaic
Python API and implement validation code. They also provide a few
convenience functions implemented in terms of the raw API.

Concrete implementations subclass the abstract base classes and
implement all the abstract properties.

"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

from abc import ABCMeta, abstractproperty
import collections
import itertools as IT
import re

import numpy as N
import numpy.linalg as LA

import mosaic.utility


# Mosaic version number (major, minor)
# An increase in the minor number indicates a superset of the preceding
# versions. An increase in the major number indicated an incompatibility
# with preceding versions.

MOSAIC_VERSION = (1, 0)


# Base class for all classes that represent top-level data items,
# i.e. items that can be stored in files, retrieved, etc.

class MosaicDataItem(object):

    """Base class for top-level data items

    Instances of subclasses of MosaicDataItem can be stored in files.
    """

    __metaclass__ = ABCMeta


# An atom is defined by the following characteristics:
#
#  - a type, which is
#    - 'element' for a standard atom that has a chemical element.
#      The element symbol is the value of the name attribute.
#    - 'cgparticle' , for particles in coarse-grained models
#      representing multiple atoms.
#    - 'dummy' for pseudo-atoms that don't physically exist,
#      such as virtual interaction sites.
#    - '' for anything else
#
# - a name, which is the chemical element symbol for atoms of type
#   'element', and any suitable identifier for the other types
#
# - a label, which identifies the atom uniquely inside its fragment
#
# - the number of sites, equal to the number of Cartesian coordinate
#   sets required by the atom. It is > 1 for path integral, wave
#   functions, atoms with alternate positions in crystal structures, etc.

class MosaicAtom(object):

    """Atom inside a :py:class:`MosaicUniverse`

    See the :ref:`data model documentation<mosaic-atom>` for atoms.

    """

    __metaclass__ = ABCMeta

    # API properties

    @abstractproperty
    def label(self):
        """An ASCII string not containing dots and identifying
        the atom uniquely inside its parent fragment.
        See the :ref:`data model documentation<mosaic-atom-label>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def name(self):
        """An ASCII string describing the type of the atom. For 'real'
        atoms in the chemical sense, this must be the chemical element
        symbol.
        See the :ref:`data model documentation<mosaic-atom-name>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def type(self):
        """A string identifying the type of the atom.
        See the :ref:`data model documentation<mosaic-atom-type>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def number_of_sites(self):
        """The number of sites associated with the atom.
        See the :ref:`data model documentation<mosaic-atom-nsites>`.
        """
        raise NotImplementedError()

    # Property shared by all atoms

    @property
    def number_of_atoms(self):
        """The number of atoms associated with the atom (always 1)
        """
        return 1

    # Equivalence test

    def validate_equivalence(self, other):
        """Verify the equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :raises ValueError: if other is not equivalent to self
        """
        if not isinstance(other, MosaicAtom):
            raise TypeError("%s is not an atom" % str(type(other)))
        if self.label != other.label:
            raise ValueError("labels differ: %s != %s"
                             % (repr(self.label), repr(other.label)))
        if self.name != other.name:
            raise ValueError("names differ: %s != %s"
                             % (repr(self.name), repr(other.name)))
        if self.type != other.type:
            raise ValueError("types differ: %s != %s"
                             % (repr(self.type), repr(other.type)))
        if self.number_of_sites != other.number_of_sites:
            raise ValueError("site numbers differ: %d != %d"
                             % (self.number_of_sites, other.number_of_sites))

    def is_equivalent(self, other):
        """Check for equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :returns: True if other is equivalent to self
        :rtype: bool
        """
        try:
            self.validate_equivalence(other)
            return True
        except:
            return False

    # Validation

    _allowed_types = ('element', 'cgparticle', 'dummy', '')

    _elements = ('Ac', 'Ag', 'Al', 'Am', 'Ar', 'As', 'At', 'Au',
                 'B', 'Ba', 'Be', 'Bh', 'Bi', 'Bk', 'Br',
                 'C', 'Ca', 'Cd', 'Ce', 'Cf', 'Cl', 'Cm',
                 'Co', 'Cn', 'Cr', 'Cs', 'Cu',
                 'D', 'Db', 'Ds', 'Dy',
                 'Er', 'Es', 'Eu',
                 'F', 'Fe', 'Fm', 'Fr',
                 'Ga', 'Gd', 'Ge',
                 'H', 'He', 'Hf', 'Hg', 'Ho', 'Hs',
                 'I', 'In', 'Ir',
                 'K', 'Kr',
                 'La', 'Li', 'Lr', 'Lu',
                 'Md', 'Mg', 'Mn', 'Mo', 'Mt',
                 'N', 'Na', 'Nb', 'Nd', 'Ne', 'Ni', 'No', 'Np',
                 'O', 'Os',
                 'P', 'Pa', 'Pb', 'Pd', 'Pm', 'Po', 'Pr', 'Pt', 'Pu',
                 'Ra', 'Rb', 'Re', 'Rf', 'Rg', 'Rh', 'Rn', 'Ru',
                 'S', 'Sb', 'Sc', 'Se', 'Sg', 'Si', 'Sm', 'Sn', 'Sr',
                 'Ta', 'Tb', 'Tc', 'Te', 'Th', 'Ti', 'Tl', 'Tm',
                 'U', 'V', 'W', 'Xe', 'Y', 'Yb', 'Zn', 'Zr')

    @classmethod
    def validate_element_name(self, name):
        validate_value(name, self._elements, "name")

    def validate(self):
        """Verify that the object satisfies the constraints of the
        Mosaic data model.

        :raises ValueError: if the object is not valid Mosaic data
        """
        validate_label(self.label, "Atom.label")
        validate_value(self.type, self._allowed_types, "Atom type")
        validate_label(self.name, "Atom.name")
        if self.type == 'element':
            validate_value(self.name, self._elements, "Atom element")
        if not isinstance(self.number_of_sites, int):
            raise ValueError("Atom.number_of_sites must be an integer")
        if self.number_of_sites <= 0:
            raise ValueError("Atom.number_of_sites must be positive")
        if self.number_of_atoms != 1:
            raise ValueError("Atom.number_of_atoms must be 1")


# A fragment describes a node in the tree defining the chemical
# structure of the molecules. It can represent a molecule,
# a supermolecule, or part of a molecule. A fragment is defined
# by
#
#  - a species, which is a text string describing the chemical
#    nature of the fragment
#
#  - a label, which describes the role of the fragment inside its
#    parent structure, and must be unique within that parent structure
#
#  - a list of sub-fragments
#
#  - a list of atoms
#
#  - a list of bonds
#
#  - the boolean flag is_polymer. Polymer fragments have an emtpy
#    atom list and an additional attribute 'polymer_type' whose
#    allowed values are
#      - 'polypeptide', for a peptide chain
#      - 'polyribonucleotide', for an RNA chain
#      - 'polydeoxyribonucleotide', for a DNA chain
#      - 'polynucleotide', for a chain that can contain
#        nucleotides with either type of sugar
#      - the empty string, for any other polymer

class MosaicFragment(collections.Mapping):

    """Fragment inside a :py:class:`MosaicUniverse`

    See the :ref:`data model documentation<mosaic-fragment>` for fragments.

    Fragments implement the Mapping interface. Valid keys are strings
    identifying atoms or sub-fragments, using dots to separate subsequent
    labels. For polymer fragments, integer keys are valid as well to refer
    to a specific sub-fragment in the chain.
    """

    # API properties

    @abstractproperty
    def label(self):
        """An ASCII string not containing dots and identifying
        the fragment uniquely inside its parent fragment.
        See the :ref:`data model documentation<mosaic-fragment-label>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def species(self):
        """An ASCII string describing the species of the fragment
        (e.g. what kind of molecule or moiety it represents).
        See the :ref:`data model documentation<mosaic-fragment-species>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def fragments(self):
        """Sequence of sub-fragments, may be empty.
        See the :ref:`data model documentation<mosaic-fragment-fragments>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def atoms(self):
        """Sequence of atoms, may be empty.
        See the :ref:`data model documentation<mosaic-fragment-atoms>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def bonds(self):
        """Sequence of bonds, may be empty. Each bond is
        repreented by a tuple (atom_ref_1, atom_ref_2, bond_order).
        See the :ref:`data model documentation<mosaic-fragment-bonds>`
        and the :ref:`bond reference documentation<mosaic-bonds>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def is_polymer(self):
        """True if the fragment is a polymer.
        See the :ref:`data model documentation<mosaic-fragment-polymer>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def polymer_type(self):
        """String identifying the polymer type if is_polymer is true.
        See the :ref:`data model documentation<mosaic-fragment-polymer>`.
        """
        raise NotImplementedError()

    # Properties that can be computed in terms of the API properties

    @property
    def number_of_atoms(self):
        """The number of atoms in the fragment (including the atoms
        in sub-fragments).
        """
        return sum(f.number_of_atoms for f in self.fragments) \
               + len(self.atoms)

    @property
    def number_of_sites(self):
        """The number of sites associated with the fragment, i.e.
        the sum of the numbers of sites of all the fragment's atoms,
        including those of sub-fragments.
        """
        return sum(f.number_of_sites for f in self.fragments) \
               + sum(a.number_of_sites for a in self.atoms)

    @property
    def number_of_bonds(self):
        """The number of bonds associated with the fragment,
        including the bonds inside sub-fragments.
        """
        return sum(ff.number_of_bonds for ff in self.fragments) \
               + len(self.bonds)

    # Equivalence test

    def validate_equivalence(self, other):
        """Verify the equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :raises ValueError: if other is not equivalent to self
        """
        if not isinstance(other, MosaicFragment):
            raise TypeError("%s is not a fragment" % str(type(other)))
        if self.label != other.label:
            raise ValueError("labels differ: %s != %s"
                             % (repr(self.label), repr(other.label)))
        if self.species != other.species:
            raise ValueError("species differ: %s != %s"
                             % (repr(self.species), repr(other.species)))
        if self.is_polymer:
            if not other.is_polymer:
                raise ValueError("%s is not a polymer" % str(other))
            if self.polymer_type != other.polymer_type:
                raise ValueError("polymer types differ: %s != %s"
                                 % (repr(self.polymer_type),
                                    repr(other.polymer_type)))
        for s, o in zip(self.fragments, other.fragments):
            s.validate_equivalence(o)
        for s, o in zip(self.atoms, other.atoms):
            s.validate_equivalence(o)
        if self._bond_set() != other._bond_set():
            raise ValueError("bonds differ: %s != %s"
                             % (repr(self._bond_set), repr(other._bond_set)))

    def _bond_set(self):
        return frozenset((frozenset((a1, a2)), order)
                         for a1, a2, order in self.bonds)

    def is_equivalent(self, other):
        """Check for equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :returns: True if other is equivalent to self
        :rtype: bool
        """
        try:
            self.validate_equivalence(other)
            return True
        except:
            return False

    # Validation

    _polymer_types = ['',
                      'polypeptide',
                      'polyribonucleotide',
                      'polydeoxyribonucleotide',
                      'polynucleotide']

    _bond_orders = ['', 'single', 'double', 'triple', 'quadruple', 'aromatic']

    def validate(self):
        """Verify that the object satisfies the constraints of the
        Mosaic data model.

        :raises ValueError: if the object is not valid Mosaic data
        """
        validate_label(self.label, "Fragment.label")
        validate_label(self.species, "Fragment.species")
        validate_value(self.is_polymer, [True, False], "Fragment.is_polymer")
        if self.is_polymer:
            validate_value(self.polymer_type, self._polymer_types,
                           "Fragment.polymer_type")
        validate_sequence(self.fragments, MosaicFragment, "Fragment.fragments")
        labels = set()
        for f in self.fragments:
            f.validate()
            if f.label in labels:
                raise ValueError("Label %s occurs more than once" % f.label)
            labels.add(f.label)
        validate_sequence(self.atoms, MosaicAtom, "Fragment.atoms")
        for a in self.atoms:
            a.validate()
            if a.label in labels:
                raise ValueError("Label %s occurs more than once" % a.label)
            labels.add(a.label)
        for property in ['number_of_atoms',
                         'number_of_sites',
                         'number_of_bonds']:
            value = getattr(self, property)
            reference = getattr(MosaicFragment, property).fget(self)
            if value != reference:
                raise ValueError("Fragment.%s is %s, should be %s"
                                 % (property, str(value), str(reference)))

    # Methods based on API properties

    def recursive_atom_iterator(self):
        """An iterator over the atoms in the fragment, including the
        atoms in sub-fragments.
        """
        return IT.chain(IT.chain.from_iterable(f.recursive_atom_iterator()
                                               for f in self.fragments),
                        iter(self.atoms))

    def recursive_atom_path_iterator(self):
        """An iterator over the atom paths in the fragment, including the
        atoms in sub-fragments.
        """
        for f in self.fragments:
            l = f.label
            for ap in f.recursive_atom_path_iterator():
                yield l + '.' + ap
        for a in self.atoms:
            yield a.label

    def recursive_bond_iterator(self):
        """An iterator over the bonds in the fragment, including the
        bonds in sub-fragments.
        """
        for f in self.fragments:
            l = f.label
            for a1, a2, order in f.recursive_bond_iterator():
                yield (l + '.' + a1, l + '.' + a2, order)
        for b in self.bonds:
            yield b

    def site_to_atom_index_mapping(self):
        """
        :returns: an array whose element [s] is the atom index
                  corresponding to site index s.
        :rtype: numpy.ndarray
        """
        ns = [a.number_of_sites for a in self.recursive_atom_iterator()]
        return N.repeat(N.arange(len(ns)), ns)

    def _atom_to_site_index_mapping(self):
        ns = [a.number_of_sites for a in self.recursive_atom_iterator()]
        return N.add.accumulate(ns)

    def atom_to_site_index_mapping(self):
        """
        :returns: an array whose element [s] is the site index
                  corresponding to the first site of the atom with
                  atom index s.
        :rtype: numpy.ndarray
        """
        m = self._atom_to_site_index_mapping()
        m[1:] = m[:-1]
        m[0] = 0
        return m

    # Mapping interface

    def __getitem__(self, item):
        if isinstance(item, int) and self.is_polymer:
            return self.fragments[item]
        # A rather inefficient implementation of substructure selection
        # using only the public API elements. Concrete implementations
        # can do better.
        assert isinstance(item, str)
        path = item.split('.')
        item = None
        for f in self.fragments:
            if f.label == path[0]:
                item = f
                break
        if item is None:
            for a in self.atoms:
                if a.label == path[0]:
                    item = a
                    break
        if item is None:
            raise KeyError(path[0])
        if len(path) > 1:
            return item['.'.join(path[1:])]
        else:
            return item

    def __len__(self):
        """
        :returns: the number of sub-elements, i.e. atoms and sub-fragments
        :rtype: int
        """
        return len(self.fragments) + len(self.atoms)

    def __iter__(self):
        """Iterate over sub-fragments first, then over the fragment's atoms.
        """
        for f in self.fragments:
            yield f.label
        for a in self.atoms:
            yield a.label

    # Override some methods from collections.Mapping for efficiency

    def values(self):
        return self.fragments + self.atoms

    def itervalues(self):
        return IT.chain(iter(self.fragments), iter(self.atoms))


# A universe description consists of
#
# - the cell shape
#
# - a list of symmetry transformations
#
# - the chemical structure of its contents
#
# - a label indicating additional conventions that this universe
#   conforms to, in particular naming conventions for atoms and molecules

class MosaicUniverse(MosaicDataItem):

    """Universe

    See the :ref:`data model documentation<mosaic-universe>` for universes.

    """

    # API properties

    @abstractproperty
    def cell_shape(self):
        """A string identifying the cell shape.
        See the :ref:`data model documentation<mosaic-universe-cell-shape>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def symmetry_transformations(self):
        """A sequence of symmetry transformations, possibly empty.
        Each symmetry tranformation is defined by a two-element tuple,
        whose first item is a rotation matrix and whose second item
        is a translation vector.
        See the :ref:`data model documentation<mosaic-universe-symmetry>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def convention(self):
        """An ASCII string naming the conventions used inside the
        definitions of fragements and atoms.
        See the :ref:`data model documentation<mosaic-universe-convention>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def molecules(self):
        """A sequence of molecule specifications. Each element is a
        two-element tuple, whose first element is a :py:class:`MosaicFragment`
        and whose second element is an integer specifying the number of
        copies of the molecule.
        See the :ref:`data model documentation<mosaic-universe-molecules>`.
        """
        raise NotImplementedError()

    # Properties that can be computed in terms of the API properties

    _cell_parameter_array_shapes = {
        'infinite': (0,),
        'cube': (),
        'cuboid': (3,),
        'parallelepiped': (3, 3)}

    @property
    def cell_parameter_array_shape(self):
        """The shape of a valid cell_parameters array
        in a :py:class:`MosaicConfiguration`.
        """
        return self._cell_parameter_array_shapes[self.cell_shape]

    @property
    def number_of_molecules(self):
        "The number of molecules in the universe."
        return sum(n for f, n in self.molecules)

    @property
    def number_of_atoms(self):
        "The number of atoms in the universe."
        return sum(n * f.number_of_atoms for f, n in self.molecules)

    @property
    def number_of_sites(self):
        "The number of sites in the universe."
        return sum(n * f.number_of_sites for f, n in self.molecules)

    @property
    def number_of_bonds(self):
        "The number of bonds in the universe."
        return sum(n * f.number_of_bonds for f, n in self.molecules)

    @property
    def number_of_template_atoms(self):
        """The number of template atoms in the universe, i.e. the
        total number of atoms in all fragment definitions. It is
        equal to the number of atoms iff all molecule repetition
        counts are 1.
        """
        return sum(f.number_of_atoms for f, n in self.molecules)

    @property
    def number_of_template_sites(self):
        """The number of template sites in the universe, i.e. the
        total number of sites in all fragment definitions. It is
        equal to the number of sites iff all molecule repetition
        counts are 1.
        """
        return sum(f.number_of_sites for f, n in self.molecules)

    # Equivalence test

    def validate_equivalence(self, other):
        """Verify the equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :raises ValueError: if other is not equivalent to self
        """
        if not isinstance(other, MosaicUniverse):
            raise TypeError("%s is not a universe" % str(type(other)))
        if self.cell_shape != other.cell_shape:
            raise ValueError("cell shapes differ: %s != %s"
                             % (repr(self.cell_shape), repr(other.cell_shape)))
        if self.convention != other.convention:
            raise ValueError("naming conventions differ: %s != %s"
                             % (repr(self.convention),
                                repr(other.convention)))
        for (r1, t1), (r2, t2) in zip(self.symmetry_transformations,
                                      other.symmetry_transformations):
            if (r1 != r2).any():
                raise ValueError("rotation matrices differ: %s != %s"
                                 % (str(r1), str(r2)))
            if (t1 != t2).any():
                raise ValueError("translation vectors differ: %s != %s"
                                 % (str(t1), str(t2)))
        for (sf, sc), (of, oc) in zip(self.molecules, other.molecules):
            sf.validate_equivalence(of)
            if sc != oc:
                raise ValueError("molecule counts differ: %d != %d" % (sc, oc))

    def is_equivalent(self, other):
        """Check for equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :returns: True if other is equivalent to self
        :rtype: bool
        """
        try:
            self.validate_equivalence(other)
            return True
        except:
            return False

    # Validation

    def validate(self):
        """Verify that the object satisfies the constraints of the
        Mosaic data model.

        :raises ValueError: if the object is not valid Mosaic data
        """
        validate_value(self.cell_shape,
                       self._cell_parameter_array_shapes.keys(),
                       "Universe.cell_shape")
        validate_label(self.convention, "Universe.convention")
        validate_sequence(self.symmetry_transformations,
                          collections.Sequence,
                          "Universe.symmetry_transformations",
                          ((lambda t: len(t) == 2,
                            "must have length 2"),
                           (lambda rot, trans:
                            getattr(rot, "shape", None) == (3, 3),
                            "rotation matrix shape is not (3, 3)"),
                           (lambda rot, trans:
                            getattr(trans, "shape", None) == (3,),
                            "translation vector shape is not (3,)"),
                           (lambda rot, trans:
                            rot.dtype == N.float64
                            and trans.dtype == N.float64,
                            "rotation and translation must be float64")))
        if self.cell_shape == "infinite" \
           and len(self.symmetry_transformations) > 0:
            raise ValueError("Symmetry transformations are allowed "
                             "only in periodic universes")
        for f, n in self.molecules:
            f.validate()
            if not isinstance(n, int):
                raise ValueError("Molecule count must be an integer")
            if n <= 0:
                raise ValueError("Molecule count must be positive")
        for property in ['cell_parameter_array_shape',
                         'number_of_molecules', 'number_of_atoms',
                         'number_of_sites', 'number_of_bonds',
                         'number_of_template_atoms',
                         'number_of_template_sites']:
            value = getattr(self, property)
            reference = getattr(MosaicUniverse, property).fget(self)
            if value != reference:
                raise ValueError("Universe.%s is %s, should be %s"
                                 % (property, str(value), str(reference)))

    # Methods based on API properties

    def recursive_atom_iterator(self):
        """An iterator over the atoms in the universe.
        """
        for fragment, count in self.molecules:
            for _ in range(count):
                for a in fragment.recursive_atom_iterator():
                    yield a

    def bond_index_array(self):
        """Returns an integer array of shape (N, 2), where N
        is the total number of bonds in the universe. The entries
        [i, 0] and [i, 1] refer to the two atoms that are connected
        by bond i. The entry [i, 0] is smaller than the entry [i, 1].

        :returns: the bond index array
        :rtype: numpy.ndarray
        """
        natoms = 0
        bonds = []
        for fragment, count in self.molecules:
            f_paths = list(fragment.recursive_atom_path_iterator())
            f_bonds = [sorted((f_paths.index(a1), f_paths.index(a2)))
                       for a1, a2, order in fragment.recursive_bond_iterator()]
            for _ in range(count):
                bonds.extend([(natoms+a1, natoms+a2) for a1, a2 in f_bonds])
                natoms += len(f_paths)
        return N.array(bonds)

    def site_to_atom_index_mapping(self):
        """
        :returns: an array whose element [s] is the atom index
                  corresponding to site index s.
        :rtype: numpy.ndarray
        """
        natoms = 0
        total = []
        for fragment, count in self.molecules:
            per_fragment = fragment.site_to_atom_index_mapping()
            f_natoms = fragment.number_of_atoms
            for _ in range(count):
                total.append(per_fragment + natoms)
                natoms += f_natoms
        return N.concatenate(total)

    def atom_to_site_index_mapping(self):
        """
        :returns: an array whose element [s] is the site index
                  corresponding to the first site of the atom with
                  atom index s.
        :rtype: numpy.ndarray
        """
        nsites = 0
        total = []
        for fragment, count in self.molecules:
            per_fragment = fragment._atom_to_site_index_mapping()
            for _ in range(count):
                total.append(per_fragment + nsites)
                nsites += per_fragment[-1]
        m = N.concatenate(total)
        m[1:] = m[:-1]
        m[0] = 0
        return m

    def template_site_to_template_atom_index_mapping(self):
        """
        :returns: an array whose element [s] is the template atom index
                  corresponding to template site index s.
        :rtype: numpy.ndarray
        """
        natoms = 0
        total = []
        for fragment, count in self.molecules:
            per_fragment = fragment.site_to_atom_index_mapping()
            f_natoms = fragment.number_of_atoms
            total.append(per_fragment + natoms)
            natoms += f_natoms
        return N.concatenate(total)

    def site_to_template_index_mapping(self):
        """
        :returns: an array whose element [s] is the template site index
                  corresponding to site index s.
        :rtype: numpy.ndarray
        """
        ntsites = 0
        total = []
        for fragment, count in self.molecules:
            f_nsites = fragment.number_of_sites
            per_fragment = N.arange(f_nsites)
            for _ in range(count):
                total.append(per_fragment + ntsites)
            ntsites += f_nsites
        return N.concatenate(total)

    def atom_to_template_index_mapping(self):
        """
        :returns: an array whose element [s] is the template atom index
                  corresponding to atom index s.
        :rtype: numpy.ndarray
        """
        ntatoms = 0
        total = []
        for fragment, count in self.molecules:
            f_natoms = fragment.number_of_atoms
            per_fragment = N.arange(f_natoms)
            for _ in range(count):
                total.append(per_fragment + ntatoms)
            ntatoms += f_natoms
        return N.concatenate(total)


# Properties associate a value with each atom or site in a
# universe. The value can be an array of any shape and any element
# type, but shape and element type are the same for all atoms.
# Properties also have a name that describes the quantitity they
# store, and the units of this quantity.
#
# Properties whose type starts with "template" are defined only for
# each atom or site in the molecule templates, not for each individual
# molecular instance. They are used for properties that are the same
# for all molecules of the same type (e.g. atomic mass).

class MosaicProperty(MosaicDataItem):

    """Property

    See the :ref:`data model documentation<mosaic-property>` for properties.

    """

    # API properties

    @abstractproperty
    def type(self):
        """A string identifying the type of the property.
        See the :ref:`data model documentation<mosaic-property-type>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def name(self):
        """An ASCII string describing the property.
        See the :ref:`data model documentation<mosaic-property-name>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def units(self):
        """A string identifying the physical units of the property.
        See the :ref:`data model documentation<mosaic-property-units>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def universe(self):
        "The :py:class:`MosaicUniverse` for which the property is defined."
        raise NotImplementedError()

    @abstractproperty
    def data(self):
        """An array containing the property's values.
        See the :ref:`data model documentation<mosaic-property-data>`.
        """
        raise NotImplementedError()

    # Properties that can be computed in terms of the API properties

    @property
    def element_shape(self):
        """The shape of the sub-array containing the propery for one
        atom or site.
        """
        return self.data.shape[1:]

    # Equivalence test

    def validate_equivalence(self, other):
        """Verify the equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :raises ValueError: if other is not equivalent to self
        """
        if not isinstance(other, MosaicProperty):
            raise TypeError("%s is not a property" % str(type(other)))
        self.universe.validate_equivalence(other.universe)
        if self.type != other.type:
            raise ValueError("types differ: %s != %s"
                             % (repr(self.type), repr(other.type)))
        if self.name != other.name:
            raise ValueError("names differ: %s != %s"
                             % (repr(self.name), repr(other.name)))
        if self.units != other.units:
            raise ValueError("units differ: %s != %s"
                             % (repr(self.units), repr(other.units)))
        if self.element_shape != other.element_shape:
            raise ValueError("element shapes differ: %s != %s"
                             % (repr(self.element_shape),
                                repr(other.element_shape)))
        if self.data.dtype != other.data.dtype:
            raise ValueError("data dtypes differ: %s != %s"
                             % (repr(self.data.dtype),
                                repr(other.data.dtype)))
        if (self.data != other.data).any():
            raise ValueError("data arrays differ")

    def is_equivalent(self, other):
        """Check for equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :returns: True if other is equivalent to self
        :rtype: bool
        """
        try:
            self.validate_equivalence(other)
            return True
        except:
            return False

    # Validation

    _allowed_dtypes = [N.int8, N.int16, N.int32, N.int64,
                       N.uint8, N.uint16, N.uint32, N.uint64,
                       N.float32, N.float64,
                       N.bool]

    _allowed_types = ["atom", "site", "template_atom", "template_site"]

    def validate(self):
        """Verify that the object satisfies the constraints of the
        Mosaic data model.

        :raises ValueError: if the object is not valid Mosaic data
        """
        validate_value(self.type, self._allowed_types, "Property.type")
        validate_label(self.name, "Property.name")
        validate_units(self.units, "Property.units")
        validate_type(self.universe, MosaicUniverse, "Property.universe")
        self.universe.validate()
        el_shape = self.element_shape
        data_shape = ({"atom": self.universe.number_of_atoms,
                       "site": self.universe.number_of_sites,
                       "template_atom": self.universe.number_of_template_atoms,
                       "template_site": self.universe.number_of_template_sites,
                       }[self.type],) + el_shape
        validate_array(self.data,
                       data_shape,
                       self._allowed_dtypes,
                       "Property.data")


# Labels associate a text string with each atom or site in a
# universe. They work much like properties, except for having a string
# value. Labels are a separate data item because the differences
# compared to numerical properties (no element shape, no unit) would
# make validation of a common data item type too complicated.
#
# Labels whose type starts with "template" are defined only for
# each atom or site in the molecule templates, not for each individual
# molecular instance.

class MosaicLabel(MosaicDataItem):

    """Label

    See the :ref:`data model documentation<mosaic-label>` for labels.

    """

    # API properties

    @abstractproperty
    def type(self):
        """A string identifying the type of the label.
        See the :ref:`data model documentation<mosaic-label-type>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def name(self):
        """An ASCII string describing the label.
        See the :ref:`data model documentation<mosaic-label-name>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def universe(self):
        "The :py:class:`MosaicUniverse` for which the property is defined."
        raise NotImplementedError()

    @abstractproperty
    def strings(self):
        """A sequence of strings representing a label for each atom or site.
        See the :ref:`data model documentation<mosaic-label-strings>`.
        """
        raise NotImplementedError()

    # Equivalence test

    def validate_equivalence(self, other):
        """Verify the equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :raises ValueError: if other is not equivalent to self
        """
        if not isinstance(other, MosaicLabel):
            raise TypeError("%s is not a label" % str(type(other)))
        self.universe.validate_equivalence(other.universe)
        if self.type != other.type:
            raise ValueError("types differ: %s != %s"
                             % (repr(self.type), repr(other.type)))
        if self.name != other.name:
            raise ValueError("names differ: %s != %s"
                             % (repr(self.name), repr(other.name)))
        for s1, s2 in zip(self.strings, other.strings):
            if s1 != s2:
                raise ValueError("labels differ: %s != %s"
                                 % (repr(s1), repr(s2)))

    def is_equivalent(self, other):
        """Check for equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :returns: True if other is equivalent to self
        :rtype: bool
        """
        try:
            self.validate_equivalence(other)
            return True
        except:
            return False

    # Validation

    _allowed_types = ["atom", "site", "template_atom", "template_site"]

    def validate(self):
        """Verify that the object satisfies the constraints of the
        Mosaic data model.

        :raises ValueError: if the object is not valid Mosaic data
        """
        validate_value(self.type, self._allowed_types, "Label.type")
        validate_label(self.name, "Label.name")
        validate_type(self.universe, MosaicUniverse, "Label.universe")
        self.universe.validate()
        validate_sequence(self.strings, str, "Label.strings")
        nstrings = {"atom": self.universe.number_of_atoms,
                    "site": self.universe.number_of_sites,
                    "template_atom": self.universe.number_of_template_atoms,
                    "template_site": self.universe.number_of_template_sites,
                    }[self.type]
        if len(self.strings) != nstrings:
            raise ValueError("incorrect number of strings")
        for s in self.strings:
            validate_ascii_string(s, "label")


# A configuration specifies a coordinate for each site and
# the shape and size of the cell.

class MosaicConfiguration(MosaicDataItem):

    """Configuration

    See the :ref:`data model documentation<mosaic-configuration>`
    for configurations.

    """

    # API properties

    @abstractproperty
    def universe(self):
        "The :py:class:`MosaicUniverse` for which the property is defined."
        raise NotImplementedError()

    @abstractproperty
    def cell_parameters(self):
        """An array containing the parameters defining the shape and size
        of the cell.
        See the :ref:`data model documentation<mosaic-configuration-cp>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def positions(self):
        """An array containing an (x, y, z) position for each site.
        See the :ref:`data model documentation<mosaic-configuration-pos>`.
        """
        raise NotImplementedError()

    # Methods based on API properties

    def lattice_vectors(self):
        """
        :returns: a sequence of arrays of shape (3,) containing the lattice
                  vectors for the simulation cell.
        :rtype: tuple
        """
        return {'infinite': lambda p: (),
                'cube': lambda p: (N.array([p, 0., 0.], dtype=p.dtype),
                                   N.array([0., p, 0.], dtype=p.dtype),
                                   N.array([0., 0., p], dtype=p.dtype)),
                'cuboid': lambda p: (N.array([p[0], 0., 0.], dtype=p.dtype),
                                     N.array([0., p[1], 0.], dtype=p.dtype),
                                     N.array([0., 0., p[2]], dtype=p.dtype)),
                'parallelepiped': lambda p: tuple(N.array(v) for v in p),
                }[self.universe.cell_shape](self.cell_parameters)

    def cell_volume(self):
        """
        :returns: the volume the simulation cell, or None for
                  infinite universes
        :rtype: float
        """
        return {'infinite': lambda p: None,
                'cube': lambda p: float(p*p*p),
                'cuboid': lambda p: p[0]*p[1]*p[2],
                'parallelepiped': lambda p: abs(LA.det(p)),
                }[self.universe.cell_shape](self.cell_parameters)

    # Equivalence test

    def validate_equivalence(self, other):
        """Verify the equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :raises ValueError: if other is not equivalent to self
        """
        if not isinstance(other, MosaicConfiguration):
            raise TypeError("%s is not a configuration" % str(type(other)))
        self.universe.validate_equivalence(other.universe)
        if (self.cell_parameters != other.cell_parameters).any():
            raise ValueError("cell parameters differ: %s != %s"
                             % (repr(self.cell_parameters),
                                repr(other.cell_parameters)))
        if self.cell_parameters.dtype != other.cell_parameters.dtype:
            raise ValueError("cell parameter dtypes differ: %s != %s"
                             % (repr(self.cell_parameters.dtype),
                                repr(other.cell_parameters.dtype)))
        if self.positions.dtype != other.positions.dtype:
            raise ValueError("position dtypes differ: %s != %s"
                             % (repr(self.positions.dtype),
                                repr(other.positions.dtype)))
        if (self.positions != other.positions).any():
            raise ValueError("position arrays differ")

    def is_equivalent(self, other):
        """Check for equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :returns: True if other is equivalent to self
        :rtype: bool
        """
        try:
            self.validate_equivalence(other)
            return True
        except:
            return False

    # Validation

    _allowed_dtypes = [N.float32, N.float64]

    def validate(self):
        """Verify that the object satisfies the constraints of the
        Mosaic data model.

        :raises ValueError: if the object is not valid Mosaic data
        """
        validate_type(self.universe, MosaicUniverse, "Configuration.universe")
        self.universe.validate()
        validate_array(self.cell_parameters,
                       self.universe.cell_parameter_array_shape,
                       self._allowed_dtypes,
                       "Configuration.cell_parameters")
        validate_array(self.positions,
                       (self.universe.number_of_sites, 3),
                       self._allowed_dtypes,
                       "Configuration.positions")
        if self.cell_parameters.dtype != self.positions.dtype:
            raise ValueError("Configuration.cell_parameters and "
                             "Configuration.positions must have same dtypes")


# Selections specify a subset of atoms or sites, either in the templates
# or in the molecules of a universe.

class MosaicSelection(MosaicDataItem):

    """Selection

    See the :ref:`data model documentation<mosaic-selection>` for selections.

    """

    # API properties

    @abstractproperty
    def type(self):
        """A string identifying the type of the selection.
        See the :ref:`data model documentation<mosaic-selection-type>`.
        """
        raise NotImplementedError()

    @abstractproperty
    def universe(self):
        "The :py:class:`MosaicUniverse` for which the selection is defined."
        raise NotImplementedError()

    @abstractproperty
    def indices(self):
        """An array containing the indices of the contained atoms or sites.
        See the :ref:`data model documentation<mosaic-selection-indices>`.
        """
        raise NotImplementedError()

    # Properties that can be computed in terms of the API properties

    @property
    def number_of_atoms(self):
        """The number of atoms in the selection.
        """
        indices = self.indices
        if self.type == "atom":
            return len(indices)
        elif self.type == "template_atom":
            natoms = 0
            ntatoms = 0
            for fragment, count in self.universe.molecules:
                nta = fragment.number_of_atoms
                natoms += count * N.sum((indices >= ntatoms)
                                        & (indices < ntatoms+nta))
                ntatoms += nta
            return natoms
        else:
            raise TypeError("number of atoms undefined in site selection")

    @property
    def number_of_sites(self):
        """The number of sites in the selection.
        """
        indices = self.indices
        if self.type == "atom":
            sites_per_atom = [a.number_of_sites
                              for a in self.universe.recursive_atom_iterator()]
            # If 'indices' is an immutable array, N.take crashes due to
            # a numpy bug. Conversion to a plain array prevents this.
            # See Github issue #3758 for numpy/numpy.
            return N.sum(N.take(sites_per_atom, N.array(indices)))
        elif self.type == "template_atom":
            sites_per_atom = [a.number_of_sites
                              for a in self.universe.recursive_atom_iterator()]
            nsites = 0
            ntatoms = 0
            for fragment, count in self.universe.molecules:
                nta = fragment.number_of_atoms
                mask = (indices >= ntatoms) & (indices < ntatoms+nta)
                ntatoms += nta
                # Conversion to plain arrays is required for the same
                # reason as above.
                nsites += count * N.sum(N.take(sites_per_atom,
                                               N.repeat(N.array(indices),
                                                        N.array(mask))))
            return nsites
        elif self.type == "site":
            return len(indices)
        elif self.type == "template_site":
            nsites = 0
            ntsites = 0
            for fragment, count in self.universe.molecules:
                nts = fragment.number_of_sites
                nsites += count * N.sum((indices >= ntsites)
                                        & (indices < ntsites+nts))
                ntsites += nts
            return nsites

    # Equivalence test

    def validate_equivalence(self, other):
        """Verify the equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :raises ValueError: if other is not equivalent to self
        """
        if not isinstance(other, MosaicSelection):
            raise TypeError("%s is not a selection" % str(type(other)))
        self.universe.validate_equivalence(other.universe)
        if self.type != other.type:
            raise ValueError("types differ: %s != %s"
                             % (repr(self.type), repr(other.type)))
        if (self.indices != other.indices).any():
            raise ValueError("indices differ")

    def is_equivalent(self, other):
        """Check for equivalence of Python objects representing
        Mosaic data. The two objects can belong to different models;
        only the common API functionality is used for the test.

        :parameter other: an arbitrary Python object
        :returns: True if other is equivalent to self
        :rtype: bool
        """
        try:
            self.validate_equivalence(other)
            return True
        except:
            return False

    # Validation

    _allowed_types = ["atom", "site", "template_atom", "template_site"]

    def validate(self):
        """Verify that the object satisfies the constraints of the
        Mosaic data model.

        :raises ValueError: if the object is not valid Mosaic data
        """
        validate_value(self.type, self._allowed_types, "Selection.type")
        validate_type(self.universe, MosaicUniverse, "Selection.universe")
        max_index = {"atom": self.universe.number_of_atoms,
                     "site": self.universe.number_of_sites,
                     "template_atom": self.universe.number_of_template_atoms,
                     "template_site": self.universe.number_of_template_sites,
                    }[self.type]
        validate_indices(self.indices, max_index, "Selection.indices")

# Validation functions

def validate_type(obj, cls, text):
    if isinstance(obj, cls):
        return
    raise TypeError("%s must be of type %s (is %s)"
                    % (text, cls.__name__, str(type(obj))))


def validate_value(value, allowed_values, name):
    if not value in allowed_values:
        raise ValueError("%s must be one of %s"
                         % (name, ", ".join([str(v) for v in allowed_values])))


def validate_ascii_string(s, text):
    if not mosaic.utility.isascii(s):
        raise ValueError("non-ASCII string in %s" % text)


def validate_label(label, text):
    validate_type(label, str, text)
    if len(label) > 32767:
        raise ValueError("%s too long, must be <= 32767 characters" % text)
    for c in label:
        if not c in _allowed_in_labels:
            raise ValueError("illegal character '%s' in %s"
                             % (c, text))

_allowed_in_labels = 'abcdefghijklmnopqrstuvwxyz' + \
                     'ABCDEFGHIJKLMNOPQRSTUVWXYZ' + \
                     '0123456789!#$%&?@^_~+-*/=,()[]' + "'"


def validate_array(array, shape, allowed_dtypes, text):
    # Don't check for N.ndarray in order to allow on-disk
    # arrays (HDF5, netCDF)
    try:
        a_shape = array.shape
        a_dtype = array.dtype
    except AttributeError:
        raise TypeError("%s must be an array" % text)
    if shape is not None and a_shape != shape:
        raise ValueError("%s must have shape %s"
                         % (text, str(shape)))
    if a_dtype not in allowed_dtypes:
        raise TypeError(" %s must have element type %s"
                        % (text,
                           " or ".join(str(t) for t in allowed_dtypes)))


def validate_sequence(obj, el_cls, text, additional_tests=()):
    if not (isinstance(obj, collections.Sequence) and
            all(isinstance(item, el_cls) for item in obj)):
        raise TypeError("%s must be a sequence of %s elements"
                        % (text, el_cls.__name__))
    for test_fn, text in additional_tests:
        for el in obj:
            if not test_fn(el):
                raise ValueError("%s: %s" % (str(el), text))

def validate_indices(indices, max_index, text):
    validate_array(indices, None,
                   [N.uint8, N.uint16, N.uint32, N.uint64],
                   text)
    if len(indices.shape) != 1:
        raise ValueError("index array not 1d")
    if (indices >= max_index).any():
        raise ValueError("index too large")
    if len(indices) > 1:
        d = indices[1:]-indices[:-1]
        if (d <= 0).any():
            raise ValueError("indices not sorted")

# The unit validator accepts a very limited syntax, which
# should be sufficient: a unit is defined by a string of
# unit names with an optional numeric suffix indication a power.
# The unit list can include integers or decimal fractions.
#
# Examples: "nm ps-1" (velocity), "nm2" (area), "kJ mol-1" (energy)
#           "0.1 nm" (length), "60 s" (time)

def validate_units(unit_spec, text):
    validate_type(unit_spec, str, text)
    first = True
    for unit in unit_spec.split():
        # number (integer or decimal fraction)
        if first:
            first = False
            if (re.match("^[1-9][0-9]*$", unit) \
                or re.match("^0\.[0-9]+$", unit)):
                continue
        # symbol+exponent
        m = re.match("^(?P<symbol>[a-zA-Z]+)([-]?[0-9]+)?$", unit)
        if m and m.group('symbol') in _allowed_units:
            continue
        raise ValueError("invalid unit '%s' in %s" % (unit, unit_spec))

_allowed_units = \
    ["pm",   "Ang",  "nm",   "um",   "mm",   "m",    
     "fs",   "ps",   "ns",   "us",   "ms",   "s",    
     "amu",  "g",    "kg",   
     "mol",  
     "J",    "kJ",   "cal",  "kcal", "eV",   
     "K",    
     "Pa",   "kPa",  "MPa",  "GPa",  "atm",  "bar",  "kbar", 
     "e",    "C",    "A",    "V",    
     "rad",  
     "c",    "h",    "me"]

# Validation as a test function rather than raising exceptions

def is_valid(obj):
    try:
        obj.validate()
        return True
    except:
        return False
