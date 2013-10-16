.. Written by Konrad Hinsen
.. License: CC-BY 3.0

.. index::
   single: HDF5

Mosaic in HDF5 files
####################

`HDF5 <http://www.hdfgroup.org/HDF5/>`_ files contain a tree structure
whose leaves are datasets. Non-leaf nodes are called groups and work
much like a directory in a file system. Each dataset is an array whose
elements can be numbers or characters, but also compound data types
(similar to record types in various programming languages) or
fixed-size arrays of numbers or characters. Groups and datasets can
have metadata tags called attributes. Each group or dataset can be
identified by a path specifying how to reach it from the root group of
a file. However, a path is not necessarily unique, because HDF5
provides links that effectively put a single node in several places in
the tree. Moreover, HDF5 provides a data type "reference" that
allows to refer to a node or a subset of a dataset.

The design criteria for the HDF5 representation of Mosaic data were
efficiency of storage and ease of use from low-level languages such as
C or Fortran. As much as possible, Mosaic data is stored as arrays of
numbers.

HDF5 has two string layouts: fixed-size strings (character arrays) and
variable-length strings. The two layouts are not interchangeable.
In order to facilitate software development, Mosaic uses only
variable-length strings.


Mosaic data items in a HDF5 file
--------------------------------

A Mosaic data item in an HDF5 file can be a dataset or a group
containing multiple datasets. It is identified by four attributes, all
of which are required:

  DATA_MODEL
     a variable-length string with the value "MOSAIC"

  DATA_MODEL_MAJOR_VERSION
     an integer

  DATA_MODEL_MINOR_VERSION
     an integer

  MOSAIC_DATA_TYPE
     a variable-length string

References between data items are stored as attributes whose value
is an HDF5 object reference.


Universes
.........

A universe is stored as a group containing several datasets. The
datasets :ref:`convention<mosaic-universe-convention>` and
:ref:`cell_shape<mosaic-universe-cell-shape>` are variable-length
strings. The :ref:`symmetry transformation<mosaic-universe-symmetry>`
list is stored in dataset ``symmetry_transformations`` as a
one-dimensional array, possibly empty, whose elements are of a
compound data type with fields

  rotation
    A 3x3 array of float64 numbers.

  translation
    An array of 3 float64 numbers.

The fragment tree is stored in several arrays. All string values are
stored in a one-dimensional array dataset ``symbols`` whose elements
are variable-length strings. In the fragment tree, the strings can
then be replaced by integers, which are indices into this symbol list.
Ideally, each string is stored only once in the symbol array, though
this is not a requirement. The tree data structure is stored in two
integer arrays, ``fragments`` and ``atoms``. Each node (fragment or
atom) has one array entry, which is a compound data type whose fields
are unsigned integers. Any size of unsigned integer can be used,
but the same type must be used everywhere for a given universe group.

The entries of the ``fragments`` array have the fields

  parent_index
    The index of the parent node in the ``fragments`` array.
    A value of 0 indicates a root node, which has no parent.

  label_symbol_index
    The index of the :ref:`label<mosaic-fragment-label>` in the
    ``symbols`` array.

  species_symbol_index
    The index of the :ref:`species<mosaic-fragment-species>` in the
    ``symbols`` array.

  number_of_fragments
    The number of sub-fragments. This is redundant information,
    provided to facilitate reading fragment-related information
    without analyzing the whole fragment tree.

The first entry (index 0) of the ``fragments`` array is unused,
in order to allow an index value of 0 to stand for "no parent".
The ``fragments`` array thus has one more entry than the number
of fragments in the universe.

The entries of the ``atoms`` array have the fields

  parent_index
    The index of the parent node in the ``fragments`` array.

  label_symbol_index
    The index of the :ref:`label<mosaic-atom-label>` in the
    ``symbols`` array.

  type_symbol_index
    The index of the :ref:`type<mosaic-atom-type>` in the
    ``symbols`` array.

  name_symbol_index
    The index of the :ref:`name<mosaic-atom-name>` in the
    ``symbols`` array.

  number_of_sites
    The :ref:`number of sites<mosaic-atom-nsites>`.
 
The entries of the ``bonds`` array have the fields

  atom_index_1
    The index of the first atom in the ``atoms`` array.

  atom_index_2
    The index of the second atom in the ``atoms`` array.

  bond_order_symbol_index
    The index of the bond-order label in the ``symbols`` array.

The entries of the ``molecules`` array have a large number of
redundant fields (all but the first two) that are provided to allow
atoms be attributed to molecules without analyzing the full fragment
tree.

  fragment_index
    The index of the fragment node in the ``fragments`` array.

  number_of_copies
    The number of copies of the molecule in the universe.

  first_atom_index
    The index of the first atom in the ``atoms`` array.

  number_of_atoms
    The number of atoms in the molecule.

  first_bond_index
    The index of the first bond in the ``bonds`` array.

  number_of_bonds
    The number of bonds in the molecule.

  first_site_index
    The index of the first site of the molecule.

  number_of_sites
    The number of sites in the molecule.

Since the atoms, sites, and bonds of a molecule have consecutive
indices, the redundant "first_index" and "number_of" values are
sufficient to locate atoms, sites, and bonds for each molecule.
For many applications this is sufficient, making it unnecessary
to use the ``fragments`` array.

Finally, the array ``polymers`` has one entry for each polymer
fragment in the universe. Its fields are

  fragment_index
    The index of the fragment node in the ``fragments`` array.

  polymer_type_symbol_index
    The index of the polymer-type label in the ``symbols`` array.

If the universe has no polymer fragments, the dataset ``polymers``
may be omitted.

Configurations
..............

A configuration is stored as a group containing two datasets:
``positions`` (required) and ``cell_parameters`` (required if the
universe's cell shape is not "infinite"). The reference to the
universe is stored in the attribute ``universe`` of the group.

The dataset ``positions`` is a one-dimensional array whose length is
equal to the number of sites in the universe. Its elements are
one-dimensional arrays of length 3 whose elements are of type
"float32" or "float64".

The dataset ``cell_parameters`` is an array whose elements are of
type "float32" or "float64", and whose shape is defined in the
:ref:`specification<mosaic-configuration-cp>`.

Properties
..........

A property data item is stored as a dataset that is a one-dimensional
array whose length is equal to the number of atoms or sites in the
universe or the universe's fragment list. Each element of this array
is an array whose shape and element type is defined by the property's
:ref:`data<mosaic-property-data>`.
The :ref:`reference<mosaic-property-universe>` to the
universe is stored in the attribute ``universe`` of the group.
The property's :ref:`name<mosaic-property-name>` and
:ref:`units<mosaic-property-units>` are stored in attributes
of the same name as variable-length strings.
The property's :ref:`type<mosaic-property-type>` is stored in the
attribute ``property_type``, also as a variable-length string.

Labels
......

A label data item is stored as a dataset that is a one-dimensional
array whose length is equal to the number of atoms or sites in the
universe or the universe's fragment list. Each element of this array
is a variable-length string.
The :ref:`reference<mosaic-label-universe>` to the
universe is stored in the attribute ``universe`` of the group.
The label's :ref:`name<mosaic-label-name>` is stored in the attribute
``name`` as a variable-length string.
The label's :ref:`type<mosaic-label-type>` is stored in the
attribute ``label_type``, also as a variable-length string.


Selections
..........

A selection is stored as a dataset that is a one-dimensional array
of integers.
The :ref:`reference<mosaic-selection-universe>` to the
universe is stored in the attribute ``universe`` of the group.
The selection's :ref:`type<mosaic-selection-type>` is stored in the
attribute ``selection_type`` as a variable-length string.
