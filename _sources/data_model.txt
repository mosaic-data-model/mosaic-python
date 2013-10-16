.. Written by Konrad Hinsen
.. License: CC-BY 3.0


The Mosaic data model
#####################

.. index::
   single: data item

Mosaic data items
-----------------

Mosaic data items are the smallest units of information that can be
stored in files. A Mosaic file can contain any number of data items,
each of which has a unique identifier. Identifiers are text strings
whose exact specification may vary between file formats. ASCII encoded
identifiers are allowed in all file formats and are therefore
preferred.


Some definitions used in the following
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   int8, int16, int32, int64
      A signed integer occupying 8, 16, 32, or 64 bits in memory.

   uint8, uint16, uint32, uint64
      An unsigned integer occupying 8, 16, 32, or 64 bits in memory.

   bool
      A truth value that is either True or False.

   float32
      An IEEE single-precision floating point number, occupying
      32 bits in memory.

   float64
      An IEEE double-precision floating point number, occupying
      64 bits in memory.

   label
      A text string in ASCII encoding containing at most 32767
      characters which are each an upper or lower case letter, a
      digit, or one of the punctuation characters in the string
      "\!\#\$\%\&\?\@\^\_\~\+\-\*\/\=\,\(\)\[\]\'".  This includes all
      the ASCII punctuation characters except for the dot. Spaces and
      control characters are not allowed.

   list
      An ordered collection of items. Corresponding data structures in
      programming languages are typically called list, vector, or
      array.

   array
      An n-dimensional ordered collection whose elements are of
      identical type.

   atom
      A point-like particle in a molecular simulation.  May represent
      a real atom, a united-atom particle in a coarse-grained model,
      or a dummy interaction point with no physical properties.

   site
      A point in space related to an atom. An atom has at least
      one site. Atoms with multiple sites can be used for representing
      quantum models (path integrals, wave functions, ...),
      atoms with multi-modal position distributions (as used in
      crystallography), etc.


.. index::
   pair: data item; universe

.. _mosaic-universe:

Data item type "universe"
-------------------------

The universe is the most central Mosaic data item because all the
other ones require a reference to a universe, because their data
contents can only be interpreted meaningfully in the context of their
universe.

A universe describes a molecular system, defining

  1. the chemical structure of the molecules it contains

  2. the topology of the whole system (periodic boundary conditions etc.)

  3. symmetries, if required


A Mosaic universe contains:

  .. _mosaic-universe-cell-shape:

  - a cell shape field, whose value is "infinite", "cube",
    "cuboid", or "parallelepiped"

  .. _mosaic-universe-convention:

  - a convention field, whose value is an ASCII-encoded text string
    naming a convention for atom names, decomposition of standard
    groups (e.g. amino acid residues) into subgroups and atoms, etc.

  .. _mosaic-universe-symmetry:

  - a (possibly empty) list of symmetry transformation. Each symmetry
    transformation is defined by a rotation matrix and a translation
    vector. The full system consists of the explicitly represented
    atoms and molecules plus their images obtained by applying all the
    symmetry transformations. Symmetry transformations are defined in
    fractional coordinates and therefore allowed only for periodic
    universes.

  .. _mosaic-universe-molecules:

  - a list of molecules, each molecule being defined by a
    (fragment, count) pair, where count is a positive integer.
    Fragments are defined below.

.. _mosaic-fragment:

A Mosaic fragment is not a data item, because it cannot be written to
a file in isolation. Fragments exist only as part of a universe
definition. A fragment loosely corresponds to the concept of a
functional group or moiety in chemistry, but its definition covers a
much wider range of chemical structures: the extreme use cases for
fragments in Mosaic are single atoms and whole molecules.

A fragment contains the following information:

  .. _mosaic-fragment-label:

  - a label field, whose value is a label that identifies the
    fragment uniquely inside its parent fragment (if any parent
    fragment exists)

  .. _mosaic-fragment-species:

  - a species field, whose value is a label that describes the
    chemical entity described by the fragment

  .. _mosaic-fragment-polymer:

  - a boolean field is_polymer, whose value is "true" if the fragment
    describes a polymer. A polymer fragment has an empty atom list,
    i.e. it contains only sub-fragments.  A polymer fragment also has
    an additional polymer_type field, whose value is "",
    "polypeptide", "polyribonucleotide", "polydeoxyribonucleotide", or
    "polynucleotide". The empty string is used for polymers that are
    not of any other type, or for polymers of unknown type.

  .. _mosaic-fragment-fragments:

  - a (possibly empty) list of sub-fragments

  .. _mosaic-fragment-atoms:

  - a (possibly empty) list of atoms

  .. _mosaic-fragment-bonds:

  - a (possibly empty) list of bonds


.. _mosaic-atom:

An atom is described by:

  .. _mosaic-atom-label:

  - a label field, whose value is a label that identifies the
    atom uniquely inside its parent fragment. Each label inside
    a parent fragment can name an atom //or// a sub-fragment,
    but not both.

  .. _mosaic-atom-type:

  - a type field, whose value is "element", "cgparticle", "dummy", or
    "". The empty string is used for any type other then the explicitly
    named ones, and for atoms of unknown type. The type "element" refers
    to a physical atom with a well-defined chemical element. The type
    "cgparticle" refers to coarse-grained particles that represent several
    physical atoms. The type "dummy" refers to interaction sites that
    have no physical reality.

  .. _mosaic-atom-name:

  - a name field, whose value is a label that describes the chemical
    nature of the atom. For atoms of type "element", it must be the
    chemical element symbol, with the first letter upper-case and
    the second letter, if one exists, in lower-case.

  .. _mosaic-atom-nsites:

  - a number of sites field, whose value is a positive integer.


.. _mosaic-bonds:

A bond is described by two atom references and a bond order
specification, whose value is "", "single", "double", "triple",
"quadruple", or "aromatic".  The empty string is used for bonds of any
other order, or for bonds of unknown order. Bonds must be defined at
the level of the smallest possible fragment that includes both atoms
implied in the bond. In other words, it must be possible to check if
two atoms in a fragment are linked by a bond without looking at parent
fragments.

An atom reference is an ASCII-encoded text string naming an atom
relative to the current fragment by the sequence of labels that define
the path to the atom. The labels in the sequence are separated by a
dot.


.. index::
   pair: data item; configuration

.. _mosaic-configuration:

Data item type "configuration"
------------------------------

A configuration contains:

  .. _mosaic-configuration-universe:

  - a reference to a universe

  .. _mosaic-configuration-pos:

  - one position vector for each site in the universe

  .. _mosaic-configuration-cp:

  - for universes with a bounded cell, the parameters of the cell,
    stored as an array whose shape is determined by the universe's
    cell shape: an empty shape vector (i.e. the array is a scalar)
    for "cube", shape (3) for "cuboid", and (3,3) for "parallelepiped".

The elements of the position vectors and the cell parameters are
floats of the same precision, either float32 or float64.


.. index::
   pair: data item; property

.. _mosaic-property:

Data item type "property"
-------------------------

A property contains:

  .. _mosaic-property-type:

  - a type field, whose value is "atom", "site", "template_atom",
    or "template_site"

  .. _mosaic-property-universe:

  - a reference to a universe

  .. _mosaic-property-data:

  - one array (see details below) for each

     * atom in the universe, if the type field is "atom"

     * site in the universe, if the type field is "site"

     * atom in the fragment templates, if the type field is "template_atom"

     * site in the fragment templates, if the type field is "template_site"

  .. _mosaic-property-name:

  - a name field, whose value is a label

  .. _mosaic-property-units:
    
  - a units field, see details below

The arrays for each atom or site have identical shapes and their
elements identical types. The type can be int8, int16, int32, int64,
uint8, uint16, uint32, uint32, uint64, float32, float64, or bool.

The value of the units field is a text string in ASCII encoding.  It
contains a sequence of unit factors separated by a space.  A unit
factor is either a number (an integer or a decimal fraction) or a unit
symbol optionally followed by a non-zero integer which indicates the
power to which this factor is taken.  Examples:

  - "nm3" stands for cubic nanometers

  - "nm ps-1" stands for nanometers per picosecond

  - "60 s" stands for a minute

Each unit symbol may occur only once in the units field. There may also
be at most one numeric factor, which must be the first one.

The following unit symbols may be used:

   +-------------+------+-----------------+
   | Length      | pm   | picometer       |
   +             +------+-----------------+
   |             | Ang  | Ångström        |
   +             +------+-----------------+
   |             | nm   | nanometer       |
   +             +------+-----------------+
   |             | um   | micrometer      |
   +             +------+-----------------+
   |             | mm   | millimeter      |
   +             +------+-----------------+
   |             | m    | meter           |
   +-------------+------+-----------------+
   | Time        | fs   | femtosecond     |
   +             +------+-----------------+
   |             | ps   | picosecond      |
   +             +------+-----------------+
   |             | ns   | nanosecond      |
   +             +------+-----------------+
   |             | us   | microsecond     |
   +             +------+-----------------+
   |             | ms   | millisecond     |
   +             +------+-----------------+
   |             | s    | second          |
   +-------------+------+-----------------+
   | Mass        | amu  | gram/mole       |
   +             +------+-----------------+
   |             | g    | gram            |
   +             +------+-----------------+
   |             | kg   | kilogram        |
   +-------------+------+-----------------+
   | Quantity    | mol  | mole            |
   +-------------+------+-----------------+
   | Energy      | J    | joule           |
   +             +------+-----------------+
   |             | kJ   | kilojoule       |
   +             +------+-----------------+
   |             | cal  | calorie         |
   +             +------+-----------------+
   |             | kcal | kilocalorie     |
   +             +------+-----------------+
   |             | eV   | electron-volt   |
   +-------------+------+-----------------+
   | Temperature | K    | Kelvin          |
   +-------------+------+-----------------+
   | Pressure    | Pa   | pascal          |
   +             +------+-----------------+
   |             | kPa  | kilopascal      |
   +             +------+-----------------+
   |             | MPa  | megapascal      |
   +             +------+-----------------+
   |             | GPa  | giggapascal     |
   +             +------+-----------------+
   |             | atm  | atmosphere      |
   +             +------+-----------------+
   |             | bar  | bar             |
   +             +------+-----------------+
   |             | kbar | kilobar         |
   +-------------+------+-----------------+
   | Electrical  | e    | proton charge   |
   + units       +------+-----------------+
   |             | C    | coulomb         |
   +             +------+-----------------+
   |             | A    | ampere          |
   +             +------+-----------------+
   |             | V    | volt            |
   +-------------+------+-----------------+
   | Angles      | rad  | radian          |
   +             +------+-----------------+
   |             | deg  | degree          |
   +-------------+------+-----------------+
   | Constants   | c    | speed of light  |
   +             +------+-----------------+
   |             | h    | Planck constant |
   +             +------+-----------------+
   |             | Nav  | Avogadro number |
   +             +------+-----------------+
   |             | me   | electron mass   |
   +-------------+------+-----------------+


.. index::
   pair: data item; label

.. _mosaic-label:

Data item type "label"
----------------------

A label contains:

  .. _mosaic-label-type:

  - a type field, whose value is "atom", "site", "template_atom",
    or "template_site"

  .. _mosaic-label-universe:

  - a reference to a universe

  .. _mosaic-label-strings:

  - one text string in ASCII encoding for each

     * atom in the universe, if the type field is "atom"

     * site in the universe, if the type field is "site"

     * atom in the fragment templates, if the type field is "template_atom"

     * site in the fragment templates, if the type field is "template_site"

  .. _mosaic-label-name:

  - a name field, whose value is a label


.. index::
   pair: data item; selection

.. _mosaic-selection:

Data item type "selection"
--------------------------

A selection contains:

  .. _mosaic-selection-type:

  - a type field, whose value is "atom", "site", "template_atom",
    or "template_site"

  .. _mosaic-selection-universe:

  - a reference to a universe

  .. _mosaic-selection-indices:

  - an array whose values are the indices of the atoms
    or sites that are part of the selection.

The index array is one-dimensional and the type of its elements is one
of the unsigned integer types: uint8, uint16, uint32, uint32, uint64.
The indices are stored in monotonously increasing order with no index
being listed more than once.
