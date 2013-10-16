.. Written by Konrad Hinsen
.. License: CC-BY 3.0

.. index::
   single: PDB

Mosaic PDB convention
#####################

Mosaic can be used to store molecular models from the
`Protein Data Bank <http://www.wwpdb.org/>`_ (PDB). The main
application is to use such models as the starting point for
molecular simulations. The following conventions describe
how a PDB structure is stored in terms of Mosaic data items.
Note that only the structure itself can be stored, but not
experimental data (structure factors etc.) or metadata
describing the experiment or the refinement process.

The PDB's official data format is called
`PDBx/mmCIF <http://mmcif.pdb.org/>`_. In the conversion from
PDBx/mmCIF to Mosaic, as much information as possible is
transposed without modification. In particular, residue and atom
names are the same.


.. index::
   single: crystallographic structures

Crystallographic structures
---------------------------

A crystallographic structure is represented by two required
data items:

 - A universe defining the molecular structures and, in the
   case of crystals, the symmetries. The atoms in the universe
   have multiple sites if the PDB structure contains alternate
   locations.

 - A configuration providing the positions for all sites and
   the shape of the unit cell in the case of crystals.

Additional information from the PDB entry can be provided by
optional data items:

 - The occupancy of each site can be provided as a
   :ref:`property<mosaic-property>` of
   :ref:`type<mosaic-property-type>` "site" or "template_site"
   with an empty :ref:`units<mosaic-property-units>` string.
   Each value is a scalar of type "float32" or "float64" in
   the interval [0..1]. If no occupancy values are provided,
   the occupancy of all sites is assumed to be 1.

 - An anisotropic displacement parameter for each site can be provided
   as a :ref:`property<mosaic-property>` of
   :ref:`type<mosaic-property-type>` "site" or "template_site".  A
   valid :ref:`units<mosaic-property-units>` string must be provided,
   the preferred units are "nm2".  Each value is an array of shape "6"
   and of type "float32" or "float64", the order of the elements is
   [1][1], [2][2], [3][3], [2][3], [1][3], [1][2]. For the precise
   definition of the anisotropic displacement parameters, see the PDB
   documentation for items ``_atom_site.aniso_U[1][1]`` to
   ``_atom_site.aniso_U[3][3]``.

 - An isotropic displacement parameter for each site can be provided
   as a :ref:`property<mosaic-property>` of
   :ref:`type<mosaic-property-type>` "site" or "template_site".  A
   valid :ref:`units<mosaic-property-units>` string must be provided,
   the preferred units are "nm2".  Each value is a scalar of type
   "float32" or "float64". An isotropic displacement parameter of
   value ``x`` is equivalent to an anisotropic displacement parameter
   of value ``[x x x 0 0 0]``.

If anisotropic displacement parameters are provided, then no isotropic
displacement parameters may be given, in order to prevent
incoherencies in the data.


.. index::
   single: NMR structures

NMR structures
--------------

An NMR structure is represented by the following data items:

 - A universe defining the molecular structures.

 - One configuration per model contained in the PDB entry.
