.. Written by Konrad Hinsen
.. License: CC-BY 3.0

The file converter ``mosaic``
=============================

The command-line script ``mosaic`` can convert between the
XML and HDF5 file formats and also import entries from the
Protein Data Bank (PDB). The command-line syntax is

  ``mosaic`` <source> <destination>

where <source> and <destination> are filenames whose extensions
are ``.xml`` or ``.h5``. The extension of <source> can also
be ``.cif`` for an mmCIF input file. Finally, <source> can
be a four-letter PDB code, in which case the corresponding
entry is downloaded automatically.



