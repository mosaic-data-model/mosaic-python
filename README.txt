Mosaic is a modular data model for molecular simulation applications.

The Mosaic Python library provides an in-memory representation of
Mosaic data items plus input/output in HDF5 and XML formats. A
command-line tool provides conversion between the HDF5 and XML
representations as well as data import from the Protein Data Bank
(PDB).

For the definition of the Mosaic data model, and a description
of the Python API, see the documentation in the directory 'doc'.

The Mosaic Python library requires Python 2.7 or Python 3.2 or 3.3
plus the following libraries:

  - NumPy 1.6 or later (http://numpy.scipy.org/)
  - HDF5 1.8.7 or later (http://www.hdfgroup.org/HDF5/)
  - h5py 2.1 or later (http://www.h5py.org/)
  - ImmutablePy 0.1 (http://bitbucket.org/khinsen/immutablepy/)

Runnning the tests:

   python setup.py test

Installation:

   python setup.py install

