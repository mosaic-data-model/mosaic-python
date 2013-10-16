.. Written by Konrad Hinsen
.. License: CC-BY 3.0

.. index::
   single: XML

Mosaic XML files
################

The Mosaic XML format is defined by a `Relax NG <http://relaxng.org/>`_
schema. The following explanations document additional constraints
that cannot be encoded in a schema.

Mosaic data items in an XML file
--------------------------------

A Mosaic XML file contains a top-level element with tag ``mosaic``
whose children are Mosaic data items. Each of these data items has a
required attribute ``id`` that gives it a name through which it can be
identified uniquely inside the ``mosaic`` element. No two data items
in the same top-level element may have the same ``id`` value.

References to a data item take the form of an empty element with the
same tag as used for the definition of the data item itself. The
empty element contains a single attribute ``ref`` whose value is the
unique identifier of the data item that is referenced.


Floating-point numbers
----------------------

:ref:`Configurations<mosaic-configuration>` and real-valued
:ref:`properties<mosaic-property>` contain floating-point data. The
Mosaic specification allows two floating-point data types, "float32"
and "float64", which are based on IEEE binary floating point
formats. In XML, numbers are usually stored in decimal form.  Since an
exact conversion between binary and decimal floating-point numbers is
not possible, a compromise must be chosen. Since a main reason for
using XML is its human-readable layout, Mosaic's XML format uses a
standard XML decimal representation, at the cost of non-exact
conversion from and to binary data layouts.

Programs writing Mosaic XML files should convert floating point values
using the largest possible number of digits and then remove trailing
zeros from the mantissa. They should not attempt any rounding. The
symbols "NaN", "+inf" and "-inf" should be used for non-numbers, although
such values rarely make sense in the context of molecular simulations.


Array data
----------

:ref:`Property<mosaic-property>` data items contain one array value
per atom or site. The N-dimensional array values in a property are
stored as a single N+1-dimensional array, whose leading dimension is
the number of atoms or sites. Likewise, the positions in a
:ref:`configuration<mosaic-configuration>` data item are stored as an
array of shape (M, 3), where M is the number of sites.

The elements of an array are stored as a single list in row-major
order.
