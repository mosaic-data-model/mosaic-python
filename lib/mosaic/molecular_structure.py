# -*- coding: utf-8 -*-
"""Molecular structure utilities

.. moduleauthor:: Konrad Hinsen

Utilities for working with hierarchical molecular structures.
"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2014 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

class FragmentIterator(object):

    """Iterator over fragments that match some predicate

    Usage example for iterating over peptide chains:

        def peptide_predicate(fragment):
            return fragment.is_polymer and fragment.polymer_type == 'polypeptide'

       for chain, atom_index, site_index in FragmentIterator(universe, peptide_predicate):
           ...
    """

    class PseudoFragment(object):

        def __init__(self, universe, template=False):
            if template:
                self.fragments = [f for f, c in universe.molecules]
            else:
                self.fragments = sum([c*[f]
                                      for f, c in universe.molecules],
                                     [])
            self.atoms = ()

        is_polymer = False

    def __init__(self, universe, predicate,
                 descend_on_match=False, template=False):
        """
        :parameter universe: a universe
        :type universe: :class:`mosaic.api.MosaicUniverse`
        :parameter predicate: an object called with the fragment as its only
                              argument and expected to return True or False
        :type predicate: callable
        :parameter descend_on_match: if True, iterate over the subfragments
                                     of fragments that match the predicate
        :type descend_on_match: bool
        :parameter template: if True, iterate over template fragments,
                             otherwise iterate over molecules
        :type template: bool
        """
        self.predicate = predicate
        self.descend_on_match = descend_on_match
        self.stack = [FragmentIterator.PseudoFragment(universe, template)]
        self.atom_index = 0
        self.site_index = 0

    def __iter__(self):
        return self

    def __next__(self):
        while self.stack:
            item = self.stack.pop(0)
            if isinstance(item, tuple):
                a, s = item
                self.atom_index += a
                self.site_index += s
            else:
                match = None
                push = []
                if not isinstance(item, FragmentIterator.PseudoFragment) \
                   and self.predicate(item):
                    match = (item, self.atom_index, self.site_index)
                    if self.descend_on_match:
                        push.extend(item.fragments)
                    else:
                        for f in item.fragments:
                            push.append((f.number_of_atoms, f.number_of_sites))
                else:
                    push.extend(item.fragments)
                push.append((len(item.atoms),
                             sum(a.number_of_sites
                                 for a in item.atoms)))
                for new_item in reversed(push):
                    self.stack.insert(0, new_item)
                if match is not None:
                    return match
        raise StopIteration

    next = __next__
