# -*- coding: utf-8 -*-

"""PDB structure import

This module converts structure files from the PDB in mmCIF format to
Mosaic data. The information about the chemical structure of the
residues is taken from the PDB Chemical Component Dictionary.

The code can be considered alpha quality, with standard proteins and
nucleic acid chains being converted correctly. The main difficulty
with handling the unusual cases is ensuring the correct interpretation
of the PDB data.

.. moduleauthor:: Konrad Hinsen <konrad.hinsen@cnrs-orleans.fr>

"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

from collections import OrderedDict
import copy
import itertools as it
import os
import sys
import traceback

import numpy as np

from mosaic_pdb.mmcif import MMCIFParser

from mosaic_pdb.space_groups import space_groups

from mosaic_pdb.pdb_chem_comp import components, variants
from mosaic.mutable_model import \
        Element, Atom, Fragment, Polymer, Universe, \
        Configuration, SiteProperty
from mosaic import api

from mosaic.utility import isstring

# Unit conversion factors

class Units:
    Ang = 0.1
    deg = np.pi / 180.


# Crystallographic unit cell

class UnitCell(object):

    def __init__(self, a, b, c, alpha, beta, gamma):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        e1 = np.array([self.a, 0, 0])
        e2 = self.b * np.array([np.cos(self.gamma), np.sin(self.gamma), 0.])
        e3_x = np.cos(self.beta)
        e3_y = (np.cos(self.alpha) - np.cos(self.beta) * np.cos(self.gamma)) \
               / np.sin(self.gamma)
        e3_z = np.sqrt(1. - e3_x ** 2 - e3_y ** 2)
        e3 = self.c * np.array([e3_x, e3_y, e3_z])
        self.basis = (e1, e2, e3)


# An MMCIFStructure accumulates all the information extracted
# from an mmCIF file.

class MMCIFStructure(object):

    def __init__(self, structure_file=None, pdb_code=None,
                 with_auth_spec=False):
        """
        Specify the data to be loaded. The following combinations
        are valid:

         - pdb_code only: the data is taken from a public or local
           PDB repository

         - structure_file only: the data is taken from the structure
           file

        @param structure_file: the name of the structure mmCIF file
        @type structure_file: C{str}
        @param pdb_code: a four-letter PDB code
        @type pdb_code: C{str}
        @param with_auth_spec: flag for reading and storing the
                               author's alternative atom site
                               specifications
        @type with_auth_spec: C{bool}
        """
        if pdb_code is not None:
            self.pdb_code = pdb_code
            assert structure_file is None
            self.structure_file = None
            parser = MMCIFParser(pdb_code=self.pdb_code)
        elif structure_file is not None:
            self.pdb_code = None
            if isstring(structure_file):
                structure_file = open(structure_file)
            self.structure_file = structure_file
            parser = MMCIFParser(file_object=self.structure_file)
        else:
            raise MMCIFError("No input structure given")

        self.experiment = {}
        self.cell_parameters = {}
        self.symmetry = {}
        self.entities = OrderedDict()
        self.components = OrderedDict()
        self.molecules = OrderedDict()
        self.ss_bridges = []
        self.extra_bonds = []
        self.models = []
        self.u = {}
        if with_auth_spec:
            self.auth_spec = {}
        else:
            self.auth_spec = None
        self.modified = False
        parser.parseToObjects(cell=self.cell_parameters,
                              symmetry=self.symmetry,
                              exptl=self.experiment,
                              entity=self.addEntity,
                              entity_name_com=self.addEntityName,
                              entity_poly=self.addPolymerInfo,
                              entity_poly_seq=self.buildSequence,
                              chem_comp=self.addComponent,
                              struct_conn=self.addConnection,
                              struct_asym=self.addMolecule,
                              atom_site=self.addSite,
                              atom_sites_alt=self.recordAltLabel,
                              atom_site_anisotrop=self.addAnisoU)
        if self.structure_file is not None:
            self.structure_file.close()
        if 'entry_id' in self.experiment:
            if self.pdb_code is None:
                self.pdb_code = self.experiment['entry_id']
            else:
                assert self.pdb_code == self.experiment['entry_id']
        if len(self.cell_parameters) == 0 \
               or (self.cell_parameters['length_a'] ==
                   self.cell_parameters['length_b'] ==
                   self.cell_parameters['length_a'] == '1'
                   and self.cell_parameters['angle_alpha'] ==
                       self.cell_parameters['angle_beta'] ==
                       self.cell_parameters['angle_gamma'] == '90') \
               or (self.cell_parameters['length_a'] ==
                   self.cell_parameters['length_b'] ==
                   self.cell_parameters['length_a'] == '?'
                   and self.cell_parameters['angle_alpha'] ==
                       self.cell_parameters['angle_beta'] ==
                       self.cell_parameters['angle_gamma'] == '?'):
            # NMR or other non-crystal structure
            self.cell = None
        else:
            self.cell = \
                UnitCell(
                    float(self.cell_parameters['length_a']) * Units.Ang,
                    float(self.cell_parameters['length_b']) * Units.Ang,
                    float(self.cell_parameters['length_c']) * Units.Ang,
                    float(self.cell_parameters['angle_alpha']) * Units.deg,
                    float(self.cell_parameters['angle_beta']) * Units.deg,
                    float(self.cell_parameters['angle_gamma']) * Units.deg)
        del self.cell_parameters
        self.space_group = None
        sg = self.symmetry.get('Int_Tables_number', None)
        if sg not in [None, '?', '.']:
            self.space_group = space_groups[int(sg)]
        else:
            for field in ['pdbx_full_space_group_name_H-M',
                          'space_group_name_H-M']:
                sg = self.symmetry.get(field, None)
                if sg not in [None, '?', '.']:
                    self.space_group = space_groups[sg]
                    break

    def getField(self, label, indices, data):
        if label not in indices:
            return None
        value = data[indices[label]]
        if value == '?':
            return None
        else:
            return value

    def getYesNoField(self, label, indices, data):
        value = self.getField(label, indices, data)
        if value in ['yes', 'y']:
            return True
        elif value in ['no', 'n']:
            return False
        assert value in [None, '.']
        return None

    def addEntity(self, indices, data):
        data = dict((label, data[index])
                    for label, index in indices.items()
                    if data[index] != '?')
        self.entities[data['id']] = data

    def addEntityName(self, indices, data):
        id = data[indices['entity_id']]
        name = data[indices['name']]
        if name != '?':
            self.entities[id]['common_name'] = name

    def addPolymerInfo(self, indices, data):
        entity_id = data[indices['entity_id']]
        polymer_type = self.getField('type', indices, data)
        if polymer_type is not None:
            self.entities[entity_id]['polymer_type'] = polymer_type
            self.entities[entity_id]['polymer_nstd_linkage'] = \
                         self.getYesNoField('nstd_linkage', indices, data)

    def buildSequence(self, indices, data):
        entity_id = data[indices['entity_id']]
        seq = self.entities[entity_id].get('sequence', [])
        self.entities[entity_id]['sequence'] = seq
        res_number = int(data[indices['num']])
        res_type = data[indices['mon_id']]
        if self.getYesNoField('hetero', indices, data):
            if seq and res_number == seq[-1][0]:
                # Consecutive monomer at the same position as the
                # preceding one.
                seq[-1] = (res_number, seq[-1][1] + (res_type,))
            else:
                # First monomer at this position.
                seq.append((res_number, (res_type,)))
        else:
            if seq and res_number - seq[-1][0] != 1:
                raise ValueError("non-consecutive residue numbers %d-%d"
                                 "in entity_poly_seq" % (res_number, seq[-1][0]))
            seq.append((res_number, res_type))

    def addMolecule(self, indices, data):
        data = dict((label, data[index])
                    for label, index in indices.items()
                    if data[index] != '?')
        self.molecules[data['id']] = data

    def addComponent(self, indices, data):
        data = dict((label, data[index])
                    for label, index in indices.items()
                    if data[index] != '?')
        self.components[data['id']] = data

    def addConnection(self, indices, data):
        if self.getField('conn_type_id', indices, data) == 'disulf':
            asym1, comp1, seq1, asym2, comp2, seq2 = \
                   [self.getField(f, indices, data)
                    for f in ['ptnr1_label_asym_id', 'ptnr1_label_comp_id',
                              'ptnr1_label_seq_id', 'ptnr2_label_asym_id',
                              'ptnr2_label_comp_id', 'ptnr2_label_seq_id']]
            self.ss_bridges.append(((asym1, comp1, int(seq1)),
                                    (asym2, comp2, int(seq2))))
        elif self.getField('conn_type_id', indices, data) == 'covale':
            asym1, comp1, seq1, atom1, asym2, comp2, seq2, atom2 = \
                   [self.getField(f, indices, data)
                    for f in ['ptnr1_label_asym_id', 'ptnr1_label_comp_id',
                              'ptnr1_label_seq_id', 'ptnr1_label_atom_id',
                              'ptnr2_label_asym_id', 'ptnr2_label_comp_id',
                              'ptnr2_label_seq_id', 'ptnr2_label_atom_id']]
            self.extra_bonds.append(((asym1, comp1, seq1, atom1),
                                     (asym2, comp2, seq2, atom2)))

    def recordAltLabel(self, indices, data):
        data = dict((label, data[index])
                    for label, index in indices.items()
                    if data[index] != '?')
        # this information is in principle required for any
        # entry with alternate locations, but in practice it is
        # very rare. This method exists solely to spot such
        # information in order to see if it can be exploited.
        raise ValueError("encountered atom_sites_alt")

    def addSite(self, indices, data):
        asym_id = data[indices['label_asym_id']]
        unique_id = data[indices['id']]
        atom_id = data[indices['label_atom_id']]
        alt_id = data[indices['label_alt_id']]
        comp_id = data[indices['label_comp_id']]

        entity = self.entities[self.molecules[asym_id]['entity_id']]
        is_polymer = entity['type']  == 'polymer' \
                     and 'polymer_type' in entity
        if is_polymer:
            seq_id_string = data[indices['label_seq_id']] 
            try:
                seq_id = int(data[indices['label_seq_id']])
            except ValueError:
                raise ValueError('non-numeric seq_id %s for atom id %s'
                                 % (data[indices['label_seq_id']], unique_id))
            if seq_id <= 0:
                raise ValueError('invalid seq_id %d for atom id %s'
                                 % (seq_id, unique_id))
            ins_code = self.getField('pdbx_PDB_ins_code', indices, data)
        else:
            seq_id = None
            ins_code = None

        try:
            model_number = int(data[indices['pdbx_PDB_model_num']])
        except (KeyError, ValueError):
            model_number = 1
        while len(self.models) < model_number:
            self.models.append(OrderedDict())
        model = self.models[model_number - 1]
        if is_polymer:
            molecule = model.get(asym_id, [[]])
            monomer = molecule[-1]
            if len(monomer) > 0 and \
                   (monomer[0][4] != seq_id or monomer[0][5] != ins_code):
                monomer = []
                molecule.append(monomer)
        else:
            molecule = model.get(asym_id, [])
            monomer = molecule
        model[asym_id] = molecule

        atom_type = data[indices['type_symbol']]
        x = float(data[indices['Cartn_x']]) * Units.Ang
        y = float(data[indices['Cartn_y']]) * Units.Ang
        z = float(data[indices['Cartn_z']]) * Units.Ang
        occupancy = float(data[indices['occupancy']])
        b = self.getField('B_iso_or_equiv', indices, data)
        if b is None:
            u_iso = None
        else:
            u_iso = float(b) * Units.Ang ** 2 / (8. * np.pi ** 2)
        monomer.append((unique_id, atom_id, alt_id, comp_id, seq_id, ins_code,
                        atom_type, x, y, z, occupancy, u_iso))

        if self.auth_spec is not None:
            auth_asym_id = self.getField('auth_asym_id', indices, data)
            auth_atom_id = self.getField('auth_atom_id', indices, data)
            auth_comp_id = self.getField('auth_comp_id', indices, data)
            auth_seq_id = self.getField('auth_seq_id', indices, data)
            self.auth_spec[unique_id] =  dict(asym_id = auth_asym_id,
                                              atom_id = auth_atom_id,
                                              comp_id = auth_comp_id,
                                              seq_id  = int(auth_seq_id))

        stop = False
        for field in ['adp_type', 'aniso_B[1][1]', 'aniso_U[1][1]']:
            value = self.getField(field, indices, data)
            if value is not None:
                sys.stderr.write("%s = %s\n" % (field, str(value)))
                stop = True
        if stop:
            raise ValueError("found aniso data in atom_site")

    def addAnisoU(self, indices, data):
        unique_id = data[indices['id']]
        factor = 1.
        letter = 'U'
        u11 = self.getField('%s[1][1]' % letter, indices, data)
        if u11 is None:
            factor = 1. / (8. * np.pi ** 2)
            letter = 'B'
            if self.getField('%s[1][1]' % letter, indices, data) is None:
                raise ValueError("atom_site_anisotrop with neither U nor B")
        u = np.array([self.getField(letter + ij, indices, data)
                         for ij in ['[1][1]', '[2][2]', '[3][3]',
                                    '[2][3]', '[1][3]', '[1][2]']],
                        dtype=np.float32) * factor * Units.Ang ** 2
        self.u[unique_id] = u


# Make a fragment for a PDB residue/monomer, using the
# chemical component dictionary.

def partition_sites(sites):
    group = []
    last_group_id = None
    for s in sites:
        unique_id, atom_id, alt_id, comp_id, seq_id, ins_code, \
            atom_type, x, y, z, occupancy, u_iso = s
        group_id = (comp_id, seq_id)
        if last_group_id is not None and last_group_id != group_id:
            yield comp_id, group
            group = [s]
        else:
            group.append(s)
        last_group_id = group_id
    if group:
        yield comp_id, group

def make_monomer(sites, site_properties):
    comp_ids = [s[3] for s in sites]
    if len(set(comp_ids)) > 1:
        # Make a heterogeneous monomer
        ms = [_make_monomer(comp_id, list(comp_sites), site_properties)
              for comp_id, comp_sites in it.groupby(sites, lambda s: s[3])]
        fs = [f for f, u in ms]
        us = [u for f, u in ms]
        species = '/'.join(f.label.split('_')[0] for f in fs)
        seq_id = fs[0].label.split('_')[1]
        return Fragment(species+'_'+seq_id, species, fs, (), ()), sum(us, ())
    else:
        # Make a standard monomer
        return _make_monomer(comp_ids[0], sites, site_properties)

def _make_monomer(comp_id, sites, site_properties):
    atom_ids = set(s[1] for s in sites)
    comp = components[comp_id]
    comp_atom_ids = set(comp.atoms['atom_id'])
    if hasattr(comp, 'bonds'):
        bond_orders = {'sing': 'single',
                       'doub': 'double',
                       'trip': 'triple',
                       'quad': 'qudruple',
                       'arom': 'aromatic'}
        bonds = tuple((x['atom_id_1'], x['atom_id_2'],
                       bond_orders.get(x['value_order'].lower(), ''))
                      for x in comp.bonds
                      if x['atom_id_1'] in atom_ids
                         and x['atom_id_2'] in atom_ids)
    else:
        bonds = ()
    if not comp_atom_ids.issuperset(atom_ids):
        if atom_ids.difference(comp_atom_ids) == set(["HO5'"]) \
           and "O5'" in atom_ids:
            # Provide a quick fix for a common case in the PDB, to avoid
            # having lots of DNA structures fail for a trivial reason
            comp_atom_ids.add("HO5'")
            bonds = bonds + (("HO5'", "O5'", 'single'),)
        else:
            # Variants exist only for amino acids. Find those that match
            # the atom_ids in the current monomer.
            vs = [v for v in variants.get(comp_id, [])
                  if set(components[v].atoms['atom_id']).issuperset(atom_ids)]
            if not vs:
                missing = atom_ids.difference(comp_atom_ids)
                raise ValueError("Atom(s) %s not in component %s"
                                 % (', '.join(missing), comp_id))
            variant = None
            if len(vs) == 1:
                variant = vs[0]
            else:
                diff = [len(set(components[v].atoms['atom_id'])
                            .difference(atom_ids))
                        for v in vs]
                if sum(d == 0 for d in diff) == 1:
                    variant = vs[diff.index(0)]
                else:
                    seq_id = sites[0][4]
                    if seq_id == 1:
                        # N-terminal variants have LSN3 in their name
                        vs = [v for v in vs if 'LSN3' in v]
                        if len(vs) == 1:
                            variant = vs[0]
                        else:
                            diff = [len(set(components[v].atoms['atom_id'])
                                        .difference(atom_ids))
                                    for v in vs]
                            if sum(d == 0 for d in diff) == 1:
                                variant = vs[diff.index(0)]
            if variant is None:
                raise ValueError("Monomer variant not unique: %s" % str(vs))
            comp = components[variant]
            comp_atom_ids = set(comp.atoms['atom_id'])
    atom_sites = {}
    atom_ids = []
    atom_types = []
    for unique_id, atom_id, alt_id, comp_id, seq_id, ins_code, \
            atom_type, x, y, z, occupancy, u_iso in sites:
        if atom_id not in comp_atom_ids:
            raise ValueError("No atom named %s in component %s"
                             % (atom_id, comp.id))
        if unique_id in site_properties:
            raise ValueError("Atom id %s not unique" % unique_id)
        if atom_id not in atom_ids:
            atom_ids.append(atom_id)
            atom_types.append(atom_type)
        atom_sites[atom_id] = atom_sites.get(atom_id, [])
        atom_sites[atom_id].append(unique_id)
        site_properties[unique_id] = (x, y, z, occupancy, u_iso)
    atoms = tuple(Atom(atom_id, Element(atom_type.capitalize()),
                       len(atom_sites[atom_id]))
                  for atom_id, atom_type in zip(atom_ids, atom_types))
    unique_ids = tuple(sum((atom_sites[atom_id] for atom_id in atom_ids), []))
    if seq_id is None:
        label = comp.three_letter_code
    else:
        label = comp.three_letter_code + '_' + str(seq_id)
        #if ins_code is not None:
        #    label = label + '_' + ins_code
    return Fragment(label, comp.id, (), atoms, bonds), unique_ids


def ss_bridges(structure, asym_id_1, asym_id_2):
    """
    Iterate over the disulfide bridges in the structure.
    """
    for (asym1, comp1, seq1), (asym2, comp2, seq2) \
            in structure.ss_bridges:
        if asym1 == asym_id_1 and asym2 == asym_id_2:
            yield seq1, seq2
        elif asym1 == asym_id_2 and asym2 == asym_id_1:
            yield seq2, seq1


def distinct_names(fragments):
    """
    Make sure that all the fragments have distinct names,
    adding distinctive integer suffixes if necessary.
    """
    name_counts = {}
    for f in fragments:
        name_counts[f.label] = name_counts.get(f.label, 0) + 1
    new_names = {}
    for name, count in name_counts.items():
        if count == 1:
            new_names[name] = [name]
        else:
            new_names[name] = [name + '_' + str(n)
                               for n in range(1, count + 1)]
    for f in fragments:
        label = new_names[f.label][0]
        del new_names[f.label][0]
        f.label = label


class MoleculeList(list):

    """
    Turn a sequence of molecules into a compressed molecule list
    with repetition counts.
    """

    def __init__(self, iterable):
        list.__init__(self, [])
        for item in iterable:
            self.append(item)

    def append(self, item):
        assert isinstance(item, Fragment)
        if len(self) == 0:
            list.append(self, (item, 1))
        else:
            last, count = self[-1]
            if last.eq_structure(item):
                if last.label == item.label:
                    f = last
                else:
                    f = copy.copy(last)
                    f.label = last.label + ',' + item.label
                self[-1] = (f, count + 1)
            else:
                list.append(self, (item, 1))


def symmetry_transformations(space_group):
    """
    Convert the symmetry transformations for the crystal's space
    group to the right format.
    """
    transformations = []
    assert (space_group.transformations[0][0] ==
            np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])).all()
    assert (space_group.transformations[0][1] == np.array([0, 0, 0])).all()
    assert (space_group.transformations[0][2] == np.array([1, 1, 1])).all()
    for rot, tr_num, tr_den in space_group.transformations[1:]:
        rot = np.array(rot, dtype=np.float64)
        tr = np.array(tr_num, dtype=np.float64) / \
             np.array(tr_den, dtype=np.float64)
        transformations.append((rot, tr))
    return transformations

def chain_links(residues, atom1, atom2, bond_order):
    def link_atoms(residue, atom):
        if residue.atoms:
            # standard monomer
            atom_ref = residue.get(atom)
            if atom_ref is None:
                return []
            else:
                return [atom_ref]
        else:
            # heterogeneous monomer
            fs = residue.fragments
            return [f[atom] for f in fs if atom in f]
    bonds = []
    for r1, r2 in zip(residues[:-1], residues[1:]):
        for a1 in link_atoms(r1, atom1):
            for a2 in link_atoms(r2, atom2):
                bonds.append((a1, a2, bond_order))
    return tuple(bonds)

def make_crystal(structure, model_number, universes):
    """
    Make a Mosaic representation of a single model for a crystal.
    Add the universe to the passed-in universe list if it's not
    already in it. This takes care of the common situation that
    multiple models share a common universe.
    """
    model = structure.models[model_number]
    site_properties = {}
    monomer_fragments = {}
    molecules = {}
    unique_ids = {}
    for asym_id, sites in model.items():
        unique_ids[asym_id] = []
        entity_id = structure.molecules[asym_id]['entity_id']
        entity = structure.entities[entity_id]
        entity_type = entity['type']
        if entity_type == 'polymer' and 'polymer_type' in entity:
            # Sugars are classified as polymers in the PDB, but have
            # no polymer type and a single residue with the
            # non-numeric seq_id '.'. This looks like a bug in the PDB
            # data, but since it's so frequent, we have to deal with
            # it. We do so by treating polymers without a polymer_type
            # as non-polymers.
            polymer_type = entity['polymer_type']
            if entity['polymer_nstd_linkage']:
                raise ValueError("Non-standard linkage not yet implemented")
            if 'sequence' not in entity:
                raise ValueError('Polymer entity has no sequence information')
            name = entity.get('pdbx_description', None)
            if name is None:
                name = entity.get('common_name', None)
            if name is None:
                name = "Entity_" + entity_id
            name = name.replace(' ', '_')
            name = name.replace('.', '_')
            # Start with empty residues
            residues = OrderedDict()
            for seq_id, comp_id in entity['sequence']:
                if isinstance(comp_id, tuple):
                    # heterogeneous monomer
                    species = '/'.join(comp_id)
                else:
                    species = comp_id
                f = Fragment(species + '_' + str(seq_id), species,
                             (), (), ())
                residues[seq_id] = f
                monomer_fragments[(asym_id, seq_id)] = f
            # For residues that have sites, make non-empty fragments
            for monomer_sites in sites:
                assert structure.modified or len(monomer_sites) > 0
                if not monomer_sites:
                    continue
                seq_id = monomer_sites[0][4]
                f, u = make_monomer(monomer_sites, site_properties)
                species = monomer_fragments[(asym_id, seq_id)].species
                assert f.species.startswith(species)
                monomer_fragments[(asym_id, seq_id)] = f
                if seq_id not in residues:
                    raise ValueError('seq_id %d in site list but not '
                                     'in sequence' % seq_id)
                residues[seq_id] = f
                unique_ids[asym_id].extend(u)
            residues = list(residues.values())
            # Add inter-monomer links
            if polymer_type in ['polypeptide(L)']:
                peptide_bonds = chain_links(residues, 'C', 'N', 'single')
                ss_bonds = tuple((monomer_fragments[(asym_id, seq1)].get('SG'),
                                  monomer_fragments[(asym_id, seq2)].get('SG'),
                                  'single')
                                 for seq1, seq2 in
                                     ss_bridges(structure, asym_id, asym_id))
                bonds = peptide_bonds + ss_bonds
                polymer_type = 'polypeptide'
            elif polymer_type in ['polyribonucleotide',
                                  'polydeoxyribonucleotide',
                                  'polydeoxyribonucleotide/'
                                      'polyribonucleotide hybrid']:
                bonds = chain_links(residues, 'O3\'', 'P', 'single')
            elif polymer_type in ['polysaccharide(D)', 'polysaccharide(L)']:
                # need to figure out how monosaccharides are linked
                # example: 1AGA
                bonds = tuple()
            else:
                raise ValueError("Polymer type %s not yet implemented"
                                 % polymer_type)
            bonds = tuple((a1, a2, order)
                          for a1, a2, order in bonds
                          if a1 is not None and a2 is not None)
            polymer_types = {'polydeoxyribonucleotide/'
                                 'polyribonucleotide hybrid':
                             'polynucleotide'}
            # Cleanup name to make it a valid species name
            name = ''.join(c for c in name if c in api._allowed_in_labels)
            if len(name) > 32767:
                name = name[:32767]
            molecules[asym_id] = [Polymer(asym_id, name, residues, bonds,
                                          polymer_types.get(polymer_type,
                                                            polymer_type))]
        elif entity_type == 'water':
            ms = []
            for atom in sites:
                if atom[6] == 'O':
                    ms.append([atom])
                else:
                    assert atom[6] == 'H'
                    ms[-1].append(atom)
            molecules[asym_id] = []
            for m in ms:
                f, u = make_monomer(m, site_properties)
                molecules[asym_id].append(f)
                unique_ids[asym_id].extend(u)
        else:
            f, u = make_monomer(sites, site_properties)
            molecules[asym_id] = [f]
            unique_ids[asym_id].extend(u)

    # Treat inter-chain disulfide bridges
    order_kept = not structure.modified
    clusters = []
    mol_to_cluster = {}
    for a in model.keys():
        mol_to_cluster[a] = len(clusters)
        clusters.append((set([a]), []))
    for asym1, asym2 in it.combinations(model.keys(), 2):
        ss = list(ss_bridges(structure, asym1, asym2))
        if ss:
            c1 = mol_to_cluster[asym1]
            c2 = mol_to_cluster[asym2]
            if c1 != c2:
                clusters[c1][0].update(clusters[c2][0])
                clusters[c1][1].extend(clusters[c2][1])
                clusters[c2] = (None, None)
                mol_to_cluster[asym2] = c1
            clusters[c1][1].extend([(asym1, s1, asym2, s2) for s1, s2 in ss])
    for asym_ids, bridges in clusters:
        if asym_ids is not None and len(asym_ids) > 1:
            asym_ids = list(asym_ids)
            order_kept = False
            # This could be improved to maintain the atom order in
            # specific and frequent special cases. For later...
            fragments = sum((molecules[a] for a in asym_ids), [])
            distinct_names(fragments)
            uids = sum((unique_ids[a] for a in asym_ids), [])
            bonds = tuple((monomer_fragments[(asym1, seq1)].get('SG'),
                           monomer_fragments[(asym2, seq2)].get('SG'),
                           'single')
                          for asym1, seq1, asym2, seq2 in bridges)
            bonds = tuple((a1, a2, order)
                          for a1, a2, order in bonds
                          if a1 is not None and a2 is not None)
            mol = Fragment('+'.join(f.label for f in fragments),
                           '+'.join(f.species for f in fragments),
                           fragments, (), bonds)
            molecules[asym_ids[0]] = [mol]
            unique_ids[asym_ids[0]] = uids
            for a in asym_ids[1:]:
                del molecules[a]
                del unique_ids[a]

    if structure.cell is None:
        cell_shape = 'infinite'
        st = []
    else:
        pdb_cell_shape = structure.symmetry.get('cell_setting', None)
        cell_shape = {'cubic': 'cube',
                      'orthorhombic': 'cuboid',
                      'tetragonal': 'cuboid',
                      }.get(pdb_cell_shape, 'parallelepiped')
        if structure.space_group is None:
            st = []
        else:
            st = symmetry_transformations(structure.space_group)
    universe = Universe(cell_shape,
                        MoleculeList(m
                                     for m in it.chain.from_iterable(
                                         molecules.get(asym_id, [])
                                         for asym_id in model.keys())),
                        symmetry_transformations=st,
                        convention='PDB')
    assert universe.number_of_sites == len(site_properties)

    # Check if an equal universe already exists
    if universe in universes:
        universe = universes[universes.index(universe)]
    else:
        universes.append(universe)

    # Configuration and occupancy
    positions = np.zeros((universe.number_of_sites, 3), np.float32)
    occupancy = np.zeros((universe.number_of_sites,), np.float32)
    u_iso = np.zeros((universe.number_of_sites,), np.float32)
    for i, uid in enumerate(sum((unique_ids.get(asym_id, [])
                                 for asym_id in model.keys()), [])):
        # The conversion process should maintain the order of the sites
        # in the mmCIF file, as long as multiple sites corresponding to the
        # same atom (different alt_ids) occur consecutively in the list.
        # This is not explicitly promised by the documentation, but seems
        # to be the case in most if not all real-life mmCIF files.
        # Some specific features (e.g. disulfide bonds between chains)
        # can cause the site order to be changed.
        if order_kept and model_number == 0:
            assert i == int(uid) - 1
        positions[i, :] = site_properties[uid][0:3]
        occupancy[i] = site_properties[uid][3]
        if u_iso is not None and site_properties[uid][4] is not None:
            u_iso[i] = site_properties[uid][4]
    if universe.cell_shape == 'infinite':
        cell_parameters = None
    else:
        cell_parameters = np.zeros((3, 3), np.float32)
        for i, v in enumerate(structure.cell.basis):
            cell_parameters[i, :] = v
        if universe.cell_shape == 'cube':
            edge = cell_parameters[0, 0]
            assert np.minimum.reduce(
                        np.fabs((cell_parameters
                                    - edge * np.eye(3)).flat)) < 1.e-7
            cell_parameters = edge
        elif universe.cell_shape == 'cuboid':
            edge1 = cell_parameters[0, 0]
            edge2 = cell_parameters[1, 1]
            edge3 = cell_parameters[2, 2]
            assert np.minimum.reduce(
                np.fabs((cell_parameters -
                            np.array([[edge1, 0., 0.],
                                         [0., edge2, 0.],
                                         [0., 0., edge3]])).flat)) < 1.e-7
            cell_parameters = np.array([edge1, edge2, edge3],
                                          dtype=np.float32)

    data = OrderedDict()
    data['configuration'] = Configuration(universe, positions, cell_parameters)
    if not (occupancy == 1.).all():
        data['occupancy'] = SiteProperty(universe, "occupancy", "", occupancy)

    # Anisotropic or isotropic U
    if len(structure.models) == 1 and len(structure.u) > 0:
        u_aniso = np.zeros((universe.number_of_sites, 6), np.float32)
        for i, uid in enumerate(sum((unique_ids.get(asym_id, [])
                                     for asym_id in model.keys()), [])):
            try:
                u_aniso[i] = structure.u[uid]
            except KeyError:
                u_aniso[i, :3] = u_iso[i]
                u_aniso[i, 3:] = 0.
        data['anisotropic_displacement_parameter'] = \
            SiteProperty(universe, 'anisotropic_displacement_parameter',
                         'nm2', u_aniso)
    elif not (u_iso == 0.).all():
        data['isotropic_displacement_parameter'] = \
            SiteProperty(universe, 'isotropic_displacement_parameter',
                         'nm2', u_iso)

    return data


def make_models(structure_file=None, pdb_code=None):
    """
    Return an iterator over all the data items needed to
    represent the structure file's contents. The first
    items are the universes, then everything else.
    """
    s = MMCIFStructure(structure_file=structure_file,
                       pdb_code=pdb_code)
    universes = []
    models = [make_crystal(s, i, universes)
              for i in range(len(s.models))]

    if len(universes) == 1:
        yield "universe", universes[0]
    else:
        for i in range(len(universes)):
            yield "universe_%d" % (i + 1), universes[i]

    if len(models) == 1:
        for name, value in models[0].items():
            yield name, value
    else:
        for i, model in enumerate(models):
            for name, value in model.items():
                yield "%s_model_%d" % (name, i + 1), value
