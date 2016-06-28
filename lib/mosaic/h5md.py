# -*- coding: utf-8 -*-
"""H5MD trajectory I/O

.. moduleauthor:: Konrad Hinsen


"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2014 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import h5py

import numpy as np

import mosaic.api as api
import mosaic.array_model as AM
from mosaic.hdf5 import HDF5Store
from mosaic.version import version as mosaic_version
from mosaic.utility import isstring

class Trajectory(object):

    """H5MD trajectory
    """

    def __init__(self, hdf5_group_or_filename, universe=None,
                 author_name='unknown', author_email=None,
                 file_mode=None,
                 float_type=np.float32):
        """
        :param hdf5_group_or_filename: a group in an open HDF5 file,
                                       or a string interpreted as a file name
        :type hdf5_group_or_filename: str or ``h5py.Group``
        :param file_mode: ``'r'`` or ``'w'``, used only with a filename argument
        :type file_mode:  str
        :param universe: the universe associated with the trajectory
        :type universe:  :class:`mosaic.api.MosaicUniverse`
        """
        if isstring(hdf5_group_or_filename):
            if file_mode is None:
                file_mode = 'r'
            self.root = h5py.File(hdf5_group_or_filename, file_mode)
            self._close = True
        else:
            api.validate_value(file_mode, [None], "file_mode")
            self.root = hdf5_group_or_filename
            self._close = False
        self.float_type = float_type

        self._mosaic_group(universe)
        self._create_h5md_group(author_name, author_email)
        self._create_particles_group()

    def _mosaic_group(self, universe):
        if universe is None:
            self.mosaic = self.root['mosaic']
            self.store = HDF5Store(self.mosaic)
            self.universe = self.store.retrieve('universe')
        else:
            self.mosaic = self.root.create_group('mosaic')
            self.store = HDF5Store(self.mosaic)
            self.universe = universe
            self.store.store('universe', universe)

    def _create_h5md_group(self, author_name, author_email=None):
        if 'h5md' not in self.root:
            h5md = self.root.create_group('h5md')
            h5md.attrs['version'] = np.array([1, 0], np.int16)
            author = h5md.create_group('author')
            author.attrs['name'] = np.string_(author_name)
            if author_email is not None:
                author.attrs['email'] = np.string_(author_email)
            creator = h5md.create_group('creator')
            creator.attrs['name'] = np.string_('pyMosaic')
            creator.attrs['version'] = np.string_(mosaic_version)
            modules = h5md.create_group('modules')
            mmosaic = modules.create_group('mosaic')
            mmosaic.attrs['version'] = np.array([0, 1], np.int16)
            munits = modules.create_group('units')
            munits.attrs['version'] = np.array([1, 0], np.int16)
            munits.attrs['system'] = np.string_('SI')

    def _create_particles_group(self):
        if 'particles' not in self.root:

            particles = self.root.create_group('particles')
            universe = particles.create_group('universe')
            position = universe.create_group('position')
            step = position.create_dataset('step',
                                           shape=(0,),
                                           maxshape=(None,),
                                           dtype=np.int32)
            time = position.create_dataset('time',
                                           shape=(0,),
                                           maxshape=(None,),
                                           dtype=self.float_type)
            time.attrs['unit'] = 'ps'
            nsites = self.universe.number_of_sites
            position_values = position.create_dataset('value',
                                                      shape=(0,nsites,3),
                                                      maxshape=(None,nsites,3),
                                                      dtype=self.float_type)
            position_values.attrs['unit'] = 'nm'
            box = universe.create_group('box')
            box.attrs['dimension'] = np.int16(3)
            if self.universe.cell_shape == 'infinite':
                boundary = 'none'
                edge_values = None
            else:
                boundary = 'periodic'
                edges = box.create_group('edges')
                if self.universe.cell_shape == 'parallelepiped':
                    e_shape = (0, 3, 3)
                    e_maxshape = (None, 3, 3)
                else:
                    e_shape = (0, 3,)
                    e_maxshape = (None, 3)
                edge_values = edges.create_dataset('value',
                                                   shape=e_shape,
                                                   maxshape=e_maxshape,
                                                   dtype=self.float_type)
                edge_values.attrs['unit'] = 'nm'
                edges['step'] = step
                edges['time'] = time
            box.attrs['boundary'] = np.array(3*[np.string_(boundary)])

            self.time_dataset = time
            self.step_dataset = step
            self.position_dataset = position_values
            self.edge_dataset = edge_values

        else:

            universe = self.root['particles/universe']
            self.time_dataset = universe['position/time']
            self.step_dataset = universe['position/step']
            self.position_dataset = universe['position/value']
            if 'edges' in universe['box']:
                self.edge_dataset = universe['box/edges/value']
            else:
                self.edge_dataset = None

    def close(self):
        self.store.close()
        if self._close:
            self.root.file.close()
        self.root = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def __len__(self):
        return self.step_dataset.shape[0]

    def flush(self):
        self.root.file.flush()

    def write_step(self, step, time, configuration):
        assert configuration.universe is self.universe
        n = self.step_dataset.shape[0]
        self.position_dataset.resize(n+1, axis=0)
        self.position_dataset[n, :, :] = configuration.positions
        cp = configuration.cell_parameters
        if self.universe.cell_shape == 'cube':
            self.edge_dataset.resize(n+1, axis=0)
            self.edge_dataset[n, :] = np.array(3*[cp], self.float_type)
        elif self.universe.cell_shape in ['cuboid', 'parallelepiped']:
            self.edge_dataset.resize(n+1, axis=0)
            self.edge_dataset[n, :] = np.array(cp, self.float_type)
        self.time_dataset.resize(n+1, axis=0)
        self.time_dataset[n] = self.float_type(time)
        self.step_dataset.resize(n+1, axis=0)
        self.step_dataset[n] = np.int32(step)

    def read_step(self, index):
        if index < 0 or index >= self.step_dataset.shape[0]:
            raise IndexError
        pos = self.position_dataset[index, ...]
        if self.edge_dataset is None:
            cp = np.empty(self.universe.cell_parameter_array_shape,
                          pos.dtype)
        else:
            cp = self.edge_dataset[index, ...]
            if self.universe.cell_shape == 'cube':
                cp = cp[0]
        return self.step_dataset[index], \
               self.time_dataset[index], \
               AM.Configuration(self.universe, pos, cp)

    def __getitem__(self, index):
        return self.read_step(index)

    def read_time(self, first=0, last=None, skip=None):
        return self.time_dataset[first:last:skip]

    def read_atom(self, index, first=0, last=None, skip=None):
        if index < 0 or index >= self.position_dataset.shape[1]:
            raise IndexError
        xyz = self.position_dataset[first:last:skip, i, :]
        if self.edge_dataset is None:
            return xyz
        edges = self.edge_dataset[first:last:skip, :]
        # convert to box (crystallographic) coordinates
        b = self._cartesian_to_box(xyz, edges)
        # remove jumps
        d = b[1:]-b[:-1]
        d = np.fmod(d, 1.)
        d -= d >= 0.5
        d += d < -0.5
        b[1:] = d
        b = np.add.accumulate(b, axis=0)
        # back to Cartesian coordinates
        return self._box_to_cartesian(b, edges)
