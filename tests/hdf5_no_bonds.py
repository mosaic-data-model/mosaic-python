import mosaic.immutable_model
from mosaic.hdf5 import HDF5Store

M = mosaic.immutable_model
C = M.element('C')
atom = M.fragment('carbon', (), (('C', M.atom(C)),), ())
cube_universe = M.universe('cube', ((atom, 'carbons', 20),),
                           convention='my_own')

with HDF5Store("test.h5", "w") as store:
    store.store("universe", cube_universe)

with HDF5Store("test.h5") as store:
    r_universe = store.retrieve("universe")
