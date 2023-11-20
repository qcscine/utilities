__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import scine_utilities as su
import numpy as np
import pytest

def test_geometry():
    structure = su.AtomCollection(3)
    structure.elements = [su.ElementType.O, su.ElementType.H, su.ElementType.H]                                                               
    structure.positions = [[-0.0066000000, -0.4027000000, 0.0000000000], 
                           [-0.8046000000, 0.2038000000, 0.0000000000],       
                           [0.8112000000, 0.1989000000, 0.0000000000]]
    structure.positions = structure.positions * su.BOHR_PER_ANGSTROM

    mode_disp = np.random.rand(3,3)
    mode = su.normal_modes.mode(0,mode_disp)
    pos_disp = su.geometry.displace_along_modes(structure.positions, [mode], [0.02])
    pos_redisp = su.geometry.displace_along_modes(pos_disp, [mode], [-0.02])
    assert np.all(np.abs(structure.positions[:,:] - pos_redisp[:,:]) < 1e-12)
