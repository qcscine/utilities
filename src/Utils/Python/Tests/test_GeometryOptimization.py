__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import scine_utilities as su
import os
import pytest

WATER_XYZ = """3

O     -0.0066000000   -0.4027000000    0.0000000000
H     -0.8046000000    0.2038000000    0.0000000000
H      0.8112000000    0.1989000000    0.0000000000
"""


def test_geometry_optimization():
    # Set-up optimization
    manager = su.core.ModuleManager.get_instance()
    calc = manager.get("calculator", "test")
    fname = "su_test_water.xyz"

    with open(fname, "w") as water_file:
        water_file.write(WATER_XYZ)

    ac, _ = su.io.read(fname)
    os.unlink(fname)
    calc.structure = ac

    logger = su.core.Log()
    # Check that geometry optimize call does not raise if no observer is provided
    su.geometry_optimize(calc, logger)
    # Check that geometry optimize call performs sanity check correctly
    calc.settings.update(su.ValueCollection({'self_consistence_criterion': 1.0e-6}))
    with pytest.raises(RuntimeError, match="The given calculator settings are too inaccurate"):
        su.geometry_optimize(calc, logger)
