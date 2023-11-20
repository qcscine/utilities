__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import numpy
import scine_utilities as scine


def test_Constants():
    # test a few symbolically
    assert scine.ELEMENTARY_CHARGE == 1.6021766208e-19
    assert scine.AVOGADRO_NUMBER == 6.022140857e23
