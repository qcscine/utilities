__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import numpy
import scine_utilities as scine


def test_ElementInfo():
    # test a few instances
    assert scine.ElementInfo.element_from_symbol("F") == scine.ElementType.F
    assert scine.ElementInfo.symbol(scine.ElementType.Ca) == "Ca"
    assert scine.ElementInfo.mass(scine.ElementType.Ru) == 101.07

    assert scine.ElementInfo.Z(scine.ElementType.F) == 9

    assert scine.ElementInfo.all_implemented_elements()
    assert [str(e) for e in scine.ElementInfo.all_implemented_elements()]
