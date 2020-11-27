__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as scine


def test_ElementTypes():
    # Test the string magic function
    assert str(scine.ElementType.K) == "K"
    assert str(scine.ElementType.Cl) == "Cl"
    # Test the arithmetic properties of the enum
    assert scine.ElementType.Ar == scine.ElementType.Ar
    assert scine.ElementType.Ar != scine.ElementType.Ti
