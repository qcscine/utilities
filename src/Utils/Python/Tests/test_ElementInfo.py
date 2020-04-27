import pytest
import numpy
import scine_utilities as scine


def test_ElementInfo():
    # test a few instances
    assert scine.ElementInfo.element_from_symbol("F") == scine.ElementType.F
    assert scine.ElementInfo.symbol(scine.ElementType.Ca) == "Ca"
    assert scine.ElementInfo.mass(scine.ElementType.Ru) == 101.07
