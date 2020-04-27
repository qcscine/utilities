import pytest
import numpy
import scine_utilities as scine


def test_Constants():
    # test a few symbolically
    assert scine.ELEMENTARY_CHARGE == 1.6021766208e-19
    assert scine.AVOGADRO_NUMBER == 6.022140857e23
