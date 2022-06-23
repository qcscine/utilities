__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import pytest

import scine_utilities as scine
import numpy


def test_size_constructor():
    a = scine.AtomCollection(4)
    assert a.size() == 4


def test_element_access():
    a = scine.AtomCollection(2)
    a.elements = [scine.ElementType.H, scine.ElementType.F]
    a.set_element(0, scine.ElementType.Ca)
    assert a.get_element(0) == scine.ElementType.Ca
    assert a.get_element(1) == scine.ElementType.F


def test_position_access():
    pos_a = numpy.array([0.0, 1.0, 2.0])
    pos_b = numpy.array([1.0, 2.0, 3.0])
    a = scine.AtomCollection(2)
    a.positions = [pos_a, pos_b]
    assert numpy.array_equal(a.get_position(0), pos_a)
    assert numpy.array_equal(a.get_position(1), pos_b)
    assert numpy.array_equal(a[1].position, pos_b)
    assert numpy.array_equal(a[-1].position, pos_b)
    assert numpy.array_equal(a[0].position, pos_a)
    assert numpy.array_equal(a[-2].position, pos_a)


def test_vector_functions():
    a = scine.AtomCollection()
    assert a.size() == 0
    a.push_back(
        scine.Atom(scine.ElementType.P)
    )
    a.push_back(
        scine.Atom(scine.ElementType.N)
    )
    assert a.size() == 2
    a.clear()
    assert a.size() == 0
    a.resize(2)
    assert a.size() == 2


def test_sequence_functions():
    a = scine.AtomCollection(3)
    a.elements = [
        scine.ElementType.F,
        scine.ElementType.Ca,
        scine.ElementType.Br
    ]

    element_symbol_strings = ",".join([str(atom.element) for atom in a])
    assert element_symbol_strings == "F,Ca,Br"

def test_out_of_range():
    a = scine.AtomCollection(3)
    with pytest.raises(IndexError) as excinfo:
        _ = a[10]
    assert "out of range" in str(excinfo.value)
    with pytest.raises(IndexError) as excinfo:
        _ = a[-10]
    assert "out of range" in str(excinfo.value)
