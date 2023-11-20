__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest

import scine_utilities as scine
import numpy


def test_size():
    elements = [scine.ElementType.H, scine.ElementType.F]
    pos_a = numpy.array([[0.0, 1.0, 2.0], [0.5, 1.5, 2.5]])
    pos_b = numpy.array([[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]])
    traj = scine.MolecularTrajectory()
    traj.elements = elements
    assert traj.empty() 
    assert traj.size() == 0
    traj.push_back(pos_a)
    traj.push_back(pos_b)
    assert not traj.empty() 
    assert traj.size() == 2


def test_access():
    elements = [scine.ElementType.H, scine.ElementType.F]
    pos_a = numpy.array([[0.0, 1.0, 2.0], [0.5, 1.5, 2.5]])
    pos_b = numpy.array([[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]])
    traj = scine.MolecularTrajectory(elements)
    assert traj.elements[0] == scine.ElementType.H
    assert traj.elements[1] == scine.ElementType.F
    traj.push_back(pos_a)
    traj.push_back(pos_b)
    traj.push_back(pos_b)
    traj.push_back(pos_b)
    traj.push_back(pos_a)
    assert traj.size() == 5
    assert numpy.array_equal(traj[0], pos_a)
    assert numpy.array_equal(traj[1], pos_b)
    assert numpy.array_equal(traj[2], pos_b)
    assert numpy.array_equal(traj[3], pos_b)
    assert numpy.array_equal(traj[4], pos_a)
    assert numpy.array_equal(traj[-1], pos_a)
    assert numpy.array_equal(traj[-2], pos_b)
    assert numpy.array_equal(traj[-3], pos_b)
    assert numpy.array_equal(traj[-4], pos_b)
    assert numpy.array_equal(traj[-5], pos_a)
    assert numpy.array_equal(traj[:3][0], pos_a)
    assert numpy.array_equal(traj[1:3][0], pos_b)
    with pytest.raises(IndexError) as excinfo:
        _ = traj[10]
    assert "out of range" in str(excinfo.value)
    with pytest.raises(IndexError) as excinfo:
        _ = traj[-10]
    assert "out of range" in str(excinfo.value)


def test_min_addition():
    elements = [scine.ElementType.H, scine.ElementType.F]
    pos_a = numpy.array([[0.0, 1.0, 2.0], [0.5, 1.5, 2.5]])
    pos_b = numpy.array([[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]])
    pos_c = numpy.array([[10.0, 2.0, 3.0], [1.5, 2.5, 3.5]])
    traj = scine.MolecularTrajectory(elements, 1.5)
    assert traj.elements[0] == scine.ElementType.H
    assert traj.elements[1] == scine.ElementType.F
    traj.push_back(pos_a)
    traj.push_back(pos_b)
    traj.push_back(pos_b)
    traj.push_back(pos_b)
    traj.push_back(pos_a)
    traj.push_back(pos_c)
    assert traj.size() == 4
    stricter_traj = scine.MolecularTrajectory(5.0)
    stricter_traj.elements = elements
    assert len(stricter_traj.elements) == 2
    stricter_traj.push_back(pos_a)
    stricter_traj.push_back(pos_b)
    stricter_traj.push_back(pos_b)
    stricter_traj.push_back(pos_b)
    stricter_traj.push_back(pos_a)
    assert stricter_traj.size() == 1
    stricter_traj.push_back(pos_c)
    assert stricter_traj.size() == 2
    stricter_traj.push_back(pos_a)
    assert stricter_traj.size() == 3
    stricter_traj.push_back(pos_a)
    assert stricter_traj.size() == 3

