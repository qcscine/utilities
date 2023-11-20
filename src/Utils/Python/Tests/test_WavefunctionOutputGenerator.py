__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as scine
import os
import os.path

def read_file(filename: str):
    with open(filename, "r") as f:
        text = ""
        for line in f:
            text = text + line
    return text

def test_WavefunctionOutputGeneratorGeneratesFile():
    manager = scine.core.ModuleManager.get_instance()
    calc = manager.get("calculator", "TEST")
    wf_gen = scine.core.to_wf_generator(calc)

    assert wf_gen is not None

    wf_gen.wavefunction2file("test.out")

    os.path.isfile("test.out")

    os.remove("test.out")

def test_WavefunctionOutputGeneratorFileIsCorrect():
    manager = scine.core.ModuleManager.get_instance()
    calc = manager.get("calculator", "TEST")
    wf_gen = scine.core.to_wf_generator(calc)

    assert wf_gen is not None

    wf_gen.wavefunction2file("test.out")

    assert read_file("test.out") == scine.core.to_wf_generator(calc).output_wavefunction()

    os.remove("test.out")

def test_NoneTypeIfNotCastable():
    manager = scine.core.ModuleManager.get_instance()
    calc = manager.get("calculator", "LENNARDJONES")
    assert scine.core.to_wf_generator(calc) is None
