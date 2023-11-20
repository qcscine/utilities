__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as scine


def test_FormulaGenerator():
    elements = [scine.ElementType.C, scine.ElementType.F, scine.ElementType.H]
    assert scine.generate_chemical_formula(elements) == "CHF"
