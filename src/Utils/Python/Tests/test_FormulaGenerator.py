import pytest
import scine_utils_os as scine

def test_FormulaGenerator():
  elements = [scine.ElementType.C, scine.ElementType.F, scine.ElementType.H]
  assert scine.generate_chemical_formula(elements) == "CHF"
