import pytest
import scine_utils_os as scine

def test_ElementTypes():
  # Test the string magic function
  assert str(scine.ElementType.K) == "K"
  assert str(scine.ElementType.Cl) == "Cl"
  # Test the arithmetic properties of the enum
  assert scine.ElementType.Ar == scine.ElementType.Ar
  assert scine.ElementType.Ar != scine.ElementType.Ti
