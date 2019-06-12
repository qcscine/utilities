import pytest
import scine_utils_os as scine

def test_ElementData():
  a = scine.ElementDataSingleton.instance()

  # Access via element type and symbol string
  element_data_h = a[scine.ElementType.H]
  element_data_f = a["F"]

  assert element_data_h.symbol == "H"
  assert element_data_f.Z == 9
