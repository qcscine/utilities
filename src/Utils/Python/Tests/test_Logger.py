import pytest
import scine_utils_os as scine

def test_log():
  # test out the severity levels
  debug_severity = scine.Log.severity_level.debug
  warning_severity = scine.Log.severity_level.warning
  assert debug_severity != warning_severity

  scine.Log.enable_logging()
  scine.Log.start_console_logging(debug_severity)
  scine.Log.debug("Test debug level severity string")
  scine.Log.stop_console_logging()
  scine.Log.warning("Not shown warning severity string")
  scine.Log.disable_logging()
