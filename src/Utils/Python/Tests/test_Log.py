__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as scine


def compose():
    return "I composed something for you"


def test_log():
    log = scine.core.Log()
    log.debug.line("Hi")
    log.debug.lazy(compose)  # Nullary compose fn
    assert not log.debug.has_sinks()

    log.debug.add("cout", scine.core.Log.cout_sink())
    assert log.debug.has_sinks()
    log.debug.remove("cout")
    assert not log.debug.has_sinks()

    log.debug.clear()
