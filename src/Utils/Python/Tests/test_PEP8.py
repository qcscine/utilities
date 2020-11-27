__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as su


def test_module_pep8():
    # Ensure all submodules are lowercase
    def is_module(x): return isinstance(x, type(su))
    assert all([x.islower() for x in dir(su) if is_module(getattr(su, x))])

    # Ensure all constants are uppercase
    def is_fp_or_int(x): return isinstance(
        x, type(3.0)) or isinstance(x, type(3))
    assert all([x.isupper() for x in dir(su) if is_fp_or_int(getattr(su, x))])
