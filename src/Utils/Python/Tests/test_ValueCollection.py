__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import pytest
import pickle
import scine_utilities as su


def test_value_collection():
    param_opt_value = su.ParametrizedOptionValue(
        "hello",
        su.ValueCollection()
    )
    generic_values = {
        "boolean": True,
        "int": 4,
        "float": 1.5,
        "string": "hello",
        "ValueCollection": su.ValueCollection(),
        "ParametrizedOptionValue": param_opt_value,
        "IntList": [1, 2, 3],
        "DoubleList": [1.1, 2.2, 3.3],
        "StringList": ["yes", "no", "maybe so"],
        "CollectionList": [
            su.ValueCollection(),
            su.ValueCollection()
        ]
    }

    collection = su.ValueCollection(generic_values)

    for key in generic_values.keys():
        assert key in collection
        assert isinstance(generic_values[key], type(collection[key]))
        assert generic_values[key] == collection[key]

    with pytest.raises(RuntimeError):
        collection.update({"int": "not an int"})


def test_pickle():
    param_opt_value = su.ParametrizedOptionValue(
        "hello",
        su.ValueCollection()
    )
    generic_values = {
        "boolean": True,
        "int": 4,
        "float": 1.5,
        "string": "hello",
        "ValueCollection": su.ValueCollection(),
        "ParametrizedOptionValue": param_opt_value,
        "IntList": [1, 2, 3],
        "DoubleList": [1.1, 2.2, 3.3],
        "StringList": ["yes", "no", "maybe so"],
        "CollectionList": [
            su.ValueCollection(),
            su.ValueCollection()
        ]
    }

    collection = su.ValueCollection(generic_values)

    file_name = ".test_valuecollection.pkl"
    with open(file_name, "wb") as f:
        pickle.dump(collection, f)

    with open(file_name, "rb") as f:
        read_vc = pickle.load(f)

    assert collection == read_vc
    if os.path.exists(file_name):
        os.remove(file_name)
