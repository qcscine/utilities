__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as su

setting_descriptors = [
    su.BoolDescriptor,
    su.IntDescriptor,
    su.DoubleDescriptor,
    su.StringDescriptor,
    su.FileDescriptor,
    su.DirectoryDescriptor,
    su.OptionListDescriptor,
    su.ParametrizedOptionListDescriptor,
    su.IntListDescriptor,
    su.DoubleListDescriptor,
    su.StringListDescriptor,
    su.CollectionListDescriptor,
    su.DescriptorCollection
]


def test_setting_descriptors():
    # Constructible with a string, exception is CollectionListDescriptor
    expect_string_constructible = filter(
        lambda t: t != su.CollectionListDescriptor,
        setting_descriptors
    )

    def string_constructible(t):
        description = "hello"
        if not isinstance(t(description), t):
            return False

        return t(description).property_description == description

    assert all(map(string_constructible, expect_string_constructible))

    # Most should have the default_value property, there are some exceptions
    default_value_exceptions = [
        su.OptionListDescriptor,
        su.ParametrizedOptionListDescriptor,
        su.CollectionListDescriptor,
        su.DescriptorCollection
    ]
    expect_default_value = filter(
        lambda t: t not in default_value_exceptions,
        setting_descriptors
    )

    def default_value_noraise(t):
        if "default_value" not in dir(t):
            print("{} does not have a default_value member".format(t))
            return False

        instance = t("description")
        _ = instance.default_value
        return True

    assert all(map(default_value_noraise, expect_default_value))

    # All should have a default generic value property whose access shouldn't
    # raise an exception
    generic_noraise_exceptions = [
        su.OptionListDescriptor,
        su.ParametrizedOptionListDescriptor,
        su.CollectionListDescriptor  # Can't construct with str only
    ]
    expect_generic_noraise = filter(
        lambda t: t not in generic_noraise_exceptions,
        setting_descriptors
    )

    def default_generic_noraise(t):
        instance = t("description")
        _ = instance.default_generic_value
        return True

    assert all(map(default_generic_noraise, expect_generic_noraise))
