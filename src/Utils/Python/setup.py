__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import setuptools

# Read README.rst for the long description
with open("README.rst", "r") as fh:
    long_description = fh.read()


class EmptyListWithLength(list):
    """ Makes the wheel a binary distribution and platlib compliant. """

    def __len__(self):
        return 1


# Define the setup
setuptools.setup(
    name="scine_utilities",
    version="@PROJECT_VERSION@",
    author="ETH Zurich, Laboratory of Physical Chemistry, Reiher Group",
    author_email="scine@phys.chem.ethz.ch",
    description="Utilities used in all SCINE modules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://www.scine.ethz.ch",
    packages=["scine_utilities"],
    package_data={"scine_utilities": ["scine_utilities.*", "*.txt", "utilsos.module.*"@utils_PY_DEPS@]},
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    install_requires=["numpy", "scipy"],
    zip_safe=False,
    test_suite='pytest',
    tests_require=['pytest'],
    ext_modules=EmptyListWithLength()
)
