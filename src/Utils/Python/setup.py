import setuptools

# Read README.md for the long description
with open("README.md", "r") as fh:
  long_description = fh.read()

# Define the setup
setuptools.setup(
  name="scine_utils_os",
  version="0.1.0",
  author="ETH Zurich, Laboratory for Physical Chemistry, Reiher Group",
  author_email="scine@phys.chem.ethz.ch",
  description="Open source utilities libraries",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://www.scine.ethz.ch",
  packages=[""],
  package_dir={"": "."},
  package_data={"": ["scine_utils_os.so"]},
  classifiers=[
    "Programming Language :: Python",
    "Programming Language :: C++",
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Topic :: Scientific/Engineering :: Chemistry"
  ],
  zip_safe=False,
  test_suite='pytest',
  tests_require=['pytest']
)
