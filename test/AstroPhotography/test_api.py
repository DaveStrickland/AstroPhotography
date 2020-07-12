""" Test suite for the api module.

The script can be executed on its own or incorporated into a larger test suite.
However the tests are run, be aware of which version of the module is actually
being tested. If the library is installed in site-packages, that version takes
precedence over the version in this project directory. Use a virtualenv test
environment or setuptools develop mode to test against the development version.

"""
import os
import os.path
import pytest
from AstroPhotography.api import *  # tests __all__

# Make the script executable.
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
