""" Test suite for the api module.

The script can be executed on its own or incorporated into a larger test suite.
However the tests are run, be aware of which version of the module is actually
being tested. If the library is installed in site-packages, that version takes
precedence over the version in this project directory. Use a virtualenv test
environment or setuptools develop mode to test against the development version.

"""
import pytest
from AstroPhotography.api import *  # tests __all__


def test_hello():
    """ Test the hello() function.

    """
    assert hello() == "Hello, World!"
    return


def test_hello_name():
    """ Test the hello() function with a name.

    """
    assert hello("foo") == "Hello, foo!"
    return

class RawConvTest(object):
    """Tests RawConv methods"""
    
    def test_split(self):
        # TODO work out how to read fixed files in tests.
        # best info seem so far is https://pybit.es/pytest-fixtures.html
        ##rawfile = 
        ##rawconv = RawConv(rawfile)
        
        # Check split without black-level subtraction
        black = False
        ##r_im, g1_im, b_im, g2_im = rawconv.split(subtract_black=black)
        assert 3 == 5
        
        # Check split with black-level subtraction.
        black = True
        ##r_im, g1_im, b_im, g2_im = rawconv.split(subtract_black=black)
        return

# Make the script executable.

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
