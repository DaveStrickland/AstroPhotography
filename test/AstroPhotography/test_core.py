""" Test suite for the core module.

The script can be executed on its own or incorporated into a larger test suite.
However the tests are run, be aware of which version of the module is actually
being tested. If the library is installed in site-packages, that version takes
precedence over the version in this project directory. Use a virtualenv test
environment or setuptools develop mode to test against the development version.

"""
import os
import os.path
import pytest
from AstroPhotography.core.RawConv import RawConv

def rawconv_tfile():
    """Returns a rawconv object using the standard test file
    
    Assumes tests are run from the distribution base directory, i..e.
    the tests directory is a subdirectory of ps.getcwd().
    """
    
    cwd = os.getcwd()
    print(dir())
    test_data_dir = 'test/AstroPhotography/data/'
    test_data_file = 'capture000003.cr2'
    rawfile = os.path.join(cwd, test_data_dir, test_data_file)
    if not os.path.isfile(rawfile):
        pytest.fail('Could not find expected test RAW file {}'.format(rawfile))
    rawconv = RawConv(rawfile)
    return rawconv
    

class RawConvTest(object):
    """Tests RawConv methods"""
        
    def test_split(self):
        rawconv = rawconv_tfile()
                
        # Check split without black-level subtraction
        black = False
        r_im, g1_im, b_im, g2_im = rawconv.split(subtract_black=black)
        assert 3 == 5
        
        # Check split with black-level subtraction.
        black = True
        r_im, g1_im, b_im, g2_im = rawconv.split(subtract_black=black)
        return

# Make the script executable.

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
