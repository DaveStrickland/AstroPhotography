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
    
    def _get_split_data(self, channel, black=False):
        """Returns expected postage-stamp data for a given channel.
        
        :param channel: Bayer channel, must be one on 'R', 'G1', 'B', 'G2'
        :param black: Boolean whether black-level subtraction should have
           been performed.
        """
        # Expected values obtained using split.m octave script.
        
        if not black:
            if 'R' in channel:
                # ../../../capture000003_noblack_r.tiff 451 464 2851 2864
                oct_ind = [451, 464, 2851, 2864]
                # Array 14x14 has 196 pixels
                meanval=120
                stdval=213
                minval=0
                maxval=621
                sumval=23579
                imgdata = [[483, 0, 446, 0, 420, 0, 389, 0, 360, 0, 337, 0, 321, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [516, 0, 489, 0, 451, 0, 412, 0, 382, 0, 357, 0, 339, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [554, 0, 528, 0, 481, 0, 445, 0, 417, 0, 380, 0, 361, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [591, 0, 566, 0, 517, 0, 481, 0, 443, 0, 401, 0, 379, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [607, 0, 586, 0, 553, 0, 524, 0, 475, 0, 442, 0, 409, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [618, 0, 611, 0, 588, 0, 556, 0, 517, 0, 473, 0, 433, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [621, 0, 615, 0, 610, 0, 588, 0, 549, 0, 501, 0, 457, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
            elif 'G1' in channel:
                # ../../../capture000003_noblack_g1.tiff 451 464 2851 2864
                oct_ind = [451, 464, 2851, 2864]
                # Array 14x14 has 196 pixels
                meanval=168
                stdval=304
                minval=0
                maxval=971
                sumval=32968
                imgdata = [[0, 678, 0, 618, 0, 552, 0, 501, 0, 451, 0, 406, 0, 367],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 740, 0, 678, 0, 609, 0, 545, 0, 493, 0, 443, 0, 404],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 818, 0, 748, 0, 652, 0, 597, 0, 544, 0, 486, 0, 437],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 875, 0, 812, 0, 733, 0, 662, 0, 591, 0, 525, 0, 480],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 919, 0, 884, 0, 813, 0, 735, 0, 663, 0, 597, 0, 530],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 953, 0, 933, 0, 862, 0, 803, 0, 726, 0, 634, 0, 581],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 960, 0, 971, 0, 947, 0, 862, 0, 802, 0, 714, 0, 634],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
                
            elif 'G2' in channel:
                # ../../../capture000003_noblack_b.tiff 451 464 2851 2864
                oct_ind = [451, 464, 2851, 2864]
                # Array 14x14 has 196 pixels
                meanval=127
                stdval=226
                minval=0
                maxval=658
                sumval=24930
                imgdata = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 521, 0, 471, 0, 441, 0, 414, 0, 390, 0, 367, 0, 336],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 551, 0, 511, 0, 474, 0, 445, 0, 408, 0, 383, 0, 346],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 593, 0, 553, 0, 508, 0, 469, 0, 435, 0, 401, 0, 379],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 612, 0, 588, 0, 550, 0, 507, 0, 467, 0, 435, 0, 406],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 634, 0, 627, 0, 595, 0, 547, 0, 509, 0, 458, 0, 435],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 658, 0, 647, 0, 618, 0, 589, 0, 550, 0, 492, 0, 455],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 643, 0, 642, 0, 639, 0, 624, 0, 592, 0, 530, 0, 485]]
                    
            elif 'B' in channel:
                # ../../../capture000003_noblack_g2.tiff 451 464 2851 2864
                oct_ind = [451, 464, 2851, 2864]
                # Array 14x14 has 196 pixels
                meanval=182
                stdval=327
                minval=0
                maxval=976
                sumval=35671
                imgdata = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [757, 0, 661, 0, 613, 0, 543, 0, 492, 0, 447, 0, 404, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [833, 0, 742, 0, 662, 0, 602, 0, 543, 0, 479, 0, 449, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [871, 0, 808, 0, 753, 0, 666, 0, 597, 0, 532, 0, 491, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [931, 0, 875, 0, 814, 0, 743, 0, 664, 0, 592, 0, 534, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [943, 0, 921, 0, 878, 0, 811, 0, 729, 0, 652, 0, 576, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [970, 0, 970, 0, 929, 0, 873, 0, 810, 0, 713, 0, 632, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [942, 0, 961, 0, 976, 0, 929, 0, 859, 0, 795, 0, 704, 0]]
            else:
                # TODO pytest error
                pytest.fail('Unexpected channel {}. Expecting on of R, G1, B. or G2'.format(channel))
        else:
            # black subtraction HAS been performed
            pytest.fail('Have not implemented test datat for black-level subtracted data')
            
        return oct_ind, meanval, stdval, minval, maxval, sumval, imgdata
            
    def _compare_channel_data(self, channel, channel_img, black):
        
        # Get expected data
        oct_ind, meanval, stdval, minval, maxval, sumval, imgdata = self._get_split_data(channel, black)
        
        # Extract sub-image from current data
        
        # Compare
        pytest.fail('Test not complete')
        return
        
    def test_split(self):
        """Tests RawConv.split(), which ultimately tests a lot of the most
        basic processing used by any call to RawConv.
        """
        
        rawconv = rawconv_tfile()
                
        # Check split without black-level subtraction
        black = False
        r_im, g1_im, b_im, g2_im = rawconv.split(subtract_black=black)
        
        
        self._compare_channel_data('R', r_im, black)
        
        # Check split with black-level subtraction.
        black = True
        r_im, g1_im, b_im, g2_im = rawconv.split(subtract_black=black)
        return

# Make the script executable.

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
