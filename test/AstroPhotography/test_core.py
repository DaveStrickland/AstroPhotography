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
import numpy as np
from AstroPhotography.core.RawConv import RawConv

@pytest.fixture
def rawconv_tfile():
    """Returns a rawconv object using the standard test file
    
    Assumes tests are run from the distribution base directory, i..e.
    the tests directory is a subdirectory of ps.getcwd().
    """
    
    cwd = os.getcwd()
    test_data_dir = 'test/AstroPhotography/data/'
    test_data_file = 'capture000003.cr2'
    rawfile = os.path.join(cwd, test_data_dir, test_data_file)
    if not os.path.isfile(rawfile):
        pytest.fail('Could not find expected test RAW file {}'.format(rawfile))
    rawconv = RawConv(rawfile)
    return rawconv

@pytest.fixture(params=("R", "G1", "B", "G2"))  
def channel(request):
    """Returns a parameterized channel for the split command tests"""
    return request.param
    
@pytest.fixture(params=("daylight", "camera"))  
def wb_method(request):
    """Returns a parameterized whitebalance method"""
    return request.param


class RawConvTest(object):
    """Tests RawConv methods"""
    
    def _get_split_data(self, channel, black=False):
        """Returns expected postage-stamp data for a given channel.
        
        This test data set is limited in the sense it is data file
        specific, and requires some by-hand work to get data from the
        octave scripts into the python. It also presumes that some
        aspects of the split command were correct in order to generate
        the test channel files read in octave.
        
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
                
            elif 'B' in channel:
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
                    
            elif 'G2' in channel:
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
            # black-level subtraction has been performed.
            if 'R' in channel:
                # black subtraction HAS been performed
                # ../../../capture000003_r.tiff 451 464 2851 2864
                oct_ind = [451, 464, 2851, 2864]
                # Array 14x14 has 196 pixels
                meanval=56
                stdval=107
                minval=0
                maxval=365
                sumval=11035
                imgdata = [[227, 0, 190, 0, 164, 0, 133, 0, 104, 0, 81, 0, 65, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [260, 0, 233, 0, 195, 0, 156, 0, 126, 0, 101, 0, 83, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [298, 0, 272, 0, 225, 0, 189, 0, 161, 0, 124, 0, 105, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [335, 0, 310, 0, 261, 0, 225, 0, 187, 0, 145, 0, 123, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [351, 0, 330, 0, 297, 0, 268, 0, 219, 0, 186, 0, 153, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [362, 0, 355, 0, 332, 0, 300, 0, 261, 0, 217, 0, 177, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [365, 0, 359, 0, 354, 0, 332, 0, 293, 0, 245, 0, 201, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
            elif 'G1' in channel:
                # ../../../capture000003_g1.tiff 451 464 2851 2864
                oct_ind = [451, 464, 2851, 2864]
                # Array 14x14 has 196 pixels
                meanval=104
                stdval=200
                minval=0
                maxval=715
                sumval=20424
                imgdata = [[0, 422, 0, 362, 0, 296, 0, 245, 0, 195, 0, 150, 0, 111],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 484, 0, 422, 0, 353, 0, 289, 0, 237, 0, 187, 0, 148],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 562, 0, 492, 0, 396, 0, 341, 0, 288, 0, 230, 0, 181],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 619, 0, 556, 0, 477, 0, 406, 0, 335, 0, 269, 0, 224],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 663, 0, 628, 0, 557, 0, 479, 0, 407, 0, 341, 0, 274],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 697, 0, 677, 0, 606, 0, 547, 0, 470, 0, 378, 0, 325],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 704, 0, 715, 0, 691, 0, 606, 0, 546, 0, 458, 0, 378],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
            elif 'B' in channel:
                # ../../../capture000003_b.tiff 451 464 2851 2864
                oct_ind = [451, 464, 2851, 2864]
                # Array 14x14 has 196 pixels
                meanval=63
                stdval=119
                minval=0
                maxval=402
                sumval=12386
                imgdata = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 265, 0, 215, 0, 185, 0, 158, 0, 134, 0, 111, 0, 80],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 295, 0, 255, 0, 218, 0, 189, 0, 152, 0, 127, 0, 90],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 337, 0, 297, 0, 252, 0, 213, 0, 179, 0, 145, 0, 123],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 356, 0, 332, 0, 294, 0, 251, 0, 211, 0, 179, 0, 150],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 378, 0, 371, 0, 339, 0, 291, 0, 253, 0, 202, 0, 179],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 402, 0, 391, 0, 362, 0, 333, 0, 294, 0, 236, 0, 199],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 387, 0, 386, 0, 383, 0, 368, 0, 336, 0, 274, 0, 229]]
            elif 'G2' in channel:
                # ../../../capture000003_g2.tiff 451 464 2851 2864
                oct_ind = [451, 464, 2851, 2864]
                # Array 14x14 has 196 pixels
                meanval=118
                stdval=221
                minval=0
                maxval=720
                sumval=23127
                imgdata = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [501, 0, 405, 0, 357, 0, 287, 0, 236, 0, 191, 0, 148, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [577, 0, 486, 0, 406, 0, 346, 0, 287, 0, 223, 0, 193, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [615, 0, 552, 0, 497, 0, 410, 0, 341, 0, 276, 0, 235, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [675, 0, 619, 0, 558, 0, 487, 0, 408, 0, 336, 0, 278, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [687, 0, 665, 0, 622, 0, 555, 0, 473, 0, 396, 0, 320, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [714, 0, 714, 0, 673, 0, 617, 0, 554, 0, 457, 0, 376, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [686, 0, 705, 0, 720, 0, 673, 0, 603, 0, 539, 0, 448, 0]]
            else:
                # TODO pytest error
                pytest.fail('Unexpected channel {}. Expecting on of R, G1, B. or G2'.format(channel))
                        
        return oct_ind, meanval, stdval, minval, maxval, sumval, np.array(imgdata)
    
    def _get_split_wb(self, black=False):
        """Returns the expected white-balance multipliers for the test postage stamps
 
        This test data set is limited in the sense it is data file
        specific, and requires some by-hand work to get data from the
        octave scripts into the python. It also presumes that some
        aspects of the split command were correct in order to generate
        the test channel files read in octave.
        
        :param black: If true, then black-level subtraction has been
          applied. Otherwise the black-level subtraction has not been
          applied.
        """
        
        # Raw image sums for the postage stamps
        # TODO: Should really read expected postage stamps and calculate sums.
        if black:
            wb_in = [11035, 20424, 12386, 23127]
        else:
            wb_in = [23579, 32968, 24930, 35671]
        
        # Convert numbers to expected conversion factors, which are of
        # the form: camera_wb [1.849074074074074, 1.0, 2.160185185185185, 1.0]
        max_val = max(wb_in)
        wb_list = []
        for idx, val in enumerate(wb_in):
            wb_list.append(max_val / val)
        
        return wb_list        
    
    def _compare_channel_data(self, channel, channel_img, black):
        
        # Get expected data
        oct_ind, meanval, stdval, minval, maxval, sumval, imgdata = self._get_split_data(channel, black)
        
        # numpy indices are 0-based, but ranges are exclusive on the upper end
        # so for an octave range a:b the numpy equivalent is a-1:b.
        r1 = oct_ind[0] - 1
        r2 = oct_ind[1]
        c1 = oct_ind[2] - 1
        c2 = oct_ind[3]
        
        # Extract sub-image from current data in channel image
        sub_img = channel_img[r1:r2,c1:c2]
        
        # Compare raw sub-image shapes, then arrays themselves
        assert imgdata.shape == sub_img.shape
        assert np.array_equal(imgdata, sub_img) == True
        
        # Compute and compare statistics (which should also match if the
        # arrays match.
        
        return
        
    def test_split(self, rawconv_tfile, channel):
        """Tests RawConv.split(), which ultimately tests a lot of the most
        basic processing used by any call to RawConv.
        
        :param rawconv_tfile: RawConv object loaded with the test file
        :param channel: A string denoting the
          channel to process.
        """
        
        # Check split without black-level subtraction
        black = False
        r_im, g1_im, b_im, g2_im = rawconv_tfile.split(subtract_black=black)
        
        if 'R' in channel:
            img = r_im
        elif 'G1' in channel:
            img = g1_im
        elif 'B' in channel:
            img = b_im
        elif 'G2' in channel:
            img = g2_im
            
        self._compare_channel_data(channel, img, black)
        
        # Check split with black-level subtraction.
        black = True
        r_im, g1_im, b_im, g2_im = rawconv_tfile.split(subtract_black=black)
        
        if 'R' in channel:
            
            img = r_im
        elif 'G1' in channel:
            img = g1_im
        elif 'B' in channel:
            img = b_im
        elif 'G2' in channel:
            img = g2_im
            
        self._compare_channel_data(channel, img, black)
        return

    def test_get_whitebalance(self, rawconv_tfile, wb_method):
        """Tests that get_whitebalance returns the expected white balance multipliers
        
        :param rawconv_tfile: RawConv object loaded with the test file
        :param wb_method: String denoting the whitebalance method to use.
        """
        
        # Allowed floating point tolerance, better than the dynamic range
        # of any realistic camera.
        toler = 1.0-6
        
        if wb_method == 'daylight':
            expected_wb = [2.63078458003413, 1.0, 1.2493076831610352, 1.0]
        elif wb_method == 'camera':
            expected_wb = [1.849074074074074, 1.0, 2.160185185185185, 1.0]
        else:
            pytest.fail('Error: unexpected whitebalance method {} tested.'.format(wb_method))
            
        wb = rawconv_tfile.get_whitebalance(wb_method)

        # Check that whitebalance lists have the same number of elements,
        # then check that each value matches with the tolerance allowed.
        assert len(wb) == len(expected_wb)
        for idx, val in enumerate(wb):
            assert val == pytest.approx(expected_wb[idx],
                abs=toler)
        
        return

# Make the script executable.

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
