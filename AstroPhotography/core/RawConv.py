""" Implements the RawConv, which encapsulates raw data manipulations.
"""

from .logger import logger
import rawpy
import numpy as np
import os.path
import ast
import exifread     # Read camera RAW exif data

from .. import __version__

class RawConv:
    """ Higher level conversions and processing of RAW image files to graphics formats 
        and/or FITS image format.
    """
    
    def __init__(self, rawfile):
        """
        RawConv constructor loads the raw file and performs initial processing.
        """
        
        self._supported_colors = ['RGBG']
        self._makernote        = False      # Don't read MakerNote
        
        # Percentiles. Note first and last must be 0.01 and 99.99 as these
        # values are used in image renormalization
        self._ipctls = [0.01, 0.1, 0.5, 1.0, 5.0, 10, 25, 50, 75, 90, 95, 99, 99.5, 99.9, 99.99]
        self._load(rawfile)
        return
        
    def __del__(self):
        """
        RawConv class destructor
        """
        
        if self._rawpy is not None:
            self._rawpy.close()
        return

    def _unsupported_colors(self):
        """
        Raises an exception if the image color description is not 
        supported.
        """
        
        err_str = 'Image colors ({}) not in supported list of {}.'.format(self._color_desc, 
            self._supported_colors)
        raise NotImplementedError(err_str)        
        return

    def _image_stats(self, numpy_array, pctiles=None):
        """
        Calculates basic statistics on an input array, including percentiles
        """
        if pctiles is None:
            ipctls = [0.01, 0.1, 0.5, 1.0, 5.0, 10, 25, 50, 75, 90, 95, 99, 99.5, 99.9, 99.99]
        else:
            ipctls = pctiles
        
        minval    = np.nanmin(numpy_array)
        maxval    = np.nanmax(numpy_array)
        meanval   = np.nanmean(numpy_array)
        stdval    = np.nanstd(numpy_array)
        medianval = np.nanmedian(numpy_array)
        
        # percentiles
        opctls = np.nanpercentile(numpy_array, ipctls)    
        return minval, maxval, meanval, stdval, medianval, opctls

    def _load(self, rawfile):
        """
        Reads the RAW file nad performs some preliminary processing on it.
        
        :param rawfile: Input RAW file recognized by rawpy.
        """
        self._rawfile    = rawfile
        self._rawpy      = None
        if os.path.isfile(rawfile):
            logger.debug('RawConv loading {}'.format(rawfile))
            self._rawpy      = rawpy.imread(rawfile)
        else:
            err_msg = 'RawConv cannot load {}. Not a valid path or file.'.format(rawfile)
            logger.error(err_msg)
            raise RuntimeError(err_msg)
        
        # Read EXIF tags
        self._read_exif(rawfile, use_makernote=self._makernote)
        
        # Common image characteristics
        self._nrows        = self._rawpy.raw_image_visible.shape[0]
        self._ncols        = self._rawpy.raw_image_visible.shape[1]
        self._color_desc   = str(self._rawpy.color_desc)

        # TODO have this based on color description
        # Know that for RGBG the color indices are 0, 1, 2, 3
        self.R  = 0
        self.G1 = 1
        self.B  = 2
        self.G2 = 3

        if self._color_desc not in self._supported_colors:
            self._unsupported_colors

        self._black_levels = self._rawpy.black_level_per_channel
        self._black_subtracted = False
        self._default_whitebalances()
        return
        
    def _build_raw_channel_images(self):
        """
        Create color map, masks, and per-channel raw images.
        
        Used only by the split command.
        """
        self._color_map  = self._rawpy.raw_colors_visible.copy()
        self._mask_r     = self._color_map == self.R
        self._mask_g1    = self._color_map == self.G1
        self._mask_b     = self._color_map == self.B
        self._mask_g2    = self._color_map == self.G2
        
        self._rawim_r    = np.where(self._mask_r,  self._rawpy.raw_image_visible, 0)
        self._rawim_g1   = np.where(self._mask_g1, self._rawpy.raw_image_visible, 0)
        self._rawim_b    = np.where(self._mask_b,  self._rawpy.raw_image_visible, 0)
        self._rawim_g2   = np.where(self._mask_g2, self._rawpy.raw_image_visible, 0)

        return

    def _default_whitebalances(self):
        """
        Convert default camera and daylight whitebalances into
        a more consistent form
           
        For the Canon Digital Rebel XTi, the format of the camera
        whitebalance and the daylight whitebalance format are
        inconsistent:
        - camera:   1997.0, 1080.0, 2333.0, 1080.0
        - daylight: 2.4238, 0.9213, 1.1510, 0.0

        This function converts them into a uniform set of R G B G
        factors where the lowest value is 1.0.
        """
        self._wb_camera   = []
        self._wb_daylight = []

        tmp_cam = self._rawpy.camera_whitebalance
        tmp_min = min(tmp_cam)
        for idx, val in enumerate(tmp_cam):
            self._wb_camera.append( val/tmp_min )
        
        # If last item is zero and corresponds to G, then use other G.
        tmp_day = self._rawpy.daylight_whitebalance
        if tmp_day[3] == 0 and self._color_desc[3] == 'G':
            tmp_day[3] = tmp_day[1]
        tmp_min = min(tmp_day)
        for idx, val in enumerate(tmp_day):
            self._wb_daylight.append( val/tmp_min )
        
        return

    def _subtract_black_levels(self):
        """
        Subtract black levels from sub-bands.
        """
        if self._black_subtracted:
            # Already black subtracted
            logger.warning('Camera black levels already subtracted. Not repeating.')
            return
            
        r_bg  = int( self._black_levels[self.R]  )
        g1_bg = int( self._black_levels[self.G1] )
        b_bg  = int( self._black_levels[self.B]  )
        g2_bg = int( self._black_levels[self.G2] )
        logger.info(f'Subtracting camera black levels: R={r_bg} G1={g1_bg} B={b_bg} G2={g2_bg}')

        self._rawim_r  = self._safe_subtract(self._rawim_r,
            self._mask_r, 
            r_bg)
        self._rawim_g1 = self._safe_subtract(self._rawim_g1,
            self._mask_g1, 
            g1_bg)
        self._rawim_b  = self._safe_subtract(self._rawim_b,
            self._mask_b, 
            b_bg)
        self._rawim_g2 = self._safe_subtract(self._rawim_g2,
            self._mask_g2, 
            g2_bg)
        self._black_subtracted = True

    def _read_exif(self, input_file, use_makernote=None):
        """
        Read EXIF metadata from the specified input file, storing it
        in self._exif_dict
        
        This function uses ExifRead, because the exif package cannot
        read RAW files. Only deciphered tags are stored, becuase we lack
        the capability to correctly rewrite them in the output file.
        
        :param input_file: Input RAW file path/name.
        :param use_makernote: In most cases RAW file MakerNote tags are not
          recognized by the exif writing code (the exif package) and
          hence cannot be re-written even if ExifRead can read them.
          However, if use_makernote is True they will be read and stored
          in the returned exif_dict. 
        """
        
        store_makernote = False
        if use_makernote is not None:
            store_makernote = use_makernote
        
        self._exif_dict = None
        # Open image file for reading (must be in binary mode)
        with open(input_file, 'rb') as f:
            # Return Exif tags
            tags = exifread.process_file(f)
        
            self._exif_dict = {}
            for tag in tags.keys():
                # Skip MakerNote?
                if not store_makernote:
                    if 'MakerNote' in tag:
                        continue
                
                if tag in ('JPEGThumbnail', 'TIFFThumbnail', 'EXIF UserComment', 'EXIF MakerNote'):
                    continue
                elif ' Tag ' in tag:
                    continue
                elif 'IFD 2 ' in tag:
                    newtag = tag.replace('IFD 2 ', 'IFD2 ')
                    self._exif_dict[newtag] = tags[tag]
                elif 'IFD 3 ' in tag:
                    newtag = tag.replace('IFD 3 ', 'IFD3 ')
                    self._exif_dict[newtag] = tags[tag]
                else:
                    self._exif_dict[tag] = tags[tag]
        
        if self._exif_dict is not None:
            logger.info(f'Input file {input_file} has {len(self._exif_dict)} items of EXIF metadata.')
            logger.debug(f'The EXIF metadata items from the file are: {self._exif_dict.keys()}')
            for key, val in self._exif_dict.items():
                logger.debug(f'  key={key:<30s}  val={val}')
            
        else:
            logger.warning(f'Input file {input_file} does not contain any EXIF metadata.')
        
        return

    def _safe_subtract(self, data_arr, data_mask, val_to_subtract):
        """Safe subtraction of one array or scalar from another array
           avoiding unsigned integer wrap-around.
        
        It is possible for a value within a data array to be less
        than the value to subtract, so we need to carefully reset those
        pixels to the black level to avoid them wrapping around to
        a very large number.
        
        :param data_array: 2-D uint array to subract values from. Only
          those pixels within data_mask will be modified. 
        :param data_mask:  2-D boolean mask specifying where data_array
          is valid.
        :param val_to_subtract: 2-D uint array of values to subtract from
          data_array, or a scalar that will be broadcast to match the shape
          of data_array.
        :returns: 2-D uint array.
        """
        odd_mask = data_arr < val_to_subtract
        num_odd  = np.sum(odd_mask)
        pct_odd  = 100 * num_odd / data_arr.size
        
        # Note: this gives a misleading answer when raw individual color
        # channels are used, because e.g. an R image will have the G1, B,
        # and G2 pixel values all set to zero.
        logger.warning('Input array has {} pixels ({:.2f}%) less than array we will subract from it.'.format(num_odd, pct_odd))
        if num_odd > 0:
            # Reset those pixels to the value we're going to subtract.
            tmp_raw = np.where(odd_mask, val_to_subtract, data_arr)
        else:
            # Just use normal rawim
            tmp_raw = data_arr.copy()
        
        # np.subtract will leave original output array values (if out 
        # supplied and not None) where the where mask is False.
        # Elsewhere, the pixels will be the first input - second input.
        return np.subtract(tmp_raw, 
            val_to_subtract, 
            where=data_mask).copy()
            
    def _get_whitebalance_from_region(self, wb_method):
        """Determine the whitebalance from a rectangular patch or the entire image
        
        :param wb_method: The method can be 'auto' to select the entire 
          visible image, or 'region[rowmin, rowmax, colmin, colmax]'
        """
        
        if wb_method == 'auto':
            pixel_list=[0, self._nrows, 0, self._ncols]
        elif 'region' in wb_method:
            # Split off region from the region spec and parse as list.
            str_list = wb_method.replace('region', '')
            pixel_list = ast.literal_eval(str_list)
            
        # Extract data sums within the regions. Masks are also necessary as there
        # may well be a different number of valid pixels per band.
        sum_r,  nvalid_r  = self._get_sum_in_region(self._rawim_r,  self._mask_r,  pixel_list)
        sum_g1, nvalid_g1 = self._get_sum_in_region(self._rawim_g1, self._mask_g1, pixel_list)
        sum_b,  nvalid_b  = self._get_sum_in_region(self._rawim_b,  self._mask_b,  pixel_list)
        sum_g2, nvalid_g2 = self._get_sum_in_region(self._rawim_g2, self._mask_g2, pixel_list)
             
        # Need to check that the number of pixels is greater than zero and the
        # data sum is greater than zero too. Otherwise, the 
        wb_avg = [] # Avg DN per pixel.
        for itr in zip([sum_r, sum_g1, sum_b, sum_g2], 
            [nvalid_r, nvalid_g1, nvalid_b, nvalid_g2], 
            ['R', 'G1', 'B', 'G2']):
            if itr[0] < 0:
                logger.warning('For band {}, the sum of data within the region is {}'.format(itr[2], itr[0]))
            elif itr[1] < 1:
                logger.error('For band {}, the number of valid pixels within the region is {}'.format(itr[2], itr[1]))
                
            wb_avg.append( itr[0]/itr[1] )
            
        max_val = max(wb_avg)
        wb_list = []
        for idx, val in enumerate(wb_avg):
            wb_list.append(max_val / val)
                
        return wb_list
        
    def _get_sum_in_region(self, data_arr, mask_arr, rpix_list):
        """Find the data sum in a rectangular region described by a list
        and optionally within a mask
        
        :param data_arr: 2-dimensional array of data values.
        :param mask_arr: An optional mask array, where only data where the mask
          is valid is selected. If not mask is to be used supply a None instead.
        :param rpix_list: A list of pixel indices correspond to rowmin,
          rowmax, colmin and colmax of the region with which the data sum is
          required. Note that these are inclusive and zero based.
        :returns: [sum, nvalid] The data sum within the specified region that
          is also valid within the optional mask, and the number of valid pixels
          within that region.
        """
        
        sumpix = 0
        numpix = 0
        
        if mask_arr is not None:
            masked = np.where(mask_arr, data_arr, 0)
        else:
            masked = data_arr
            
        # Note that numpy indices are (start,end] not (start, end)
        subregion = masked[rpix_list[0]:rpix_list[1]+1, rpix_list[2]:rpix_list[3]+1]
        sumpix = np.sum(subregion)
        
        if mask_arr is not None:
            numpix = np.sum( mask_arr[rpix_list[0]:rpix_list[1]+1, rpix_list[2]:rpix_list[3]+1] )
        else:
            numpix = subregion.size

        return sumpix, numpix
        
    def get_whitebalance(self, wb_method):
        """Return the whitebalance multiplies for a given white-balance
           calculation method.
        """
        
        # Check whitebalance method supplied.
        allowed_methods = ['daylight', 'camera', 'auto', 'region', 'user']
        # Split off method[specifier] to just check the method part.
        method = wb_method.split('[')[0]
        
        if method not in allowed_methods:
            err_str = f'Unexpected white balance method "{method}" not one of the allowed method: {allowed_methods}'
            logger.error(err_str)
            raise RuntimeError(err_str)
        
        wb_list = [1, 1, 1, 1]
        if wb_method == 'daylight':   
            wb_list = self._wb_daylight
        elif wb_method == 'camera':
            wb_list = self._wb_camera
        elif wb_method == "auto":
            wb_list = self._get_whitebalance_from_region(wb_method)
        elif 'region' in wb_method:
            wb_list = self._get_whitebalance_from_region(wb_method)
        elif 'user' in wb_method:
            # TODO parse user spec, hand off to function
            logger.warning(f'Whitebalance method {wb_method} not yet implemented')
            
            
        logger.debug(f'White balance values using method {method} adopted: {wb_list}')
        return wb_list

    def rgb(self, luminance_method='linear', subtract_black=True, 
        wb_method='auto', print_stats=None, renorm=None):
        """
        Create a 3-channel RGB image using the supplied white-balance
        values, with or without black-level subtraction.

        :param luminance_method: Luminance method to use. At present only the
           following method(s) are supported:
           - linear: Applies Bayer correction with a gamma=1 value to
             generate an RGB image, then applies the CCIR 601 coefficients
             to convert that to greyscale.
        :param subtract_black: If true the camera black levels will be 
          subtracted from the channel data.
        :param wb_method: Whitebalance method to user to determine white 
          balances. See get_whitebalance() for details of the allowed
          methods.
        :param print_stats: Output image statistics to the logger if True.
        :param renorm: Renormalize the image counts to fill the full
          16-bit dynamic range if True. Specifically it will linearly
          stretch the dynamic range from the 0.01th to 99.99th percentiles
          to fill the range 0 to (2^16)-1 ADU.
          
        Returns a numpy array of the image along with a dictionary of EXIF
        tags.
        """
        
        min_adu = 0
        max_adu = 2**16 - 1
        
        if print_stats is None:
            print_stats = False
        if renorm is None:
            renorm = False
        
        allowed_methods=['linear']
        if luminance_method not in allowed_methods:
            err_str = 'Unexpected luminance calculate method supplied to RawConv.rgb: {}. Allowed methods are: {}'.format(luminance_method, allowed_methods)
            logger.error(err_str)
        else:
            logger.info(f'Generating rgb image using {luminance_method} method.')
        
        # TODO this isn't used with linear method
        ##if subtract_black:
        ##    self._subtract_black_levels()
            
        wb_list = self.get_whitebalance(wb_method)
        
        if 'linear' in luminance_method:
            linear_gamma  = (1,1)
            default_gamma = (2.222, 4.5)
            rgb_im = self._rawpy.postprocess(gamma=linear_gamma, 
                no_auto_bright=True, no_auto_scale=True,
                output_bps=16, user_wb=wb_list)
                        
        # Calculate initial image statistics.
        minv, maxv, meanv, stdv, medianv, pctiles = self._image_stats(rgb_im, 
            self._ipctls)

        # TODO Renormalization?
        if renorm:
            # Use 0.01 and 99.99 percentiles            
            raw_min = pctiles[0]
            raw_max = pctiles[ len(pctiles)-1 ]
            logger.info(f'Renormalizing image from range [{raw_min:.2f}, {raw_max:.2f}] to [{min_adu}, {max_adu}] ADU.')
            logger.debug(f'Original image statistics before renormalization: min={minv} max={maxv} mean={meanv:.2f}+/-{stdv:.2f} median={medianv} ADU.')
            dx      = raw_max - raw_min
            dy      = max_adu - min_adu
            tmp     = rgb_im - raw_min
            rgb_im  = tmp * (dy/dx) + min_adu
            
            # recompute statistics.
            minv, maxv, meanv, stdv, medianv, pctiles = self._image_stats(rgb_im, 
                self._ipctls)

        if print_stats:
            logger.info(f'Image statistics: min={minv} max={maxv} mean={meanv:.2f}+/-{stdv:.2f} median={medianv} ADU.')
            logger.debug('Image percentiles follow:')
            for idx, pct in enumerate(self._ipctls):
                ovalue = pctiles[idx]
                logger.debug(f'  {pct:6.3f} percentile: {ovalue:8.3f} ADU')
        
        np.clip(rgb_im, min_adu, max_adu, out=rgb_im)
        logger.debug(f'RGB image of shape {rgb_im.shape} generated.')
        return rgb_im.astype(np.uint16), self._exif_dict

    def grey(self, luminance_method='linear', subtract_black=True, 
        wb_method='auto', print_stats=None, renorm=None):
        """
        Create a luminance image using the supplied white-balance
        values, with or without black-level subtraction.

        :param luminance_method: Luminance method to use. At present only the
           following method(s) are supported:
           - linear: Applies Bayer correction with a gamma=1 value to
             generate an RGB image, then applies the CCIR 601 coefficients
             to convert that to greyscale.
           - direct: Each channel is multiplied by its white-balance factor
             and the added to the final image. There is no interpolation.
        :param subtract_black: If true the camera black levels will be 
          subtracted from the channel data.
        :param wb_method: Whitebalance method to user to determine white 
          balances. See get_whitebalance() for details of the allowed
          methods.
        :param print_stats: Output image statistics to the logger if True.
        :param renorm: Renormalize the image counts to fill the full
          16-bit dynamic range if True. Specifically it will linearly
          stretch the dynamic range from the 0.01th to 99.99th percentiles
          to fill the range 0 to (2^16)-1 ADU.
          
        Returns a numpy array of the image along with a dictionary of EXIF
        tags.
        """
        
        min_adu = 0
        max_adu = 2**16 - 1
        
        if print_stats is None:
            print_stats = False
        if renorm is None:
            renorm = False
        
        allowed_methods=['linear', 'direct']
        if luminance_method not in allowed_methods:
            err_str = 'Unexpected luminance calculate method supplied to RawConv.grey: {}. Allowed methods are: {}'.format(luminance_method, allowed_methods)
            logger.error(err_str)
        else:
            logger.info(f'Generating greyscale image using {luminance_method} method.')
                
        if 'direct' in luminance_method:
            # In this case we do need the raw images to be built.
            r_im, g1_im, b_im, g2_im, exif_dict = self.split(subtract_black)
            
            if subtract_black:
                self._subtract_black_levels(subtract_black)
            
            wb_list = self.get_whitebalance(wb_method)
            
            # Perform calculations as double
            grey_im = np.zeros(self._rawim_r.shape[0:2], dtype=np.float64)
            grey_im += wb_list[self.R]  * r_im
            grey_im += wb_list[self.G1] * g1_im
            grey_im += wb_list[self.B]  * b_im
            grey_im += wb_list[self.G2] * g2_im
            
        elif 'linear' in luminance_method:
            rgb_coeff         = [0.299, 0.587, 0.114] # CCIR 601
            rgb_im, exif_dict = self.rgb(luminance_method, subtract_black,
                wb_method, False, False)
                
            grey_im = np.zeros(rgb_im.shape[0:2], dtype=np.float64)
            for idx in range(3):
                grey_im += rgb_im[:,:,idx] * rgb_coeff[idx]
        
        # Calculate initial image statistics.
        minv, maxv, meanv, stdv, medianv, pctiles = self._image_stats(rgb_im, 
            self._ipctls)

        # TODO Renormalization?
        if renorm:
            # Use 0.01 and 99.99 percentiles            
            raw_min = pctiles[0]
            raw_max = pctiles[ len(pctiles)-1 ]
            logger.info(f'Renormalizing image from range [{raw_min:.2f}, {raw_max:.2f}] to [{min_adu}, {max_adu}] ADU.')
            logger.debug(f'Original image statistics before renormalization: min={minv} max={maxv} mean={meanv:.2f}+/-{stdv:.2f} median={medianv} ADU.')
            dx      = raw_max - raw_min
            dy      = max_adu - min_adu
            tmp     = grey_im - raw_min
            grey_im = tmp * (dy/dx) + min_adu
            
            # recompute statistics.
            minv, maxv, meanv, stdv, medianv, pctiles = self._image_stats(grey_im, 
                self._ipctls)

        if print_stats:
            logger.info(f'Image statistics: min={minv} max={maxv} mean={meanv:.2f}+/-{stdv:.2f} median={medianv} ADU.')
            logger.debug('Image percentiles follow:')
            for idx, pct in enumerate(self._ipctls):
                ovalue = pctiles[idx]
                logger.debug(f'  {pct:6.3f} percentile: {ovalue:8.3f} ADU')
        
        np.clip(grey_im, min_adu, max_adu, out=grey_im)
        logger.debug(f'Greyscale image of shape {grey_im.shape} generated.')
        return grey_im.astype(np.uint16), self._exif_dict

    def split(self, subtract_black=True):
        """Exports the raw, unprocessed, bayer RGBG as four separate 
        uint16 numpy arrays.
        
        For each band (R, G etc) the only non-zero pixels will be that 
        band, and pixels associated with other bands will be set to zero.
        
        :param subtract_black: If true the camera black levels will be 
          subtracted from the channel data.
        
        Returns four numpy arrasy of the image, along with a dictionary of EXIF
        tags.
        """
        
        logger.debug('split: black_levels {}'.format( self._black_levels ))
        logger.debug('split: camera_wb    {}'.format( self._wb_camera ))
        logger.debug('split: daylight_wb  {}'.format( self._wb_daylight ))
    
        self._build_raw_channel_images()
    
        if subtract_black:
            self._subtract_black_levels()
           
        r_im  = np.where(self._mask_r,  self._rawim_r,  0)
        g1_im = np.where(self._mask_g1, self._rawim_g1, 0)
        b_im  = np.where(self._mask_b,  self._rawim_b,  0)
        g2_im = np.where(self._mask_g2, self._rawim_g2, 0)
        
        return r_im, g1_im, b_im, g2_im, self._exif_dict