"""
Function based interface to FileWriter
"""

from ..core.logger import logger
import imageio
import time
import os.path
import yaml
from pathlib import Path
from datetime import datetime, timezone
from astropy.io import fits
from typing import Any
import numpy as np
import numpy.typing as npt

# AstroPhotography includes
from .. import __version__


def _does_file_exist(filename: str, verbose: bool = False) -> bool:
    """
    Returns True if the file name or path exists, false otherwise

    Parameters
    ----------
    filename : str
        Name (optionally including path) of the file we want to check
        the existence of.
    verbose : bool, optional, default=False
        If True then writes to stdout.

    Returns
    -------
    exists : bool
        Returns True if the file path exists, False otherwise.
    """
    if verbose and not Path(filename).exists():
        print(f"Cannot find {filename}. Not a valid path or file.")
        return False
    else:
        print(f"Found {filename}.")
    return True


def read_fits(
    image_filename: str, image_extension: str | int, a_logger: Any, remove_pedestal: str = "always"
) -> tuple[npt.NDArray, Any]:
    """
    Read a single extension's data and header from a FITS file

    Parameters
    ----------
    image_filename : str
        Name/path of the FITS file to be opended.
    image_extension : str or int
        FITS extension number or name from which to read the header
        and data.
    a_logger : logger instance
        A logger instance must be provided.
    remove_pedestal : {'always', 'positive', 'never'}
        FITS ``PEDESTAL`` keyword handling. If the specified file/extension
        has the ``PEDESTAL`` keyword then ``always`` will remove the pedestal
        from the data using the value of the keyword, ``never`` will not
        modify the raw data array values irrespective of the value of the
        pedestal keyword, and ``positive`` will only remove the pedestal
        if doing so would leave a positive median pixel value. The latter
        option can be necessary in cases where other tools have modified
        the data but **not** removed the ``PEDSTAL`` keyword or set its
        value to zero.

    Returns
    -------
    ext_data : ndarray
        The data array associated with the specified file/extension
    ext_hdr : astropy.io.fits.Header
        The FITS header associated with the specified file/extension
    """

    if not _does_file_exist(image_filename, False):
        err_msg = ""
        raise RuntimeError(err_msg)

    a_logger.info("Loading extension {} of FITS file {}".format(image_extension, image_filename))

    # open() parameters that can be important.
    # Default values used here.
    # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
    uint_handling = True
    image_scaling = False

    with fits.open(
        image_filename, uint=uint_handling, do_not_scale_image_data=image_scaling
    ) as hdu_list:
        ext_hdr = hdu_list[image_extension].header
        ext_data = hdu_list[image_extension].data

    ndim = ext_hdr["NAXIS"]
    cols = ext_hdr["NAXIS1"]
    rows = ext_hdr["NAXIS2"]
    bitpix = ext_hdr["BITPIX"]
    info_str = f"{ndim}-D BITPIX={bitpix} image with {cols} columns, {rows} rows"

    if ndim == 3:
        layers = ext_hdr["NAXIS3"]
        info_str += f", {layers} layers"

    if "BSCALE" in ext_hdr:
        bscale = ext_hdr["BSCALE"]
        info_str += f", BSCALE={bscale}"

    if "BZERO" in ext_hdr:
        bzero = ext_hdr["BZERO"]
        info_str += f", BZERO={bzero}"

    a_logger.debug(info_str)

    if ndim == 3:
        err_msg = "Error, 3-D handling has not been implemented yet."
        a_logger.error(err_msg)
        raise RuntimeError(err_msg)

    # Get data absolute limits.
    minval = np.amin(ext_data)
    maxval = np.amax(ext_data)
    medval = np.median(ext_data)
    a_logger.debug(
        f"Raw data statistics are min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}"
    )

    # Is there a PEDESTAL value? MaximDL likes to add an offset, and
    # the PEDESTAL value is the value to ADD to the data to remove the
    # pedestal.
    if "PEDESTAL" in ext_hdr:
        if remove_pedestal in ["always", "positive"]:
            pedestal = float(ext_hdr["PEDESTAL"])
            condition_met = True
            if pedestal == 0:
                condition_met = False
            elif "positive" in remove_pedestal:
                # values would be on average negative if pedestal corrected for
                bg_negative = medval < pedestal
                if bg_negative:
                    condition_met = False
            if condition_met:
                a_logger.debug(f"Removing a PEDESTAL value of {pedestal} ADU.")
                ext_data += pedestal
                minval = np.amin(ext_data)
                maxval = np.amax(ext_data)
                medval = np.median(ext_data)
                a_logger.debug(
                    (
                        "After PEDESTAL removal, "
                        f"min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}"
                    )
                )

    return ext_data, ext_hdr


def read_yaml_into_metadatadict(yamlfile: str) -> dict[str, tuple[Any, str]]:
    """
    Read a YaML file of FITS header metadata and converts that into a format
    that can easily be used with FITS files

    It is convenient to interact with FITS header records as a dictionary
    of keys (str) with assoiated tuples of value and comment,
    ``key_str: (value, comment_str)``. When output to YaML the tuples become
    lists, but astropy FITS I/O does not accept lists as the value of keyword
    keys. This function reads a YaML file (assumed to be in the correct format)
    and converts the lists back into tuples.

    For example, this metadata dictionary::

        meta_dict = {'OBJNAME':  ('M101 Supernova 2023ixf', 'Target'),
             'OBJECT':   ('M101 Supernova 2023ixf', 'Target'),
             'TELESCOP': ('iTelescope 24', 'Telescope name'),
             'INSTRUME': ('FLI-PL09000', 'Instrument name'),
             'XPIXSZ':   (12.0, 'Pixel Width in microns (after binning)'),
             'YPIXSZ':   (12.0, 'Pixel Height in microns (after binning)'),
             'EGAIN':    (1.53, 'Electronic gain in e-/ADU')}

    would appear as the following YaML::

        ---
        EGAIN:
        - 1.53
        - Electronic gain in e-/ADU
        INSTRUME:
        - FLI-PL09000
        - Instrument name
        OBJECT: &id001
        - M101 Supernova 2023ixf
        - Target
        OBJNAME: *id001
        TELESCOP:
        - iTelescope 24
        - Telescope name
        XPIXSZ:
        - 12.0
        - Pixel Width in microns (after binning)
        YPIXSZ:
        - 12.0
        - Pixel Height in microns (after binning)
        ...

    The normal ``yaml.safe_load`` converts that back into::

        # Reloaded data, back in python
        {'EGAIN': [1.53, 'Electronic gain in e-/ADU'],
            'INSTRUME': ['FLI-PL09000', 'Instrument name'],
            'OBJECT': ['M101 Supernova 2023ixf', 'Target'],
            'OBJNAME': ['M101 Supernova 2023ixf', 'Target'],
            'TELESCOP': ['iTelescope 24', 'Telescope name'],
            'XPIXSZ': [12.0, 'Pixel Width in microns (after binning)'],
            'YPIXSZ': [12.0, 'Pixel Height in microns (after binning)']
        }

    while this function aims to convert it back into the original format.

    Parameters
    ----------
    yamlfile : str
        Name of the yaml file to read

    Returns
    -------
    metadata_dict : dict[str, tuple[Any, str]]
        Metadata dictionary, of the form expected by (for example)
        :func:`ApProcess.process_all`
    """

    raw_dict = {}
    with open(yamlfile, "r") as f2:
        raw_dict = yaml.safe_load(f2)

    metadata_dict = {key: tuple(val) for key, val in raw_dict.items()}
    return metadata_dict


def file_writer(out_file: str, data_array: Any, exif_dict: Any):
    """
    Abstracts away details of writing an image to a graphics
    file format or astronomical FITS file.

    Parameters
    ----------
    out_file: str
        File path/name for output file.
    data_array: ndarray
        2-D (monochrome) or 3-D (rgb) ndarray of
        data values. Usually but not always of dtype.uint16
    exif_dict: dict
        Dictionary of EXIF tags and values in the format
        produced by ExifRead. This will be processed to provide a format
        suitable for the exif package or astropy FITS to write.
    """

    ftype = determine_file_type(out_file)

    # Image properties.
    ndim = data_array.ndim
    nrows = data_array.shape[0]
    ncols = data_array.shape[1]
    nlayers = 1
    if ndim == 3:
        nlayers = data_array.shape[2]
    elif ndim == 1 or ndim > 3:
        err_str = f"Error: file_writer cannot handle {ndim}-dimension arrays."
        logger.error(err_str)
        return  # maybe throw?

    dtype = data_array.dtype

    if ndim == 2:
        info_str = f"Writing {nrows}x{ncols} {dtype} image to {out_file}"
    elif ndim == 3:
        info_str = f"Writing {nrows}x{ncols}x{nlayers} {ndim}-D {dtype} image to {out_file}"
    logger.info(info_str)

    # TODO make sure write is succesful.
    t_start = time.perf_counter()
    if ftype == "graphics":
        # Special handling of format?
        our_format = None
        if "jp2" in out_file:
            our_format = "JPEG2000-PIL"
        imageio.imsave(uri=out_file, im=data_array, format=our_format)
    elif ftype == "fits":
        if ndim == 2:
            hdu1 = fits.PrimaryHDU(data_array)
        elif ndim == 3:
            # Store RGB image as three separate 2-D fits images
            # Zero'th plane as first image
            hdu1 = fits.PrimaryHDU(data_array[:, :, 0])

        hdu_list = fits.HDUList([hdu1])
        hdr = update_fits_header_with_exif(hdu_list[0].header, exif_dict)
        hdu_list[0].header = hdr

        if ndim == 3:
            hdu_list[0].header["FILTER"] = ("Red", "Red channel of RGB image")

            # Green
            hdr2 = fits.Header(hdr, copy=True)
            hdr2["FILTER"] = ("Green", "Green channel of RGB image")
            hdu2 = fits.ImageHDU(data_array[:, :, 1], header=hdr2)
            hdu_list.append(hdu2)

            # Blue
            hdr3 = fits.Header(hdr, copy=True)
            hdr3["FILTER"] = ("Blue", "Blue channel of RGB image")
            hdu3 = fits.ImageHDU(data_array[:, :, 2], header=hdr3)
            hdu_list.append(hdu3)

            logger.debug("3-D data array written to FITS file as three 2-D image extensions.")

        hdu_list.writeto(out_file, output_verify="warn", overwrite=True)
    else:
        raise RuntimeError("Could not determine type for output file {}".format(out_file))
    t_end = time.perf_counter()

    io_size = os.path.getsize(out_file) / 1e6  # size in MB
    io_time = t_end - t_start  # seconds
    io_rate = io_size / io_time  # MB/s
    logger.debug(
        "Wrote {:.3f} MB file successfully in {:.3f} seconds ({:.3f} MB/s)".format(
            io_size, io_time, io_rate
        )
    )

    # This should be moved to a function with error/exception checking.
    if ftype == "graphics":
        logger.warning("Writing EXIF metadata not yet supported.")
    return


def update_fits_header_with_exif(hdr, exif_dict):
    """
    Update a FITS header with items from an EXIF dictionary, returning
    the modified header.

    This will attempt to add the following metadata in all cases:

    - DATE      (now)                   Date/time this file was generated
    - DATE-OBS  Image DateTime          Corresponds to date/time RAW was taken
    - INSTRUME  Image Model             Camera make and model
    - EXPOSURE  EXIF ExposureTime       Exposure time as rational fraction (s)?
    - EXPTIME                           Exposure as real number (s)
    - FNUMBER   EXIF FNumber            F number
    - ISONUM    EXIF ISOSpeedRatings    ISO number
    - FOCALLEN  EXIF FocalLength        Focal length (mm)

    :param hdr: Input FITS header
    :param exif_dict: EXIF dictionary in the format supplied by ExifRead
    """

    # file_writer is only called by RawConv
    name = "RawConv"

    # Metadata coming from the exif data.
    metadata_dict = {
        "DATE-OBS": {"exifkw": "Image DateTime", "cmt": "Date/time RAW image captured"},
        "INSTRUME": {"exifkw": "Image Model", "cmt": "Camera make and model"},
        "EXPOSURE": {
            "exifkw": "EXIF ExposureTime",
            "cmt": "[sec] Exposure time as rational number",
        },
        "FNUMBER": {"exifkw": "EXIF FNumber", "cmt": "f-number, F/D"},
        "ISONUM": {"exifkw": "EXIF ISOSpeedRatings", "cmt": "ISO number"},
        "FOCALLEN": {"exifkw": "EXIF FocalLength", "cmt": "[mm] Focal length F"},
    }

    logger.info("Updating FITS header with RAW file EXIF information.")
    new_hdr = hdr

    tnow = datetime.now().isoformat(timespec="milliseconds")
    new_hdr["HISTORY"] = f"Processed by {name} {__version__} at {tnow}"
    new_hdr["DATE"] = (tnow, "Date/time FITS file  generated.")

    exptime = None
    for kw in metadata_dict:
        val = None
        exif_kw = metadata_dict[kw]["exifkw"]
        if exif_kw in exif_dict:
            # Must convert out of Ifd format into a string before anything else.
            val = str(exif_dict[exif_kw])
            if "EXPOSURE" in kw:
                exptime = float(eval(val))  # convert string to float
            if "DATE-OBS" in kw:
                dtime = datetime.strptime(val, "%Y:%m:%d %H:%M:%S")
                val = dtime.isoformat()
            if kw in ["FNUMBER", "FOCALLEN"]:
                val = float(eval(val))
        new_hdr[kw] = (val, metadata_dict[kw]["cmt"])

    new_hdr["EXPTIME"] = (exptime, "[sec] Exposure time as a floating point number.")

    return new_hdr


def CapitalCase_to_snake_case(input_str: str) -> str:  # noqa: N802
    """
    Converts CapitalCase or camelCase strings to snake_case

    Based on: https://www.geeksforgeeks.org/python-program-to-convert-camel-case-string-to-snake-case/

    Limitations:

    - Messes up consecutive capital letters, e.g. ISO become i_s_o, or
      FocalLengthMM becomes focal_length_m_m

    :param input_str: Input CapitalCase or camelCase string
    """

    output_str = "".join(["_" + i.lower() if i.isupper() else i for i in input_str]).lstrip("_")

    return output_str


def determine_file_type(fname: str) -> str:
    """
    Determine if file is a common graphics format or FITS.
    """

    # Graphics file formats, delegated to imageio
    graphics_list = [".tif", ".tiff", ".jpg", ".jp2", "jpeg", ".png", ".gif"]

    # FITS, delegate to astropy
    fits_list = [".fits", ".ftz", ".fit", ".fits.gz"]

    # want file extension, e.g .png, .fits but also .fits.gz
    root, ext = os.path.splitext(fname.lower())
    if ".gz" in ext:
        root, ext2 = os.path.splitext(root)
        ext = ext2 + ext
    ext = ext.lower()

    if ext in fits_list:
        ftype = "fits"
    elif ext in graphics_list:
        ftype = "graphics"
    else:
        ftype = "unknown"

    return ftype
