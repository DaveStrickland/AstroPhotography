#!/usr/bin/env python3
#
#  ap_process.py
#
#  Copyright 2024 Dave Strickland <dave.strickland@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import sys
import argparse
from astropy.table import Table
import AstroPhotography as ap  # noqa: N813


def command_line_opts(argv):
    """Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(
        prog="ap_process",
        description=(
            "Performs FITS image astrometry and stacking for a set of calibrated"
            " FITS images. Optional preprocessing to update FITS header"
            " metadata, bad-pixel and cosmic ray correction is also supported."
        ),
    )

    # Required
    parser.add_argument(
        "data_dir",
        metavar="DATA_DIR",
        help=("Base-level data directory associated with the input" " data files, e.g. './'"),
    )

    # Optional default parameters
    defaults = {
        "loglevel": "INFO",
        "include_pattern": None,
        "exclude_pattern": None,
        "input_file_list": None,
        "input_rootname": None,
        "input_suffix": ".fits",
        "final_quality_file": None,
        "clean_star_detection": False,
        "clean_astrometry": False,
        "stop_on_error": False,
        "extnum": 0,
        "search_fwhm": 3.0,
        "search_nsigma": 7.0,
        "detector_bitdepth": 16,
        "max_sources": 200,
        "nosatmask": True,
        "sat_frac": 0.8,
        "quiet": True,
        "srclist_extname": "AP_XYPOS",
        "astnet_key": None,
        "use_sip": False,
        "user_scale": None,
        "scale_err_ratio": None,
        "target_wcs_file": None,
        "resampled_file_prefix": "resampled_",
        "resampled_file_suffix": "_resamp_weighted.fits",
        "resampled_dir": "./",
        "filter_list": None,
        "resampled_summary_plot": None,
        "composite_summary_plot": None,
        "preprocess_replace": None,
        "preprocess_with": None,
        "find_exposure_time": False,
        "metadata_yaml": None,
        "replace_keywords": False,
        "badpixelfile": None,
        "deltapix": 2,
        "fix_cosmic_rays": False,
    }

    def_value = defaults["loglevel"]
    parser.add_argument(
        "--loglevel",
        default=defaults["loglevel"],
        help=(
            "Logging message level. 'DEBUG' is recommended when processing"
            " a new set of images for the first time."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["include_pattern"]
    parser.add_argument(
        "--include_pattern",
        metavar="INCLUDE_PATTERN",
        default=def_value,
        help=(
            "Globbing pattern for input files we want included, specified relative"
            " to data_dir. If not specified all FITS files will be included. This"
            " parameter is ignored if file_list is not None."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["exclude_pattern"]
    parser.add_argument(
        "--exclude_pattern",
        metavar="EXCLUDE_PATTERN",
        default=def_value,
        help=(
            "Globbing pattern of files we want excluded, specified relative"
            " to data_dir. If not specified no FITS files will be excluded."
            " This parameter is ignored if file_list is not None."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["input_file_list"]
    parser.add_argument(
        "--input_file_list",
        metavar="QUOTED_FILE_LIST",
        default=def_value,
        help=(
            "An explicit space-separated list of files encapsulted in double"
            " quotes, paths relative to data_dir, may be specified. If specified"
            " then only those files in input_file_list are looked for, and the"
            " include and exclude patterns are ignored."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["input_rootname"]
    parser.add_argument(
        "--input_rootname",
        default=def_value,
        help=(
            "The part of the input files names that are shared and that"
            " designate them being the calibrated files. If not specified it"
            " is assumed that the files are from iTelescope or were created"
            " by the AstroPhotography module itself, and will have input_rootnames"
            " of either 'Calibrated-iTelescope', 'calibrated', or just 'cal' ."
            " The output files from this function replace the rootname with a"
            " file-type specific prefix, as described in the ApProcess documentation"
            " and the ApProcess.get_file_names function documentation."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["input_suffix"]
    parser.add_argument(
        "--input_suffix",
        default=def_value,
        help=(
            "String denoting the file type suffix of the input image file. For"
            " example, '.fits' or '.fits' or '.ftz' or 'fits.gz' or '.fits.bz2'"
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["final_quality_file"]
    parser.add_argument(
        "--final_quality_file",
        default=def_value,
        help=(
            "If specified, this is the name for a quality summary file in CSV"
            " format generated using ApQualitySummarizer, which summarizes all"
            " YaML quality files found in the standard quality summary file"
            " directory (typically ./MetaData/)."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["clean_star_detection"]
    parser.add_argument(
        "--clean_star_detection",
        default=False,
        action="store_true",
        help=(
            "If specified then existing star detection outputs (srclist, regile,"
            " plotfile, qualfile, and fwhmplot files) will be deleted and regenerated."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["clean_astrometry"]
    parser.add_argument(
        "--clean_astrometry",
        default=False,
        action="store_true",
        help=(
            " If specified then existing astrometry outputs ('navfile' files)"
            " will be deleted and regenerated."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["stop_on_error"]
    parser.add_argument(
        "--stop_on_error",
        default=False,
        action="store_true",
        help=(
            "By default ApProcess will continue attempting to process all the"
            " input files even if a given stage fails on one of the inputs."
            " That failure will be recorded in the output nav_status_table."
            " If you want it to instead stop processing when a failure is"
            " detected, specify stop_on_error."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["extnum"]
    parser.add_argument(
        "--extnum",
        default=def_value,
        help=(
            "Extension number or name for the extension holding the image"
            " data. Usually this is 0, for the PrimaryHDU. "
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["search_fwhm"]
    parser.add_argument(
        "--search_fwhm",
        default=def_value,
        type=float,
        help=(
            "Initial guess or estimate of the stellar PSF FWHM in pixels in this image."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["search_nsigma"]
    parser.add_argument(
        "--search_nsigma",
        default=def_value,
        type=float,
        help=(
            "Minimumn number of sigma above background for a detection."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["detector_bitdepth"]
    parser.add_argument(
        "--detector_bitdepth",
        default=def_value,
        type=int,
        help=(
            "Detector bit-depth, used in estimating which pixels are saturated."
            " (16 for most CCDs, ?? for CMOS, type dependent for cameras)."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["max_sources"]
    parser.add_argument(
        "--max_sources",
        default=def_value,
        help=(
            "Maximum number of sources to output or None (to output all sources)."
            " If not None then only the max_sources brightest sources will be"
            " output. Typically there is no advantage to having very large"
            " numbers of detected sources when performing astrometry, and may"
            " well slow it down. I normally use max_sources=200."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["nosatmask"]
    parser.add_argument(
        "--nosatmask",
        default=True,
        action="store_false",
        help=(
            " If specified then possibly saturated stars will be removed from"
            " the detected source lists and PSF-fitting. Mildly saturated stars"
            " typically do not compromise astrometry, and have little effect on the"
            " averaged PSF fit results, so leaving them in is usually OK."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["sat_frac"]
    parser.add_argument(
        "--sat_frac",
        default=def_value,
        type=float,
        help=(
            "Fraction of full well at which we assume star saturated."
            " (The detector bit depth is used to determine the number of"
            " electrons associated with full well.)"
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["quiet"]
    parser.add_argument(
        "--quiet",
        default=False,
        action="store_true",
        help=(
            "If specified this suppresses the runtime source list printing to STDOUT."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["srclist_extname"]
    parser.add_argument(
        "--srclist_extname",
        default=def_value,
        type=str,
        help=("FITS extension name for star X,Y position data." f" Default value: {def_value}"),
    )

    def_value = defaults["astnet_key"]
    parser.add_argument(
        "--astnet_key",
        default=def_value,
        type=str,
        help=(
            "Your personal Astrometry.net API key, if you have not already"
            " added it to your ~/.astropy/config/astroquery.cfg config file."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["use_sip"]
    parser.add_argument(
        "--use_sip",
        default=False,
        action="store_true",
        help=(
            "Allow astrometry.net to fit SIP polynomial distortion terms."
            " This may be necessary for very large fields of view (>10 deg),"
            " but SIP is not treated correctly by swarp (and possibly other software)."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["user_scale"]
    parser.add_argument(
        "--user_scale",
        default=def_value,
        help=(
            " If specified, override the estimate plate scale in the source"
            " list file (if any) and instead use a user spacified estimate of"
            " the plate scale. The units are arcseconds/pixel."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["scale_err_ratio"]
    parser.add_argument(
        "--scale_err_ratio",
        default=def_value,
        help=(
            "The relative uncertainty in the estimated plate scale, expressed"
            " as a ratio. This applies to either the default estimate from the"
            " source list, or a user-supplied plate scale. For example, if the"
            " estimated plate scale is 2.0 arcsec/pix and the scale_err_ratio=1.5,"
            " then the plate scale range that will be search by Astrometry.net"
            " is 2/1.5 (=4/3) to 2*1.5 (=3) arcseconds. If not specified ApAstrometry"
            " will use a value of 1.3. Using a larger value can help in cases where"
            " astrometric solutions fail, for example if incorrect telescope"
            " metadata leads to inaccurate estimated plate scales."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["target_wcs_file"]
    parser.add_argument(
        "--target_wcs_file",
        default=def_value,
        help=(
            "File name, including path, to the FITS file than contains the WCS"
            " that all navigated images should be resampled to match. If set"
            " to None then the first navigated file is used."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["resampled_file_prefix"]
    parser.add_argument(
        "--resampled_file_prefix",
        default=def_value,
        help=(
            "The part of the resampled image file name in front of the filter"
            " name. For example, for a Green filter output image name of"
            " M101_Supernova_2023ixf_Green_resamp_weighted.fits then"
            " resampled_file_prefix='M101_Supernova_2023ixf_'."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["resampled_file_suffix"]
    parser.add_argument(
        "--resampled_file_suffix",
        default=def_value,
        help=(
            "The part of the resampled image file name after the filter name."
            " Note that this should include the file suffix, which should include"
            " '.fits' (or further, e.g. '.fits.gz') For example, for a Green"
            " filter output image name of M101_Supernova_2023ixf_Green_resamp_weighted.fits"
            " then resampled_file_suffix='_resamp_weighted.fits'."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["resampled_dir"]
    parser.add_argument(
        "--resampled_dir",
        default=def_value,
        help=(
            " Directory in which resampled output images will be written. If"
            " None then the current directory is used. If not None then the"
            " named directory will be created if it does not exist."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["filter_list"]
    parser.add_argument(
        "--filter_list",
        metavar="QUOTED_FILTER_LIST",
        default=def_value,
        help=(
            "If specified then only images with FILTER header keywords matching"
            " one of the input list strings will be resampled. For example, to"
            " only resample the H-alpha and luminance images the filter list"
            " would be 'Lum Ha'. Note this must be a quoted space separated list."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["resampled_summary_plot"]
    parser.add_argument(
        "--resampled_summary_plot",
        default=def_value,
        help=(
            "If not None then a PNG file with the specified name containing"
            " plots of each output resampled image will be generated. The plots"
            " are generated using AstroPhotography.util.load_three_images_and_plot()."
            " Default intensity scaling and WCS plotting options will be used."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["composite_summary_plot"]
    parser.add_argument(
        "--composite_summary_plot",
        default=def_value,
        help=(
            " If there are resampled images in all three Red, Green, and Blue"
            " filters and/or all three Ha, SII, and OIII filters and this argument"
            " is not None then a PNG file with summary plots of three color"
            " composities will be generated. If Red, Green, and Blue filters are"
            " present then an RGB three color composite plot will be generated."
            " If the SII, H-alpha, and OIII filters are present then an SHO"
            " color-scheme composite plot will be generated. The plots are"
            " generated using AstroPhotography.util.plot_lupton_threecolor()."
            " Default intensity scaling and WCS plotting options will be used. "
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["preprocess_replace"]
    parser.add_argument(
        "--preprocess_replace",
        default=def_value,
        help=(
            "Substring common to all the input calibrated files that will be"
            " replaced with preprocess_with. For example 'calibrated' or"
            "'Calibrated' would be typical values when using iTelescope files."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["preprocess_with"]
    parser.add_argument(
        "--preprocess_with",
        default=def_value,
        help=(
            "String that will replace preprocess_replace in all input file"
            " names. For example if preprocess_replace='calibrated' I might"
            " use preprocess_with='recalibrated'."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["find_exposure_time"]
    parser.add_argument(
        "--find_exposure_time",
        default=False,
        action="store_true",
        help=(
            "If True then ApUtil.get_exposure_time() will be used to extract"
            " the exposure time from the input file headers and set the"
            " EXPOSURE keyword if it is not already present."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["metadata_yaml"]
    parser.add_argument(
        "--metadata_yaml",
        default=def_value,
        help=(
            "The name of an optional YaML file containing FITS header metadata"
            " values to be added to the preprocessed output files. The data in"
            " this file will be used as the keyword_dict parameter used by"
            " ApProcess.process_all. See the AstroPhototography.util.read_yaml_into_metadatadict"
            " documentation on the format of the YaML files."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["replace_keywords"]
    parser.add_argument(
        "--replace_keywords",
        default=False,
        action="store_true",
        help=(
            " If specified then existing FITS header keys that match the keys"
            " in the parsed metadata_yaml file will be over-written with the new values."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["badpixelfile"]
    parser.add_argument(
        "--badpixelfile",
        default=def_value,
        help=(
            " File path and name to the bad pixel file to apply. This file"
            " should conform to the format generated by ApFindBadPixels"
            " and used by ApFixBadPixels."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["deltapix"]
    parser.add_argument(
        "--deltapix",
        type=int,
        default=def_value,
        help=(
            "Linear distance away from a bad pixel from which the median value"
            " of the good pixels will be drawn. If 1 then the median value of"
            " good pixels within the surrounding 8 pixels will be used. If 2"
            " then the median of the good pixels within the surrounding 24"
            " pixels will be used. Values above 2 are not recommended."
            f" Default value: {def_value}"
        ),
    )

    def_value = defaults["fix_cosmic_rays"]
    parser.add_argument(
        "--fix_cosmic_rays",
        default=False,
        action="store_true",
        help=(
            "If specified then perform Cosmic Ray rejection on the images."
            f" Default value: {def_value}"
        ),
    )

    def_value = True
    parser.add_argument(
        "--no_summary_info",
        default=def_value,
        action="store_false",
        help=(
            "If specified then do not print the summary tables and lists output by"
            " ApProcess.process_all to STDOUT at the end of the run."
            f" Default value: {def_value}"
        ),
    )

    args = parser.parse_args(argv)

    # Correct any problematic values where the user may have input the string
    # None or its variant, instead of omitting the option.
    for key in defaults.keys():
        if hasattr(args, key):
            if "none" in str(getattr(args, key)).lower():
                setattr(args, key, None)

    # parse input file list into a list, should be a single space separated string
    if args.input_file_list is not None:
        args.input_file_list = (
            args.input_file_list.split()
        )  # handles multiple white space on its own

    # parse input filter list into a list, should be single space separated string
    if args.filter_list is not None:
        args.filter_list = args.filter_list.split()

    return args


def main(args):
    p_args = command_line_opts(args)

    # read yaml
    meta_dict = None
    if p_args.metadata_yaml is not None:
        meta_dict = ap.util.read_yaml_into_metadatadict(p_args.metadata_yaml)

    processor = ap.ApProcess(p_args.loglevel)

    quality_summary_file = p_args.final_quality_file

    nav_stat_tbl, stacked_imlist, stacked_im_tbl, preproc_flist = processor.process_all(
        p_args.data_dir,
        include_pattern=p_args.include_pattern,
        exclude_pattern=p_args.exclude_pattern,
        input_file_list=p_args.input_file_list,
        input_rootname=p_args.input_rootname,
        input_suffix=p_args.input_suffix,
        final_quality_file=quality_summary_file,
        clean_star_detection=p_args.clean_star_detection,
        clean_astrometry=p_args.clean_astrometry,
        stop_on_error=p_args.stop_on_error,
        extnum=p_args.extnum,
        search_fwhm=p_args.search_fwhm,
        search_nsigma=p_args.search_nsigma,
        detector_bitdepth=p_args.detector_bitdepth,
        max_sources=p_args.max_sources,
        nosatmask=p_args.nosatmask,
        sat_frac=p_args.sat_frac,
        quiet=p_args.quiet,
        srclist_extname=p_args.srclist_extname,
        astnet_key=p_args.astnet_key,
        use_sip=p_args.use_sip,
        user_scale=p_args.user_scale,
        scale_err_ratio=p_args.scale_err_ratio,
        target_wcs_file=p_args.target_wcs_file,
        resampled_file_prefix=p_args.resampled_file_prefix,
        resampled_file_suffix=p_args.resampled_file_suffix,
        resampled_dir=p_args.resampled_dir,
        filter_list=p_args.filter_list,
        resampled_summary_plot=p_args.resampled_summary_plot,
        composite_summary_plot=p_args.composite_summary_plot,
        preprocess_replace=p_args.preprocess_replace,
        preprocess_with=p_args.preprocess_with,
        find_exposure_time=p_args.find_exposure_time,
        keyword_dict=meta_dict,
        replace_keywords=p_args.replace_keywords,
        badpixelfile=p_args.badpixelfile,
        deltapix=p_args.deltapix,
        fix_cosmic_rays=p_args.fix_cosmic_rays,
    )

    # Output files generated by preprocessing
    if preproc_flist is not None:
        print("--- Files generated by preprocessing: ---")
        for im in preproc_flist:
            print(f"  {im}")
    else:
        print("--- No preprocessing applied. ---")

    # show status table
    print("--- Navigation status table ---")
    nav_stat_tbl.pprint(max_width=200)

    # # Load and display the quality summary table (NOTE this is a very wide table)
    if quality_summary_file is not None:
        print("--- Source detection image quality summary table ---")
        try:
            qual_summ_table = Table.read(quality_summary_file, format="ascii.csv")
            qual_summ_table[
                "file",
                "ncols",
                "nrows",
                "filter",
                "median",
                "stddev",
                "num_detected",
                "fwhm_val_pix",
                "fwhm_err_pix",
                "fwhm_val_arcs",
                "fwhm_err_arcs",
                "num_data_pts",
            ].pprint_all()
        except Exception as e:
            print(f"Quality summary table was not generated, but was expected. Error={e}")
            # continue as we're just reporting things at this stage.
    else:
        print("Quality summary table was not requested.")

    # Stacked images etc
    print("--- List of stacked/resampled images: ---")
    for im in stacked_imlist.items():
        print(f"  {im}")
    print("--- Resampled Image Info Table ---")
    stacked_im_tbl.pprint_all()

    if p_args.resampled_summary_plot is not None:
        print(f"Summary plot of resampled images: {p_args.resampled_summary_plot}")

    if p_args.composite_summary_plot is not None:
        print(f"Summary plot of 3-color composite images: {p_args.composite_summary_plot}")

    return


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
