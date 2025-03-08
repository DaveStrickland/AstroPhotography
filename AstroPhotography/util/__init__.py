"""
AstroPhotography Utilities __init__ file

Warnings
--------
Functions within the ``utils`` package should try to avoid dependencies on
the ``core`` package, in particular **not** importing from ``.core``. Any
``core`` package that is absolutely needed should be passed in via dependency
injection.
"""

# This seems overly verbose, but it works.
from .check import (
    does_file_exist,
    _does_file_exist,
    CapitalCase_to_snake_case,
    determine_file_type,
)
from .datautils import (
    regionprops_to_astropy_table,
    img_stats,
)
from .fitsutils import (
    load_wcs_from_file,
    summarize_wcs,
    make_dothead_from_keywords,
    make_dothead_from_file,
    make_dothead_from_wcs,
    get_exposure_time,
)
from .names import namefn_calibrated_input, namefn_getdir, _get_name_conv_dict
from .plotting import (
    load_image_and_plot,
    load_imlist_and_plot,
    plot_lupton_threecolor,
    plot_stiff_threecolor,
    make_lupton_threecolor_plots,
    make_stiff_threecolor_plots,
)
from .read import read_fits, read_yaml_into_metadatadict
from .write import file_writer, update_fits_header_with_exif

__all__ = [
    "CapitalCase_to_snake_case",
    "determine_file_type",
    "_does_file_exist",
    "does_file_exist",
    "file_writer",
    "get_exposure_time",
    "_get_name_conv_dict",
    "img_stats",
    "load_image_and_plot",
    "load_imlist_and_plot",
    "load_wcs_from_file",
    "make_dothead_from_file",
    "make_dothead_from_keywords",
    "make_dothead_from_wcs",
    "make_lupton_threecolor_plots",
    "make_stiff_threecolor_plots",
    "namefn_calibrated_input",
    "namefn_getdir",
    "plot_lupton_threecolor",
    "plot_stiff_threecolor",
    "read_fits",
    "read_yaml_into_metadatadict",
    "regionprops_to_astropy_table",
    "summarize_wcs",
    "update_fits_header_with_exif",
]
