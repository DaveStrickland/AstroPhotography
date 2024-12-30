"""
AstroPhotography Utilities
"""

# This seems overly verbose, but it works.
from .ApUtil import (
    regionprops_to_astropy_table,
    img_stats,
    does_file_exist,
    load_image_and_plot,
    load_imlist_and_plot,
    plot_lupton_threecolor,
    make_lupton_threecolor_plots,
    namefn_calibrated_input,
    namefn_getdir,
    load_wcs_from_file,
    summarize_wcs,
    get_exposure_time,
    make_dothead_from_file,
)
from .file_writer import (
    file_writer,
    determine_file_type,
    CapitalCase_to_snake_case,
    update_fits_header_with_exif,
    read_yaml_into_metadatadict,
)

__all__ = [
    "regionprops_to_astropy_table",
    "img_stats",
    "does_file_exist",
    "load_image_and_plot",
    "load_imlist_and_plot",
    "plot_lupton_threecolor",
    "make_lupton_threecolor_plots",
    "namefn_calibrated_input",
    "namefn_getdir",
    "load_wcs_from_file",
    "summarize_wcs",
    "get_exposure_time",
    "make_dothead_from_file",
    "file_writer",
    "determine_file_type",
    "CapitalCase_to_snake_case",
    "update_fits_header_with_exif",
    "read_yaml_into_metadatadict",
]
