# -*- coding: utf-8 -*-
#
#  Contains the implementation of various Astrophotography utilities
#  associated with FITS file handling and associated concepts, e.g.
#  WCS.
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

# 2024-04-07 dks : issue-002 Initial coding
# 2025-03-01 dks : Split off from ApUtil.py

##import logging
##import os.path
from typing import Any
from pathlib import Path
import numpy as np
import numpy.typing as npt
import skimage.measure as skm

from astropy.table import Table

# AstroPhotography includes
##from .. import __version__


def regionprops_to_astropy_table(
    label_image,
    intensity_image=None,
    properties=("label", "bbox"),
    *,
    cache=True,
    separator="-",
    extra_properties=None,
    spacing=None,
):
    """
    Compute image properties and return them as an AstroPy Table.

    This function is a wrapper around ``skimage.measure.regionprops_table`` that
    converts the output dictionary into an AstroPy Table.

    Parameters
    ----------
    label_image : (M, N[, P]) ndarray
        Labeled input image. Labels with value 0 are ignored.
    intensity_image : (M, N[, P][, C]) ndarray, optional
        Intensity (i.e., input) image with same size as labeled image, plus
        optionally an extra dimension for multichannel data. The channel dimension,
        if present, must be the last axis. Default is None.
    properties : tuple or list of str, optional
        Properties that will be included in the resulting dictionary
        For a list of available properties, please see :func:`regionprops`.
        Users should remember to add "label" to keep track of region
        identities.
    cache : bool, optional
        Determine whether to cache calculated properties. The computation is
        much faster for cached properties, whereas the memory consumption
        increases.
    separator : str, optional
        For non-scalar properties not listed in OBJECT_COLUMNS, each element
        will appear in its own column, with the index of that element separated
        from the property name by this separator. For example, the inertia
        tensor of a 2D region will appear in four columns:
        ``inertia_tensor-0-0``, ``inertia_tensor-0-1``, ``inertia_tensor-1-0``,
        and ``inertia_tensor-1-1`` (where the separator is ``-``).

        Object columns are those that cannot be split in this way because the
        number of columns would change depending on the object. For example,
        ``image`` and ``coords``.
    extra_properties : Iterable of callables
        Add extra property computation functions that are not included with
        skimage. The name of the property is derived from the function name,
        the dtype is inferred by calling the function on a small sample.
        If the name of an extra property clashes with the name of an existing
        property the extra property will not be visible and a UserWarning is
        issued. A property computation function must take a region mask as its
        first argument. If the property requires an intensity image, it must
        accept the intensity image as the second argument.
    spacing: tuple of float, shape (ndim,)
        The pixel spacing along each axis of the image.

    Returns
    -------
    out_table : astropy.Table
        An AstroPy containing one column for each dictionary item.

    Notes
    -----
    See skimage.measure.regionprops and skimage.measure.regionprops_table
    for more information.
    """

    if "label" not in properties:
        err_msg = 'Error, "label" must be one of the properties'
        raise RuntimeError(err_msg)

    props_tbl = skm.regionprops_table(
        label_image,
        intensity_image=intensity_image,
        properties=properties,
        cache=cache,
        separator=separator,
        extra_properties=extra_properties,
        spacing=spacing,
    )
    ##for key, val in props_tbl.items():
    ##    print(key, val)

    # properties currently handled. Note that non-scalar properties like the
    # bbox will be split into multiple keys and not caught by these lists.
    # They are handled separated below.
    int_props = ["label", "num_pixels"]
    float_props = [
        "area",
        "area_bbox",
        "area_convex",
        "area_filled",
        "axis_major_length",
        "axis_minor_length",
        "eccentricity",
        "intensity_max",
        "intensity_min",
        "intensity_mean",
        "intensity_std",
        "orientation",
    ]

    dtype_list = []
    for prop_key in props_tbl.keys():
        # determine format to use for the column based on the name of the data
        table_format = "f8"
        if prop_key in int_props:
            table_format = "i4"
        elif prop_key in float_props:
            table_format = "f8"
        elif "bbox" in prop_key:
            table_format = "i4"
        else:
            print(
                f"Warning: unexpected property {prop_key} encountered."
                " This will not be added to the Table."
            )
            continue

        dtype_list.append((prop_key, table_format))

    nrows = len(props_tbl["label"])
    ##print(f'dtypes={dtype_list}')
    ##print(f'number of rows = {nrows}')
    out_table = Table(data=np.zeros(nrows, dtype=dtype_list))

    for prop_key in props_tbl.keys():
        out_table[prop_key] = props_tbl[prop_key]

    rename_col_dict = {
        "bbox-0": "bbox_rmin",
        "bbox-1": "bbox_cmin",
        "bbox-2": "bbox_rmax",
        "bbox-3": "bbox_cmax",
    }
    for key, val in rename_col_dict.items():
        if key in out_table.colnames:
            out_table.rename_column(key, val)

    return out_table


def img_stats(data, label, verbose=False):
    """
    Calculate and optionally display some image statistics, returning
    a list of the minimum, maximum, mean, standard deviation, and
    median values.

    If verbose=True then the computed statistics are also written
    to the log at INFO level along with the specified informative
    label.

    Parameters
    ----------
    data: ndarray
        Data array to compute statistics of
    label : str
        Descriptive label, only used if verbose=True
    verbose : bool, optional, default=False
        If true the the computed statistics are also written
        to the log at INFO level

    Returns
    -------
    minval :float
        Minimum value in NaN-filtered input data
    maxval : float
        Maximum value in NaN-filtered input data
    meanval : float
        Mean value in NaN-filtered input data
    stdval : float
        Standard deviation in NaN-filtered input data
    medval : float
        Median value in NaN-filtered input data
    """

    minval = np.nanmin(data)
    maxval = np.nanmax(data)
    meanval = np.nanmean(data)
    stdval = np.nanstd(data)

    # percentiles, 50th percentile is the median
    #         0    1    2    3   4   5   6   7   8   9   10
    ipctls = [0.1, 1.0, 5.0, 10, 25, 50, 75, 90, 95, 99, 99.9]
    opctls = np.nanpercentile(data, ipctls)
    medval = opctls[5]

    if verbose:
        msg1: str = (
            f"{label} data min={minval:.2f}, max={maxval:.2f},"
            f" mean={meanval:.2f} +/- {stdval:.2f}, median={medval:.2f} ADU."
        )
        msg2: str = (
            f"  90% of data between {opctls[2]:.2f} and {opctls[8]:.2f} ADU (5-95 precentiles)"
        )
        msg3: str = (
            f"  98% of data between {opctls[1]:.2f} and {opctls[9]:.2f} ADU (1-99 precentiles)"
        )
        print(msg1)
        print(msg2)
        print(msg3)
    return [minval, maxval, meanval, stdval, medval]
