"""
AstroPhotography Core Python Interface Classes
----------------------------------------------

The `Astrophotography.core` module contains classes for FITS
image reduction and processing.
"""

# This seems overly verbose, but it works.
from .RawConv import RawConv as RawConv
from .ApAddMetadata import ApAddMetadata as ApAddMetadata
from .ApAstrometry import ApAstrometry as ApAstrometry
from .ApAutoBadcols import ApAutoBadcols as ApAutoBadcols
from .ApCalibrate import ApCalibrate as ApCalibrate
from .ApFindBadPixels import ApFindBadPixels as ApFindBadPixels
from .ApFindStars import ApFindStars as ApFindStars
from .ApFitsRasterizer import ApFitsRasterizer as ApFitsRasterizer
from .ApFixBadPixels import ApFixBadPixels as ApFixBadPixels
from .ApFixCosmicRays import ApFixCosmicRays as ApFixCosmicRays
from .ApFixHoles import ApFixHoles as ApFixHoles
from .ApImArith import ApImArith as ApImArith
from .ApMasterCal import ApMasterCal as ApMasterCal
from .ApMeasureBackground import ApMeasureBackground as ApMeasureBackground
from .ApMeasureStars import ApMeasureStars as ApMeasureStars
from .ApProcess import ApProcess as ApProcess
from .ApQualitySummarizer import ApQualitySummarizer as ApQualitySummarizer


__all__ = [
    "RawConv",
    "ApCalibrate",
    "ApFindBadPixels",
    "ApFitsRasterizer",
    "ApFixHoles",
    "ApFixBadPixels",
    "ApFixCosmicRays",
    "ApAddMetadata",
    "ApAutoBadcols",
    "ApFindStars",
    "ApMeasureStars",
    "ApQualitySummarizer",
    "ApAstrometry",
    "ApImArith",
    "ApMeasureBackground",
    "ApProcess",
    "ApMasterCal",
]
