"""
AstroPhotography Core Python Interface Classes
----------------------------------------------

The `Astrophotography.core` module contains classes for FITS
image reduction and processing.
"""

# This seems overly verbose, but it works.
from .RawConv import RawConv as RawConv
from .ApFixBadPixels import ApFixBadPixels as ApFixBadPixels
from .ApFindBadPixels import ApFindBadPixels as ApFindBadPixels
from .ApAddMetadata import ApAddMetadata as ApAddMetadata
from .ApAutoBadcols import ApAutoBadcols as ApAutoBadcols
from .ApCalibrate import ApCalibrate as ApCalibrate
from .ApFixCosmicRays import ApFixCosmicRays as ApFixCosmicRays
from .ApFindStars import ApFindStars as ApFindStars
from .ApMeasureStars import ApMeasureStars as ApMeasureStars
from .ApQualitySummarizer import ApQualitySummarizer as ApQualitySummarizer
from .ApAstrometry import ApAstrometry as ApAstrometry
from .ApImArith import ApImArith as ApImArith
from .ApMeasureBackground import ApMeasureBackground as ApMeasureBackground 
from .ApProcess import ApProcess as ApProcess

__all__ = ["RawConv", 
    "ApCalibrate",
    "ApFindBadPixels", 
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
    "ApProcess"
    ]
