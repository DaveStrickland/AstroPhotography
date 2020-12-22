# Processing images from iTelescope.net

## Processing Pipelines

The following table outlines conceptual pipeline processing stages
for astronomical FITS imagery, primarily those obtained using iTelescope.

| Stage Name          | Description                     | Simple Proc | Package Used      | Implemented?           |
| ------------------- | ------------------------------- | ----------- | ----------------- | ---------------------- |
| Build Master Cal    | Build master calibration files  | N/A         | astropy.ccdproc   | ap_combine_darks.py    |
|                     |                                 |             |                   | ap_calc_read_noise.py  |
| Find Bad Pixels     | Generate a bad pixel map        | Yes         |                   | ap_find_badpix.py      |
| Apply Master Cal    | Apply bias/dark/flat correction | N/A         |                   | ap_calibrate.py        |
| Add Metadata        | Add metadata to calibrated FITS | N/A         | astroplan         | ap_add_metadata.py     |
| Fix Bad Pixels      | Correct bad pixels              | Yes         |                   | ap_fix_badpix.py       |
| Quality Map         | Build data quality map          | N/A         | ?                 | No, not needed?        |
| Cosmic Ray Reject   | Find and correct cosmic rays    | Yes         | astropy.ccdproc   | Not yet                |
| Find Stars And PSF  | Find stars, measure PSF         | Yes         | astropy.photutils | ap_find_stars.py       |
| Astrometry          | Astrometric solution per image  | Yes         | astrometry.net    | ap_astrometry.py       |
| Quality Filter      | Identify "poor quality" images  | Yes         |                   | ap_find_stars.py       |
|                     | and summarize set of images     | Yes         |                   | ap_quality_summary.py  |
| Combine By Filter   | Resample images                 | Yes         | astromatic swarp? | Not yet                |
| Continuum Subtract  | Continuum scaling/subtraction   | Yes         | ?                 | Not yet                |
| RGB Color Composite | Make 3-color composites         | Yes         | astromatic stiff? | Not yet                | 
