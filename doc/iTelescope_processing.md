# Processing images from iTelescope.net

## Processing Pipelines

The following table outlines a conceptual pipeline for processing and 
combining multiple FITS images. It lists a simplified, less ideal,
process based on the calibrated FITS images supplied by iTelescope,
suitable for "good enough" final imagery. Finally it notes which
packages are used in the stage and which parts of the pipeline 
are implemented in AstroPhotography.

| Stage Name        | Description                     | Simple Proc | Package Used      | Implemented? |
| ----------------- | ------------------------------- | ----------- | ----------------- | ------------ |
| BuildMasterCal    | Build master calibration files  | N/A         | astropy.ccdproc   | Not yet      |
| ApplyMasterCal    | Apply bias/dark/flat correction | N/A         | astropy.ccdproc   | Not yet      |
| BadPix            | Find/correct cold/hot pixels    | N/A         | astropy.ccdproc   | Not yet      |
| QualityMap        | Build data quality map          | N/A         | ?                 | Not yet      |
| CosmicRayReject   | Find and correct cosmic rays    | Yes         | astropy.ccdproc   | Not yet      |
| FindStarsAndPSF   | Find stars, measure PSF         | Yes         | astropy.photutils | ap_star_find, partially |
| Astrometry        | Astrometric solution per image  | Yes         | astrometry.net    | by hand      |
| CombineByFilter   | Resample images                 | Yes         | astromatic swarp  | by hand      |
| ContinuumSubtract | Continuum scaling/subtraction   | Yes         | ?                 | Not yet      |
| RGBColorComposite | Make 3-color composites         | Yes         | astromatic stiff  | Not yet      | 
