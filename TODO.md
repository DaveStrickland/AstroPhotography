# AstroPhotography TODO

## General

- [X] Package layout
- [X] Work out how to add a command line executable.
- [ ] Work out how pytest works with the command line executables 
- [ ] Clean up or remove the baselevel __main__.py.

## RAW Conversion

- [ ] whitebalance: options to select camera, daylight, auto, region, user-specified
- [ ] spit: Option to output unprocessed RGBG as separate 16-bit greyscal TIFF
- [ ] grey, rgb: Output options for FITS, TIFF (?)
- [ ] hist: Calculate histograms (rgb or greyscale).
- [ ] hist: Plot histograms (matplotlib)
- [ ] hist: Write histograms to CSV or other ascii format.
- [ ] general: percentiles, used for preparing for image segmentation.
- [ ] features: Image segmentation and feature finding.
- [ ] features: Indexed FITS or PNG image of features
- [ ] features: CSV or other ascii file of feature locations.
- [ ] general: black-level subtraction
- [ ] general: Camera metadata conversion to FITS kw or image metadata
- [ ] general: Optional flipping
