# AstroPhotography TODO

## General Tool Features

- [X] Package layout
- [X] Work out how to add a command line executable.
- [ ] Clean up or remove the baselevel __main__.py.
- [ ] Check it works with DNG format files.
- [ ] Write basic FITS files (Correct metadata is separate item)

## RAW Conversion Sub Tools

- [X] whitebalance: options to select camera, daylight, auto, region, user-specified
- [X] whitebalance: Implement auto
- [X] whitebalance: Implement region
- [ ] whitebalance: Implement user
- [X] spit: Option to output unprocessed RGBG as separate 16-bit greyscale TIFFs
- [ ] grey, rgb: Output options for FITS, TIFF (?)
- [ ] hist: Calculate histograms (rgb or greyscale).
- [ ] hist: Plot histograms (matplotlib)
- [ ] hist: Write histograms to CSV or other ascii format.
- [ ] general: Percentiles, used for preparing for image segmentation.
- [ ] general: Image stats output in verbose mode.
- [ ] features: Image segmentation and feature finding.
- [ ] features: Indexed FITS or PNG image of features
- [ ] features: CSV or other ascii file of feature locations.
- [X] general: black-level subtraction
- [ ] general: Camera metadata conversion to FITS kw or image metadata
- [ ] general: Optional flipping
- [X] general: Check rawfile exists and can be read using os.path
- [ ] general: Optional normalization to 16 bit dynamic range. 
- [ ] general: Optional asinh strech.
- [ ] general: Try dng file type.
- [ ] Switch TODO file to list task backlog by planned release.

## Metadata

- [ ] info: Extract and print image metadata to stdout.
- [ ] Investigate what constitutes valid FITS metadata.
- [ ] Investigate valid PNG metadata.
- [X] See what metadata dcraw writes when creating tiff with -T
- [ ] Plate scale based on camera?
- [X] Use imageio image class that does metadata. 
      This doesn't seem to read metadata for dcraw-created TIFF files.
- [ ] Investigate switching to pillow as imageio cannot handle metadata,
      png output is slow, and it doesn't support 16-bit jp2.

## Configuration

- [ ] Work out how to get config values through to commands.

## Testing

- [X] Work out how pytest works with the command line executables 
- [X] Get test fixtures working? DONE (within single file)
- [ ] Move basic test machinary into a separate file.
- [ ] Separate tests of individual commands into separate files.
- [X] Work how to do test coverage in python.
- [X] RawConv.split() unit test.
- [ ] RawConv.get_whitebalance() unit test. IN PROGRESS
- [ ] RawConv.grey() unit test.

## Bugs

- [X] Canon CR2 has minimum values below the black levels causing uint16 wrapping.
