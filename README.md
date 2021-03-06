# AstroPhotography

The aim of the AstroPhotography python package is to provides python 
classes and command line applications for amateur astronomy, specifically:

- quick inspection and conversion of RAW digital camera format to common 
  graphical image formats (e.g. PNG) and astronomical FITS format.
- reduction and combination of multiple FITS images, including calibration,
  artifact removal, star detection, and astrometry.

FITS images can then be viewed with the powerful features SAO `ds9`
provides, and/or within python using `astropy`.

*Current status:* This is still very much a work in progress. RAW to 
image/FITS conversion is partially implemented, but development stalled
when I became dissatified with the clunky nature of the unit tests. I've
been putting more effort into the FITS data reduction side of the project,
but still have a long way to go.

## Command Line Functionality

The command line python program `dksraw` will provide a simple method of
quickly converting RAW files into useful images or FITS files without 
the user having to mess around with `dcraw`, `gimp`, `photoshop` or 
the equivalent.

A series of python scripts, all beginning with the prefix `ap_`, will 
perform separate stages of traditional astronomical image data reduction
given FITS files generated either by `dksraw` or obtained directly
from some other sources (e.g. iTelescope, an archive, etc).

### dksraw

**Note** Planned, only some are partially implemented at this stage. 
See TODO.md for plan and implementation status.

The command line `dksraw` application will provide the following subcommands:
- grey: Convert a RAW file into a single channel (greyscale) 16-bit PNG or 
        FITS file using one of several methods. **Partially implemented.**
- rgb: Convert a RAW file into an RGB PNG image. **Not yet implemented.**
- hist: Calculate and plot RGB or greyscale levels from a RAW image.  **Not yet implemented.** 
- split: Splits the input RAW file into separate 16-bit PNG files for each
         of the R, G, B and G channel in the Bayer mask. **Implemented.**
- whitebalance: Perform whitebalance calculations on the input RAW file in one
                of several ways. **Partially implemented as part of `grey`.**
- info: Print metadata about the input RAW file to stdout.  **Not yet implemented.**
- features: Performs image segmentation on the input RAW image and outputs a
            list of the coordinates and sizes of features thus found, along
            with an indexed image that can be used to view them.  **Not yet implemented.**

### ap_ scripts

(To be described)

# Installation Instructions

## Minimum Requirements

See `requirements.txt` for full dependency list. This is a Python 3 
project, with no intention to support Python 2.

## Optional Requirements

- `pytest` http://pytest.org (for running the test suite)
- `Sphinx` http://sphinx-doc.org (for generating documentation)


## Basic Setup

Install:

```bash
# Install for user (if not using a virtual environment).
$ python3 -m pip install . -r requirements.txt --user
# or install for the system (as root, or as a user in a virtual environment).
$ pip3 install . -r requirements.txt
```

To install in developer mode replace the last line with 
`pip3 install -e . -r requirements.txt`.

Get general command line application help:

```bash
$ dksraw --help
# Or for a specific subcommand, e.g. split
$ dksraw split --help
```

Run the test suite. For reasons I haven't resolved running pytest as a stand-alone
runs into the path issues related to https://docs.pytest.org/en/latest/pythonpath.html

```bash
    # Run tests capturing stdout
    $ python3 -m pytest -rfsP test/
    
    # Runs tests with a short summary of each test run
    $ python3 -m pytest -rfsp test/
```

To generate test coverage:
```bash
# Generates html files in the directory ./htmlcov
python3 -m pytest --cov-report html --cov=AstroPhotography test/

# Open htmlcov/index.html with a browser...
```

Build documentation:

```bash
    $ sphinx-build -b html doc doc/_build/html
```
