# AstroPhotography

The AstroPhotography python package provides python classes and command line 
applications for quick and convenient processing amateur astronomical images
from RAW digital camera format to common graphical image formats (e.g. PNG)
and astronomical FITS format.

The provided command line application `dksraw` provides a simple method of
quickly converting RAW files into useful images without the user having 
to mess around with `dcraw`, `gimp`, `photoshop` or the equivalent.

FITS images can then be viewed with the powerful features SAO `ds9`
provides.

## Command Line Functionality

**Note** Planned, not yet fully implemented. See TODO.md for plan and
implementation status.

The command line `dksraw` application will provide the following subcommands:
- grey: Convert a RAW file into a single channel (greyscale) 16-bit PNG or 
        FITS file using one of several methods.
- rgb: Convert a RAW file into an RGB PNG image.
- hist: Calculate and plot RGB or greyscale levels from a RAW image. 
- split: Splits the input RAW file into separate 16-bit PNG files for each
         of the R, G, B and G channel in the Bayer mask.
- whitebalance: Perform whitebalance calculations on the input RAW file in one
                of several ways.
- info: Print metadata about the input RAW file to stdout.
- features: Performs image segmentation on the input RAW image and outputs a
            list of the coordinates and sizes of features thus found, along
            with an indexed image that can be used to view them.

# Installation Instructions

## Minimum Requirements

See `requirements.txt` for full dependency list.
- Python 3.5+
- `numpy` and `matplotlib`
- `PyYAML`
- `rawpy` (python binding and interface to libraw)


## Optional Requirements

- `pytest` http://pytest.org (for running the test suite)
- `Sphinx` http://sphinx-doc.org (for generating documentation)


## Basic Setup

Install:

```bash
# Install for user
$ python3 -m pip install . -r requirements.txt --user
# or install for the system (as root)
$ pip3 install . -r requirements.txt
```

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

Build documentation:

```bash
    $ sphinx-build -b html doc doc/_build/html
```
