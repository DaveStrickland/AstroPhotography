# AstroPhotography

The aim of the AstroPhotography python package is to provide Python 
classes and command line applications for amateur astronomy, specifically:

- Reduction and combination of multiple FITS images, including calibration,
  artifact removal, star detection, and astrometry.
- Quick inspection and conversion of RAW digital camera formats to common 
  graphical image formats (e.g., PNG) and astronomical FITS format.

FITS images can then be viewed with the powerful features SAO `ds9`
provides, and/or within Python using `astropy`.

*Current status:* This is still very much a work in progress. RAW to 
image/FITS conversion is partially implemented, but development stalled
when I became dissatisfied with the clunky nature of the unit tests. I've
been putting more effort into the FITS data reduction side of the project,
but still have a long way to go.

See [RELEASE_NOTES.md](RELEASE_NOTES.md) for a summary of what has 
changed recently.

## Command Line Functionality

A series of Python scripts, all beginning with the prefix `ap_`, 
perform separate command-line stages of astronomical 
image data reduction using FITS files generated either by `dksraw` or obtained directly
from other sources (e.g., iTelescope, an archive).

These scripts use Python classes and functions (typically starting with `Ap`)
that provide the main functionality. These can also be called from Python or Jupyter notebooks.
A set of notebooks provides example walk-throughs, particularly for processing calibrated FITS images from iTelescope.

The command line Python program `dksraw` provides a simple method of
quickly converting RAW files into images or FITS files without the user needing `dcraw`, `gimp`, `photoshop`, etc.

### ap_ scripts

A series of Python classes for FITS data processing (with names beginning 
with `Ap`) can be used from the command line using scripts (names 
beginning with `ap_`).

(To be described, but see [doc/iTelescope_processing.md](doc/iTelescope_processing.md)
for a high-level summary. The notebooks in `AstroPhotography/notebooks` are more up-to-date.)

### dksraw

**Note:** Some features are partially implemented. 

The command line `dksraw` application provides the following subcommands:

- `grey`: Convert a RAW file into a single-channel (greyscale) 16-bit PNG/JPG/TIFF or FITS file. **Working with limited options.**
- `rgb`: Convert a RAW file into an RGB PNG/JPG/TIFF image or FITS file. **Working with limited options.**
- `split`: Splits the input RAW file into separate 16-bit PNG/JPG/TIFF/FITS 
  files for each of the R, G, B, and G channels in the Bayer mask. **Implemented.**
- `whitebalance`: Perform whitebalance calculations on the input RAW file. **Partially implemented as part of `grey`.**
- `info`: Print metadata about the input RAW file. **Not yet implemented.**

Development of RAW fie processing features has been limited as I am between
cameras and lack access to a broad range of good datasets with different good 
RAW file formats

## Installation Instructions

### Minimum Requirements

This is a Python 3 project (>=3.6), with no intention to support Python 2.  
Dependencies are managed in `pyproject.toml`. The runtime dependencies include:

- PyYAML
- matplotlib
- numpy
- rawpy
- imageio
- astropy
- regions
- astroplan
- astroquery
- astroscrappy
- ccdproc
- photutils
- ExifRead

Optional dependencies for testing and documentation are defined as "extras" (`test`, `docs`) in `pyproject.toml`.

---

### Basic Setup (Windows / Linux / macOS)

> **Note:** Windows users may need to upgrade pip for editable installs.

```powershell
# Upgrade pip (required for editable installs)
python -m pip install --upgrade pip

# Install build tools
pip install --upgrade build setuptools wheel

# Optional: clean old build artifacts
Remove-Item -Recurse -Force build, dist, *.egg-info -ErrorAction Ignore
```

### Install the package
```powershell
# Normal install
pip install .

# Install with test and docs extras
pip install .[test,docs]

# Edqitable install (for development)
pip install -e .[test,docs]
```

### Using conda (optional)
```powershell
conda env create -f ap-env.yml
conda activate ap-env

# Then install the package
pip install -e .[test,docs]
```

### Verifying Installation
```powershell
# CLI help
dksraw --help
# Specific subcommand
dksraw split --help
```

### Developers Only

Formatting and linting follow Ruff. [TBA]

Run the test suite:
'''powershell
# Basic test run
python -m pytest -rfsP test/

# Short summary of each test run
python -m pytest -rfsp test/
'''

Generate test coverage:
```powershell
python -m pytest --cov-report html --cov=AstroPhotography test/
# Open htmlcov/index.html in a browser
```

To run a subset of tests::
```bash
python -m pytest tests.test_astrophotography
```

### Documentation
```powershell
cd doc
make html
# Open doc/_build/html/index.html in a browser
```
