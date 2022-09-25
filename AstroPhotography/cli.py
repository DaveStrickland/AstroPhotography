"""
Implementation of the command line interface for the dksraw
command within the AstroPhotography python module.
    
The dksraw command is intended to provide the following subcommands:

- grey: Convert a RAW file into a single channel (greyscale) 16-bit PNG or 
  FITS file using one of several methods.
- split: Splits the input RAW file into separate 16-bit PNG files for each
  of the R, G, B and G channel in the Bayer mask.
- rgb: Convert a RAW file into an RGB PNG image.
- whitebalance: Perform whitebalance calculations on the input RAW file in one
  of several ways.
- info: Print metadata about the input RAW file to stdout.

Only the first threee (rg, grey and split) are currently implemented, although
parts of the whitebalance calculations are implemented.

Typical use-cases might be:

- Taking a RAW image using a digital camera attached to a telescope,
  and generating a quick look image using `dksraw grey --renomalize --loglevel DEBUG`
  to check both framing and raw image count percentiles.
- Later, converting all the RAW images taken in that session to FITS
  files, using the default options of `dksraw grey` or `dksraw rgb`,
  and then processing those with the ap_ scripts to build calibrated
  and resampled composite image.  
- Looking at basic information from a set of RAW files using 
  `dksraw info`. (Still to be implemented.)

"""
from argparse import ArgumentParser
from inspect import getfullargspec

from . import __version__
from .api import split
from .api import grey
from .api import rgb
from .core.config import config
from .core.logger import logger
import sys

__all__ = "main",


def main(argv=None) -> int:
    """
    Execute the application CLI.

    :param argv: argument list to parse (sys.argv by default)
    :return: exit status
    """
    args = _args(argv)
    logger.start(args.loglevel or "DEBUG")  # can't use default from config yet
    logger.debug("Starting execution")
    config.load(args.config)
    config.core.config = args.config
    if args.loglevel:
        config.core.logging = args.loglevel
    logger.stop()  # clear handlers to prevent duplicate records
    logger.start(config.core.logging)
    command = args.command
    args = vars(args)
    spec = getfullargspec(command)
    if not spec.varkw:
        # No kwargs, remove unexpected arguments.
        args = {key: args[key] for key in args if key in spec.args}
    try:
        command(**args)
    except RuntimeError as err:
        logger.critical(err)
        return 1
    logger.debug("Successful completion")
    return 0
 

def _args(argv):
    """
    Parse command line arguments.

    :param argv: argument list to parse
    """
    p_loglevel = "INFO"
    parser = ArgumentParser(prog='dksraw',
        description='Dave\'s RAW file conversion tool for AstroPhotography')
    parser.add_argument("-c", "--config", action="append",
            help="config file [etc/config.yml]")
            
    # Common options for all command parsers
    common = ArgumentParser(add_help=False)
    common.add_argument("rawfile",
        help="RAW file to process.") 
    common.add_argument("--output", "-o",
        default=None,
        help="Root name for output files." + 
        " If omitted then the name of the input RAW file " +
        "(with extension removed) will be used as the root name.")
    common.add_argument("-v", "--version", action="version",
            version="AstroPhotography {:s}".format(__version__),
            help="print version and exit")
    common.add_argument("-l", "--loglevel", default=p_loglevel,
            help=("Logger informational level, one of NOTSET, DEBUG, INFO"
            ", WARNING, ERROR, or CRITICAL in order of increasing severity."
            f" Default: {p_loglevel}"))
    
    # Command parsers
    subparsers = parser.add_subparsers(title="Commands", 
        description="RAW file processing commands." +
        " One command MUST be selected by the user.",
        help="Command-specfic help.")
    _rgb(subparsers, common)
    _grey(subparsers, common)
    _split(subparsers, common)
    
    # Process
    args = parser.parse_args(argv)
    if not args.config:
        # Don't specify this as an argument default or else it will always be
        # included in the list.
        args.config = "etc/config.yml"
    
    # Perform any additional command line argument checking or processing
    _check_args(parser, args)    
        
    return args
 
def _check_args(parser, args):
    """Runs any additional command line argument validation or processing
    """
    if not hasattr(args, 'command'):
        # Command not defined.
        parser.print_help()
        print('Error: no command was specified.')
        sys.exit(1)
    
    if ( args.command == split ):
        _check_output_args(args)
    elif ( (args.command == grey) or (args.command == rgb)):
        default_suffix = '.fits'
        if args.output is None:
            _check_output_args(args)
            args.output += default_suffix
    return
    
def _check_output_args(args):
    """Additional processing for the output argument.
    """
    # If args.output is not defined we need to craft a suitable output
    # prefix to which split will add vaious band-specific suffixes.
    if args.output is None:
        # For split command output is root name for output PNG images.
        idx = args.rawfile.rfind('.')
        if ( idx != -1 ):
            slash_idx = args.rawfile.rfind('/')
            if (slash_idx == -1):
                start_idx = 0
            else:
                start_idx = slash_idx + 1
            args.output = args.rawfile[start_idx:idx]
        else:
            raise RuntimeError('Could not determine root name from {}'.format(args.rawfile))
    return

def _rgb(subparsers, common):
    """Defines command line options for the rgb command.

    :param subparsers: subcommand parsers
    :param common: parser for common subcommand arguments
    """
    allowed_wb = ['daylight', 'camera', 'auto', 'region[regspec]', 'user[userspec]']
    default_wb = 'camera'
    allowed_method = ['linear']
    default_method = 'linear'
    
    parser = subparsers.add_parser("rgb", 
        description="Creates a 3-channel RGB output image using the specified method and white-balance.",
        parents=[common], 
        help='Creates a RGB output image using the specified method and white-balance.')
    parser.add_argument('-m', '--method',
        default=default_method,
        help=('Method used to assemble the monochrome luminance image'
            f' from the Bayer sub channels. Available options are: {allowed_method}'
            ' linear: Performs the rawpy linear postprocess to an RBG image,'
            ' then applies CCIR 601 luma coefficients to obtain greyscale.'
            ' This de-Bayers the image, but unfortunately renormalizes it'
            ' to fill the full dynamic range available in the output.'
            f' Default: {default_method}'))
    parser.add_argument('-w', '--whitebalance',
        default=default_wb,
        help='Whitebalance to use when convert R, G and B channels.' + 
            ' Allowed whitebalance methods are: {}'.format(allowed_wb) +
            ' To calculate the whitebalance from the entire image use "auto".' +
            ' To calculate the whitebalance from part of an image use "region"' +
            ' with a region specifier of the form [minrow, maxrow, mincol, maxcol],' +
            ' where the pixel indices are zero-based and inclusive. For example:' +
            ' "region[450, 463, 2850, 2863]".' +
            ' To specify a user selected whitebalance use "user" with a user specifier' +
            ' of form [Rmult, G1mult, Bmult, G2mult].' +
            ' For example "user[1.85, 1.0, 2.01, 1.0]".' +
            ' The region and user options should be enclosed in quotes to prevent shell expansion.' +
            ' Default: {}'.format(default_wb))
    parser.add_argument('--keepblack',
        default=False,
        action='store_true',
        help=('Retain the camera band-specific black levels in the data.'
            ' These are roughly equivalent to a CCD bias level.'
            ' Default: False'))
    parser.add_argument('--renormalize',
        default=False,
        action='store_true',
        help=('If specified the output image will be linearly renormalized'
            ' to span the dynamic range 0 to (2^16)-1 ADU, similar to the'
            ' ImageMagick -normalize option.'
            ' This option is useful for quick look purposes, but as it'
            ' alters the output data values it is not suitable for later'
            ' image combination or processing.'))
    parser.set_defaults(command=rgb)
    return

def _grey(subparsers, common):
    """
    Defines command line options for the grey command.

    :param subparsers: subcommand parsers
    :param common: parser for common subcommand arguments
    """
    allowed_wb = ['daylight', 'camera', 'auto', 'region[regspec]', 'user[userspec]']
    default_wb = 'camera'
    allowed_method = ['linear', 'direct']
    default_method = 'linear'
    
    parser = subparsers.add_parser("grey", 
        description="Creates a monochrome output image using the specified method and white-balance.",
        parents=[common], 
        help='Creates a monochrome output image using the specified method and white-balance.')
    parser.add_argument('-m', '--method',
        default=default_method,
        help=('Method used to assemble the monochrome luminance image'
            f' from the Bayer sub channels. Available options are: {allowed_method}'
            ' linear: Performs the rawpy linear postprocess to an RBG image,'
            ' then applies CCIR 601 luma coefficients to obtain greyscale.'
            ' This de-Bayers the image, but unfortunately renormalizes it'
            ' to fill the full dynamic range available in the output.'
            ' direct: Direct does not perform any de-Bayer calculation, instead'
            ' just setting each pixel to its whitebalance-scaled value from'
            ' which ever RGBG subband it can from. This will look spotty'
            ' unless the correct whitebalance is chosen.'
            f' Default: {default_method}'))
    parser.add_argument('-w', '--whitebalance',
        default=default_wb,
        help='Whitebalance to use when convert R, G and B channels.' + 
            ' Allowed whitebalance methods are: {}'.format(allowed_wb) +
            ' To calculate the whitebalance from the entire image use "auto".' +
            ' To calculate the whitebalance from part of an image use "region"' +
            ' with a region specifier of the form [minrow, maxrow, mincol, maxcol],' +
            ' where the pixel indices are zero-based and inclusive. For example:' +
            ' "region[450, 463, 2850, 2863]".' +
            ' To specify a user selected whitebalance use "user" with a user specifier' +
            ' of form [Rmult, G1mult, Bmult, G2mult].' +
            ' For example "user[1.85, 1.0, 2.01, 1.0]".' +
            ' The region and user options should be enclosed in quotes to prevent shell expansion.' +
            ' Default: {}'.format(default_wb))
    parser.add_argument('--keepblack',
        default=False,
        action='store_true',
        help=('Retain the camera band-specific black levels in the data.'
            ' These are roughly equivalent to a CCD bias level.'
            ' Default: False'))
    parser.add_argument('--renormalize',
        default=False,
        action='store_true',
        help=('If specified the output image will be linearly renormalized'
            ' to span the dynamic range 0 to (2^16)-1 ADU, similar to the'
            ' ImageMagick -normalize option.'
            ' This option is useful for quick look purposes, but as it'
            ' alters the output data values it is not suitable for later'
            ' image combination or processing.'))
    parser.set_defaults(command=grey)
    return
    
def _split(subparsers, common):
    """
    Defines command line options for the split command.

    :param subparsers: subcommand parsers
    :param common: parser for common subcommand arguments
    """
    
    # png is very slow with imageio, and imageio doesn't
    # support 16-bit jp2 or jpeg.
    default_ext = 'tiff'
    
    parser = subparsers.add_parser("split", 
        description="Exports raw Bayer map as separate images with suffix _r.tiff," +
        " g1.tiff, _b.tiff and _g2.tiff.",
        parents=[common], 
        help='Outputs raw, unmodified, R, G, B and G as separate TIFF files.' +
            'Output files use the specified or default output root name followed' +
            ' by a channel-specific suffix, and specified or default file name extension.')
    parser.add_argument('--keepblack',
        default=False,
        action='store_true',
        help=('Retain the camera band-specific black levels in the data.'
            ' These are roughly equivalent to a CCD bias level.'
            ' Default: False'))
    parser.add_argument('--extension',
        default=default_ext,
        help='File name extension (i.e. type) for output files.' +
            ' Default: {}'.format(default_ext))
    parser.set_defaults(command=split)
    return

# Make the module executable.
if __name__ == "__main__":
    try:
        status = main()
    except:
        logger.critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
