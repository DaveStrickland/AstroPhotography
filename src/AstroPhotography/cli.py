""" Implementation of the command line interface for the AstroPhotography
    module.
    
This provides the following subcommands:
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

Typical workflow might consist of:
- Taking a RAW image using a digital camera attached to a telescope, 
  processing it with default options using `dksraw grey` in a FITS image, 
  viewing the FITS image with ds9 with asinh scaling to see both low and 
  high brightness features, and then adjusting the exposure time of further
  images.
- Looking at the levels using `dksraw hist`.
- Finding the location of specific objects within an image (e.g. Jupiter)
  using `dksraw features` to find the pixel locations to use in whitebalance
  calculations using `dksraw whitebalance`.

"""
from argparse import ArgumentParser
from inspect import getfullargspec

from . import __version__
from .api import hello
from .api import split
from .api import grey
from .core.config import config
from .core.logger import logger


__all__ = "main",


def main(argv=None) -> int:
    """ Execute the application CLI.

    :param argv: argument list to parse (sys.argv by default)
    :return: exit status
    """
    args = _args(argv)
    logger.start(args.warn or "DEBUG")  # can't use default from config yet
    logger.debug("Starting execution")
    config.load(args.config)
    config.core.config = args.config
    if args.warn:
        config.core.logging = args.warn
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
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", action="append",
            help="config file [etc/config.yml]")
    parser.add_argument("-v", "--version", action="version",
            version="AstroPhotography {:s}".format(__version__),
            help="print version and exit")
    parser.add_argument("-w", "--warn", default="WARN",
            help="logger warning level [WARN]")
            
    # Common options for all command parsers
    common = ArgumentParser(add_help=False)
    common.add_argument("rawfile",
        help="RAW file to process.") 
    common.add_argument("--output", "-o",
        default=None,
        help="Root name for output files." + 
        " If omitted then the name of the input RAW file " +
        "(with extension removed) will be used as the root name.")
    
    # Command parsers
    subparsers = parser.add_subparsers(title="Commands", 
        description="RAW file processing commands." +
        " One command must be selected by the user.",
        help="Command-specfic help.")
    _hello(subparsers, common)
    _grey(subparsers, common)
    _split(subparsers, common)
    
    # Process
    args = parser.parse_args(argv)
    if not args.config:
        # Don't specify this as an argument default or else it will always be
        # included in the list.
        args.config = "etc/config.yml"
    
    # Perform any additional command line argument checking or processing
    _check_args(args)    
        
    return args
 
def _check_args(args):
    """Runs any additional command line argument validation or processing
    """
    if ( args.command == split ):
        _check_output_args(args)
    elif ( args.command == grey ):
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

def _grey(subparsers, common):
    """Defines command line options for the grey command.

    :param subparsers: subcommand parsers
    :param common: parser for common subcommand arguments
    """
    default_wb = 'daylight'
    
    parser = subparsers.add_parser("grey", 
        description="Creates a monochrome output image using the specified white-balance.",
        parents=[common], 
        help='Creates a monochrome output image using the specified white-balance.')
    parser.add_argument('-w', '--whitebalance',
        default=default_wb,
        help='Whitebalance to use when convert R, G and B channels.' + 
            ' Default: {}'.format(default_wb))
    parser.add_argument('-b', '--black',
        default=False,
        action='store_true',
        help='Subtract camera band-specific black levels from the data. Default: False')
    parser.set_defaults(command=grey)
    return
    
def _split(subparsers, common):
    """Defines command line options for the split command.

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
    parser.add_argument('-b', '--black',
        default=False,
        action='store_true',
        help='Subtract camera band-specific black levels from the data. Default: False')
    parser.add_argument('--extension',
        default=default_ext,
        help='File name extension (i.e. type) for output files.' +
            ' Default: {}'.format(default_ext))
    parser.set_defaults(command=split)
    return

     
def _hello(subparsers, common):
    """Defines command line options for the hello command.

    :param subparsers: subcommand parsers
    :param common: parser for common subcommand arguments
    """
    parser = subparsers.add_parser("hello", 
        description="Hello world command used only for testing/debugging.",
        parents=[common], help='Test command for debugging purposes.')
    parser.add_argument("--name", "-n", 
        default="World", help="Greeting name")
    parser.set_defaults(command=hello)
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
