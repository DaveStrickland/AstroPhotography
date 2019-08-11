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
        help="Name or root name for output files." + 
        " If omitted then the name of the input RAW file " +
        "(with extension removed) will be used as the root name.")
    
    # Command parsers
    subparsers = parser.add_subparsers(title="Commands", 
        description="RAW file processing commands." +
        "One command must be selected by the user.",
        help="Command-specfic help.")
    _hello(subparsers, common)
    _split(subparsers, common)
    
    # Process
    args = parser.parse_args(argv)
    if not args.config:
        # Don't specify this as an argument default or else it will always be
        # included in the list.
        args.config = "etc/config.yml"
    return args
 

def _split(subparsers, common):
    """ CLI adaptor for the api.split command.

    :param subparsers: subcommand parsers
    :param common: parser for common subcommand arguments
    """
    parser = subparsers.add_parser("split", 
        description="Exports raw Bayer map as separate images with suffix _r.tiff," +
        " g1.tiff, _b.tiff and _g2.tiff.",
        parents=[common], help='Outputs raw, unmodified, R, G, B and G as TIFF files.')
    #parser.add_argument("--name", "-n", 
    #    default="World", help="Greeting name")
    parser.set_defaults(command=split)
    return
     
def _hello(subparsers, common):
    """ CLI adaptor for the api.hello command.

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
