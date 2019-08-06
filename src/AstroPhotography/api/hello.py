""" Implements the hello command.

"""
from ..core.logger import logger


def main(name="World") -> str:
    """ Execute the command.
    
    :param name: name to use in greeting
    """
    logger.debug("Executing hello command")
    string =  "Hello, {:s}!".format(name)  # TODO: use f-string for Python 3.6+
    print(string)
    return string
