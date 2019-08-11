""" Main application entry point.

    Actually this is unneccesary, as setup creates and installs dksraw 
    from cli.py, so this entire  file may disappear.
    
    python -m AstroPhotography  ...

"""

import sys

def warn():
    """ Inform the user that they should use dksraw directlly. 
    """
    print("To run AstroPhotography's command line tools call dksraw directly.")
    sys.exit(1)



# Make the script executable.
if __name__ == "__main__":
    warn()
