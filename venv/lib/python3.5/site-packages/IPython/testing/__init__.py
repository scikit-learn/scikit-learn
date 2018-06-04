"""Testing support (tools to test IPython itself).
"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2009-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

# User-level entry point for testing
def test(**kwargs):
    """Run the entire IPython test suite.

    Any of the options for run_iptestall() may be passed as keyword arguments.

    For example::
    
        IPython.test(testgroups=['lib', 'config', 'utils'], fast=2)
    
    will run those three sections of the test suite, using two processes.
    """

    # Do the import internally, so that this function doesn't increase total
    # import time
    from .iptestcontroller import run_iptestall, default_options
    options = default_options()
    for name, val in kwargs.items():
        setattr(options, name, val)
    run_iptestall(options)

# So nose doesn't try to run this as a test itself and we end up with an
# infinite test loop
test.__test__ = False
