"""
Test the logger module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org> 
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

import shutil
import os
from tempfile import mkdtemp

from ..logger import PrintTime

################################################################################
# Test fixtures
env = dict()

def setup():
    """ Test setup.
    """
    cachedir = mkdtemp()
    if os.path.exists(cachedir):
        shutil.rmtree(cachedir)
    env['dir'] = cachedir
    

def teardown():
    """ Test teardown.
    """
    #return True
    shutil.rmtree(env['dir'])


################################################################################
# Tests
def test_print_time():
    """ A simple smoke test for PrintTime.
    """
    print_time = PrintTime(logfile=os.path.join(env['dir'], 'test.log'))
    print_time('Foo')
    # Create a second time, to smoke test log rotation.
    print_time = PrintTime(logfile=os.path.join(env['dir'], 'test.log'))
    print_time('Foo')
    # And a third time 
    print_time = PrintTime(logfile=os.path.join(env['dir'], 'test.log'))
    print_time('Foo')

