# -*- coding: utf-8 -*-
"""Test IPython.core.logger"""

import os.path

import nose.tools as nt
from IPython.utils.tempdir import TemporaryDirectory

_ip = get_ipython()

def test_logstart_inaccessible_file():
    try:
        _ip.logger.logstart(logfname="/")   # Opening that filename will fail.
    except IOError:
        pass
    else:
        nt.assert_true(False)           # The try block should never pass.
    
    try:
        _ip.run_cell("a=1")                 # Check it doesn't try to log this
    finally:
        _ip.logger.log_active = False  # If this fails, don't let later tests fail

def test_logstart_unicode():
    with TemporaryDirectory() as tdir:
        logfname = os.path.join(tdir, "test_unicode.log")
        _ip.run_cell("'abc€'")
        try:
            _ip.magic("logstart -to %s" % logfname)
            _ip.run_cell("'abc€'")
        finally:
            _ip.logger.logstop()
