#! /usr/bin/env python
# Last Change: Mon Aug 20 02:00 PM 2007 J

try:
    from functools import partial
    print "Using built-in partial"
except ImportError:
    from myfunctools import partial
