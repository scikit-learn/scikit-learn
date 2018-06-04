from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import gc
import os
from nose.plugins import Plugin


class PerformGC(Plugin):
    """This plugin adds option to call ``gc.collect`` after each test"""
    enabled = False

    def options(self, parser, env=os.environ):
        env_opt = 'PERFORM_GC'
        parser.add_option('--perform-gc', action='store_true',
                          dest='performGC', default=env.get(env_opt, False),
                          help='Call gc.collect() after each test')

    def configure(self, options, conf):
        if not self.can_configure:
            return

        self.enabled = getattr(options, 'performGC', False)

    def afterTest(self, test):
        gc.collect()
