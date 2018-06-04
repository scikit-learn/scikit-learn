from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
from nose.plugins.errorclass import ErrorClass, ErrorClassPlugin
from ..exceptions import KnownFailureTest


class KnownFailure(ErrorClassPlugin):
    '''Plugin that installs a KNOWNFAIL error class for the
    KnownFailureClass exception.  When KnownFailureTest is raised,
    the exception will be logged in the knownfail attribute of the
    result, 'K' or 'KNOWNFAIL' (verbose) will be output, and the
    exception will not be counted as an error or failure.

    This is based on numpy.testing.noseclasses.KnownFailure.
    '''
    enabled = True
    knownfail = ErrorClass(KnownFailureTest,
                           label='KNOWNFAIL',
                           isfailure=False)

    def options(self, parser, env=os.environ):
        env_opt = 'NOSE_WITHOUT_KNOWNFAIL'
        parser.add_option('--no-knownfail', action='store_true',
                          dest='noKnownFail', default=env.get(env_opt, False),
                          help='Disable special handling of KnownFailureTest '
                               'exceptions')

    def configure(self, options, conf):
        if not self.can_configure:
            return
        self.conf = conf
        disable = getattr(options, 'noKnownFail', False)
        if disable:
            self.enabled = False

    def addError(self, test, err, *zero_nine_capt_args):
        # Fixme (Really weird): if I don't leave empty method here,
        # nose gets confused and KnownFails become testing errors when
        # using the MplNosePlugin and MplTestCase.

        # The *zero_nine_capt_args captures an extra argument. There
        # seems to be a bug in
        # nose.testing.manager.ZeroNinePlugin.addError() in which a
        # 3rd positional argument ("capt") is passed to the plugin's
        # addError() method, even if one is not explicitly using the
        # ZeroNinePlugin.
        pass
