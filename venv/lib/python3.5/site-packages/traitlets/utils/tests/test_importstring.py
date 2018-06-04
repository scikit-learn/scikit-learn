# encoding: utf-8
# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
#
# Adapted from enthought.traits, Copyright (c) Enthought, Inc.,
# also under the terms of the Modified BSD License.
"""Tests for traitlets.utils.importstring."""

import os
from unittest import TestCase

from ..importstring import import_item


class TestImportItem(TestCase):

    def test_import_unicode(self):
        self.assertIs(os, import_item(u'os'))
        self.assertIs(os.path, import_item(u'os.path'))
        self.assertIs(os.path.join, import_item(u'os.path.join'))

    def test_bad_input(self):
        class NotAString(object):
            pass
        msg = (
            "import_item accepts strings, "
            "not '%s'." % NotAString
        )
        with self.assertRaisesRegexp(TypeError, msg):
            import_item(NotAString())
