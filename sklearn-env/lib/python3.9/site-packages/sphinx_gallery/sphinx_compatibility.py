# -*- coding: utf-8 -*-
"""
Backwards-compatility shims for Sphinx
======================================

"""
from __future__ import division, absolute_import, print_function

import sphinx
import sphinx.util


# This gets set when the extension is initialized.
_app = None


def _app_get_logger(name):
    class SphinxLoggerAdapter:
        def _color_to_func(self, kwargs, default=''):
            return getattr(sphinx.util.console,
                           kwargs.pop('color', default),
                           None)

        def error(self, msg, *args, **kwargs):
            msg = msg % args
            colorfunc = self._color_to_func(kwargs, default='red')
            return _app.warn(colorfunc(msg), **kwargs)

        def critical(self, msg, *args, **kwargs):
            msg = msg % args
            colorfunc = self._color_to_func(kwargs, default='red')
            return _app.warn(colorfunc(msg), **kwargs)

        def warning(self, msg, *args, **kwargs):
            msg = msg % args
            colorfunc = self._color_to_func(kwargs)
            if colorfunc:
                # colorfunc is a valid kwarg in 1.5, but not older, so we just
                # apply it ourselves.
                msg = colorfunc(msg)
            return _app.warn(msg, **kwargs)

        def info(self, msg='', *args, **kwargs):
            msg = msg % args
            colorfunc = self._color_to_func(kwargs)
            if colorfunc:
                msg = colorfunc(msg)
            return _app.info(msg, **kwargs)

        def verbose(self, msg, *args, **kwargs):
            return _app.verbose(msg, *args, **kwargs)

        def debug(self, msg, *args, **kwargs):
            return _app.debug(msg, *args, **kwargs)

    return SphinxLoggerAdapter()


def _app_status_iterator(iterable, summary, **kwargs):
    global _app

    color = kwargs.pop('color', None)
    if color is not None:
        kwargs['colorfunc'] = getattr(sphinx.util.console, color)

    for item in _app.status_iterator(iterable, summary, **kwargs):
        yield item


getLogger = sphinx.util.logging.getLogger
status_iterator = sphinx.util.status_iterator
