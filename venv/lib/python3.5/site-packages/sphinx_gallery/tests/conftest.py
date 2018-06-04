# -*- coding: utf-8 -*-
"""
Pytest fixtures
"""
from __future__ import division, absolute_import, print_function

import collections
import logging

import pytest

import sphinx_gallery.docs_resolv
import sphinx_gallery.gen_gallery
import sphinx_gallery.gen_rst
from sphinx_gallery import sphinx_compatibility


Params = collections.namedtuple('Params', 'args kwargs')


class FakeSphinxApp:
    def __init__(self):
        self.calls = collections.defaultdict(list)

    def status_iterator(self, *args, **kwargs):
        self.calls['status_iterator'].append(Params(args, kwargs))
        for it in args[0]:
            yield it

    def warning(self, *args, **kwargs):
        self.calls['warning'].append(Params(args, kwargs))

    def warn(self, *args, **kwargs):
        self.calls['warn'].append(Params(args, kwargs))

    def info(self, *args, **kwargs):
        self.calls['info'].append(Params(args, kwargs))

    def verbose(self, *args, **kwargs):
        self.calls['verbose'].append(Params(args, kwargs))

    def debug(self, *args, **kwargs):
        self.calls['debug'].append(Params(args, kwargs))


@pytest.fixture
def fakesphinxapp():
    orig_app = sphinx_gallery.sphinx_compatibility._app
    sphinx_gallery.sphinx_compatibility._app = app = FakeSphinxApp()
    try:
        yield app
    finally:
        sphinx_gallery.sphinx_compatibility._app = orig_app


@pytest.fixture
def log_collector():
    orig_dr_logger = sphinx_gallery.docs_resolv.logger
    orig_gg_logger = sphinx_gallery.gen_gallery.logger
    orig_gr_logger = sphinx_gallery.gen_rst.logger
    app = FakeSphinxApp()
    sphinx_gallery.docs_resolv.logger = app
    sphinx_gallery.gen_gallery.logger = app
    sphinx_gallery.gen_rst.logger = app
    try:
        yield app
    finally:
        sphinx_gallery.docs_resolv.logger = orig_dr_logger
        sphinx_gallery.gen_gallery.logger = orig_gg_logger
        sphinx_gallery.gen_rst.logger = orig_gr_logger
