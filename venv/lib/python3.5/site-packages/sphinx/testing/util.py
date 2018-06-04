# -*- coding: utf-8 -*-
"""
    sphinx.testing.util
    ~~~~~~~~~~~~~~~~~~~

    Sphinx test suite utilities

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""
import os
import re
import sys
import warnings
from xml.etree import ElementTree

from docutils import nodes
from docutils.parsers.rst import directives, roles
from six import string_types

from sphinx import application
from sphinx.builders.latex import LaTeXBuilder
from sphinx.ext.autodoc import AutoDirective
from sphinx.pycode import ModuleAnalyzer
from sphinx.testing.path import path
from sphinx.util.osutil import relpath

if False:
    from typing import List  # NOQA


__all__ = [
    'Struct',
    'SphinxTestApp', 'SphinxTestAppWrapperForSkipBuilding',
    'remove_unicode_literals',
]


def assert_re_search(regex, text, flags=0):
    if not re.search(regex, text, flags):
        assert False, '%r did not match %r' % (regex, text)


def assert_not_re_search(regex, text, flags=0):
    if re.search(regex, text, flags):
        assert False, '%r did match %r' % (regex, text)


def assert_startswith(thing, prefix):
    if not thing.startswith(prefix):
        assert False, '%r does not start with %r' % (thing, prefix)


def assert_node(node, cls=None, xpath="", **kwargs):
    if cls:
        if isinstance(cls, list):
            assert_node(node, cls[0], xpath=xpath, **kwargs)
            if cls[1:]:
                if isinstance(cls[1], tuple):
                    assert_node(node, cls[1], xpath=xpath, **kwargs)
                else:
                    assert len(node) == 1, \
                        'The node%s has %d child nodes, not one' % (xpath, len(node))
                    assert_node(node[0], cls[1:], xpath=xpath + "[0]", **kwargs)
        elif isinstance(cls, tuple):
            assert len(node) == len(cls), \
                'The node%s has %d child nodes, not %r' % (xpath, len(node), len(cls))
            for i, nodecls in enumerate(cls):
                path = xpath + "[%d]" % i
                assert_node(node[i], nodecls, xpath=path, **kwargs)
        elif isinstance(cls, string_types):
            assert node == cls, 'The node %r is not %r: %r' % (xpath, cls, node)
        else:
            assert isinstance(node, cls), \
                'The node%s is not subclass of %r: %r' % (xpath, cls, node)

    for key, value in kwargs.items():
        assert key in node, 'The node%s does not have %r attribute: %r' % (xpath, key, node)
        assert node[key] == value, \
            'The node%s[%s] is not %r: %r' % (xpath, key, value, node[key])


def etree_parse(path):
    with warnings.catch_warnings(record=False):
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        return ElementTree.parse(path)


class Struct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class SphinxTestApp(application.Sphinx):
    """
    A subclass of :class:`Sphinx` that runs on the test root, with some
    better default values for the initialization parameters.
    """

    def __init__(self, buildername='html', srcdir=None,
                 freshenv=False, confoverrides=None, status=None, warning=None,
                 tags=None, docutilsconf=None):

        if docutilsconf is not None:
            (srcdir / 'docutils.conf').write_text(docutilsconf)

        builddir = srcdir / '_build'
#        if confdir is None:
        confdir = srcdir
#        if outdir is None:
        outdir = builddir.joinpath(buildername)
        if not outdir.isdir():
            outdir.makedirs()
#        if doctreedir is None:
        doctreedir = builddir.joinpath('doctrees')
        if not doctreedir.isdir():
            doctreedir.makedirs()
        if confoverrides is None:
            confoverrides = {}
#        if warningiserror is None:
        warningiserror = False

        self._saved_path = sys.path[:]
        self._saved_directives = directives._directives.copy()
        self._saved_roles = roles._roles.copy()

        self._saved_nodeclasses = set(v for v in dir(nodes.GenericNodeVisitor)
                                      if v.startswith('visit_'))

        try:
            application.Sphinx.__init__(self, srcdir, confdir, outdir, doctreedir,
                                        buildername, confoverrides, status, warning,
                                        freshenv, warningiserror, tags)
        except Exception:
            self.cleanup()
            raise

    def cleanup(self, doctrees=False):
        AutoDirective._registry.clear()
        ModuleAnalyzer.cache.clear()
        LaTeXBuilder.usepackages = []
        sys.path[:] = self._saved_path
        sys.modules.pop('autodoc_fodder', None)
        directives._directives = self._saved_directives
        roles._roles = self._saved_roles
        for method in dir(nodes.GenericNodeVisitor):
            if method.startswith('visit_') and \
               method not in self._saved_nodeclasses:
                delattr(nodes.GenericNodeVisitor, 'visit_' + method[6:])
                delattr(nodes.GenericNodeVisitor, 'depart_' + method[6:])

    def __repr__(self):
        return '<%s buildername=%r>' % (self.__class__.__name__, self.builder.name)


class SphinxTestAppWrapperForSkipBuilding(object):
    """
    This class is a wrapper for SphinxTestApp to speed up the test by skipping
    `app.build` process if it is already built and there is even one output
    file.
    """

    def __init__(self, app_):
        self.app = app_

    def __getattr__(self, name):
        return getattr(self.app, name)

    def build(self, *args, **kw):
        if not self.app.outdir.listdir():
            # if listdir is empty, do build.
            self.app.build(*args, **kw)
            # otherwise, we can use built cache


_unicode_literals_re = re.compile(r'u(".*?")|u(\'.*?\')')


def remove_unicode_literals(s):
    return _unicode_literals_re.sub(lambda x: x.group(1) or x.group(2), s)


def find_files(root, suffix=None):
    for dirpath, dirs, files in os.walk(root, followlinks=True):
        dirpath = path(dirpath)
        for f in [f for f in files if not suffix or f.endswith(suffix)]:
            fpath = dirpath / f
            yield relpath(fpath, root)


def strip_escseq(text):
    return re.sub('\x1b.*?m', '', text)
