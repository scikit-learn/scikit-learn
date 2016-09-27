"""
=========
traitsdoc
=========

Sphinx extension that handles docstrings in the Numpy standard format, [1]
and support Traits [2].

This extension can be used as a replacement for ``numpydoc`` when support
for Traits is required.

.. [1] http://projects.scipy.org/numpy/wiki/CodingStyleGuidelines#docstring-standard
.. [2] http://code.enthought.com/projects/traits/

"""
from __future__ import division, absolute_import, print_function

import inspect
import os
import pydoc
import collections

from . import docscrape
from . import docscrape_sphinx
from .docscrape_sphinx import SphinxClassDoc, SphinxFunctionDoc, SphinxDocString

from . import numpydoc

from . import comment_eater

class SphinxTraitsDoc(SphinxClassDoc):
    def __init__(self, cls, modulename='', func_doc=SphinxFunctionDoc):
        if not inspect.isclass(cls):
            raise ValueError("Initialise using a class. Got %r" % cls)
        self._cls = cls

        if modulename and not modulename.endswith('.'):
            modulename += '.'
        self._mod = modulename
        self._name = cls.__name__
        self._func_doc = func_doc

        docstring = pydoc.getdoc(cls)
        docstring = docstring.split('\n')

        # De-indent paragraph
        try:
            indent = min(len(s) - len(s.lstrip()) for s in docstring
                         if s.strip())
        except ValueError:
            indent = 0

        for n,line in enumerate(docstring):
            docstring[n] = docstring[n][indent:]

        self._doc = docscrape.Reader(docstring)
        self._parsed_data = {
            'Signature': '',
            'Summary': '',
            'Description': [],
            'Extended Summary': [],
            'Parameters': [],
            'Returns': [],
            'Yields': [],
            'Raises': [],
            'Warns': [],
            'Other Parameters': [],
            'Traits': [],
            'Methods': [],
            'See Also': [],
            'Notes': [],
            'References': '',
            'Example': '',
            'Examples': '',
            'index': {}
            }

        self._parse()

    def _str_summary(self):
        return self['Summary'] + ['']

    def _str_extended_summary(self):
        return self['Description'] + self['Extended Summary'] + ['']

    def __str__(self, indent=0, func_role="func"):
        out = []
        out += self._str_signature()
        out += self._str_index() + ['']
        out += self._str_summary()
        out += self._str_extended_summary()
        for param_list in ('Parameters', 'Traits', 'Methods',
                           'Returns', 'Yields', 'Raises'):
            out += self._str_param_list(param_list)
        out += self._str_see_also("obj")
        out += self._str_section('Notes')
        out += self._str_references()
        out += self._str_section('Example')
        out += self._str_section('Examples')
        out = self._str_indent(out,indent)
        return '\n'.join(out)

def looks_like_issubclass(obj, classname):
    """ Return True if the object has a class or superclass with the given class
    name.

    Ignores old-style classes.
    """
    t = obj
    if t.__name__ == classname:
        return True
    for klass in t.__mro__:
        if klass.__name__ == classname:
            return True
    return False

def get_doc_object(obj, what=None, config=None):
    if what is None:
        if inspect.isclass(obj):
            what = 'class'
        elif inspect.ismodule(obj):
            what = 'module'
        elif isinstance(obj, collections.Callable):
            what = 'function'
        else:
            what = 'object'
    if what == 'class':
        doc = SphinxTraitsDoc(obj, '', func_doc=SphinxFunctionDoc, config=config)
        if looks_like_issubclass(obj, 'HasTraits'):
            for name, trait, comment in comment_eater.get_class_traits(obj):
                # Exclude private traits.
                if not name.startswith('_'):
                    doc['Traits'].append((name, trait, comment.splitlines()))
        return doc
    elif what in ('function', 'method'):
        return SphinxFunctionDoc(obj, '', config=config)
    else:
        return SphinxDocString(pydoc.getdoc(obj), config=config)

def setup(app):
    # init numpydoc
    numpydoc.setup(app, get_doc_object)

