"""
========
numpydoc
========

Sphinx extension that handles docstrings in the Numpy standard format. [1]

It will:

- Convert Parameters etc. sections to field lists.
- Convert See Also section to a See also entry.
- Renumber references.
- Extract the signature from the docstring, if it can't be determined
  otherwise.

.. [1] http://projects.scipy.org/numpy/wiki/CodingStyleGuidelines#docstring-standard

"""

from __future__ import unicode_literals

import sys # Only needed to check Python version
import os
import re
import pydoc
from .docscrape_sphinx import get_doc_object
from .docscrape_sphinx import SphinxDocString
import inspect


def mangle_docstrings(app, what, name, obj, options, lines,
                      reference_offset=[0]):

    cfg = dict(use_plots=app.config.numpydoc_use_plots,
               show_class_members=app.config.numpydoc_show_class_members)

    if what == 'module':
        # Strip top title
        title_re = re.compile(r'^\s*[#*=]{4,}\n[a-z0-9 -]+\n[#*=]{4,}\s*',
                              re.I | re.S)
        lines[:] = title_re.sub('', "\n".join(lines)).split("\n")
    else:
        doc = get_doc_object(obj, what, "\n".join(lines), config=cfg)
        if sys.version_info[0] < 3:
            lines[:] = unicode(doc).splitlines()
        else:
            lines[:] = str(doc).splitlines()

    if app.config.numpydoc_edit_link and hasattr(obj, '__name__') and \
           obj.__name__:
        if hasattr(obj, '__module__'):
            v = dict(full_name="%s.%s" % (obj.__module__, obj.__name__))
        else:
            v = dict(full_name=obj.__name__)
        lines += [u'', u'.. htmlonly::', '']
        lines += [u'    %s' % x for x in
                  (app.config.numpydoc_edit_link % v).split("\n")]

    # replace reference numbers so that there are no duplicates
    references = []
    for line in lines:
        line = line.strip()
        m = re.match(r'^.. \[([a-z0-9_.-])\]', line, re.I)
        if m:
            references.append(m.group(1))

    # start renaming from the longest string, to avoid overwriting parts
    references.sort(key=lambda x: -len(x))
    if references:
        for i, line in enumerate(lines):
            for r in references:
                if re.match(r'^\d+$', r):
                    new_r = "R%d" % (reference_offset[0] + int(r))
                else:
                    new_r = u"%s%d" % (r, reference_offset[0])
                lines[i] = lines[i].replace(u'[%s]_' % r,
                                            u'[%s]_' % new_r)
                lines[i] = lines[i].replace(u'.. [%s]' % r,
                                            u'.. [%s]' % new_r)

    reference_offset[0] += len(references)


def mangle_signature(app, what, name, obj,
                     options, sig, retann):
    # Do not try to inspect classes that don't define `__init__`
    if (inspect.isclass(obj) and
        (not hasattr(obj, '__init__') or
        'initializes x; see ' in pydoc.getdoc(obj.__init__))):
        return '', ''

    if not (callable(obj) or hasattr(obj, '__argspec_is_invalid_')):
        return
    if not hasattr(obj, '__doc__'):
        return

    doc = SphinxDocString(pydoc.getdoc(obj))
    if doc['Signature']:
        sig = re.sub("^[^(]*", "", doc['Signature'])
        return sig, ''


def setup(app, get_doc_object_=get_doc_object):
    global get_doc_object
    get_doc_object = get_doc_object_

    if sys.version_info[0] < 3:
        app.connect(b'autodoc-process-docstring', mangle_docstrings)
        app.connect(b'autodoc-process-signature', mangle_signature)
    else:
        app.connect('autodoc-process-docstring', mangle_docstrings)
        app.connect('autodoc-process-signature', mangle_signature)
    app.add_config_value('numpydoc_edit_link', None, False)
    app.add_config_value('numpydoc_use_plots', None, False)
    app.add_config_value('numpydoc_show_class_members', True, True)

    # Extra mangling domains
    app.add_domain(NumpyPythonDomain)
    app.add_domain(NumpyCDomain)

#-----------------------------------------------------------------------------
# Docstring-mangling domains
#-----------------------------------------------------------------------------

try:
    import sphinx  # lazy to avoid test dependency
except ImportError:
    CDomain = PythonDomain = object
else:
    from sphinx.domains.c import CDomain
    from sphinx.domains.python import PythonDomain


class ManglingDomainBase(object):
    directive_mangling_map = {}

    def __init__(self, *a, **kw):
        super(ManglingDomainBase, self).__init__(*a, **kw)
        self.wrap_mangling_directives()

    def wrap_mangling_directives(self):
        for name, objtype in self.directive_mangling_map.items():
            self.directives[name] = wrap_mangling_directive(
                self.directives[name], objtype)


class NumpyPythonDomain(ManglingDomainBase, PythonDomain):
    name = 'np'
    directive_mangling_map = {
        'function': 'function',
        'class': 'class',
        'exception': 'class',
        'method': 'function',
        'classmethod': 'function',
        'staticmethod': 'function',
        'attribute': 'attribute',
    }


class NumpyCDomain(ManglingDomainBase, CDomain):
    name = 'np-c'
    directive_mangling_map = {
        'function': 'function',
        'member': 'attribute',
        'macro': 'function',
        'type': 'class',
        'var': 'object',
    }


def wrap_mangling_directive(base_directive, objtype):
    class directive(base_directive):
        def run(self):
            env = self.state.document.settings.env

            name = None
            if self.arguments:
                m = re.match(r'^(.*\s+)?(.*?)(\(.*)?', self.arguments[0])
                name = m.group(2).strip()

            if not name:
                name = self.arguments[0]

            lines = list(self.content)
            mangle_docstrings(env.app, objtype, name, None, None, lines)
            # local import to avoid testing dependency
            from docutils.statemachine import ViewList
            self.content = ViewList(lines, self.content.parent)

            return base_directive.run(self)

    return directive
