"""
========
numpydoc
========

Sphinx extension that handles docstrings in the Numpy standard format. [1]

It will:

- Convert Parameters etc. sections to field lists.
- Convert See Also section to a See also entry.
- Renumber references.
- Extract the signature from the docstring, if it can't be determined otherwise.

.. [1] http://projects.scipy.org/scipy/numpy/wiki/CodingStyleGuidelines#docstring-standard

"""

import os, re, pydoc
from docscrape_sphinx import get_doc_object, SphinxDocString
import inspect

def mangle_docstrings(app, what, name, obj, options, lines,
                      reference_offset=[0]):
    if what == 'module':
        # Strip top title
        title_re = re.compile(r'^\s*[#*=]{4,}\n[a-z0-9 -]+\n[#*=]{4,}\s*',
                              re.I|re.S)
        lines[:] = title_re.sub('', "\n".join(lines)).split("\n")
    else:
        doc = get_doc_object(obj, what)
        lines[:] = str(doc).split("\n")

    if app.config.numpydoc_edit_link and hasattr(obj, '__name__') and \
           obj.__name__:
        v = dict(full_name=obj.__name__)
        lines += [''] + (app.config.numpydoc_edit_link % v).split("\n")

    # replace reference numbers so that there are no duplicates
    references = []
    for l in lines:
        l = l.strip()
        if l.startswith('.. ['):
            try:
                references.append(int(l[len('.. ['):l.index(']')]))
            except ValueError:
                print "WARNING: invalid reference in %s docstring" % name

    # Start renaming from the biggest number, otherwise we may
    # overwrite references.
    references.sort()
    if references:
        for i, line in enumerate(lines):
            for r in references:
                new_r = reference_offset[0] + r
                lines[i] = lines[i].replace('[%d]_' % r,
                                            '[%d]_' % new_r)
                lines[i] = lines[i].replace('.. [%d]' % r,
                                            '.. [%d]' % new_r)

    reference_offset[0] += len(references)

def mangle_signature(app, what, name, obj, options, sig, retann):
    # Do not try to inspect classes that don't define `__init__`
    if (inspect.isclass(obj) and
        'initializes x; see ' in pydoc.getdoc(obj.__init__)):
        return '', ''

    if not (callable(obj) or hasattr(obj, '__argspec_is_invalid_')): return
    if not hasattr(obj, '__doc__'): return

    doc = SphinxDocString(pydoc.getdoc(obj))
    if doc['Signature']:
        sig = re.sub("^[^(]*", "", doc['Signature'])
        return sig, ''

def initialize(app):
    try:
        app.connect('autodoc-process-signature', mangle_signature)
    except:
        monkeypatch_sphinx_ext_autodoc()

def setup(app, get_doc_object_=get_doc_object):
    global get_doc_object
    get_doc_object = get_doc_object_

    app.connect('autodoc-process-docstring', mangle_docstrings)
    app.connect('builder-inited', initialize)
    app.add_config_value('numpydoc_edit_link', None, True)

#------------------------------------------------------------------------------
# Monkeypatch sphinx.ext.autodoc to accept argspecless autodocs (Sphinx < 0.5)
#------------------------------------------------------------------------------

def monkeypatch_sphinx_ext_autodoc():
    global _original_format_signature
    import sphinx.ext.autodoc

    if sphinx.ext.autodoc.format_signature is our_format_signature:
        return

    print "[numpydoc] Monkeypatching sphinx.ext.autodoc ..."
    _original_format_signature = sphinx.ext.autodoc.format_signature
    sphinx.ext.autodoc.format_signature = our_format_signature

def our_format_signature(what, obj):
    r = mangle_signature(None, what, None, obj, None, None, None)
    if r is not None:
        return r[0]
    else:
        return _original_format_signature(what, obj)
