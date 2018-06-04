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

.. [1] https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

"""
from __future__ import division, absolute_import, print_function

import sys
import re
import pydoc
import inspect
import collections
import hashlib

from docutils.nodes import citation, Text
import sphinx
from sphinx.addnodes import pending_xref, desc_content

if sphinx.__version__ < '1.0.1':
    raise RuntimeError("Sphinx 1.0.1 or newer is required")

from .docscrape_sphinx import get_doc_object, SphinxDocString
from . import __version__

if sys.version_info[0] >= 3:
    sixu = lambda s: s
else:
    sixu = lambda s: unicode(s, 'unicode_escape')


HASH_LEN = 12


def rename_references(app, what, name, obj, options, lines):
    # decorate reference numbers so that there are no duplicates
    # these are later undecorated in the doctree, in relabel_references
    references = set()
    for line in lines:
        line = line.strip()
        m = re.match(sixu('^.. \\[(%s)\\]') % app.config.numpydoc_citation_re,
                     line, re.I)
        if m:
            references.add(m.group(1))

    if references:
        # we use a hash to mangle the reference name to avoid invalid names
        sha = hashlib.sha256()
        sha.update(name.encode('utf8'))
        prefix = 'R' + sha.hexdigest()[:HASH_LEN]

        for r in references:
            new_r = prefix + '-' + r
            for i, line in enumerate(lines):
                lines[i] = lines[i].replace(sixu('[%s]_') % r,
                                            sixu('[%s]_') % new_r)
                lines[i] = lines[i].replace(sixu('.. [%s]') % r,
                                            sixu('.. [%s]') % new_r)


def _ascend(node, cls):
    while node and not isinstance(node, cls):
        node = node.parent
    return node


def relabel_references(app, doc):
    # Change 'hash-ref' to 'ref' in label text
    for citation_node in doc.traverse(citation):
        if _ascend(citation_node, desc_content) is None:
            # no desc node in ancestry -> not in a docstring
            # XXX: should we also somehow check it's in a References section?
            continue
        label_node = citation_node[0]
        prefix, _, new_label = label_node[0].astext().partition('-')
        assert len(prefix) == HASH_LEN + 1
        new_text = Text(new_label)
        label_node.replace(label_node[0], new_text)

        for id in citation_node['backrefs']:
            ref = doc.ids[id]
            ref_text = ref[0]

            # Sphinx has created pending_xref nodes with [reftext] text.
            def matching_pending_xref(node):
                return (isinstance(node, pending_xref) and
                        node[0].astext() == '[%s]' % ref_text)

            for xref_node in ref.parent.traverse(matching_pending_xref):
                xref_node.replace(xref_node[0], Text('[%s]' % new_text))
            ref.replace(ref_text, new_text.copy())


DEDUPLICATION_TAG = '    !! processed by numpydoc !!'


def mangle_docstrings(app, what, name, obj, options, lines):
    if DEDUPLICATION_TAG in lines:
        return

    cfg = {'use_plots': app.config.numpydoc_use_plots,
           'use_blockquotes': app.config.numpydoc_use_blockquotes,
           'show_class_members': app.config.numpydoc_show_class_members,
           'show_inherited_class_members':
           app.config.numpydoc_show_inherited_class_members,
           'class_members_toctree': app.config.numpydoc_class_members_toctree}

    u_NL = sixu('\n')
    if what == 'module':
        # Strip top title
        pattern = '^\\s*[#*=]{4,}\\n[a-z0-9 -]+\\n[#*=]{4,}\\s*'
        title_re = re.compile(sixu(pattern), re.I | re.S)
        lines[:] = title_re.sub(sixu(''), u_NL.join(lines)).split(u_NL)
    else:
        doc = get_doc_object(obj, what, u_NL.join(lines), config=cfg,
                             builder=app.builder)
        if sys.version_info[0] >= 3:
            doc = str(doc)
        else:
            doc = unicode(doc)
        lines[:] = doc.split(u_NL)

    if (app.config.numpydoc_edit_link and hasattr(obj, '__name__') and
            obj.__name__):
        if hasattr(obj, '__module__'):
            v = dict(full_name=sixu("%s.%s") % (obj.__module__, obj.__name__))
        else:
            v = dict(full_name=obj.__name__)
        lines += [sixu(''), sixu('.. htmlonly::'), sixu('')]
        lines += [sixu('    %s') % x for x in
                  (app.config.numpydoc_edit_link % v).split("\n")]

    # call function to replace reference numbers so that there are no
    # duplicates
    rename_references(app, what, name, obj, options, lines)

    lines += ['..', DEDUPLICATION_TAG]


def mangle_signature(app, what, name, obj, options, sig, retann):
    # Do not try to inspect classes that don't define `__init__`
    if (inspect.isclass(obj) and
        (not hasattr(obj, '__init__') or
            'initializes x; see ' in pydoc.getdoc(obj.__init__))):
        return '', ''

    if not (isinstance(obj, collections.Callable) or
            hasattr(obj, '__argspec_is_invalid_')):
        return

    if not hasattr(obj, '__doc__'):
        return
    doc = SphinxDocString(pydoc.getdoc(obj))
    sig = doc['Signature'] or getattr(obj, '__text_signature__', None)
    if sig:
        sig = re.sub(sixu("^[^(]*"), sixu(""), sig)
        return sig, sixu('')


def setup(app, get_doc_object_=get_doc_object):
    if not hasattr(app, 'add_config_value'):
        return  # probably called by nose, better bail out

    global get_doc_object
    get_doc_object = get_doc_object_

    app.connect('autodoc-process-docstring', mangle_docstrings)
    app.connect('autodoc-process-signature', mangle_signature)
    app.connect('doctree-read', relabel_references)
    app.add_config_value('numpydoc_edit_link', None, False)
    app.add_config_value('numpydoc_use_plots', None, False)
    app.add_config_value('numpydoc_use_blockquotes', None, False)
    app.add_config_value('numpydoc_show_class_members', True, True)
    app.add_config_value('numpydoc_show_inherited_class_members', True, True)
    app.add_config_value('numpydoc_class_members_toctree', True, True)
    app.add_config_value('numpydoc_citation_re', '[a-z0-9_.-]+', True)

    # Extra mangling domains
    app.add_domain(NumpyPythonDomain)
    app.add_domain(NumpyCDomain)

    app.setup_extension('sphinx.ext.autosummary')

    metadata = {'version': __version__,
                'parallel_read_safe': True}
    return metadata

# ------------------------------------------------------------------------------
# Docstring-mangling domains
# ------------------------------------------------------------------------------

from docutils.statemachine import ViewList
from sphinx.domains.c import CDomain
from sphinx.domains.python import PythonDomain


class ManglingDomainBase(object):
    directive_mangling_map = {}

    def __init__(self, *a, **kw):
        super(ManglingDomainBase, self).__init__(*a, **kw)
        self.wrap_mangling_directives()

    def wrap_mangling_directives(self):
        for name, objtype in list(self.directive_mangling_map.items()):
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
    indices = []


class NumpyCDomain(ManglingDomainBase, CDomain):
    name = 'np-c'
    directive_mangling_map = {
        'function': 'function',
        'member': 'attribute',
        'macro': 'function',
        'type': 'class',
        'var': 'object',
    }


def match_items(lines, content_old):
    """Create items for mangled lines.

    This function tries to match the lines in ``lines`` with the items (source
    file references and line numbers) in ``content_old``. The
    ``mangle_docstrings`` function changes the actual docstrings, but doesn't
    keep track of where each line came from. The manging does many operations
    on the original lines, which are hard to track afterwards.

    Many of the line changes come from deleting or inserting blank lines. This
    function tries to match lines by ignoring blank lines. All other changes
    (such as inserting figures or changes in the references) are completely
    ignored, so the generated line numbers will be off if ``mangle_docstrings``
    does anything non-trivial.

    This is a best-effort function and the real fix would be to make
    ``mangle_docstrings`` actually keep track of the ``items`` together with
    the ``lines``.

    Examples
    --------
    >>> lines = ['', 'A', '', 'B', '   ', '', 'C', 'D']
    >>> lines_old = ['a', '', '', 'b', '', 'c']
    >>> items_old = [('file1.py', 0), ('file1.py', 1), ('file1.py', 2),
    ...              ('file2.py', 0), ('file2.py', 1), ('file2.py', 2)]
    >>> content_old = ViewList(lines_old, items=items_old)
    >>> match_items(lines, content_old) # doctest: +NORMALIZE_WHITESPACE
    [('file1.py', 0), ('file1.py', 0), ('file2.py', 0), ('file2.py', 0),
     ('file2.py', 2), ('file2.py', 2), ('file2.py', 2), ('file2.py', 2)]
    >>> # first 2 ``lines`` are matched to 'a', second 2 to 'b', rest to 'c'
    >>> # actual content is completely ignored.

    Notes
    -----
    The algorithm tries to match any line in ``lines`` with one in
    ``lines_old``.  It skips over all empty lines in ``lines_old`` and assigns
    this line number to all lines in ``lines``, unless a non-empty line is
    found in ``lines`` in which case it goes to the next line in ``lines_old``.

    """
    items_new = []
    lines_old = content_old.data
    items_old = content_old.items
    j = 0
    for i, line in enumerate(lines):
        # go to next non-empty line in old:
        # line.strip() checks whether the string is all whitespace
        while j < len(lines_old) - 1 and not lines_old[j].strip():
            j += 1
        items_new.append(items_old[j])
        if line.strip() and j < len(lines_old) - 1:
            j += 1
    assert(len(items_new) == len(lines))
    return items_new


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
            if self.content:
                items = match_items(lines, self.content)
                self.content = ViewList(lines, items=items,
                                        parent=self.content.parent)

            return base_directive.run(self)

    return directive
