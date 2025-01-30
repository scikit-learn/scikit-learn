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

.. [1] https://github.com/numpy/numpydoc

"""

import hashlib
import inspect
import itertools
import pydoc
import re
from collections.abc import Callable
from copy import deepcopy

from docutils.nodes import Text, citation, comment, inline, reference, section
from sphinx.addnodes import desc_content, pending_xref
from sphinx.util import logging

from . import __version__
from .docscrape_sphinx import get_doc_object
from .validate import get_validation_checks, validate
from .xref import DEFAULT_LINKS

logger = logging.getLogger(__name__)

HASH_LEN = 12


def _traverse_or_findall(node, condition, **kwargs):
    """Triage node.traverse (docutils <0.18.1) vs node.findall.

    TODO: This check can be removed when the minimum supported docutils version
    for numpydoc is docutils>=0.18.1
    """
    return (
        node.findall(condition, **kwargs)
        if hasattr(node, "findall")
        else node.traverse(condition, **kwargs)
    )


def rename_references(app, what, name, obj, options, lines):
    # decorate reference numbers so that there are no duplicates
    # these are later undecorated in the doctree, in relabel_references
    references = set()
    for line in lines:
        line = line.strip()
        m = re.match(
            r"^\.\. +\[(%s)\]" % app.config.numpydoc_citation_re, line, re.IGNORECASE
        )
        if m:
            references.add(m.group(1))

    if references:
        # we use a hash to mangle the reference name to avoid invalid names
        sha = hashlib.sha256()
        sha.update(name.encode("utf8"))
        prefix = "R" + sha.hexdigest()[:HASH_LEN]

        for r in references:
            new_r = prefix + "-" + r
            for i, line in enumerate(lines):
                lines[i] = lines[i].replace(f"[{r}]_", f"[{new_r}]_")
                lines[i] = lines[i].replace(f".. [{r}]", f".. [{new_r}]")


def _is_cite_in_numpydoc_docstring(citation_node):
    # Find DEDUPLICATION_TAG in comment as last node of sibling section

    # XXX: I failed to use citation_node.traverse to do this:
    section_node = citation_node.parent

    def is_docstring_section(node):
        return isinstance(node, (section, desc_content))

    while not is_docstring_section(section_node):
        section_node = section_node.parent
        if section_node is None:
            return False

    sibling_sections = itertools.chain(
        _traverse_or_findall(
            section_node,
            is_docstring_section,
            include_self=True,
            descend=False,
            siblings=True,
        )
    )
    for sibling_section in sibling_sections:
        if not sibling_section.children:
            continue

        for child in sibling_section.children[::-1]:
            if not isinstance(child, comment):
                continue

            if child.rawsource.strip() == DEDUPLICATION_TAG.strip():
                return True

    return False


def relabel_references(app, doc):
    # Change 'hash-ref' to 'ref' in label text
    for citation_node in _traverse_or_findall(doc, citation):
        if not _is_cite_in_numpydoc_docstring(citation_node):
            continue
        label_node = citation_node[0]
        prefix, _, new_label = label_node[0].astext().partition("-")
        assert len(prefix) == HASH_LEN + 1
        new_text = Text(new_label)
        label_node.replace(label_node[0], new_text)

        for id_ in citation_node["backrefs"]:
            ref = doc.ids[id_]
            ref_text = ref[0]

            # Sphinx has created pending_xref nodes with [reftext] text.
            def matching_pending_xref(node):
                return (
                    isinstance(node, pending_xref)
                    and node[0].astext() == f"[{ref_text}]"
                )

            for xref_node in _traverse_or_findall(ref.parent, matching_pending_xref):
                xref_node.replace(xref_node[0], Text(f"[{new_text}]"))
            ref.replace(ref_text, new_text.copy())


def clean_backrefs(app, doc, docname):
    # only::latex directive has resulted in citation backrefs without reference
    known_ref_ids = set()
    for ref in _traverse_or_findall(doc, reference, descend=True):
        for id_ in ref["ids"]:
            known_ref_ids.add(id_)
    # some extensions produce backrefs to inline elements
    for ref in _traverse_or_findall(doc, inline, descend=True):
        for id_ in ref["ids"]:
            known_ref_ids.add(id_)
    for citation_node in _traverse_or_findall(doc, citation, descend=True):
        # remove backrefs to non-existent refs
        citation_node["backrefs"] = [
            id_ for id_ in citation_node["backrefs"] if id_ in known_ref_ids
        ]


DEDUPLICATION_TAG = "    !! processed by numpydoc !!"


def mangle_docstrings(app, what, name, obj, options, lines):
    if DEDUPLICATION_TAG in lines:
        return
    show_inherited_class_members = app.config.numpydoc_show_inherited_class_members
    if isinstance(show_inherited_class_members, dict):
        try:
            show_inherited_class_members = show_inherited_class_members[name]
        except KeyError:
            show_inherited_class_members = True

    cfg = {
        "use_plots": app.config.numpydoc_use_plots,
        "show_class_members": app.config.numpydoc_show_class_members,
        "show_inherited_class_members": show_inherited_class_members,
        "class_members_toctree": app.config.numpydoc_class_members_toctree,
        "attributes_as_param_list": app.config.numpydoc_attributes_as_param_list,
        "xref_param_type": app.config.numpydoc_xref_param_type,
        "xref_aliases": app.config.numpydoc_xref_aliases_complete,
        "xref_ignore": app.config.numpydoc_xref_ignore,
    }

    cfg.update(options or {})
    u_NL = "\n"
    if what == "module":
        # Strip top title
        pattern = "^\\s*[#*=]{4,}\\n[a-z0-9 -]+\\n[#*=]{4,}\\s*"
        title_re = re.compile(pattern, re.IGNORECASE | re.DOTALL)
        lines[:] = title_re.sub("", u_NL.join(lines)).split(u_NL)
    else:
        try:
            doc = get_doc_object(
                obj, what, u_NL.join(lines), config=cfg, builder=app.builder
            )
            lines[:] = str(doc).split(u_NL)
        except Exception:
            logger.error("[numpydoc] While processing docstring for %r", name)
            raise

        if app.config.numpydoc_validation_checks:
            # If the user has supplied patterns to ignore via the
            # numpydoc_validation_exclude config option, skip validation for
            # any objs whose name matches any of the patterns
            excluder = app.config.numpydoc_validation_excluder
            exclude_from_validation = excluder.search(name) if excluder else False
            if not exclude_from_validation:
                # TODO: Currently, all validation checks are run and only those
                # selected via config are reported. It would be more efficient to
                # only run the selected checks.
                report = validate(doc)
                errors = [
                    err
                    for err in report["errors"]
                    if not (
                        (
                            overrides := app.config.numpydoc_validation_overrides.get(
                                err[0]
                            )
                        )
                        and re.search(overrides, report["docstring"])
                    )
                ]
                if {err[0] for err in errors} & app.config.numpydoc_validation_checks:
                    msg = (
                        f"[numpydoc] Validation warnings while processing "
                        f"docstring for {name!r}:\n"
                    )
                    for err in errors:
                        if err[0] in app.config.numpydoc_validation_checks:
                            msg += f"  {err[0]}: {err[1]}\n"
                    logger.warning(msg)

    # call function to replace reference numbers so that there are no
    # duplicates
    rename_references(app, what, name, obj, options, lines)

    lines += ["..", DEDUPLICATION_TAG]


def mangle_signature(app, what, name, obj, options, sig, retann):
    # Do not try to inspect classes that don't define `__init__`
    if inspect.isclass(obj) and (
        not hasattr(obj, "__init__")
        or "initializes x; see " in pydoc.getdoc(obj.__init__)
    ):
        return "", ""

    if not (isinstance(obj, Callable) or hasattr(obj, "__argspec_is_invalid_")):
        return None

    if not hasattr(obj, "__doc__"):
        return None
    doc = get_doc_object(obj, config={"show_class_members": False})
    sig = doc["Signature"] or _clean_text_signature(
        getattr(obj, "__text_signature__", None)
    )
    if sig:
        sig = re.sub("^[^(]*", "", sig)
        return sig, ""


def _clean_text_signature(sig):
    if sig is None:
        return None
    start_pattern = re.compile(r"^[^(]*\(")
    start, end = start_pattern.search(sig).span()
    start_sig = sig[start:end]
    sig = sig[end:-1]
    sig = re.sub(r"^\$(self|module|type)(,\s|$)", "", sig, count=1)
    sig = re.sub(r"(^|(?<=,\s))/,\s\*", "*", sig, count=1)
    return start_sig + sig + ")"


def setup(app, get_doc_object_=get_doc_object):
    if not hasattr(app, "add_config_value"):
        return None  # probably called by nose, better bail out

    global get_doc_object
    get_doc_object = get_doc_object_

    app.setup_extension("sphinx.ext.autosummary")
    app.connect("config-inited", update_config)
    app.connect("autodoc-process-docstring", mangle_docstrings)
    app.connect("autodoc-process-signature", mangle_signature)
    app.connect("doctree-read", relabel_references)
    app.connect("doctree-resolved", clean_backrefs)
    app.add_config_value("numpydoc_use_plots", None, False)
    app.add_config_value("numpydoc_show_class_members", True, True)
    app.add_config_value(
        "numpydoc_show_inherited_class_members", True, True, types=(bool, dict)
    )
    app.add_config_value("numpydoc_class_members_toctree", True, True)
    app.add_config_value("numpydoc_citation_re", "[a-z0-9_.-]+", True)
    app.add_config_value("numpydoc_attributes_as_param_list", True, True)
    app.add_config_value("numpydoc_xref_param_type", False, True)
    app.add_config_value("numpydoc_xref_aliases", dict(), True)
    app.add_config_value("numpydoc_xref_ignore", set(), True)
    app.add_config_value("numpydoc_validation_checks", set(), True)
    app.add_config_value("numpydoc_validation_exclude", set(), False)
    app.add_config_value("numpydoc_validation_overrides", dict(), False)

    # Extra mangling domains
    app.add_domain(NumpyPythonDomain)
    app.add_domain(NumpyCDomain)

    metadata = {"version": __version__, "parallel_read_safe": True}
    return metadata


def update_config(app, config=None):
    """Update the configuration with default values."""
    if config is None:  # needed for testing and old Sphinx
        config = app.config
    # Do not simply overwrite the `app.config.numpydoc_xref_aliases`
    # otherwise the next sphinx-build will compare the incoming values (without
    # our additions) to the old values (with our additions) and trigger
    # a full rebuild!
    numpydoc_xref_aliases_complete = deepcopy(config.numpydoc_xref_aliases)
    for key, value in DEFAULT_LINKS.items():
        if key not in numpydoc_xref_aliases_complete:
            numpydoc_xref_aliases_complete[key] = value
    config.numpydoc_xref_aliases_complete = numpydoc_xref_aliases_complete

    # Processing to determine whether numpydoc_validation_checks is treated
    # as a blocklist or allowlist
    config.numpydoc_validation_checks = get_validation_checks(
        config.numpydoc_validation_checks
    )

    # Generate the regexp for docstrings to ignore during validation
    if isinstance(config.numpydoc_validation_exclude, str):
        raise ValueError(
            f"numpydoc_validation_exclude must be a container of strings, "
            f"e.g. [{config.numpydoc_validation_exclude!r}]."
        )
    config.numpydoc_validation_excluder = None
    if config.numpydoc_validation_exclude:
        exclude_expr = re.compile(
            r"|".join(exp for exp in config.numpydoc_validation_exclude)
        )
        config.numpydoc_validation_excluder = exclude_expr

    for check, patterns in config.numpydoc_validation_overrides.items():
        config.numpydoc_validation_overrides[check] = re.compile(
            r"|".join(exp for exp in patterns)
        )


# ------------------------------------------------------------------------------
# Docstring-mangling domains
# ------------------------------------------------------------------------------

from docutils.statemachine import ViewList
from sphinx.domains.c import CDomain
from sphinx.domains.python import PythonDomain


class ManglingDomainBase:
    directive_mangling_map = {}

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.wrap_mangling_directives()

    def wrap_mangling_directives(self):
        for name, objtype in list(self.directive_mangling_map.items()):
            self.directives[name] = wrap_mangling_directive(
                self.directives[name], objtype
            )


class NumpyPythonDomain(ManglingDomainBase, PythonDomain):
    name = "np"
    directive_mangling_map = {
        "function": "function",
        "class": "class",
        "exception": "class",
        "method": "function",
        "classmethod": "function",
        "staticmethod": "function",
        "attribute": "attribute",
    }
    indices = []


class NumpyCDomain(ManglingDomainBase, CDomain):
    name = "np-c"
    directive_mangling_map = {
        "function": "function",
        "member": "attribute",
        "macro": "function",
        "type": "class",
        "var": "object",
    }


def match_items(lines, content_old):
    """Create items for mangled lines.

    This function tries to match the lines in ``lines`` with the items (source
    file references and line numbers) in ``content_old``. The
    ``mangle_docstrings`` function changes the actual docstrings, but doesn't
    keep track of where each line came from. The mangling does many operations
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
    >>> lines = ["", "A", "", "B", "   ", "", "C", "D"]
    >>> lines_old = ["a", "", "", "b", "", "c"]
    >>> items_old = [
    ...     ("file1.py", 0),
    ...     ("file1.py", 1),
    ...     ("file1.py", 2),
    ...     ("file2.py", 0),
    ...     ("file2.py", 1),
    ...     ("file2.py", 2),
    ... ]
    >>> content_old = ViewList(lines_old, items=items_old)
    >>> match_items(lines, content_old)  # doctest: +NORMALIZE_WHITESPACE
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
    assert len(items_new) == len(lines)
    return items_new


def wrap_mangling_directive(base_directive, objtype):
    class directive(base_directive):
        def run(self):
            env = self.state.document.settings.env

            name = None
            if self.arguments:
                m = re.match(r"^(.*\s+)?(.*?)(\(.*)?", self.arguments[0])
                name = m.group(2).strip()

            if not name:
                name = self.arguments[0]

            lines = list(self.content)
            mangle_docstrings(env.app, objtype, name, None, None, lines)
            if self.content:
                items = match_items(lines, self.content)
                self.content = ViewList(lines, items=items, parent=self.content.parent)

            return base_directive.run(self)

    return directive
