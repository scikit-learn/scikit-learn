from __future__ import annotations

import contextlib
import re
from typing import TYPE_CHECKING, ClassVar

from docutils import nodes
from docutils.parsers.rst import directives

from sphinx import addnodes
from sphinx.addnodes import desc_signature, pending_xref, pending_xref_condition
from sphinx.directives import ObjectDescription
from sphinx.domains.python._annotations import (
    _parse_annotation,
    _parse_arglist,
    _parse_type_list,
    _pseudo_parse_arglist,
    parse_reftarget,
)
from sphinx.locale import _
from sphinx.util import logging
from sphinx.util.docfields import Field, GroupedField, TypedField
from sphinx.util.nodes import (
    make_id,
)

if TYPE_CHECKING:
    from docutils.nodes import Node
    from docutils.parsers.rst.states import Inliner

    from sphinx.environment import BuildEnvironment
    from sphinx.util.typing import OptionSpec, TextlikeNode

logger = logging.getLogger(__name__)

# REs for Python signatures
py_sig_re = re.compile(
    r'''^ ([\w.]*\.)?            # class name(s)
          (\w+)  \s*             # thing name
          (?: \[\s*(.*)\s*])?    # optional: type parameters list
          (?: \(\s*(.*)\s*\)     # optional: arguments
           (?:\s* -> \s* (.*))?  #           return annotation
          )? $                   # and nothing more
          ''', re.VERBOSE)


# This override allows our inline type specifiers to behave like :class: link
# when it comes to handling "." and "~" prefixes.
class PyXrefMixin:
    def make_xref(
        self,
        rolename: str,
        domain: str,
        target: str,
        innernode: type[TextlikeNode] = nodes.emphasis,
        contnode: Node | None = None,
        env: BuildEnvironment | None = None,
        inliner: Inliner | None = None,
        location: Node | None = None,
    ) -> Node:
        # we use inliner=None to make sure we get the old behaviour with a single
        # pending_xref node
        result = super().make_xref(rolename, domain, target,  # type: ignore[misc]
                                   innernode, contnode,
                                   env, inliner=None, location=None)
        if isinstance(result, pending_xref):
            assert env is not None
            result['refspecific'] = True
            result['py:module'] = env.ref_context.get('py:module')
            result['py:class'] = env.ref_context.get('py:class')

            reftype, reftarget, reftitle, _ = parse_reftarget(target)
            if reftarget != reftitle:
                result['reftype'] = reftype
                result['reftarget'] = reftarget

                result.clear()
                result += innernode(reftitle, reftitle)  # type: ignore[call-arg]
            elif env.config.python_use_unqualified_type_names:
                children = result.children
                result.clear()

                shortname = target.split('.')[-1]
                textnode = innernode('', shortname)  # type: ignore[call-arg]
                contnodes = [pending_xref_condition('', '', textnode, condition='resolved'),
                             pending_xref_condition('', '', *children, condition='*')]
                result.extend(contnodes)

        return result

    _delimiters_re = re.compile(
        r'(\s*[\[\]\(\),](?:\s*o[rf]\s)?\s*|\s+o[rf]\s+|\s*\|\s*|\.\.\.)'
    )

    def make_xrefs(
        self,
        rolename: str,
        domain: str,
        target: str,
        innernode: type[TextlikeNode] = nodes.emphasis,
        contnode: Node | None = None,
        env: BuildEnvironment | None = None,
        inliner: Inliner | None = None,
        location: Node | None = None,
    ) -> list[Node]:
        sub_targets = self._delimiters_re.split(target)

        split_contnode = bool(contnode and contnode.astext() == target)

        in_literal = False
        results = []
        for sub_target in filter(None, sub_targets):
            if split_contnode:
                contnode = nodes.Text(sub_target)

            if in_literal or self._delimiters_re.match(sub_target):
                results.append(contnode or innernode(sub_target, sub_target))  # type: ignore[call-arg]
            else:
                results.append(self.make_xref(rolename, domain, sub_target,
                                              innernode, contnode, env, inliner, location))

            if sub_target in {'Literal', 'typing.Literal', '~typing.Literal'}:
                in_literal = True

        return results


class PyField(PyXrefMixin, Field):
    pass


class PyGroupedField(PyXrefMixin, GroupedField):
    pass


class PyTypedField(PyXrefMixin, TypedField):
    pass


class PyObject(ObjectDescription[tuple[str, str]]):
    """
    Description of a general Python object.

    :cvar allow_nesting: Class is an object that allows for nested namespaces
    :vartype allow_nesting: bool
    """

    option_spec: ClassVar[OptionSpec] = {
        'no-index': directives.flag,
        'no-index-entry': directives.flag,
        'no-contents-entry': directives.flag,
        'no-typesetting': directives.flag,
        'noindex': directives.flag,
        'noindexentry': directives.flag,
        'nocontentsentry': directives.flag,
        'single-line-parameter-list': directives.flag,
        'single-line-type-parameter-list': directives.flag,
        'module': directives.unchanged,
        'canonical': directives.unchanged,
        'annotation': directives.unchanged,
    }

    doc_field_types = [
        PyTypedField('parameter', label=_('Parameters'),
                     names=('param', 'parameter', 'arg', 'argument',
                            'keyword', 'kwarg', 'kwparam'),
                     typerolename='class', typenames=('paramtype', 'type'),
                     can_collapse=True),
        PyTypedField('variable', label=_('Variables'),
                     names=('var', 'ivar', 'cvar'),
                     typerolename='class', typenames=('vartype',),
                     can_collapse=True),
        PyGroupedField('exceptions', label=_('Raises'), rolename='exc',
                       names=('raises', 'raise', 'exception', 'except'),
                       can_collapse=True),
        Field('returnvalue', label=_('Returns'), has_arg=False,
              names=('returns', 'return')),
        PyField('returntype', label=_('Return type'), has_arg=False,
                names=('rtype',), bodyrolename='class'),
    ]

    allow_nesting = False

    def get_signature_prefix(self, sig: str) -> list[nodes.Node]:
        """May return a prefix to put before the object name in the
        signature.
        """
        return []

    def needs_arglist(self) -> bool:
        """May return true if an empty argument list is to be generated even if
        the document contains none.
        """
        return False

    def handle_signature(self, sig: str, signode: desc_signature) -> tuple[str, str]:
        """Transform a Python signature into RST nodes.

        Return (fully qualified name of the thing, classname if any).

        If inside a class, the current class name is handled intelligently:
        * it is stripped from the displayed name if present
        * it is added to the full name (return value) if not present
        """
        m = py_sig_re.match(sig)
        if m is None:
            raise ValueError
        prefix, name, tp_list, arglist, retann = m.groups()

        # determine module and class name (if applicable), as well as full name
        modname = self.options.get('module', self.env.ref_context.get('py:module'))
        classname = self.env.ref_context.get('py:class')
        if classname:
            add_module = False
            if prefix and (prefix == classname or
                           prefix.startswith(classname + ".")):
                fullname = prefix + name
                # class name is given again in the signature
                prefix = prefix[len(classname):].lstrip('.')
            elif prefix:
                # class name is given in the signature, but different
                # (shouldn't happen)
                fullname = classname + '.' + prefix + name
            else:
                # class name is not given in the signature
                fullname = classname + '.' + name
        else:
            add_module = True
            if prefix:
                classname = prefix.rstrip('.')
                fullname = prefix + name
            else:
                classname = ''
                fullname = name

        signode['module'] = modname
        signode['class'] = classname
        signode['fullname'] = fullname

        max_len = (self.env.config.python_maximum_signature_line_length
                   or self.env.config.maximum_signature_line_length
                   or 0)

        # determine if the function arguments (without its type parameters)
        # should be formatted on a multiline or not by removing the width of
        # the type parameters list (if any)
        sig_len = len(sig)
        tp_list_span = m.span(3)
        multi_line_parameter_list = (
            'single-line-parameter-list' not in self.options
            and (sig_len - (tp_list_span[1] - tp_list_span[0])) > max_len > 0
        )

        # determine whether the type parameter list must be wrapped or not
        arglist_span = m.span(4)
        multi_line_type_parameter_list = (
            'single-line-type-parameter-list' not in self.options
            and (sig_len - (arglist_span[1] - arglist_span[0])) > max_len > 0
        )

        sig_prefix = self.get_signature_prefix(sig)
        if sig_prefix:
            if type(sig_prefix) is str:
                msg = ("Python directive method get_signature_prefix()"
                       " must return a list of nodes."
                       f" Return value was '{sig_prefix}'.")
                raise TypeError(msg)
            signode += addnodes.desc_annotation(str(sig_prefix), '', *sig_prefix)

        if prefix:
            signode += addnodes.desc_addname(prefix, prefix)
        elif modname and add_module and self.env.config.add_module_names:
            nodetext = modname + '.'
            signode += addnodes.desc_addname(nodetext, nodetext)

        signode += addnodes.desc_name(name, name)

        if tp_list:
            try:
                signode += _parse_type_list(tp_list, self.env, multi_line_type_parameter_list)
            except Exception as exc:
                logger.warning("could not parse tp_list (%r): %s", tp_list, exc,
                               location=signode)

        if arglist:
            try:
                signode += _parse_arglist(arglist, self.env, multi_line_parameter_list)
            except SyntaxError:
                # fallback to parse arglist original parser
                # (this may happen if the argument list is incorrectly used
                # as a list of bases when documenting a class)
                # it supports to represent optional arguments (ex. "func(foo [, bar])")
                _pseudo_parse_arglist(signode, arglist, multi_line_parameter_list)
            except (NotImplementedError, ValueError) as exc:
                # duplicated parameter names raise ValueError and not a SyntaxError
                logger.warning("could not parse arglist (%r): %s", arglist, exc,
                               location=signode)
                _pseudo_parse_arglist(signode, arglist, multi_line_parameter_list)
        else:
            if self.needs_arglist():
                # for callables, add an empty parameter list
                signode += addnodes.desc_parameterlist()

        if retann:
            children = _parse_annotation(retann, self.env)
            signode += addnodes.desc_returns(retann, '', *children)

        anno = self.options.get('annotation')
        if anno:
            signode += addnodes.desc_annotation(' ' + anno, '',
                                                addnodes.desc_sig_space(),
                                                nodes.Text(anno))

        return fullname, prefix

    def _object_hierarchy_parts(self, sig_node: desc_signature) -> tuple[str, ...]:
        if 'fullname' not in sig_node:
            return ()
        modname = sig_node.get('module')
        fullname = sig_node['fullname']

        if modname:
            return (modname, *fullname.split('.'))
        else:
            return tuple(fullname.split('.'))

    def get_index_text(self, modname: str, name: tuple[str, str]) -> str:
        """Return the text for the index entry of the object."""
        msg = 'must be implemented in subclasses'
        raise NotImplementedError(msg)

    def add_target_and_index(self, name_cls: tuple[str, str], sig: str,
                             signode: desc_signature) -> None:
        modname = self.options.get('module', self.env.ref_context.get('py:module'))
        fullname = (modname + '.' if modname else '') + name_cls[0]
        node_id = make_id(self.env, self.state.document, '', fullname)
        signode['ids'].append(node_id)
        self.state.document.note_explicit_target(signode)

        domain = self.env.domains.python_domain
        domain.note_object(fullname, self.objtype, node_id, location=signode)

        canonical_name = self.options.get('canonical')
        if canonical_name:
            domain.note_object(canonical_name, self.objtype, node_id, aliased=True,
                               location=signode)

        if 'no-index-entry' not in self.options:
            indextext = self.get_index_text(modname, name_cls)
            if indextext:
                self.indexnode['entries'].append(('single', indextext, node_id, '', None))

    def before_content(self) -> None:
        """Handle object nesting before content

        :py:class:`PyObject` represents Python language constructs. For
        constructs that are nestable, such as a Python classes, this method will
        build up a stack of the nesting hierarchy so that it can be later
        de-nested correctly, in :py:meth:`after_content`.

        For constructs that aren't nestable, the stack is bypassed, and instead
        only the most recent object is tracked. This object prefix name will be
        removed with :py:meth:`after_content`.
        """
        prefix = None
        if self.names:
            # fullname and name_prefix come from the `handle_signature` method.
            # fullname represents the full object name that is constructed using
            # object nesting and explicit prefixes. `name_prefix` is the
            # explicit prefix given in a signature
            (fullname, name_prefix) = self.names[-1]
            if self.allow_nesting:
                prefix = fullname
            elif name_prefix:
                prefix = name_prefix.strip('.')
        if prefix:
            self.env.ref_context['py:class'] = prefix
            if self.allow_nesting:
                classes = self.env.ref_context.setdefault('py:classes', [])
                classes.append(prefix)
        if 'module' in self.options:
            modules = self.env.ref_context.setdefault('py:modules', [])
            modules.append(self.env.ref_context.get('py:module'))
            self.env.ref_context['py:module'] = self.options['module']

    def after_content(self) -> None:
        """Handle object de-nesting after content

        If this class is a nestable object, removing the last nested class prefix
        ends further nesting in the object.

        If this class is not a nestable object, the list of classes should not
        be altered as we didn't affect the nesting levels in
        :py:meth:`before_content`.
        """
        classes = self.env.ref_context.setdefault('py:classes', [])
        if self.allow_nesting:
            with contextlib.suppress(IndexError):
                classes.pop()

        self.env.ref_context['py:class'] = (classes[-1] if len(classes) > 0
                                            else None)
        if 'module' in self.options:
            modules = self.env.ref_context.setdefault('py:modules', [])
            if modules:
                self.env.ref_context['py:module'] = modules.pop()
            else:
                self.env.ref_context.pop('py:module')

    def _toc_entry_name(self, sig_node: desc_signature) -> str:
        if not sig_node.get('_toc_parts'):
            return ''

        config = self.env.app.config
        objtype = sig_node.parent.get('objtype')
        if config.add_function_parentheses and objtype in {'function', 'method'}:
            parens = '()'
        else:
            parens = ''
        *parents, name = sig_node['_toc_parts']
        if config.toc_object_entries_show_parents == 'domain':
            return sig_node.get('fullname', name) + parens
        if config.toc_object_entries_show_parents == 'hide':
            return name + parens
        if config.toc_object_entries_show_parents == 'all':
            return '.'.join([*parents, name + parens])
        return ''
