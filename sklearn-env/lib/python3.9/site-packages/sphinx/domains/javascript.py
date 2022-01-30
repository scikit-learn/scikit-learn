"""
    sphinx.domains.javascript
    ~~~~~~~~~~~~~~~~~~~~~~~~~

    The JavaScript domain.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from typing import Any, Dict, Iterator, List, Optional, Tuple, cast

from docutils import nodes
from docutils.nodes import Element, Node
from docutils.parsers.rst import directives

from sphinx import addnodes
from sphinx.addnodes import desc_signature, pending_xref
from sphinx.application import Sphinx
from sphinx.builders import Builder
from sphinx.directives import ObjectDescription
from sphinx.domains import Domain, ObjType
from sphinx.domains.python import _pseudo_parse_arglist
from sphinx.environment import BuildEnvironment
from sphinx.locale import _, __
from sphinx.roles import XRefRole
from sphinx.util import logging
from sphinx.util.docfields import Field, GroupedField, TypedField
from sphinx.util.docutils import SphinxDirective
from sphinx.util.nodes import make_id, make_refnode
from sphinx.util.typing import OptionSpec

logger = logging.getLogger(__name__)


class JSObject(ObjectDescription[Tuple[str, str]]):
    """
    Description of a JavaScript object.
    """
    #: If set to ``True`` this object is callable and a `desc_parameterlist` is
    #: added
    has_arguments = False

    #: If ``allow_nesting`` is ``True``, the object prefixes will be accumulated
    #: based on directive nesting
    allow_nesting = False

    option_spec: OptionSpec = {
        'noindex': directives.flag,
        'noindexentry': directives.flag,
    }

    def get_display_prefix(self) -> List[Node]:
        #: what is displayed right before the documentation entry
        return []

    def handle_signature(self, sig: str, signode: desc_signature) -> Tuple[str, str]:
        """Breaks down construct signatures

        Parses out prefix and argument list from construct definition. The
        namespace and class will be determined by the nesting of domain
        directives.
        """
        sig = sig.strip()
        if '(' in sig and sig[-1:] == ')':
            member, arglist = sig.split('(', 1)
            member = member.strip()
            arglist = arglist[:-1].strip()
        else:
            member = sig
            arglist = None
        # If construct is nested, prefix the current prefix
        prefix = self.env.ref_context.get('js:object', None)
        mod_name = self.env.ref_context.get('js:module')

        name = member
        try:
            member_prefix, member_name = member.rsplit('.', 1)
        except ValueError:
            member_name = name
            member_prefix = ''
        finally:
            name = member_name
            if prefix and member_prefix:
                prefix = '.'.join([prefix, member_prefix])
            elif prefix is None and member_prefix:
                prefix = member_prefix
        fullname = name
        if prefix:
            fullname = '.'.join([prefix, name])

        signode['module'] = mod_name
        signode['object'] = prefix
        signode['fullname'] = fullname

        display_prefix = self.get_display_prefix()
        if display_prefix:
            signode += addnodes.desc_annotation('', '', *display_prefix)

        actual_prefix = None
        if prefix:
            actual_prefix = prefix
        elif mod_name:
            actual_prefix = mod_name
        if actual_prefix:
            addName = addnodes.desc_addname('', '')
            for p in actual_prefix.split('.'):
                addName += addnodes.desc_sig_name(p, p)
                addName += addnodes.desc_sig_punctuation('.', '.')
            signode += addName
        signode += addnodes.desc_name('', '', addnodes.desc_sig_name(name, name))
        if self.has_arguments:
            if not arglist:
                signode += addnodes.desc_parameterlist()
            else:
                _pseudo_parse_arglist(signode, arglist)
        return fullname, prefix

    def add_target_and_index(self, name_obj: Tuple[str, str], sig: str,
                             signode: desc_signature) -> None:
        mod_name = self.env.ref_context.get('js:module')
        fullname = (mod_name + '.' if mod_name else '') + name_obj[0]
        node_id = make_id(self.env, self.state.document, '', fullname)
        signode['ids'].append(node_id)

        # Assign old styled node_id not to break old hyperlinks (if possible)
        # Note: Will be removed in Sphinx-5.0 (RemovedInSphinx50Warning)
        old_node_id = self.make_old_id(fullname)
        if old_node_id not in self.state.document.ids and old_node_id not in signode['ids']:
            signode['ids'].append(old_node_id)

        self.state.document.note_explicit_target(signode)

        domain = cast(JavaScriptDomain, self.env.get_domain('js'))
        domain.note_object(fullname, self.objtype, node_id, location=signode)

        if 'noindexentry' not in self.options:
            indextext = self.get_index_text(mod_name, name_obj)
            if indextext:
                self.indexnode['entries'].append(('single', indextext, node_id, '', None))

    def get_index_text(self, objectname: str, name_obj: Tuple[str, str]) -> str:
        name, obj = name_obj
        if self.objtype == 'function':
            if not obj:
                return _('%s() (built-in function)') % name
            return _('%s() (%s method)') % (name, obj)
        elif self.objtype == 'class':
            return _('%s() (class)') % name
        elif self.objtype == 'data':
            return _('%s (global variable or constant)') % name
        elif self.objtype == 'attribute':
            return _('%s (%s attribute)') % (name, obj)
        return ''

    def before_content(self) -> None:
        """Handle object nesting before content

        :py:class:`JSObject` represents JavaScript language constructs. For
        constructs that are nestable, this method will build up a stack of the
        nesting hierarchy so that it can be later de-nested correctly, in
        :py:meth:`after_content`.

        For constructs that aren't nestable, the stack is bypassed, and instead
        only the most recent object is tracked. This object prefix name will be
        removed with :py:meth:`after_content`.

        The following keys are used in ``self.env.ref_context``:

            js:objects
                Stores the object prefix history. With each nested element, we
                add the object prefix to this list. When we exit that object's
                nesting level, :py:meth:`after_content` is triggered and the
                prefix is removed from the end of the list.

            js:object
                Current object prefix. This should generally reflect the last
                element in the prefix history
        """
        prefix = None
        if self.names:
            (obj_name, obj_name_prefix) = self.names.pop()
            prefix = obj_name_prefix.strip('.') if obj_name_prefix else None
            if self.allow_nesting:
                prefix = obj_name
        if prefix:
            self.env.ref_context['js:object'] = prefix
            if self.allow_nesting:
                objects = self.env.ref_context.setdefault('js:objects', [])
                objects.append(prefix)

    def after_content(self) -> None:
        """Handle object de-nesting after content

        If this class is a nestable object, removing the last nested class prefix
        ends further nesting in the object.

        If this class is not a nestable object, the list of classes should not
        be altered as we didn't affect the nesting levels in
        :py:meth:`before_content`.
        """
        objects = self.env.ref_context.setdefault('js:objects', [])
        if self.allow_nesting:
            try:
                objects.pop()
            except IndexError:
                pass
        self.env.ref_context['js:object'] = (objects[-1] if len(objects) > 0
                                             else None)

    def make_old_id(self, fullname: str) -> str:
        """Generate old styled node_id for JS objects.

        .. note:: Old Styled node_id was used until Sphinx-3.0.
                  This will be removed in Sphinx-5.0.
        """
        return fullname.replace('$', '_S_')


class JSCallable(JSObject):
    """Description of a JavaScript function, method or constructor."""
    has_arguments = True

    doc_field_types = [
        TypedField('arguments', label=_('Arguments'),
                   names=('argument', 'arg', 'parameter', 'param'),
                   typerolename='func', typenames=('paramtype', 'type')),
        GroupedField('errors', label=_('Throws'), rolename='func',
                     names=('throws', ),
                     can_collapse=True),
        Field('returnvalue', label=_('Returns'), has_arg=False,
              names=('returns', 'return')),
        Field('returntype', label=_('Return type'), has_arg=False,
              names=('rtype',)),
    ]


class JSConstructor(JSCallable):
    """Like a callable but with a different prefix."""

    allow_nesting = True

    def get_display_prefix(self) -> List[Node]:
        return [addnodes.desc_sig_keyword('class', 'class'),
                addnodes.desc_sig_space()]


class JSModule(SphinxDirective):
    """
    Directive to mark description of a new JavaScript module.

    This directive specifies the module name that will be used by objects that
    follow this directive.

    Options
    -------

    noindex
        If the ``noindex`` option is specified, no linkable elements will be
        created, and the module won't be added to the global module index. This
        is useful for splitting up the module definition across multiple
        sections or files.

    :param mod_name: Module name
    """

    has_content = False
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec: OptionSpec = {
        'noindex': directives.flag
    }

    def run(self) -> List[Node]:
        mod_name = self.arguments[0].strip()
        self.env.ref_context['js:module'] = mod_name
        noindex = 'noindex' in self.options
        ret: List[Node] = []
        if not noindex:
            domain = cast(JavaScriptDomain, self.env.get_domain('js'))

            node_id = make_id(self.env, self.state.document, 'module', mod_name)
            domain.note_module(mod_name, node_id)
            # Make a duplicate entry in 'objects' to facilitate searching for
            # the module in JavaScriptDomain.find_obj()
            domain.note_object(mod_name, 'module', node_id,
                               location=(self.env.docname, self.lineno))

            target = nodes.target('', '', ids=[node_id], ismod=True)

            # Assign old styled node_id not to break old hyperlinks (if possible)
            # Note: Will be removed in Sphinx-5.0 (RemovedInSphinx50Warning)
            old_node_id = self.make_old_id(mod_name)
            if old_node_id not in self.state.document.ids and old_node_id not in target['ids']:
                target['ids'].append(old_node_id)

            self.state.document.note_explicit_target(target)
            ret.append(target)
            indextext = _('%s (module)') % mod_name
            inode = addnodes.index(entries=[('single', indextext, node_id, '', None)])
            ret.append(inode)
        return ret

    def make_old_id(self, modname: str) -> str:
        """Generate old styled node_id for JS modules.

        .. note:: Old Styled node_id was used until Sphinx-3.0.
                  This will be removed in Sphinx-5.0.
        """
        return 'module-' + modname


class JSXRefRole(XRefRole):
    def process_link(self, env: BuildEnvironment, refnode: Element,
                     has_explicit_title: bool, title: str, target: str) -> Tuple[str, str]:
        # basically what sphinx.domains.python.PyXRefRole does
        refnode['js:object'] = env.ref_context.get('js:object')
        refnode['js:module'] = env.ref_context.get('js:module')
        if not has_explicit_title:
            title = title.lstrip('.')
            target = target.lstrip('~')
            if title[0:1] == '~':
                title = title[1:]
                dot = title.rfind('.')
                if dot != -1:
                    title = title[dot + 1:]
        if target[0:1] == '.':
            target = target[1:]
            refnode['refspecific'] = True
        return title, target


class JavaScriptDomain(Domain):
    """JavaScript language domain."""
    name = 'js'
    label = 'JavaScript'
    # if you add a new object type make sure to edit JSObject.get_index_string
    object_types = {
        'function':  ObjType(_('function'),  'func'),
        'method':    ObjType(_('method'),    'meth'),
        'class':     ObjType(_('class'),     'class'),
        'data':      ObjType(_('data'),      'data'),
        'attribute': ObjType(_('attribute'), 'attr'),
        'module':    ObjType(_('module'),    'mod'),
    }
    directives = {
        'function':  JSCallable,
        'method':    JSCallable,
        'class':     JSConstructor,
        'data':      JSObject,
        'attribute': JSObject,
        'module':    JSModule,
    }
    roles = {
        'func':  JSXRefRole(fix_parens=True),
        'meth':  JSXRefRole(fix_parens=True),
        'class': JSXRefRole(fix_parens=True),
        'data':  JSXRefRole(),
        'attr':  JSXRefRole(),
        'mod':   JSXRefRole(),
    }
    initial_data: Dict[str, Dict[str, Tuple[str, str]]] = {
        'objects': {},  # fullname -> docname, node_id, objtype
        'modules': {},  # modname  -> docname, node_id
    }

    @property
    def objects(self) -> Dict[str, Tuple[str, str, str]]:
        return self.data.setdefault('objects', {})  # fullname -> docname, node_id, objtype

    def note_object(self, fullname: str, objtype: str, node_id: str,
                    location: Any = None) -> None:
        if fullname in self.objects:
            docname = self.objects[fullname][0]
            logger.warning(__('duplicate %s description of %s, other %s in %s'),
                           objtype, fullname, objtype, docname, location=location)
        self.objects[fullname] = (self.env.docname, node_id, objtype)

    @property
    def modules(self) -> Dict[str, Tuple[str, str]]:
        return self.data.setdefault('modules', {})  # modname -> docname, node_id

    def note_module(self, modname: str, node_id: str) -> None:
        self.modules[modname] = (self.env.docname, node_id)

    def clear_doc(self, docname: str) -> None:
        for fullname, (pkg_docname, _node_id, _l) in list(self.objects.items()):
            if pkg_docname == docname:
                del self.objects[fullname]
        for modname, (pkg_docname, _node_id) in list(self.modules.items()):
            if pkg_docname == docname:
                del self.modules[modname]

    def merge_domaindata(self, docnames: List[str], otherdata: Dict) -> None:
        # XXX check duplicates
        for fullname, (fn, node_id, objtype) in otherdata['objects'].items():
            if fn in docnames:
                self.objects[fullname] = (fn, node_id, objtype)
        for mod_name, (pkg_docname, node_id) in otherdata['modules'].items():
            if pkg_docname in docnames:
                self.modules[mod_name] = (pkg_docname, node_id)

    def find_obj(self, env: BuildEnvironment, mod_name: str, prefix: str, name: str,
                 typ: str, searchorder: int = 0) -> Tuple[str, Tuple[str, str, str]]:
        if name[-2:] == '()':
            name = name[:-2]

        searches = []
        if mod_name and prefix:
            searches.append('.'.join([mod_name, prefix, name]))
        if mod_name:
            searches.append('.'.join([mod_name, name]))
        if prefix:
            searches.append('.'.join([prefix, name]))
        searches.append(name)

        if searchorder == 0:
            searches.reverse()

        newname = None
        for search_name in searches:
            if search_name in self.objects:
                newname = search_name

        return newname, self.objects.get(newname)

    def resolve_xref(self, env: BuildEnvironment, fromdocname: str, builder: Builder,
                     typ: str, target: str, node: pending_xref, contnode: Element
                     ) -> Optional[Element]:
        mod_name = node.get('js:module')
        prefix = node.get('js:object')
        searchorder = 1 if node.hasattr('refspecific') else 0
        name, obj = self.find_obj(env, mod_name, prefix, target, typ, searchorder)
        if not obj:
            return None
        return make_refnode(builder, fromdocname, obj[0], obj[1], contnode, name)

    def resolve_any_xref(self, env: BuildEnvironment, fromdocname: str, builder: Builder,
                         target: str, node: pending_xref, contnode: Element
                         ) -> List[Tuple[str, Element]]:
        mod_name = node.get('js:module')
        prefix = node.get('js:object')
        name, obj = self.find_obj(env, mod_name, prefix, target, None, 1)
        if not obj:
            return []
        return [('js:' + self.role_for_objtype(obj[2]),
                 make_refnode(builder, fromdocname, obj[0], obj[1], contnode, name))]

    def get_objects(self) -> Iterator[Tuple[str, str, str, str, str, int]]:
        for refname, (docname, node_id, typ) in list(self.objects.items()):
            yield refname, refname, typ, docname, node_id, 1

    def get_full_qualified_name(self, node: Element) -> str:
        modname = node.get('js:module')
        prefix = node.get('js:object')
        target = node.get('reftarget')
        if target is None:
            return None
        else:
            return '.'.join(filter(None, [modname, prefix, target]))


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_domain(JavaScriptDomain)

    return {
        'version': 'builtin',
        'env_version': 2,
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
