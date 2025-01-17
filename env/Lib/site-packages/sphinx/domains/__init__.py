"""Support for domains.

Domains are groupings of description directives
and roles describing e.g. constructs of one programming language.
"""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING

from sphinx.domains._index import Index, IndexEntry
from sphinx.locale import _

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence, Set
    from typing import Any

    from docutils import nodes
    from docutils.nodes import Element, Node
    from docutils.parsers.rst import Directive
    from docutils.parsers.rst.states import Inliner

    from sphinx.addnodes import pending_xref
    from sphinx.builders import Builder
    from sphinx.environment import BuildEnvironment
    from sphinx.roles import XRefRole
    from sphinx.util.typing import RoleFunction, TitleGetter

__all__ = (
    'Domain',
    'Index',
    'IndexEntry',
    'ObjType',
)


class ObjType:
    """
    An ObjType is the description for a type of object that a domain can
    document.  In the object_types attribute of Domain subclasses, object type
    names are mapped to instances of this class.

    Constructor arguments:

    - *lname*: localized name of the type (do not include domain name)
    - *roles*: all the roles that can refer to an object of this type
    - *attrs*: object attributes -- currently only "searchprio" is known,
      which defines the object's priority in the full-text search index,
      see :meth:`Domain.get_objects()`.
    """

    known_attrs = {
        'searchprio': 1,
    }

    def __init__(self, lname: str, /, *roles: Any, **attrs: Any) -> None:
        self.lname: str = lname
        self.roles: tuple[Any, ...] = roles
        self.attrs: dict[str, Any] = self.known_attrs | attrs


class Domain:
    """
    A Domain is meant to be a group of "object" description directives for
    objects of a similar nature, and corresponding roles to create references to
    them.  Examples would be Python modules, classes, functions etc., elements
    of a templating language, Sphinx roles and directives, etc.

    Each domain has a separate storage for information about existing objects
    and how to reference them in `self.data`, which must be a dictionary.  It
    also must implement several functions that expose the object information in
    a uniform way to parts of Sphinx that allow the user to reference or search
    for objects in a domain-agnostic way.

    About `self.data`: since all object and cross-referencing information is
    stored on a BuildEnvironment instance, the `domain.data` object is also
    stored in the `env.domaindata` dict under the key `domain.name`.  Before the
    build process starts, every active domain is instantiated and given the
    environment object; the `domaindata` dict must then either be nonexistent or
    a dictionary whose 'version' key is equal to the domain class'
    :attr:`data_version` attribute.  Otherwise, `OSError` is raised and the
    pickled environment is discarded.
    """

    #: domain name: should be short, but unique
    name = ''
    #: domain label: longer, more descriptive (used in messages)
    label = ''
    #: type (usually directive) name -> ObjType instance
    object_types: dict[str, ObjType] = {}
    #: directive name -> directive class
    directives: dict[str, type[Directive]] = {}
    #: role name -> role callable
    roles: dict[str, RoleFunction | XRefRole] = {}
    #: a list of Index subclasses
    indices: list[type[Index]] = []
    #: role name -> a warning message if reference is missing
    dangling_warnings: dict[str, str] = {}
    #: node_class -> (enum_node_type, title_getter)
    enumerable_nodes: dict[type[Node], tuple[str, TitleGetter | None]] = {}
    #: data value for a fresh environment
    initial_data: dict = {}
    #: data value
    data: dict[str, Any]
    #: data version, bump this when the format of `self.data` changes
    data_version = 0

    def __init__(self, env: BuildEnvironment) -> None:
        self.env: BuildEnvironment = env
        self._role_cache: dict[str, Callable] = {}
        self._directive_cache: dict[str, Callable] = {}
        self._role2type: dict[str, list[str]] = {}
        self._type2role: dict[str, str] = {}

        # convert class variables to instance one (to enhance through API)
        self.object_types = dict(self.object_types)
        self.directives = dict(self.directives)
        self.roles = dict(self.roles)
        self.indices = list(self.indices)

        if self.name not in env.domaindata:
            assert isinstance(self.initial_data, dict)
            new_data = copy.deepcopy(self.initial_data)
            new_data['version'] = self.data_version
            self.data = env.domaindata[self.name] = new_data
        else:
            self.data = env.domaindata[self.name]
            if self.data['version'] != self.data_version:
                raise OSError('data of %r domain out of date' % self.label)
        for name, obj in self.object_types.items():
            for rolename in obj.roles:
                self._role2type.setdefault(rolename, []).append(name)
            self._type2role[name] = obj.roles[0] if obj.roles else ''
        self.objtypes_for_role = self._role2type.get
        self.role_for_objtype = self._type2role.get

    def setup(self) -> None:
        """Set up domain object."""
        # Add special hyperlink target for index pages (ex. py-modindex)
        std = self.env.domains.standard_domain
        for index in self.indices:
            if index.name and index.localname:
                docname = f"{self.name}-{index.name}"
                std.note_hyperlink_target(docname, docname, '', index.localname)

    def add_object_type(self, name: str, objtype: ObjType) -> None:
        """Add an object type."""
        self.object_types[name] = objtype
        if objtype.roles:
            self._type2role[name] = objtype.roles[0]
        else:
            self._type2role[name] = ''

        for role in objtype.roles:
            self._role2type.setdefault(role, []).append(name)

    def role(self, name: str) -> RoleFunction | None:
        """Return a role adapter function that always gives the registered
        role its full name ('domain:name') as the first argument.
        """
        if name in self._role_cache:
            return self._role_cache[name]
        if name not in self.roles:
            return None
        fullname = f'{self.name}:{name}'

        def role_adapter(typ: str, rawtext: str, text: str, lineno: int,
                         inliner: Inliner, options: dict | None = None,
                         content: Sequence[str] = (),
                         ) -> tuple[list[Node], list[nodes.system_message]]:
            return self.roles[name](fullname, rawtext, text, lineno,
                                    inliner, options or {}, content)
        self._role_cache[name] = role_adapter
        return role_adapter

    def directive(self, name: str) -> Callable | None:
        """Return a directive adapter class that always gives the registered
        directive its full name ('domain:name') as ``self.name``.
        """
        if name in self._directive_cache:
            return self._directive_cache[name]
        if name not in self.directives:
            return None
        fullname = f'{self.name}:{name}'
        BaseDirective = self.directives[name]

        class DirectiveAdapter(BaseDirective):  # type: ignore[valid-type,misc]
            def run(self) -> list[Node]:
                self.name = fullname
                return super().run()
        self._directive_cache[name] = DirectiveAdapter
        return DirectiveAdapter

    # methods that should be overwritten

    def clear_doc(self, docname: str) -> None:
        """Remove traces of a document in the domain-specific inventories."""
        pass

    def merge_domaindata(self, docnames: Set[str], otherdata: dict[str, Any]) -> None:
        """Merge in data regarding *docnames* from a different domaindata
        inventory (coming from a subprocess in parallel builds).
        """
        raise NotImplementedError('merge_domaindata must be implemented in %s '
                                  'to be able to do parallel builds!' %
                                  self.__class__)

    def process_doc(self, env: BuildEnvironment, docname: str,
                    document: nodes.document) -> None:
        """Process a document after it is read by the environment."""
        pass

    def check_consistency(self) -> None:
        """Do consistency checks (**experimental**)."""
        pass

    def process_field_xref(self, pnode: pending_xref) -> None:
        """Process a pending xref created in a doc field.
        For example, attach information about the current scope.
        """
        pass

    def resolve_xref(self, env: BuildEnvironment, fromdocname: str, builder: Builder,
                     typ: str, target: str, node: pending_xref, contnode: Element,
                     ) -> Element | None:
        """Resolve the pending_xref *node* with the given *typ* and *target*.

        This method should return a new node, to replace the xref node,
        containing the *contnode* which is the markup content of the
        cross-reference.

        If no resolution can be found, None can be returned; the xref node will
        then given to the :event:`missing-reference` event, and if that yields no
        resolution, replaced by *contnode*.

        The method can also raise :exc:`sphinx.environment.NoUri` to suppress
        the :event:`missing-reference` event being emitted.
        """
        pass

    def resolve_any_xref(self, env: BuildEnvironment, fromdocname: str, builder: Builder,
                         target: str, node: pending_xref, contnode: Element,
                         ) -> list[tuple[str, Element]]:
        """Resolve the pending_xref *node* with the given *target*.

        The reference comes from an "any" or similar role, which means that we
        don't know the type.  Otherwise, the arguments are the same as for
        :meth:`resolve_xref`.

        The method must return a list (potentially empty) of tuples
        ``('domain:role', newnode)``, where ``'domain:role'`` is the name of a
        role that could have created the same reference, e.g. ``'py:func'``.
        ``newnode`` is what :meth:`resolve_xref` would return.

        .. versionadded:: 1.3
        """
        raise NotImplementedError

    def get_objects(self) -> Iterable[tuple[str, str, str, str, str, int]]:
        """Return an iterable of "object descriptions".

        Object descriptions are tuples with six items:

        ``name``
          Fully qualified name.

        ``dispname``
          Name to display when searching/linking.

        ``type``
          Object type, a key in ``self.object_types``.

        ``docname``
          The document where it is to be found.

        ``anchor``
          The anchor name for the object.

        ``priority``
          How "important" the object is (determines placement in search
          results). One of:

          ``1``
            Default priority (placed before full-text matches).
          ``0``
            Object is important (placed before default-priority objects).
          ``2``
            Object is unimportant (placed after full-text matches).
          ``-1``
            Object should not show up in search at all.
        """
        return []

    def get_type_name(self, type: ObjType, primary: bool = False) -> str:
        """Return full name for given ObjType."""
        if primary:
            return type.lname
        return _('%s %s') % (self.label, type.lname)

    def get_enumerable_node_type(self, node: Node) -> str | None:
        """Get type of enumerable nodes (experimental)."""
        enum_node_type, _ = self.enumerable_nodes.get(node.__class__, (None, None))
        return enum_node_type

    def get_full_qualified_name(self, node: Element) -> str | None:
        """Return full qualified name for given node."""
        pass
