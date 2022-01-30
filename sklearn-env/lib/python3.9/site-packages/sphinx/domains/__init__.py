"""
    sphinx.domains
    ~~~~~~~~~~~~~~

    Support for domains, which are groupings of description directives
    and roles describing e.g. constructs of one programming language.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import copy
from abc import ABC, abstractmethod
from typing import (TYPE_CHECKING, Any, Callable, Dict, Iterable, List, NamedTuple, Optional,
                    Tuple, Type, Union, cast)

from docutils import nodes
from docutils.nodes import Element, Node, system_message
from docutils.parsers.rst.states import Inliner

from sphinx.addnodes import pending_xref
from sphinx.errors import SphinxError
from sphinx.locale import _
from sphinx.roles import XRefRole
from sphinx.util.typing import RoleFunction

if TYPE_CHECKING:
    from sphinx.builders import Builder
    from sphinx.environment import BuildEnvironment


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

    def __init__(self, lname: str, *roles: Any, **attrs: Any) -> None:
        self.lname = lname
        self.roles: Tuple = roles
        self.attrs: Dict = self.known_attrs.copy()
        self.attrs.update(attrs)


class IndexEntry(NamedTuple):
    name: str
    subtype: int
    docname: str
    anchor: str
    extra: str
    qualifier: str
    descr: str


class Index(ABC):
    """
    An Index is the description for a domain-specific index.  To add an index to
    a domain, subclass Index, overriding the three name attributes:

    * `name` is an identifier used for generating file names.
      It is also used for a hyperlink target for the index. Therefore, users can
      refer the index page using ``ref`` role and a string which is combined
      domain name and ``name`` attribute (ex. ``:ref:`py-modindex```).
    * `localname` is the section title for the index.
    * `shortname` is a short name for the index, for use in the relation bar in
      HTML output.  Can be empty to disable entries in the relation bar.

    and providing a :meth:`generate()` method.  Then, add the index class to
    your domain's `indices` list.  Extensions can add indices to existing
    domains using :meth:`~sphinx.application.Sphinx.add_index_to_domain()`.

    .. versionchanged:: 3.0

       Index pages can be referred by domain name and index name via
       :rst:role:`ref` role.
    """

    name: str = None
    localname: str = None
    shortname: str = None

    def __init__(self, domain: "Domain") -> None:
        if self.name is None or self.localname is None:
            raise SphinxError('Index subclass %s has no valid name or localname'
                              % self.__class__.__name__)
        self.domain = domain

    @abstractmethod
    def generate(self, docnames: Iterable[str] = None
                 ) -> Tuple[List[Tuple[str, List[IndexEntry]]], bool]:
        """Get entries for the index.

        If ``docnames`` is given, restrict to entries referring to these
        docnames.

        The return value is a tuple of ``(content, collapse)``:

        ``collapse``
          A boolean that determines if sub-entries should start collapsed (for
          output formats that support collapsing sub-entries).

        ``content``:
          A sequence of ``(letter, entries)`` tuples, where ``letter`` is the
          "heading" for the given ``entries``, usually the starting letter, and
          ``entries`` is a sequence of single entries. Each entry is a sequence
          ``[name, subtype, docname, anchor, extra, qualifier, descr]``. The
          items in this sequence have the following meaning:

          ``name``
            The name of the index entry to be displayed.

          ``subtype``
            The sub-entry related type. One of:

            ``0``
              A normal entry.
            ``1``
              An entry with sub-entries.
            ``2``
              A sub-entry.

          ``docname``
            *docname* where the entry is located.

          ``anchor``
            Anchor for the entry within ``docname``

          ``extra``
            Extra info for the entry.

          ``qualifier``
            Qualifier for the description.

          ``descr``
            Description for the entry.

        Qualifier and description are not rendered for some output formats such
        as LaTeX.
        """
        raise NotImplementedError


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
    object_types: Dict[str, ObjType] = {}
    #: directive name -> directive class
    directives: Dict[str, Any] = {}
    #: role name -> role callable
    roles: Dict[str, Union[RoleFunction, XRefRole]] = {}
    #: a list of Index subclasses
    indices: List[Type[Index]] = []
    #: role name -> a warning message if reference is missing
    dangling_warnings: Dict[str, str] = {}
    #: node_class -> (enum_node_type, title_getter)
    enumerable_nodes: Dict[Type[Node], Tuple[str, Callable]] = {}

    #: data value for a fresh environment
    initial_data: Dict = {}
    #: data value
    data: Dict
    #: data version, bump this when the format of `self.data` changes
    data_version = 0

    def __init__(self, env: "BuildEnvironment") -> None:
        self.env: BuildEnvironment = env
        self._role_cache: Dict[str, Callable] = {}
        self._directive_cache: Dict[str, Callable] = {}
        self._role2type: Dict[str, List[str]] = {}
        self._type2role: Dict[str, str] = {}

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
        self.objtypes_for_role: Callable[[str], List[str]] = self._role2type.get
        self.role_for_objtype: Callable[[str], str] = self._type2role.get

    def setup(self) -> None:
        """Set up domain object."""
        from sphinx.domains.std import StandardDomain

        # Add special hyperlink target for index pages (ex. py-modindex)
        std = cast(StandardDomain, self.env.get_domain('std'))
        for index in self.indices:
            if index.name and index.localname:
                docname = "%s-%s" % (self.name, index.name)
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

    def role(self, name: str) -> Optional[RoleFunction]:
        """Return a role adapter function that always gives the registered
        role its full name ('domain:name') as the first argument.
        """
        if name in self._role_cache:
            return self._role_cache[name]
        if name not in self.roles:
            return None
        fullname = '%s:%s' % (self.name, name)

        def role_adapter(typ: str, rawtext: str, text: str, lineno: int,
                         inliner: Inliner, options: Dict = {}, content: List[str] = []
                         ) -> Tuple[List[Node], List[system_message]]:
            return self.roles[name](fullname, rawtext, text, lineno,
                                    inliner, options, content)
        self._role_cache[name] = role_adapter
        return role_adapter

    def directive(self, name: str) -> Optional[Callable]:
        """Return a directive adapter class that always gives the registered
        directive its full name ('domain:name') as ``self.name``.
        """
        if name in self._directive_cache:
            return self._directive_cache[name]
        if name not in self.directives:
            return None
        fullname = '%s:%s' % (self.name, name)
        BaseDirective = self.directives[name]

        class DirectiveAdapter(BaseDirective):  # type: ignore
            def run(self) -> List[Node]:
                self.name = fullname
                return super().run()
        self._directive_cache[name] = DirectiveAdapter
        return DirectiveAdapter

    # methods that should be overwritten

    def clear_doc(self, docname: str) -> None:
        """Remove traces of a document in the domain-specific inventories."""
        pass

    def merge_domaindata(self, docnames: List[str], otherdata: Dict) -> None:
        """Merge in data regarding *docnames* from a different domaindata
        inventory (coming from a subprocess in parallel builds).
        """
        raise NotImplementedError('merge_domaindata must be implemented in %s '
                                  'to be able to do parallel builds!' %
                                  self.__class__)

    def process_doc(self, env: "BuildEnvironment", docname: str,
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

    def resolve_xref(self, env: "BuildEnvironment", fromdocname: str, builder: "Builder",
                     typ: str, target: str, node: pending_xref, contnode: Element
                     ) -> Optional[Element]:
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

    def resolve_any_xref(self, env: "BuildEnvironment", fromdocname: str, builder: "Builder",
                         target: str, node: pending_xref, contnode: Element
                         ) -> List[Tuple[str, Element]]:
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

    def get_objects(self) -> Iterable[Tuple[str, str, str, str, str, int]]:
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

    def get_enumerable_node_type(self, node: Node) -> Optional[str]:
        """Get type of enumerable nodes (experimental)."""
        enum_node_type, _ = self.enumerable_nodes.get(node.__class__, (None, None))
        return enum_node_type

    def get_full_qualified_name(self, node: Element) -> Optional[str]:
        """Return full qualified name for given node."""
        return None
