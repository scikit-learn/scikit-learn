"""
    sphinx.util.docutils
    ~~~~~~~~~~~~~~~~~~~~

    Utility functions for docutils.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import os
import re
from contextlib import contextmanager
from copy import copy
from os import path
from types import ModuleType
from typing import (IO, TYPE_CHECKING, Any, Callable, Dict, Generator, List, Optional, Set,
                    Tuple, Type, cast)

import docutils
from docutils import nodes
from docutils.io import FileOutput
from docutils.nodes import Element, Node, system_message
from docutils.parsers.rst import Directive, directives, roles
from docutils.parsers.rst.states import Inliner
from docutils.statemachine import State, StateMachine, StringList
from docutils.utils import Reporter, unescape
from packaging import version

from sphinx.errors import SphinxError
from sphinx.locale import _, __
from sphinx.util import logging
from sphinx.util.typing import RoleFunction

logger = logging.getLogger(__name__)
report_re = re.compile('^(.+?:(?:\\d+)?): \\((DEBUG|INFO|WARNING|ERROR|SEVERE)/(\\d+)?\\) ')

if TYPE_CHECKING:
    from sphinx.builders import Builder
    from sphinx.config import Config
    from sphinx.environment import BuildEnvironment


__version_info__ = version.parse(docutils.__version__).release
additional_nodes: Set[Type[Element]] = set()


@contextmanager
def docutils_namespace() -> Generator[None, None, None]:
    """Create namespace for reST parsers."""
    try:
        _directives = copy(directives._directives)  # type: ignore
        _roles = copy(roles._roles)  # type: ignore

        yield
    finally:
        directives._directives = _directives  # type: ignore
        roles._roles = _roles  # type: ignore

        for node in list(additional_nodes):
            unregister_node(node)
            additional_nodes.discard(node)


def is_directive_registered(name: str) -> bool:
    """Check the *name* directive is already registered."""
    return name in directives._directives  # type: ignore


def register_directive(name: str, directive: Type[Directive]) -> None:
    """Register a directive to docutils.

    This modifies global state of docutils.  So it is better to use this
    inside ``docutils_namespace()`` to prevent side-effects.
    """
    directives.register_directive(name, directive)


def is_role_registered(name: str) -> bool:
    """Check the *name* role is already registered."""
    return name in roles._roles  # type: ignore


def register_role(name: str, role: RoleFunction) -> None:
    """Register a role to docutils.

    This modifies global state of docutils.  So it is better to use this
    inside ``docutils_namespace()`` to prevent side-effects.
    """
    roles.register_local_role(name, role)


def unregister_role(name: str) -> None:
    """Unregister a role from docutils."""
    roles._roles.pop(name, None)  # type: ignore


def is_node_registered(node: Type[Element]) -> bool:
    """Check the *node* is already registered."""
    return hasattr(nodes.GenericNodeVisitor, 'visit_' + node.__name__)


def register_node(node: Type[Element]) -> None:
    """Register a node to docutils.

    This modifies global state of some visitors.  So it is better to use this
    inside ``docutils_namespace()`` to prevent side-effects.
    """
    if not hasattr(nodes.GenericNodeVisitor, 'visit_' + node.__name__):
        nodes._add_node_class_names([node.__name__])  # type: ignore
        additional_nodes.add(node)


def unregister_node(node: Type[Element]) -> None:
    """Unregister a node from docutils.

    This is inverse of ``nodes._add_nodes_class_names()``.
    """
    if hasattr(nodes.GenericNodeVisitor, 'visit_' + node.__name__):
        delattr(nodes.GenericNodeVisitor, "visit_" + node.__name__)
        delattr(nodes.GenericNodeVisitor, "depart_" + node.__name__)
        delattr(nodes.SparseNodeVisitor, 'visit_' + node.__name__)
        delattr(nodes.SparseNodeVisitor, 'depart_' + node.__name__)


@contextmanager
def patched_get_language() -> Generator[None, None, None]:
    """Patch docutils.languages.get_language() temporarily.

    This ignores the second argument ``reporter`` to suppress warnings.
    refs: https://github.com/sphinx-doc/sphinx/issues/3788
    """
    from docutils.languages import get_language

    def patched_get_language(language_code: str, reporter: Reporter = None) -> Any:
        return get_language(language_code)

    try:
        docutils.languages.get_language = patched_get_language
        yield
    finally:
        # restore original implementations
        docutils.languages.get_language = get_language


@contextmanager
def using_user_docutils_conf(confdir: Optional[str]) -> Generator[None, None, None]:
    """Let docutils know the location of ``docutils.conf`` for Sphinx."""
    try:
        docutilsconfig = os.environ.get('DOCUTILSCONFIG', None)
        if confdir:
            os.environ['DOCUTILSCONFIG'] = path.join(path.abspath(confdir), 'docutils.conf')

        yield
    finally:
        if docutilsconfig is None:
            os.environ.pop('DOCUTILSCONFIG', None)
        else:
            os.environ['DOCUTILSCONFIG'] = docutilsconfig


@contextmanager
def patch_docutils(confdir: Optional[str] = None) -> Generator[None, None, None]:
    """Patch to docutils temporarily."""
    with patched_get_language(), using_user_docutils_conf(confdir):
        yield


class CustomReSTDispatcher:
    """Custom reST's mark-up dispatcher.

    This replaces docutils's directives and roles dispatch mechanism for reST parser
    by original one temporarily.
    """

    def __init__(self) -> None:
        self.directive_func: Callable = lambda *args: (None, [])
        self.roles_func: Callable = lambda *args: (None, [])

    def __enter__(self) -> None:
        self.enable()

    def __exit__(self, exc_type: Type[Exception], exc_value: Exception, traceback: Any) -> None:  # NOQA
        self.disable()

    def enable(self) -> None:
        self.directive_func = directives.directive
        self.role_func = roles.role

        directives.directive = self.directive
        roles.role = self.role

    def disable(self) -> None:
        directives.directive = self.directive_func
        roles.role = self.role_func

    def directive(self,
                  directive_name: str, language_module: ModuleType, document: nodes.document
                  ) -> Tuple[Optional[Type[Directive]], List[system_message]]:
        return self.directive_func(directive_name, language_module, document)

    def role(self, role_name: str, language_module: ModuleType, lineno: int, reporter: Reporter
             ) -> Tuple[RoleFunction, List[system_message]]:
        return self.role_func(role_name, language_module, lineno, reporter)


class ElementLookupError(Exception):
    pass


class sphinx_domains(CustomReSTDispatcher):
    """Monkey-patch directive and role dispatch, so that domain-specific
    markup takes precedence.
    """
    def __init__(self, env: "BuildEnvironment") -> None:
        self.env = env
        super().__init__()

    def lookup_domain_element(self, type: str, name: str) -> Any:
        """Lookup a markup element (directive or role), given its name which can
        be a full name (with domain).
        """
        name = name.lower()
        # explicit domain given?
        if ':' in name:
            domain_name, name = name.split(':', 1)
            if domain_name in self.env.domains:
                domain = self.env.get_domain(domain_name)
                element = getattr(domain, type)(name)
                if element is not None:
                    return element, []
            else:
                logger.warning(_('unknown directive or role name: %s:%s'), domain_name, name)
        # else look in the default domain
        else:
            def_domain = self.env.temp_data.get('default_domain')
            if def_domain is not None:
                element = getattr(def_domain, type)(name)
                if element is not None:
                    return element, []

        # always look in the std domain
        element = getattr(self.env.get_domain('std'), type)(name)
        if element is not None:
            return element, []

        raise ElementLookupError

    def directive(self,
                  directive_name: str, language_module: ModuleType, document: nodes.document
                  ) -> Tuple[Optional[Type[Directive]], List[system_message]]:
        try:
            return self.lookup_domain_element('directive', directive_name)
        except ElementLookupError:
            return super().directive(directive_name, language_module, document)

    def role(self, role_name: str, language_module: ModuleType, lineno: int, reporter: Reporter
             ) -> Tuple[RoleFunction, List[system_message]]:
        try:
            return self.lookup_domain_element('role', role_name)
        except ElementLookupError:
            return super().role(role_name, language_module, lineno, reporter)


class WarningStream:
    def write(self, text: str) -> None:
        matched = report_re.search(text)
        if not matched:
            logger.warning(text.rstrip("\r\n"))
        else:
            location, type, level = matched.groups()
            message = report_re.sub('', text).rstrip()
            logger.log(type, message, location=location)


class LoggingReporter(Reporter):
    @classmethod
    def from_reporter(cls, reporter: Reporter) -> "LoggingReporter":
        """Create an instance of LoggingReporter from other reporter object."""
        return cls(reporter.source, reporter.report_level, reporter.halt_level,
                   reporter.debug_flag, reporter.error_handler)

    def __init__(self, source: str, report_level: int = Reporter.WARNING_LEVEL,
                 halt_level: int = Reporter.SEVERE_LEVEL, debug: bool = False,
                 error_handler: str = 'backslashreplace') -> None:
        stream = cast(IO, WarningStream())
        super().__init__(source, report_level, halt_level,
                         stream, debug, error_handler=error_handler)


class NullReporter(Reporter):
    """A dummy reporter; write nothing."""

    def __init__(self) -> None:
        super().__init__('', 999, 4)


def is_html5_writer_available() -> bool:
    return __version_info__ > (0, 13, 0)


@contextmanager
def switch_source_input(state: State, content: StringList) -> Generator[None, None, None]:
    """Switch current source input of state temporarily."""
    try:
        # remember the original ``get_source_and_line()`` method
        get_source_and_line = state.memo.reporter.get_source_and_line  # type: ignore

        # replace it by new one
        state_machine = StateMachine([], None)
        state_machine.input_lines = content
        state.memo.reporter.get_source_and_line = state_machine.get_source_and_line  # type: ignore  # NOQA

        yield
    finally:
        # restore the method
        state.memo.reporter.get_source_and_line = get_source_and_line  # type: ignore


class SphinxFileOutput(FileOutput):
    """Better FileOutput class for Sphinx."""

    def __init__(self, **kwargs: Any) -> None:
        self.overwrite_if_changed = kwargs.pop('overwrite_if_changed', False)
        super().__init__(**kwargs)

    def write(self, data: str) -> str:
        if (self.destination_path and self.autoclose and 'b' not in self.mode and
                self.overwrite_if_changed and os.path.exists(self.destination_path)):
            with open(self.destination_path, encoding=self.encoding) as f:
                # skip writing: content not changed
                if f.read() == data:
                    return data

        return super().write(data)


class SphinxDirective(Directive):
    """A base class for Sphinx directives.

    This class provides helper methods for Sphinx directives.

    .. note:: The subclasses of this class might not work with docutils.
              This class is strongly coupled with Sphinx.
    """

    @property
    def env(self) -> "BuildEnvironment":
        """Reference to the :class:`.BuildEnvironment` object."""
        return self.state.document.settings.env

    @property
    def config(self) -> "Config":
        """Reference to the :class:`.Config` object."""
        return self.env.config

    def get_source_info(self) -> Tuple[str, int]:
        """Get source and line number."""
        return self.state_machine.get_source_and_line(self.lineno)

    def set_source_info(self, node: Node) -> None:
        """Set source and line number to the node."""
        node.source, node.line = self.get_source_info()

    def get_location(self) -> str:
        """Get current location info for logging."""
        return ':'.join(str(s) for s in self.get_source_info())


class SphinxRole:
    """A base class for Sphinx roles.

    This class provides helper methods for Sphinx roles.

    .. note:: The subclasses of this class might not work with docutils.
              This class is strongly coupled with Sphinx.
    """
    name: str           #: The role name actually used in the document.
    rawtext: str        #: A string containing the entire interpreted text input.
    text: str           #: The interpreted text content.
    lineno: int         #: The line number where the interpreted text begins.
    inliner: Inliner    #: The ``docutils.parsers.rst.states.Inliner`` object.
    options: Dict       #: A dictionary of directive options for customization
                        #: (from the "role" directive).
    content: List[str]  #: A list of strings, the directive content for customization
                        #: (from the "role" directive).

    def __call__(self, name: str, rawtext: str, text: str, lineno: int,
                 inliner: Inliner, options: Dict = {}, content: List[str] = []
                 ) -> Tuple[List[Node], List[system_message]]:
        self.rawtext = rawtext
        self.text = unescape(text)
        self.lineno = lineno
        self.inliner = inliner
        self.options = options
        self.content = content

        # guess role type
        if name:
            self.name = name.lower()
        else:
            self.name = self.env.temp_data.get('default_role', '')
            if not self.name:
                self.name = self.env.config.default_role
            if not self.name:
                raise SphinxError('cannot determine default role!')

        return self.run()

    def run(self) -> Tuple[List[Node], List[system_message]]:
        raise NotImplementedError

    @property
    def env(self) -> "BuildEnvironment":
        """Reference to the :class:`.BuildEnvironment` object."""
        return self.inliner.document.settings.env

    @property
    def config(self) -> "Config":
        """Reference to the :class:`.Config` object."""
        return self.env.config

    def get_source_info(self, lineno: int = None) -> Tuple[str, int]:
        if lineno is None:
            lineno = self.lineno
        return self.inliner.reporter.get_source_and_line(lineno)  # type: ignore

    def set_source_info(self, node: Node, lineno: int = None) -> None:
        node.source, node.line = self.get_source_info(lineno)

    def get_location(self) -> str:
        """Get current location info for logging."""
        return ':'.join(str(s) for s in self.get_source_info())


class ReferenceRole(SphinxRole):
    """A base class for reference roles.

    The reference roles can accept ``link title <target>`` style as a text for
    the role.  The parsed result; link title and target will be stored to
    ``self.title`` and ``self.target``.
    """
    has_explicit_title: bool    #: A boolean indicates the role has explicit title or not.
    disabled: bool              #: A boolean indicates the reference is disabled.
    title: str                  #: The link title for the interpreted text.
    target: str                 #: The link target for the interpreted text.

    # \x00 means the "<" was backslash-escaped
    explicit_title_re = re.compile(r'^(.+?)\s*(?<!\x00)<(.*?)>$', re.DOTALL)

    def __call__(self, name: str, rawtext: str, text: str, lineno: int,
                 inliner: Inliner, options: Dict = {}, content: List[str] = []
                 ) -> Tuple[List[Node], List[system_message]]:
        # if the first character is a bang, don't cross-reference at all
        self.disabled = text.startswith('!')

        matched = self.explicit_title_re.match(text)
        if matched:
            self.has_explicit_title = True
            self.title = unescape(matched.group(1))
            self.target = unescape(matched.group(2))
        else:
            self.has_explicit_title = False
            self.title = unescape(text)
            self.target = unescape(text)

        return super().__call__(name, rawtext, text, lineno, inliner, options, content)


class SphinxTranslator(nodes.NodeVisitor):
    """A base class for Sphinx translators.

    This class adds a support for visitor/departure method for super node class
    if visitor/departure method for node class is not found.

    It also provides helper methods for Sphinx translators.

    .. note:: The subclasses of this class might not work with docutils.
              This class is strongly coupled with Sphinx.
    """

    def __init__(self, document: nodes.document, builder: "Builder") -> None:
        super().__init__(document)
        self.builder = builder
        self.config = builder.config
        self.settings = document.settings

    def dispatch_visit(self, node: Node) -> None:
        """
        Dispatch node to appropriate visitor method.
        The priority of visitor method is:

        1. ``self.visit_{node_class}()``
        2. ``self.visit_{super_node_class}()``
        3. ``self.unknown_visit()``
        """
        for node_class in node.__class__.__mro__:
            method = getattr(self, 'visit_%s' % (node_class.__name__), None)
            if method:
                method(node)
                break
        else:
            super().dispatch_visit(node)

    def dispatch_departure(self, node: Node) -> None:
        """
        Dispatch node to appropriate departure method.
        The priority of departure method is:

        1. ``self.depart_{node_class}()``
        2. ``self.depart_{super_node_class}()``
        3. ``self.unknown_departure()``
        """
        for node_class in node.__class__.__mro__:
            method = getattr(self, 'depart_%s' % (node_class.__name__), None)
            if method:
                method(node)
                break
        else:
            super().dispatch_departure(node)

    def unknown_visit(self, node: Node) -> None:
        logger.warning(__('unknown node type: %r'), node, location=node)


# Node.findall() is a new interface to traverse a doctree since docutils-0.18.
# This applies a patch docutils-0.17 or older to be available Node.findall()
# method to use it from our codebase.
if __version_info__ < (0, 18):
    def findall(self, *args, **kwargs):
        return iter(self.traverse(*args, **kwargs))

    Node.findall = findall  # type: ignore


# cache a vanilla instance of nodes.document
# Used in new_document() function
__document_cache__: Optional[nodes.document] = None


def new_document(source_path: str, settings: Any = None) -> nodes.document:
    """Return a new empty document object.  This is an alternative of docutils'.

    This is a simple wrapper for ``docutils.utils.new_document()``.  It
    caches the result of docutils' and use it on second call for instantiation.
    This makes an instantiation of document nodes much faster.
    """
    global __document_cache__
    if __document_cache__ is None:
        __document_cache__ = docutils.utils.new_document(source_path)

    if settings is None:
        # Make a copy of ``settings`` from cache to accelerate instantiation
        settings = copy(__document_cache__.settings)

    # Create a new instance of nodes.document using cached reporter
    from sphinx import addnodes
    document = addnodes.document(settings, __document_cache__.reporter, source=source_path)
    document.note_source(source_path, -1)
    return document
