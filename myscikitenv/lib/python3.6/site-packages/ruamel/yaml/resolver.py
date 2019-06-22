# coding: utf-8

from __future__ import absolute_import

import re

if False:  # MYPY
    from typing import Any, Dict, List, Union, Text, Optional  # NOQA
    from ruamel.yaml.compat import VersionType  # NOQA

from ruamel.yaml.compat import string_types, _DEFAULT_YAML_VERSION  # NOQA
from ruamel.yaml.error import *  # NOQA
from ruamel.yaml.nodes import MappingNode, ScalarNode, SequenceNode  # NOQA
from ruamel.yaml.util import RegExp  # NOQA

__all__ = ['BaseResolver', 'Resolver', 'VersionedResolver']


# fmt: off
# resolvers consist of
# - a list of applicable version
# - a tag
# - a regexp
# - a list of first characters to match
implicit_resolvers = [
    ([(1, 2)],
        u'tag:yaml.org,2002:bool',
        RegExp(u'''^(?:true|True|TRUE|false|False|FALSE)$''', re.X),
        list(u'tTfF')),
    ([(1, 1)],
        u'tag:yaml.org,2002:bool',
        RegExp(u'''^(?:y|Y|yes|Yes|YES|n|N|no|No|NO
        |true|True|TRUE|false|False|FALSE
        |on|On|ON|off|Off|OFF)$''', re.X),
        list(u'yYnNtTfFoO')),
    ([(1, 2)],
        u'tag:yaml.org,2002:float',
        RegExp(u'''^(?:
         [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |[-+]?\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.')),
    ([(1, 1)],
        u'tag:yaml.org,2002:float',
        RegExp(u'''^(?:
         [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*  # sexagesimal float
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.')),
    ([(1, 2)],
        u'tag:yaml.org,2002:int',
        RegExp(u'''^(?:[-+]?0b[0-1_]+
        |[-+]?0o?[0-7_]+
        |[-+]?[0-9_]+
        |[-+]?0x[0-9a-fA-F_]+)$''', re.X),
        list(u'-+0123456789')),
    ([(1, 1)],
        u'tag:yaml.org,2002:int',
        RegExp(u'''^(?:[-+]?0b[0-1_]+
        |[-+]?0?[0-7_]+
        |[-+]?(?:0|[1-9][0-9_]*)
        |[-+]?0x[0-9a-fA-F_]+
        |[-+]?[1-9][0-9_]*(?::[0-5]?[0-9])+)$''', re.X),  # sexagesimal int
        list(u'-+0123456789')),
    ([(1, 2), (1, 1)],
        u'tag:yaml.org,2002:merge',
        RegExp(u'^(?:<<)$'),
        [u'<']),
    ([(1, 2), (1, 1)],
        u'tag:yaml.org,2002:null',
        RegExp(u'''^(?: ~
        |null|Null|NULL
        | )$''', re.X),
        [u'~', u'n', u'N', u'']),
    ([(1, 2), (1, 1)],
        u'tag:yaml.org,2002:timestamp',
        RegExp(u'''^(?:[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]
        |[0-9][0-9][0-9][0-9] -[0-9][0-9]? -[0-9][0-9]?
        (?:[Tt]|[ \\t]+)[0-9][0-9]?
        :[0-9][0-9] :[0-9][0-9] (?:\\.[0-9]*)?
        (?:[ \\t]*(?:Z|[-+][0-9][0-9]?(?::[0-9][0-9])?))?)$''', re.X),
        list(u'0123456789')),
    ([(1, 2), (1, 1)],
        u'tag:yaml.org,2002:value',
        RegExp(u'^(?:=)$'),
        [u'=']),
    # The following resolver is only for documentation purposes. It cannot work
    # because plain scalars cannot start with '!', '&', or '*'.
    ([(1, 2), (1, 1)],
        u'tag:yaml.org,2002:yaml',
        RegExp(u'^(?:!|&|\\*)$'),
        list(u'!&*')),
]
# fmt: on


class ResolverError(YAMLError):
    pass


class BaseResolver(object):

    DEFAULT_SCALAR_TAG = u'tag:yaml.org,2002:str'
    DEFAULT_SEQUENCE_TAG = u'tag:yaml.org,2002:seq'
    DEFAULT_MAPPING_TAG = u'tag:yaml.org,2002:map'

    yaml_implicit_resolvers = {}  # type: Dict[Any, Any]
    yaml_path_resolvers = {}  # type: Dict[Any, Any]

    def __init__(self, loadumper=None):
        # type: (Any, Any) -> None
        self.loadumper = loadumper
        if self.loadumper is not None and getattr(self.loadumper, '_resolver', None) is None:
            self.loadumper._resolver = self.loadumper
        self._loader_version = None  # type: Any
        self.resolver_exact_paths = []  # type: List[Any]
        self.resolver_prefix_paths = []  # type: List[Any]

    @property
    def parser(self):
        # type: () -> Any
        if self.loadumper is not None:
            if hasattr(self.loadumper, 'typ'):
                return self.loadumper.parser
            return self.loadumper._parser
        return None

    @classmethod
    def add_implicit_resolver_base(cls, tag, regexp, first):
        # type: (Any, Any, Any) -> None
        if 'yaml_implicit_resolvers' not in cls.__dict__:
            # deepcopy doesn't work here
            cls.yaml_implicit_resolvers = dict(
                (k, cls.yaml_implicit_resolvers[k][:]) for k in cls.yaml_implicit_resolvers
            )
        if first is None:
            first = [None]
        for ch in first:
            cls.yaml_implicit_resolvers.setdefault(ch, []).append((tag, regexp))

    @classmethod
    def add_implicit_resolver(cls, tag, regexp, first):
        # type: (Any, Any, Any) -> None
        if 'yaml_implicit_resolvers' not in cls.__dict__:
            # deepcopy doesn't work here
            cls.yaml_implicit_resolvers = dict(
                (k, cls.yaml_implicit_resolvers[k][:]) for k in cls.yaml_implicit_resolvers
            )
        if first is None:
            first = [None]
        for ch in first:
            cls.yaml_implicit_resolvers.setdefault(ch, []).append((tag, regexp))
        implicit_resolvers.append(([(1, 2), (1, 1)], tag, regexp, first))

    # @classmethod
    # def add_implicit_resolver(cls, tag, regexp, first):

    @classmethod
    def add_path_resolver(cls, tag, path, kind=None):
        # type: (Any, Any, Any) -> None
        # Note: `add_path_resolver` is experimental.  The API could be changed.
        # `new_path` is a pattern that is matched against the path from the
        # root to the node that is being considered.  `node_path` elements are
        # tuples `(node_check, index_check)`.  `node_check` is a node class:
        # `ScalarNode`, `SequenceNode`, `MappingNode` or `None`.  `None`
        # matches any kind of a node.  `index_check` could be `None`, a boolean
        # value, a string value, or a number.  `None` and `False` match against
        # any _value_ of sequence and mapping nodes.  `True` matches against
        # any _key_ of a mapping node.  A string `index_check` matches against
        # a mapping value that corresponds to a scalar key which content is
        # equal to the `index_check` value.  An integer `index_check` matches
        # against a sequence value with the index equal to `index_check`.
        if 'yaml_path_resolvers' not in cls.__dict__:
            cls.yaml_path_resolvers = cls.yaml_path_resolvers.copy()
        new_path = []  # type: List[Any]
        for element in path:
            if isinstance(element, (list, tuple)):
                if len(element) == 2:
                    node_check, index_check = element
                elif len(element) == 1:
                    node_check = element[0]
                    index_check = True
                else:
                    raise ResolverError('Invalid path element: %s' % (element,))
            else:
                node_check = None
                index_check = element
            if node_check is str:
                node_check = ScalarNode
            elif node_check is list:
                node_check = SequenceNode
            elif node_check is dict:
                node_check = MappingNode
            elif (
                node_check not in [ScalarNode, SequenceNode, MappingNode]
                and not isinstance(node_check, string_types)
                and node_check is not None
            ):
                raise ResolverError('Invalid node checker: %s' % (node_check,))
            if not isinstance(index_check, (string_types, int)) and index_check is not None:
                raise ResolverError('Invalid index checker: %s' % (index_check,))
            new_path.append((node_check, index_check))
        if kind is str:
            kind = ScalarNode
        elif kind is list:
            kind = SequenceNode
        elif kind is dict:
            kind = MappingNode
        elif kind not in [ScalarNode, SequenceNode, MappingNode] and kind is not None:
            raise ResolverError('Invalid node kind: %s' % (kind,))
        cls.yaml_path_resolvers[tuple(new_path), kind] = tag

    def descend_resolver(self, current_node, current_index):
        # type: (Any, Any) -> None
        if not self.yaml_path_resolvers:
            return
        exact_paths = {}
        prefix_paths = []
        if current_node:
            depth = len(self.resolver_prefix_paths)
            for path, kind in self.resolver_prefix_paths[-1]:
                if self.check_resolver_prefix(depth, path, kind, current_node, current_index):
                    if len(path) > depth:
                        prefix_paths.append((path, kind))
                    else:
                        exact_paths[kind] = self.yaml_path_resolvers[path, kind]
        else:
            for path, kind in self.yaml_path_resolvers:
                if not path:
                    exact_paths[kind] = self.yaml_path_resolvers[path, kind]
                else:
                    prefix_paths.append((path, kind))
        self.resolver_exact_paths.append(exact_paths)
        self.resolver_prefix_paths.append(prefix_paths)

    def ascend_resolver(self):
        # type: () -> None
        if not self.yaml_path_resolvers:
            return
        self.resolver_exact_paths.pop()
        self.resolver_prefix_paths.pop()

    def check_resolver_prefix(self, depth, path, kind, current_node, current_index):
        # type: (int, Text, Any, Any, Any) -> bool
        node_check, index_check = path[depth - 1]
        if isinstance(node_check, string_types):
            if current_node.tag != node_check:
                return False
        elif node_check is not None:
            if not isinstance(current_node, node_check):
                return False
        if index_check is True and current_index is not None:
            return False
        if (index_check is False or index_check is None) and current_index is None:
            return False
        if isinstance(index_check, string_types):
            if not (
                isinstance(current_index, ScalarNode) and index_check == current_index.value
            ):
                return False
        elif isinstance(index_check, int) and not isinstance(index_check, bool):
            if index_check != current_index:
                return False
        return True

    def resolve(self, kind, value, implicit):
        # type: (Any, Any, Any) -> Any
        if kind is ScalarNode and implicit[0]:
            if value == "":
                resolvers = self.yaml_implicit_resolvers.get("", [])
            else:
                resolvers = self.yaml_implicit_resolvers.get(value[0], [])
            resolvers += self.yaml_implicit_resolvers.get(None, [])
            for tag, regexp in resolvers:
                if regexp.match(value):
                    return tag
            implicit = implicit[1]
        if bool(self.yaml_path_resolvers):
            exact_paths = self.resolver_exact_paths[-1]
            if kind in exact_paths:
                return exact_paths[kind]
            if None in exact_paths:
                return exact_paths[None]
        if kind is ScalarNode:
            return self.DEFAULT_SCALAR_TAG
        elif kind is SequenceNode:
            return self.DEFAULT_SEQUENCE_TAG
        elif kind is MappingNode:
            return self.DEFAULT_MAPPING_TAG

    @property
    def processing_version(self):
        # type: () -> Any
        return None


class Resolver(BaseResolver):
    pass


for ir in implicit_resolvers:
    if (1, 2) in ir[0]:
        Resolver.add_implicit_resolver_base(*ir[1:])


class VersionedResolver(BaseResolver):
    """
    contrary to the "normal" resolver, the smart resolver delays loading
    the pattern matching rules. That way it can decide to load 1.1 rules
    or the (default) 1.2 rules, that no longer support octal without 0o, sexagesimals
    and Yes/No/On/Off booleans.
    """

    def __init__(self, version=None, loader=None, loadumper=None):
        # type: (Optional[VersionType], Any, Any) -> None
        if loader is None and loadumper is not None:
            loader = loadumper
        BaseResolver.__init__(self, loader)
        self._loader_version = self.get_loader_version(version)
        self._version_implicit_resolver = {}  # type: Dict[Any, Any]

    def add_version_implicit_resolver(self, version, tag, regexp, first):
        # type: (VersionType, Any, Any, Any) -> None
        if first is None:
            first = [None]
        impl_resolver = self._version_implicit_resolver.setdefault(version, {})
        for ch in first:
            impl_resolver.setdefault(ch, []).append((tag, regexp))

    def get_loader_version(self, version):
        # type: (Optional[VersionType]) -> Any
        if version is None or isinstance(version, tuple):
            return version
        if isinstance(version, list):
            return tuple(version)
        # assume string
        return tuple(map(int, version.split(u'.')))

    @property
    def versioned_resolver(self):
        # type: () -> Any
        """
        select the resolver based on the version we are parsing
        """
        version = self.processing_version
        if version not in self._version_implicit_resolver:
            for x in implicit_resolvers:
                if version in x[0]:
                    self.add_version_implicit_resolver(version, x[1], x[2], x[3])
        return self._version_implicit_resolver[version]

    def resolve(self, kind, value, implicit):
        # type: (Any, Any, Any) -> Any
        if kind is ScalarNode and implicit[0]:
            if value == "":
                resolvers = self.versioned_resolver.get("", [])
            else:
                resolvers = self.versioned_resolver.get(value[0], [])
            resolvers += self.versioned_resolver.get(None, [])
            for tag, regexp in resolvers:
                if regexp.match(value):
                    return tag
            implicit = implicit[1]
        if bool(self.yaml_path_resolvers):
            exact_paths = self.resolver_exact_paths[-1]
            if kind in exact_paths:
                return exact_paths[kind]
            if None in exact_paths:
                return exact_paths[None]
        if kind is ScalarNode:
            return self.DEFAULT_SCALAR_TAG
        elif kind is SequenceNode:
            return self.DEFAULT_SEQUENCE_TAG
        elif kind is MappingNode:
            return self.DEFAULT_MAPPING_TAG

    @property
    def processing_version(self):
        # type: () -> Any
        try:
            version = self.parser.yaml_version
        except AttributeError:
            try:
                if hasattr(self.loadumper, 'typ'):
                    version = self.loadumper.version
                else:
                    version = self.loadumper._serializer.use_version  # dumping
            except AttributeError:
                version = None
        if version is None:
            version = self._loader_version
            if version is None:
                version = _DEFAULT_YAML_VERSION
        return version
