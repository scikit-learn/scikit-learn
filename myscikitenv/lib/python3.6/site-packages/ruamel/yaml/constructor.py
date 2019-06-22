# coding: utf-8

from __future__ import print_function, absolute_import, division

import datetime
import base64
import binascii
import re
import sys
import types
import warnings

# fmt: off
from ruamel.yaml.error import (MarkedYAMLError, MarkedYAMLFutureWarning,
                               MantissaNoDotYAML1_1Warning)
from ruamel.yaml.nodes import *                               # NOQA
from ruamel.yaml.nodes import (SequenceNode, MappingNode, ScalarNode)
from ruamel.yaml.compat import (utf8, builtins_module, to_str, PY2, PY3,  # NOQA
                                ordereddict, text_type, nprint, nprintf, version_tnf,
                                Hashable, MutableSequence, MutableMapping)
from ruamel.yaml.comments import *                               # NOQA
from ruamel.yaml.comments import (CommentedMap, CommentedOrderedMap, CommentedSet,
                                  CommentedKeySeq, CommentedSeq, TaggedScalar,
                                  CommentedKeyMap)
from ruamel.yaml.scalarstring import (SingleQuotedScalarString, DoubleQuotedScalarString,
                                      LiteralScalarString, FoldedScalarString,
                                      PlainScalarString, ScalarString,)
from ruamel.yaml.scalarint import ScalarInt, BinaryInt, OctalInt, HexInt, HexCapsInt
from ruamel.yaml.scalarfloat import ScalarFloat
from ruamel.yaml.scalarbool import ScalarBoolean
from ruamel.yaml.timestamp import TimeStamp
from ruamel.yaml.util import RegExp

if False:  # MYPY
    from typing import Any, Dict, List, Set, Generator, Union, Optional  # NOQA


__all__ = ['BaseConstructor', 'SafeConstructor', 'Constructor',
           'ConstructorError', 'RoundTripConstructor']
# fmt: on


class ConstructorError(MarkedYAMLError):
    pass


class DuplicateKeyFutureWarning(MarkedYAMLFutureWarning):
    pass


class DuplicateKeyError(MarkedYAMLFutureWarning):
    pass


class BaseConstructor(object):

    yaml_constructors = {}  # type: Dict[Any, Any]
    yaml_multi_constructors = {}  # type: Dict[Any, Any]

    def __init__(self, preserve_quotes=None, loader=None):
        # type: (Optional[bool], Any) -> None
        self.loader = loader
        if self.loader is not None and getattr(self.loader, '_constructor', None) is None:
            self.loader._constructor = self
        self.loader = loader
        self.yaml_base_dict_type = dict
        self.yaml_base_list_type = list
        self.constructed_objects = {}  # type: Dict[Any, Any]
        self.recursive_objects = {}  # type: Dict[Any, Any]
        self.state_generators = []  # type: List[Any]
        self.deep_construct = False
        self._preserve_quotes = preserve_quotes
        self.allow_duplicate_keys = version_tnf((0, 15, 1), (0, 16))

    @property
    def composer(self):
        # type: () -> Any
        if hasattr(self.loader, 'typ'):
            return self.loader.composer
        try:
            return self.loader._composer
        except AttributeError:
            sys.stdout.write('slt {}\n'.format(type(self)))
            sys.stdout.write('slc {}\n'.format(self.loader._composer))
            sys.stdout.write('{}\n'.format(dir(self)))
            raise

    @property
    def resolver(self):
        # type: () -> Any
        if hasattr(self.loader, 'typ'):
            return self.loader.resolver
        return self.loader._resolver

    def check_data(self):
        # type: () -> Any
        # If there are more documents available?
        return self.composer.check_node()

    def get_data(self):
        # type: () -> Any
        # Construct and return the next document.
        if self.composer.check_node():
            return self.construct_document(self.composer.get_node())

    def get_single_data(self):
        # type: () -> Any
        # Ensure that the stream contains a single document and construct it.
        node = self.composer.get_single_node()
        if node is not None:
            return self.construct_document(node)
        return None

    def construct_document(self, node):
        # type: (Any) -> Any
        data = self.construct_object(node)
        while bool(self.state_generators):
            state_generators = self.state_generators
            self.state_generators = []
            for generator in state_generators:
                for _dummy in generator:
                    pass
        self.constructed_objects = {}
        self.recursive_objects = {}
        self.deep_construct = False
        return data

    def construct_object(self, node, deep=False):
        # type: (Any, bool) -> Any
        """deep is True when creating an object/mapping recursively,
        in that case want the underlying elements available during construction
        """
        if node in self.constructed_objects:
            return self.constructed_objects[node]
        if deep:
            old_deep = self.deep_construct
            self.deep_construct = True
        if node in self.recursive_objects:
            return self.recursive_objects[node]
            # raise ConstructorError(
            #     None, None, 'found unconstructable recursive node', node.start_mark
            # )
        self.recursive_objects[node] = None
        constructor = None  # type: Any
        tag_suffix = None
        if node.tag in self.yaml_constructors:
            constructor = self.yaml_constructors[node.tag]
        else:
            for tag_prefix in self.yaml_multi_constructors:
                if node.tag.startswith(tag_prefix):
                    tag_suffix = node.tag[len(tag_prefix) :]
                    constructor = self.yaml_multi_constructors[tag_prefix]
                    break
            else:
                if None in self.yaml_multi_constructors:
                    tag_suffix = node.tag
                    constructor = self.yaml_multi_constructors[None]
                elif None in self.yaml_constructors:
                    constructor = self.yaml_constructors[None]
                elif isinstance(node, ScalarNode):
                    constructor = self.__class__.construct_scalar
                elif isinstance(node, SequenceNode):
                    constructor = self.__class__.construct_sequence
                elif isinstance(node, MappingNode):
                    constructor = self.__class__.construct_mapping
        if tag_suffix is None:
            data = constructor(self, node)
        else:
            data = constructor(self, tag_suffix, node)
        if isinstance(data, types.GeneratorType):
            generator = data
            data = next(generator)
            if self.deep_construct:
                for _dummy in generator:
                    pass
            else:
                self.state_generators.append(generator)
        self.constructed_objects[node] = data
        del self.recursive_objects[node]
        if deep:
            self.deep_construct = old_deep
        return data

    def construct_scalar(self, node):
        # type: (Any) -> Any
        if not isinstance(node, ScalarNode):
            raise ConstructorError(
                None, None, 'expected a scalar node, but found %s' % node.id, node.start_mark
            )
        return node.value

    def construct_sequence(self, node, deep=False):
        # type: (Any, bool) -> Any
        """deep is True when creating an object/mapping recursively,
        in that case want the underlying elements available during construction
        """
        if not isinstance(node, SequenceNode):
            raise ConstructorError(
                None, None, 'expected a sequence node, but found %s' % node.id, node.start_mark
            )
        return [self.construct_object(child, deep=deep) for child in node.value]

    def construct_mapping(self, node, deep=False):
        # type: (Any, bool) -> Any
        """deep is True when creating an object/mapping recursively,
        in that case want the underlying elements available during construction
        """
        if not isinstance(node, MappingNode):
            raise ConstructorError(
                None, None, 'expected a mapping node, but found %s' % node.id, node.start_mark
            )
        total_mapping = self.yaml_base_dict_type()
        if getattr(node, 'merge', None) is not None:
            todo = [(node.merge, False), (node.value, False)]
        else:
            todo = [(node.value, True)]
        for values, check in todo:
            mapping = self.yaml_base_dict_type()  # type: Dict[Any, Any]
            for key_node, value_node in values:
                # keys can be list -> deep
                key = self.construct_object(key_node, deep=True)
                # lists are not hashable, but tuples are
                if not isinstance(key, Hashable):
                    if isinstance(key, list):
                        key = tuple(key)
                if PY2:
                    try:
                        hash(key)
                    except TypeError as exc:
                        raise ConstructorError(
                            'while constructing a mapping',
                            node.start_mark,
                            'found unacceptable key (%s)' % exc,
                            key_node.start_mark,
                        )
                else:
                    if not isinstance(key, Hashable):
                        raise ConstructorError(
                            'while constructing a mapping',
                            node.start_mark,
                            'found unhashable key',
                            key_node.start_mark,
                        )

                value = self.construct_object(value_node, deep=deep)
                if check:
                    if self.check_mapping_key(node, key_node, mapping, key, value):
                        mapping[key] = value
                else:
                    mapping[key] = value
            total_mapping.update(mapping)
        return total_mapping

    def check_mapping_key(self, node, key_node, mapping, key, value):
        # type: (Any, Any, Any, Any, Any) -> bool
        """return True if key is unique"""
        if key in mapping:
            if not self.allow_duplicate_keys:
                mk = mapping.get(key)
                if PY2:
                    if isinstance(key, unicode):
                        key = key.encode('utf-8')
                    if isinstance(value, unicode):
                        value = value.encode('utf-8')
                    if isinstance(mk, unicode):
                        mk = mk.encode('utf-8')
                args = [
                    'while constructing a mapping',
                    node.start_mark,
                    'found duplicate key "{}" with value "{}" '
                    '(original value: "{}")'.format(key, value, mk),
                    key_node.start_mark,
                    """
                    To suppress this check see:
                        http://yaml.readthedocs.io/en/latest/api.html#duplicate-keys
                    """,
                    """\
                    Duplicate keys will become an error in future releases, and are errors
                    by default when using the new API.
                    """,
                ]
                if self.allow_duplicate_keys is None:
                    warnings.warn(DuplicateKeyFutureWarning(*args))
                else:
                    raise DuplicateKeyError(*args)
            return False
        return True

    def check_set_key(self, node, key_node, setting, key):
        # type: (Any, Any, Any, Any, Any) -> None
        if key in setting:
            if not self.allow_duplicate_keys:
                if PY2:
                    if isinstance(key, unicode):
                        key = key.encode('utf-8')
                args = [
                    'while constructing a set',
                    node.start_mark,
                    'found duplicate key "{}"'.format(key),
                    key_node.start_mark,
                    """
                    To suppress this check see:
                        http://yaml.readthedocs.io/en/latest/api.html#duplicate-keys
                    """,
                    """\
                    Duplicate keys will become an error in future releases, and are errors
                    by default when using the new API.
                    """,
                ]
                if self.allow_duplicate_keys is None:
                    warnings.warn(DuplicateKeyFutureWarning(*args))
                else:
                    raise DuplicateKeyError(*args)

    def construct_pairs(self, node, deep=False):
        # type: (Any, bool) -> Any
        if not isinstance(node, MappingNode):
            raise ConstructorError(
                None, None, 'expected a mapping node, but found %s' % node.id, node.start_mark
            )
        pairs = []
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            value = self.construct_object(value_node, deep=deep)
            pairs.append((key, value))
        return pairs

    @classmethod
    def add_constructor(cls, tag, constructor):
        # type: (Any, Any) -> None
        if 'yaml_constructors' not in cls.__dict__:
            cls.yaml_constructors = cls.yaml_constructors.copy()
        cls.yaml_constructors[tag] = constructor

    @classmethod
    def add_multi_constructor(cls, tag_prefix, multi_constructor):
        # type: (Any, Any) -> None
        if 'yaml_multi_constructors' not in cls.__dict__:
            cls.yaml_multi_constructors = cls.yaml_multi_constructors.copy()
        cls.yaml_multi_constructors[tag_prefix] = multi_constructor


class SafeConstructor(BaseConstructor):
    def construct_scalar(self, node):
        # type: (Any) -> Any
        if isinstance(node, MappingNode):
            for key_node, value_node in node.value:
                if key_node.tag == u'tag:yaml.org,2002:value':
                    return self.construct_scalar(value_node)
        return BaseConstructor.construct_scalar(self, node)

    def flatten_mapping(self, node):
        # type: (Any) -> Any
        """
        This implements the merge key feature http://yaml.org/type/merge.html
        by inserting keys from the merge dict/list of dicts if not yet
        available in this node
        """
        merge = []  # type: List[Any]
        index = 0
        while index < len(node.value):
            key_node, value_node = node.value[index]
            if key_node.tag == u'tag:yaml.org,2002:merge':
                if merge:  # double << key
                    if self.allow_duplicate_keys:
                        del node.value[index]
                        index += 1
                        continue
                    args = [
                        'while constructing a mapping',
                        node.start_mark,
                        'found duplicate key "{}"'.format(key_node.value),
                        key_node.start_mark,
                        """
                        To suppress this check see:
                           http://yaml.readthedocs.io/en/latest/api.html#duplicate-keys
                        """,
                        """\
                        Duplicate keys will become an error in future releases, and are errors
                        by default when using the new API.
                        """,
                    ]
                    if self.allow_duplicate_keys is None:
                        warnings.warn(DuplicateKeyFutureWarning(*args))
                    else:
                        raise DuplicateKeyError(*args)
                del node.value[index]
                if isinstance(value_node, MappingNode):
                    self.flatten_mapping(value_node)
                    merge.extend(value_node.value)
                elif isinstance(value_node, SequenceNode):
                    submerge = []
                    for subnode in value_node.value:
                        if not isinstance(subnode, MappingNode):
                            raise ConstructorError(
                                'while constructing a mapping',
                                node.start_mark,
                                'expected a mapping for merging, but found %s' % subnode.id,
                                subnode.start_mark,
                            )
                        self.flatten_mapping(subnode)
                        submerge.append(subnode.value)
                    submerge.reverse()
                    for value in submerge:
                        merge.extend(value)
                else:
                    raise ConstructorError(
                        'while constructing a mapping',
                        node.start_mark,
                        'expected a mapping or list of mappings for merging, '
                        'but found %s' % value_node.id,
                        value_node.start_mark,
                    )
            elif key_node.tag == u'tag:yaml.org,2002:value':
                key_node.tag = u'tag:yaml.org,2002:str'
                index += 1
            else:
                index += 1
        if bool(merge):
            node.merge = merge  # separate merge keys to be able to update without duplicate
            node.value = merge + node.value

    def construct_mapping(self, node, deep=False):
        # type: (Any, bool) -> Any
        """deep is True when creating an object/mapping recursively,
        in that case want the underlying elements available during construction
        """
        if isinstance(node, MappingNode):
            self.flatten_mapping(node)
        return BaseConstructor.construct_mapping(self, node, deep=deep)

    def construct_yaml_null(self, node):
        # type: (Any) -> Any
        self.construct_scalar(node)
        return None

    # YAML 1.2 spec doesn't mention yes/no etc any more, 1.1 does
    bool_values = {
        u'yes': True,
        u'no': False,
        u'y': True,
        u'n': False,
        u'true': True,
        u'false': False,
        u'on': True,
        u'off': False,
    }

    def construct_yaml_bool(self, node):
        # type: (Any) -> bool
        value = self.construct_scalar(node)
        return self.bool_values[value.lower()]

    def construct_yaml_int(self, node):
        # type: (Any) -> int
        value_s = to_str(self.construct_scalar(node))
        value_s = value_s.replace('_', "")
        sign = +1
        if value_s[0] == '-':
            sign = -1
        if value_s[0] in '+-':
            value_s = value_s[1:]
        if value_s == '0':
            return 0
        elif value_s.startswith('0b'):
            return sign * int(value_s[2:], 2)
        elif value_s.startswith('0x'):
            return sign * int(value_s[2:], 16)
        elif value_s.startswith('0o'):
            return sign * int(value_s[2:], 8)
        elif self.resolver.processing_version == (1, 1) and value_s[0] == '0':
            return sign * int(value_s, 8)
        elif self.resolver.processing_version == (1, 1) and ':' in value_s:
            digits = [int(part) for part in value_s.split(':')]
            digits.reverse()
            base = 1
            value = 0
            for digit in digits:
                value += digit * base
                base *= 60
            return sign * value
        else:
            return sign * int(value_s)

    inf_value = 1e300
    while inf_value != inf_value * inf_value:
        inf_value *= inf_value
    nan_value = -inf_value / inf_value  # Trying to make a quiet NaN (like C99).

    def construct_yaml_float(self, node):
        # type: (Any) -> float
        value_so = to_str(self.construct_scalar(node))
        value_s = value_so.replace('_', "").lower()
        sign = +1
        if value_s[0] == '-':
            sign = -1
        if value_s[0] in '+-':
            value_s = value_s[1:]
        if value_s == '.inf':
            return sign * self.inf_value
        elif value_s == '.nan':
            return self.nan_value
        elif self.resolver.processing_version != (1, 2) and ':' in value_s:
            digits = [float(part) for part in value_s.split(':')]
            digits.reverse()
            base = 1
            value = 0.0
            for digit in digits:
                value += digit * base
                base *= 60
            return sign * value
        else:
            if self.resolver.processing_version != (1, 2) and 'e' in value_s:
                # value_s is lower case independent of input
                mantissa, exponent = value_s.split('e')
                if '.' not in mantissa:
                    warnings.warn(MantissaNoDotYAML1_1Warning(node, value_so))
            return sign * float(value_s)

    if PY3:

        def construct_yaml_binary(self, node):
            # type: (Any) -> Any
            try:
                value = self.construct_scalar(node).encode('ascii')
            except UnicodeEncodeError as exc:
                raise ConstructorError(
                    None,
                    None,
                    'failed to convert base64 data into ascii: %s' % exc,
                    node.start_mark,
                )
            try:
                if hasattr(base64, 'decodebytes'):
                    return base64.decodebytes(value)
                else:
                    return base64.decodestring(value)
            except binascii.Error as exc:
                raise ConstructorError(
                    None, None, 'failed to decode base64 data: %s' % exc, node.start_mark
                )

    else:

        def construct_yaml_binary(self, node):
            # type: (Any) -> Any
            value = self.construct_scalar(node)
            try:
                return to_str(value).decode('base64')
            except (binascii.Error, UnicodeEncodeError) as exc:
                raise ConstructorError(
                    None, None, 'failed to decode base64 data: %s' % exc, node.start_mark
                )

    timestamp_regexp = RegExp(
        u"""^(?P<year>[0-9][0-9][0-9][0-9])
          -(?P<month>[0-9][0-9]?)
          -(?P<day>[0-9][0-9]?)
          (?:((?P<t>[Tt])|[ \\t]+)   # explictly not retaining extra spaces
          (?P<hour>[0-9][0-9]?)
          :(?P<minute>[0-9][0-9])
          :(?P<second>[0-9][0-9])
          (?:\\.(?P<fraction>[0-9]*))?
          (?:[ \\t]*(?P<tz>Z|(?P<tz_sign>[-+])(?P<tz_hour>[0-9][0-9]?)
          (?::(?P<tz_minute>[0-9][0-9]))?))?)?$""",
        re.X,
    )

    def construct_yaml_timestamp(self, node, values=None):
        # type: (Any, Any) -> Any
        if values is None:
            try:
                match = self.timestamp_regexp.match(node.value)
            except TypeError:
                match = None
            if match is None:
                raise ConstructorError(
                    None,
                    None,
                    'failed to construct timestamp from "{}"'.format(node.value),
                    node.start_mark,
                )
            values = match.groupdict()
        year = int(values['year'])
        month = int(values['month'])
        day = int(values['day'])
        if not values['hour']:
            return datetime.date(year, month, day)
        hour = int(values['hour'])
        minute = int(values['minute'])
        second = int(values['second'])
        fraction = 0
        if values['fraction']:
            fraction_s = values['fraction'][:6]
            while len(fraction_s) < 6:
                fraction_s += '0'
            fraction = int(fraction_s)
            if len(values['fraction']) > 6 and int(values['fraction'][6]) > 4:
                fraction += 1
        delta = None
        if values['tz_sign']:
            tz_hour = int(values['tz_hour'])
            minutes = values['tz_minute']
            tz_minute = int(minutes) if minutes else 0
            delta = datetime.timedelta(hours=tz_hour, minutes=tz_minute)
            if values['tz_sign'] == '-':
                delta = -delta
        # should do something else instead (or hook this up to the preceding if statement
        # in reverse
        #  if delta is None:
        #      return datetime.datetime(year, month, day, hour, minute, second, fraction)
        #  return datetime.datetime(year, month, day, hour, minute, second, fraction,
        #                           datetime.timezone.utc)
        # the above is not good enough though, should provide tzinfo. In Python3 that is easily
        # doable drop that kind of support for Python2 as it has not native tzinfo
        data = datetime.datetime(year, month, day, hour, minute, second, fraction)
        if delta:
            data -= delta
        return data

    def construct_yaml_omap(self, node):
        # type: (Any) -> Any
        # Note: we do now check for duplicate keys
        omap = ordereddict()
        yield omap
        if not isinstance(node, SequenceNode):
            raise ConstructorError(
                'while constructing an ordered map',
                node.start_mark,
                'expected a sequence, but found %s' % node.id,
                node.start_mark,
            )
        for subnode in node.value:
            if not isinstance(subnode, MappingNode):
                raise ConstructorError(
                    'while constructing an ordered map',
                    node.start_mark,
                    'expected a mapping of length 1, but found %s' % subnode.id,
                    subnode.start_mark,
                )
            if len(subnode.value) != 1:
                raise ConstructorError(
                    'while constructing an ordered map',
                    node.start_mark,
                    'expected a single mapping item, but found %d items' % len(subnode.value),
                    subnode.start_mark,
                )
            key_node, value_node = subnode.value[0]
            key = self.construct_object(key_node)
            assert key not in omap
            value = self.construct_object(value_node)
            omap[key] = value

    def construct_yaml_pairs(self, node):
        # type: (Any) -> Any
        # Note: the same code as `construct_yaml_omap`.
        pairs = []  # type: List[Any]
        yield pairs
        if not isinstance(node, SequenceNode):
            raise ConstructorError(
                'while constructing pairs',
                node.start_mark,
                'expected a sequence, but found %s' % node.id,
                node.start_mark,
            )
        for subnode in node.value:
            if not isinstance(subnode, MappingNode):
                raise ConstructorError(
                    'while constructing pairs',
                    node.start_mark,
                    'expected a mapping of length 1, but found %s' % subnode.id,
                    subnode.start_mark,
                )
            if len(subnode.value) != 1:
                raise ConstructorError(
                    'while constructing pairs',
                    node.start_mark,
                    'expected a single mapping item, but found %d items' % len(subnode.value),
                    subnode.start_mark,
                )
            key_node, value_node = subnode.value[0]
            key = self.construct_object(key_node)
            value = self.construct_object(value_node)
            pairs.append((key, value))

    def construct_yaml_set(self, node):
        # type: (Any) -> Any
        data = set()  # type: Set[Any]
        yield data
        value = self.construct_mapping(node)
        data.update(value)

    def construct_yaml_str(self, node):
        # type: (Any) -> Any
        value = self.construct_scalar(node)
        if PY3:
            return value
        try:
            return value.encode('ascii')
        except UnicodeEncodeError:
            return value

    def construct_yaml_seq(self, node):
        # type: (Any) -> Any
        data = self.yaml_base_list_type()  # type: List[Any]
        yield data
        data.extend(self.construct_sequence(node))

    def construct_yaml_map(self, node):
        # type: (Any) -> Any
        data = self.yaml_base_dict_type()  # type: Dict[Any, Any]
        yield data
        value = self.construct_mapping(node)
        data.update(value)

    def construct_yaml_object(self, node, cls):
        # type: (Any, Any) -> Any
        data = cls.__new__(cls)
        yield data
        if hasattr(data, '__setstate__'):
            state = self.construct_mapping(node, deep=True)
            data.__setstate__(state)
        else:
            state = self.construct_mapping(node)
            data.__dict__.update(state)

    def construct_undefined(self, node):
        # type: (Any) -> None
        raise ConstructorError(
            None,
            None,
            'could not determine a constructor for the tag %r' % utf8(node.tag),
            node.start_mark,
        )


SafeConstructor.add_constructor(u'tag:yaml.org,2002:null', SafeConstructor.construct_yaml_null)

SafeConstructor.add_constructor(u'tag:yaml.org,2002:bool', SafeConstructor.construct_yaml_bool)

SafeConstructor.add_constructor(u'tag:yaml.org,2002:int', SafeConstructor.construct_yaml_int)

SafeConstructor.add_constructor(
    u'tag:yaml.org,2002:float', SafeConstructor.construct_yaml_float
)

SafeConstructor.add_constructor(
    u'tag:yaml.org,2002:binary', SafeConstructor.construct_yaml_binary
)

SafeConstructor.add_constructor(
    u'tag:yaml.org,2002:timestamp', SafeConstructor.construct_yaml_timestamp
)

SafeConstructor.add_constructor(u'tag:yaml.org,2002:omap', SafeConstructor.construct_yaml_omap)

SafeConstructor.add_constructor(
    u'tag:yaml.org,2002:pairs', SafeConstructor.construct_yaml_pairs
)

SafeConstructor.add_constructor(u'tag:yaml.org,2002:set', SafeConstructor.construct_yaml_set)

SafeConstructor.add_constructor(u'tag:yaml.org,2002:str', SafeConstructor.construct_yaml_str)

SafeConstructor.add_constructor(u'tag:yaml.org,2002:seq', SafeConstructor.construct_yaml_seq)

SafeConstructor.add_constructor(u'tag:yaml.org,2002:map', SafeConstructor.construct_yaml_map)

SafeConstructor.add_constructor(None, SafeConstructor.construct_undefined)

if PY2:

    class classobj:
        pass


class Constructor(SafeConstructor):
    def construct_python_str(self, node):
        # type: (Any) -> Any
        return utf8(self.construct_scalar(node))

    def construct_python_unicode(self, node):
        # type: (Any) -> Any
        return self.construct_scalar(node)

    if PY3:

        def construct_python_bytes(self, node):
            # type: (Any) -> Any
            try:
                value = self.construct_scalar(node).encode('ascii')
            except UnicodeEncodeError as exc:
                raise ConstructorError(
                    None,
                    None,
                    'failed to convert base64 data into ascii: %s' % exc,
                    node.start_mark,
                )
            try:
                if hasattr(base64, 'decodebytes'):
                    return base64.decodebytes(value)
                else:
                    return base64.decodestring(value)
            except binascii.Error as exc:
                raise ConstructorError(
                    None, None, 'failed to decode base64 data: %s' % exc, node.start_mark
                )

    def construct_python_long(self, node):
        # type: (Any) -> int
        val = self.construct_yaml_int(node)
        if PY3:
            return val
        return int(val)

    def construct_python_complex(self, node):
        # type: (Any) -> Any
        return complex(self.construct_scalar(node))

    def construct_python_tuple(self, node):
        # type: (Any) -> Any
        return tuple(self.construct_sequence(node))

    def find_python_module(self, name, mark):
        # type: (Any, Any) -> Any
        if not name:
            raise ConstructorError(
                'while constructing a Python module',
                mark,
                'expected non-empty name appended to the tag',
                mark,
            )
        try:
            __import__(name)
        except ImportError as exc:
            raise ConstructorError(
                'while constructing a Python module',
                mark,
                'cannot find module %r (%s)' % (utf8(name), exc),
                mark,
            )
        return sys.modules[name]

    def find_python_name(self, name, mark):
        # type: (Any, Any) -> Any
        if not name:
            raise ConstructorError(
                'while constructing a Python object',
                mark,
                'expected non-empty name appended to the tag',
                mark,
            )
        if u'.' in name:
            lname = name.split('.')
            lmodule_name = lname
            lobject_name = []  # type: List[Any]
            while len(lmodule_name) > 1:
                lobject_name.insert(0, lmodule_name.pop())
                module_name = '.'.join(lmodule_name)
                try:
                    __import__(module_name)
                    # object_name = '.'.join(object_name)
                    break
                except ImportError:
                    continue
        else:
            module_name = builtins_module
            lobject_name = [name]
        try:
            __import__(module_name)
        except ImportError as exc:
            raise ConstructorError(
                'while constructing a Python object',
                mark,
                'cannot find module %r (%s)' % (utf8(module_name), exc),
                mark,
            )
        module = sys.modules[module_name]
        object_name = '.'.join(lobject_name)
        obj = module
        while lobject_name:
            if not hasattr(obj, lobject_name[0]):

                raise ConstructorError(
                    'while constructing a Python object',
                    mark,
                    'cannot find %r in the module %r' % (utf8(object_name), module.__name__),
                    mark,
                )
            obj = getattr(obj, lobject_name.pop(0))
        return obj

    def construct_python_name(self, suffix, node):
        # type: (Any, Any) -> Any
        value = self.construct_scalar(node)
        if value:
            raise ConstructorError(
                'while constructing a Python name',
                node.start_mark,
                'expected the empty value, but found %r' % utf8(value),
                node.start_mark,
            )
        return self.find_python_name(suffix, node.start_mark)

    def construct_python_module(self, suffix, node):
        # type: (Any, Any) -> Any
        value = self.construct_scalar(node)
        if value:
            raise ConstructorError(
                'while constructing a Python module',
                node.start_mark,
                'expected the empty value, but found %r' % utf8(value),
                node.start_mark,
            )
        return self.find_python_module(suffix, node.start_mark)

    def make_python_instance(self, suffix, node, args=None, kwds=None, newobj=False):
        # type: (Any, Any, Any, Any, bool) -> Any
        if not args:
            args = []
        if not kwds:
            kwds = {}
        cls = self.find_python_name(suffix, node.start_mark)
        if PY3:
            if newobj and isinstance(cls, type):
                return cls.__new__(cls, *args, **kwds)
            else:
                return cls(*args, **kwds)
        else:
            if newobj and isinstance(cls, type(classobj)) and not args and not kwds:
                instance = classobj()
                instance.__class__ = cls
                return instance
            elif newobj and isinstance(cls, type):
                return cls.__new__(cls, *args, **kwds)
            else:
                return cls(*args, **kwds)

    def set_python_instance_state(self, instance, state):
        # type: (Any, Any) -> None
        if hasattr(instance, '__setstate__'):
            instance.__setstate__(state)
        else:
            slotstate = {}  # type: Dict[Any, Any]
            if isinstance(state, tuple) and len(state) == 2:
                state, slotstate = state
            if hasattr(instance, '__dict__'):
                instance.__dict__.update(state)
            elif state:
                slotstate.update(state)
            for key, value in slotstate.items():
                setattr(object, key, value)

    def construct_python_object(self, suffix, node):
        # type: (Any, Any) -> Any
        # Format:
        #   !!python/object:module.name { ... state ... }
        instance = self.make_python_instance(suffix, node, newobj=True)
        self.recursive_objects[node] = instance
        yield instance
        deep = hasattr(instance, '__setstate__')
        state = self.construct_mapping(node, deep=deep)
        self.set_python_instance_state(instance, state)

    def construct_python_object_apply(self, suffix, node, newobj=False):
        # type: (Any, Any, bool) -> Any
        # Format:
        #   !!python/object/apply       # (or !!python/object/new)
        #   args: [ ... arguments ... ]
        #   kwds: { ... keywords ... }
        #   state: ... state ...
        #   listitems: [ ... listitems ... ]
        #   dictitems: { ... dictitems ... }
        # or short format:
        #   !!python/object/apply [ ... arguments ... ]
        # The difference between !!python/object/apply and !!python/object/new
        # is how an object is created, check make_python_instance for details.
        if isinstance(node, SequenceNode):
            args = self.construct_sequence(node, deep=True)
            kwds = {}  # type: Dict[Any, Any]
            state = {}  # type: Dict[Any, Any]
            listitems = []  # type: List[Any]
            dictitems = {}  # type: Dict[Any, Any]
        else:
            value = self.construct_mapping(node, deep=True)
            args = value.get('args', [])
            kwds = value.get('kwds', {})
            state = value.get('state', {})
            listitems = value.get('listitems', [])
            dictitems = value.get('dictitems', {})
        instance = self.make_python_instance(suffix, node, args, kwds, newobj)
        if bool(state):
            self.set_python_instance_state(instance, state)
        if bool(listitems):
            instance.extend(listitems)
        if bool(dictitems):
            for key in dictitems:
                instance[key] = dictitems[key]
        return instance

    def construct_python_object_new(self, suffix, node):
        # type: (Any, Any) -> Any
        return self.construct_python_object_apply(suffix, node, newobj=True)


Constructor.add_constructor(u'tag:yaml.org,2002:python/none', Constructor.construct_yaml_null)

Constructor.add_constructor(u'tag:yaml.org,2002:python/bool', Constructor.construct_yaml_bool)

Constructor.add_constructor(u'tag:yaml.org,2002:python/str', Constructor.construct_python_str)

Constructor.add_constructor(
    u'tag:yaml.org,2002:python/unicode', Constructor.construct_python_unicode
)

if PY3:
    Constructor.add_constructor(
        u'tag:yaml.org,2002:python/bytes', Constructor.construct_python_bytes
    )

Constructor.add_constructor(u'tag:yaml.org,2002:python/int', Constructor.construct_yaml_int)

Constructor.add_constructor(
    u'tag:yaml.org,2002:python/long', Constructor.construct_python_long
)

Constructor.add_constructor(
    u'tag:yaml.org,2002:python/float', Constructor.construct_yaml_float
)

Constructor.add_constructor(
    u'tag:yaml.org,2002:python/complex', Constructor.construct_python_complex
)

Constructor.add_constructor(u'tag:yaml.org,2002:python/list', Constructor.construct_yaml_seq)

Constructor.add_constructor(
    u'tag:yaml.org,2002:python/tuple', Constructor.construct_python_tuple
)

Constructor.add_constructor(u'tag:yaml.org,2002:python/dict', Constructor.construct_yaml_map)

Constructor.add_multi_constructor(
    u'tag:yaml.org,2002:python/name:', Constructor.construct_python_name
)

Constructor.add_multi_constructor(
    u'tag:yaml.org,2002:python/module:', Constructor.construct_python_module
)

Constructor.add_multi_constructor(
    u'tag:yaml.org,2002:python/object:', Constructor.construct_python_object
)

Constructor.add_multi_constructor(
    u'tag:yaml.org,2002:python/object/apply:', Constructor.construct_python_object_apply
)

Constructor.add_multi_constructor(
    u'tag:yaml.org,2002:python/object/new:', Constructor.construct_python_object_new
)


class RoundTripConstructor(SafeConstructor):
    """need to store the comments on the node itself,
    as well as on the items
    """

    def construct_scalar(self, node):
        # type: (Any) -> Any
        if not isinstance(node, ScalarNode):
            raise ConstructorError(
                None, None, 'expected a scalar node, but found %s' % node.id, node.start_mark
            )

        if node.style == '|' and isinstance(node.value, text_type):
            lss = LiteralScalarString(node.value, anchor=node.anchor)
            if node.comment and node.comment[1]:
                lss.comment = node.comment[1][0]  # type: ignore
            return lss
        if node.style == '>' and isinstance(node.value, text_type):
            fold_positions = []  # type: List[int]
            idx = -1
            while True:
                idx = node.value.find('\a', idx + 1)
                if idx < 0:
                    break
                fold_positions.append(idx - len(fold_positions))
            fss = FoldedScalarString(node.value.replace('\a', ''), anchor=node.anchor)
            if node.comment and node.comment[1]:
                fss.comment = node.comment[1][0]  # type: ignore
            if fold_positions:
                fss.fold_pos = fold_positions  # type: ignore
            return fss
        elif bool(self._preserve_quotes) and isinstance(node.value, text_type):
            if node.style == "'":
                return SingleQuotedScalarString(node.value, anchor=node.anchor)
            if node.style == '"':
                return DoubleQuotedScalarString(node.value, anchor=node.anchor)
        if node.anchor:
            return PlainScalarString(node.value, anchor=node.anchor)
        return node.value

    def construct_yaml_int(self, node):
        # type: (Any) -> Any
        width = None  # type: Any
        value_su = to_str(self.construct_scalar(node))
        try:
            sx = value_su.rstrip('_')
            underscore = [len(sx) - sx.rindex('_') - 1, False, False]  # type: Any
        except ValueError:
            underscore = None
        except IndexError:
            underscore = None
        value_s = value_su.replace('_', "")
        sign = +1
        if value_s[0] == '-':
            sign = -1
        if value_s[0] in '+-':
            value_s = value_s[1:]
        if value_s == '0':
            return 0
        elif value_s.startswith('0b'):
            if self.resolver.processing_version > (1, 1) and value_s[2] == '0':
                width = len(value_s[2:])
            if underscore is not None:
                underscore[1] = value_su[2] == '_'
                underscore[2] = len(value_su[2:]) > 1 and value_su[-1] == '_'
            return BinaryInt(
                sign * int(value_s[2:], 2),
                width=width,
                underscore=underscore,
                anchor=node.anchor,
            )
        elif value_s.startswith('0x'):
            # default to lower-case if no a-fA-F in string
            if self.resolver.processing_version > (1, 1) and value_s[2] == '0':
                width = len(value_s[2:])
            hex_fun = HexInt  # type: Any
            for ch in value_s[2:]:
                if ch in 'ABCDEF':  # first non-digit is capital
                    hex_fun = HexCapsInt
                    break
                if ch in 'abcdef':
                    break
            if underscore is not None:
                underscore[1] = value_su[2] == '_'
                underscore[2] = len(value_su[2:]) > 1 and value_su[-1] == '_'
            return hex_fun(
                sign * int(value_s[2:], 16),
                width=width,
                underscore=underscore,
                anchor=node.anchor,
            )
        elif value_s.startswith('0o'):
            if self.resolver.processing_version > (1, 1) and value_s[2] == '0':
                width = len(value_s[2:])
            if underscore is not None:
                underscore[1] = value_su[2] == '_'
                underscore[2] = len(value_su[2:]) > 1 and value_su[-1] == '_'
            return OctalInt(
                sign * int(value_s[2:], 8),
                width=width,
                underscore=underscore,
                anchor=node.anchor,
            )
        elif self.resolver.processing_version != (1, 2) and value_s[0] == '0':
            return sign * int(value_s, 8)
        elif self.resolver.processing_version != (1, 2) and ':' in value_s:
            digits = [int(part) for part in value_s.split(':')]
            digits.reverse()
            base = 1
            value = 0
            for digit in digits:
                value += digit * base
                base *= 60
            return sign * value
        elif self.resolver.processing_version > (1, 1) and value_s[0] == '0':
            # not an octal, an integer with leading zero(s)
            if underscore is not None:
                # cannot have a leading underscore
                underscore[2] = len(value_su) > 1 and value_su[-1] == '_'
            return ScalarInt(sign * int(value_s), width=len(value_s), underscore=underscore)
        elif underscore:
            # cannot have a leading underscore
            underscore[2] = len(value_su) > 1 and value_su[-1] == '_'
            return ScalarInt(
                sign * int(value_s), width=None, underscore=underscore, anchor=node.anchor
            )
        elif node.anchor:
            return ScalarInt(sign * int(value_s), width=None, anchor=node.anchor)
        else:
            return sign * int(value_s)

    def construct_yaml_float(self, node):
        # type: (Any) -> Any
        def leading_zeros(v):
            # type: (Any) -> int
            lead0 = 0
            idx = 0
            while idx < len(v) and v[idx] in '0.':
                if v[idx] == '0':
                    lead0 += 1
                idx += 1
            return lead0

        # underscore = None
        m_sign = False  # type: Any
        value_so = to_str(self.construct_scalar(node))
        value_s = value_so.replace('_', "").lower()
        sign = +1
        if value_s[0] == '-':
            sign = -1
        if value_s[0] in '+-':
            m_sign = value_s[0]
            value_s = value_s[1:]
        if value_s == '.inf':
            return sign * self.inf_value
        if value_s == '.nan':
            return self.nan_value
        if self.resolver.processing_version != (1, 2) and ':' in value_s:
            digits = [float(part) for part in value_s.split(':')]
            digits.reverse()
            base = 1
            value = 0.0
            for digit in digits:
                value += digit * base
                base *= 60
            return sign * value
        if 'e' in value_s:
            try:
                mantissa, exponent = value_so.split('e')
                exp = 'e'
            except ValueError:
                mantissa, exponent = value_so.split('E')
                exp = 'E'
            if self.resolver.processing_version != (1, 2):
                # value_s is lower case independent of input
                if '.' not in mantissa:
                    warnings.warn(MantissaNoDotYAML1_1Warning(node, value_so))
            lead0 = leading_zeros(mantissa)
            width = len(mantissa)
            prec = mantissa.find('.')
            if m_sign:
                width -= 1
            e_width = len(exponent)
            e_sign = exponent[0] in '+-'
            # nprint('sf', width, prec, m_sign, exp, e_width, e_sign)
            return ScalarFloat(
                sign * float(value_s),
                width=width,
                prec=prec,
                m_sign=m_sign,
                m_lead0=lead0,
                exp=exp,
                e_width=e_width,
                e_sign=e_sign,
                anchor=node.anchor,
            )
        width = len(value_so)
        prec = value_so.index('.')  # you can use index, this would not be float without dot
        lead0 = leading_zeros(value_so)
        return ScalarFloat(
            sign * float(value_s),
            width=width,
            prec=prec,
            m_sign=m_sign,
            m_lead0=lead0,
            anchor=node.anchor,
        )

    def construct_yaml_str(self, node):
        # type: (Any) -> Any
        value = self.construct_scalar(node)
        if isinstance(value, ScalarString):
            return value
        if PY3:
            return value
        try:
            return value.encode('ascii')
        except AttributeError:
            # in case you replace the node dynamically e.g. with a dict
            return value
        except UnicodeEncodeError:
            return value

    def construct_rt_sequence(self, node, seqtyp, deep=False):
        # type: (Any, Any, bool) -> Any
        if not isinstance(node, SequenceNode):
            raise ConstructorError(
                None, None, 'expected a sequence node, but found %s' % node.id, node.start_mark
            )
        ret_val = []
        if node.comment:
            seqtyp._yaml_add_comment(node.comment[:2])
            if len(node.comment) > 2:
                seqtyp.yaml_end_comment_extend(node.comment[2], clear=True)
        if node.anchor:
            from ruamel.yaml.serializer import templated_id

            if not templated_id(node.anchor):
                seqtyp.yaml_set_anchor(node.anchor)
        for idx, child in enumerate(node.value):
            ret_val.append(self.construct_object(child, deep=deep))
            if child.comment:
                seqtyp._yaml_add_comment(child.comment, key=idx)
            seqtyp._yaml_set_idx_line_col(
                idx, [child.start_mark.line, child.start_mark.column]
            )
        return ret_val

    def flatten_mapping(self, node):
        # type: (Any) -> Any
        """
        This implements the merge key feature http://yaml.org/type/merge.html
        by inserting keys from the merge dict/list of dicts if not yet
        available in this node
        """

        def constructed(value_node):
            # type: (Any) -> Any
            # If the contents of a merge are defined within the
            # merge marker, then they won't have been constructed
            # yet. But if they were already constructed, we need to use
            # the existing object.
            if value_node in self.constructed_objects:
                value = self.constructed_objects[value_node]
            else:
                value = self.construct_object(value_node, deep=False)
            return value

        # merge = []
        merge_map_list = []  # type: List[Any]
        index = 0
        while index < len(node.value):
            key_node, value_node = node.value[index]
            if key_node.tag == u'tag:yaml.org,2002:merge':
                if merge_map_list:  # double << key
                    if self.allow_duplicate_keys:
                        del node.value[index]
                        index += 1
                        continue
                    args = [
                        'while constructing a mapping',
                        node.start_mark,
                        'found duplicate key "{}"'.format(key_node.value),
                        key_node.start_mark,
                        """
                        To suppress this check see:
                           http://yaml.readthedocs.io/en/latest/api.html#duplicate-keys
                        """,
                        """\
                        Duplicate keys will become an error in future releases, and are errors
                        by default when using the new API.
                        """,
                    ]
                    if self.allow_duplicate_keys is None:
                        warnings.warn(DuplicateKeyFutureWarning(*args))
                    else:
                        raise DuplicateKeyError(*args)
                del node.value[index]
                if isinstance(value_node, MappingNode):
                    merge_map_list.append((index, constructed(value_node)))
                    # self.flatten_mapping(value_node)
                    # merge.extend(value_node.value)
                elif isinstance(value_node, SequenceNode):
                    # submerge = []
                    for subnode in value_node.value:
                        if not isinstance(subnode, MappingNode):
                            raise ConstructorError(
                                'while constructing a mapping',
                                node.start_mark,
                                'expected a mapping for merging, but found %s' % subnode.id,
                                subnode.start_mark,
                            )
                        merge_map_list.append((index, constructed(subnode)))
                    #     self.flatten_mapping(subnode)
                    #     submerge.append(subnode.value)
                    # submerge.reverse()
                    # for value in submerge:
                    #     merge.extend(value)
                else:
                    raise ConstructorError(
                        'while constructing a mapping',
                        node.start_mark,
                        'expected a mapping or list of mappings for merging, '
                        'but found %s' % value_node.id,
                        value_node.start_mark,
                    )
            elif key_node.tag == u'tag:yaml.org,2002:value':
                key_node.tag = u'tag:yaml.org,2002:str'
                index += 1
            else:
                index += 1
        return merge_map_list
        # if merge:
        #     node.value = merge + node.value

    def _sentinel(self):
        # type: () -> None
        pass

    def construct_mapping(self, node, maptyp, deep=False):  # type: ignore
        # type: (Any, Any, bool) -> Any
        if not isinstance(node, MappingNode):
            raise ConstructorError(
                None, None, 'expected a mapping node, but found %s' % node.id, node.start_mark
            )
        merge_map = self.flatten_mapping(node)
        # mapping = {}
        if node.comment:
            maptyp._yaml_add_comment(node.comment[:2])
            if len(node.comment) > 2:
                maptyp.yaml_end_comment_extend(node.comment[2], clear=True)
        if node.anchor:
            from ruamel.yaml.serializer import templated_id

            if not templated_id(node.anchor):
                maptyp.yaml_set_anchor(node.anchor)
        last_key, last_value = None, self._sentinel
        for key_node, value_node in node.value:
            # keys can be list -> deep
            key = self.construct_object(key_node, deep=True)
            # lists are not hashable, but tuples are
            if not isinstance(key, Hashable):
                if isinstance(key, MutableSequence):
                    key_s = CommentedKeySeq(key)
                    if key_node.flow_style is True:
                        key_s.fa.set_flow_style()
                    elif key_node.flow_style is False:
                        key_s.fa.set_block_style()
                    key = key_s
                elif isinstance(key, MutableMapping):
                    key_m = CommentedKeyMap(key)
                    if key_node.flow_style is True:
                        key_m.fa.set_flow_style()
                    elif key_node.flow_style is False:
                        key_m.fa.set_block_style()
                    key = key_m
            if PY2:
                try:
                    hash(key)
                except TypeError as exc:
                    raise ConstructorError(
                        'while constructing a mapping',
                        node.start_mark,
                        'found unacceptable key (%s)' % exc,
                        key_node.start_mark,
                    )
            else:
                if not isinstance(key, Hashable):
                    raise ConstructorError(
                        'while constructing a mapping',
                        node.start_mark,
                        'found unhashable key',
                        key_node.start_mark,
                    )
            value = self.construct_object(value_node, deep=deep)
            if self.check_mapping_key(node, key_node, maptyp, key, value):

                if key_node.comment and len(key_node.comment) > 4 and key_node.comment[4]:
                    if last_value is None:
                        key_node.comment[0] = key_node.comment.pop(4)
                        maptyp._yaml_add_comment(key_node.comment, value=last_key)
                    else:
                        key_node.comment[2] = key_node.comment.pop(4)
                        maptyp._yaml_add_comment(key_node.comment, key=key)
                    key_node.comment = None
                if key_node.comment:
                    maptyp._yaml_add_comment(key_node.comment, key=key)
                if value_node.comment:
                    maptyp._yaml_add_comment(value_node.comment, value=key)
                maptyp._yaml_set_kv_line_col(
                    key,
                    [
                        key_node.start_mark.line,
                        key_node.start_mark.column,
                        value_node.start_mark.line,
                        value_node.start_mark.column,
                    ],
                )
                maptyp[key] = value
                last_key, last_value = key, value  # could use indexing
        # do this last, or <<: before a key will prevent insertion in instances
        # of collections.OrderedDict (as they have no __contains__
        if merge_map:
            maptyp.add_yaml_merge(merge_map)

    def construct_setting(self, node, typ, deep=False):
        # type: (Any, Any, bool) -> Any
        if not isinstance(node, MappingNode):
            raise ConstructorError(
                None, None, 'expected a mapping node, but found %s' % node.id, node.start_mark
            )
        if node.comment:
            typ._yaml_add_comment(node.comment[:2])
            if len(node.comment) > 2:
                typ.yaml_end_comment_extend(node.comment[2], clear=True)
        if node.anchor:
            from ruamel.yaml.serializer import templated_id

            if not templated_id(node.anchor):
                typ.yaml_set_anchor(node.anchor)
        for key_node, value_node in node.value:
            # keys can be list -> deep
            key = self.construct_object(key_node, deep=True)
            # lists are not hashable, but tuples are
            if not isinstance(key, Hashable):
                if isinstance(key, list):
                    key = tuple(key)
            if PY2:
                try:
                    hash(key)
                except TypeError as exc:
                    raise ConstructorError(
                        'while constructing a mapping',
                        node.start_mark,
                        'found unacceptable key (%s)' % exc,
                        key_node.start_mark,
                    )
            else:
                if not isinstance(key, Hashable):
                    raise ConstructorError(
                        'while constructing a mapping',
                        node.start_mark,
                        'found unhashable key',
                        key_node.start_mark,
                    )
            # construct but should be null
            value = self.construct_object(value_node, deep=deep)  # NOQA
            self.check_set_key(node, key_node, typ, key)
            if key_node.comment:
                typ._yaml_add_comment(key_node.comment, key=key)
            if value_node.comment:
                typ._yaml_add_comment(value_node.comment, value=key)
            typ.add(key)

    def construct_yaml_seq(self, node):
        # type: (Any) -> Any
        data = CommentedSeq()
        data._yaml_set_line_col(node.start_mark.line, node.start_mark.column)
        if node.comment:
            data._yaml_add_comment(node.comment)
        yield data
        data.extend(self.construct_rt_sequence(node, data))
        self.set_collection_style(data, node)

    def construct_yaml_map(self, node):
        # type: (Any) -> Any
        data = CommentedMap()
        data._yaml_set_line_col(node.start_mark.line, node.start_mark.column)
        yield data
        self.construct_mapping(node, data, deep=True)
        self.set_collection_style(data, node)

    def set_collection_style(self, data, node):
        # type: (Any, Any) -> None
        if len(data) == 0:
            return
        if node.flow_style is True:
            data.fa.set_flow_style()
        elif node.flow_style is False:
            data.fa.set_block_style()

    def construct_yaml_object(self, node, cls):
        # type: (Any, Any) -> Any
        data = cls.__new__(cls)
        yield data
        if hasattr(data, '__setstate__'):
            state = SafeConstructor.construct_mapping(self, node, deep=True)
            data.__setstate__(state)
        else:
            state = SafeConstructor.construct_mapping(self, node)
            data.__dict__.update(state)

    def construct_yaml_omap(self, node):
        # type: (Any) -> Any
        # Note: we do now check for duplicate keys
        omap = CommentedOrderedMap()
        omap._yaml_set_line_col(node.start_mark.line, node.start_mark.column)
        if node.flow_style is True:
            omap.fa.set_flow_style()
        elif node.flow_style is False:
            omap.fa.set_block_style()
        yield omap
        if node.comment:
            omap._yaml_add_comment(node.comment[:2])
            if len(node.comment) > 2:
                omap.yaml_end_comment_extend(node.comment[2], clear=True)
        if not isinstance(node, SequenceNode):
            raise ConstructorError(
                'while constructing an ordered map',
                node.start_mark,
                'expected a sequence, but found %s' % node.id,
                node.start_mark,
            )
        for subnode in node.value:
            if not isinstance(subnode, MappingNode):
                raise ConstructorError(
                    'while constructing an ordered map',
                    node.start_mark,
                    'expected a mapping of length 1, but found %s' % subnode.id,
                    subnode.start_mark,
                )
            if len(subnode.value) != 1:
                raise ConstructorError(
                    'while constructing an ordered map',
                    node.start_mark,
                    'expected a single mapping item, but found %d items' % len(subnode.value),
                    subnode.start_mark,
                )
            key_node, value_node = subnode.value[0]
            key = self.construct_object(key_node)
            assert key not in omap
            value = self.construct_object(value_node)
            if key_node.comment:
                omap._yaml_add_comment(key_node.comment, key=key)
            if subnode.comment:
                omap._yaml_add_comment(subnode.comment, key=key)
            if value_node.comment:
                omap._yaml_add_comment(value_node.comment, value=key)
            omap[key] = value

    def construct_yaml_set(self, node):
        # type: (Any) -> Any
        data = CommentedSet()
        data._yaml_set_line_col(node.start_mark.line, node.start_mark.column)
        yield data
        self.construct_setting(node, data)

    def construct_undefined(self, node):
        # type: (Any) -> Any
        try:
            if isinstance(node, MappingNode):
                data = CommentedMap()
                data._yaml_set_line_col(node.start_mark.line, node.start_mark.column)
                if node.flow_style is True:
                    data.fa.set_flow_style()
                elif node.flow_style is False:
                    data.fa.set_block_style()
                data.yaml_set_tag(node.tag)
                yield data
                if node.anchor:
                    data.yaml_set_anchor(node.anchor)
                self.construct_mapping(node, data)
                return
            elif isinstance(node, ScalarNode):
                data2 = TaggedScalar()
                data2.value = self.construct_scalar(node)
                data2.style = node.style
                data2.yaml_set_tag(node.tag)
                yield data2
                if node.anchor:
                    data2.yaml_set_anchor(node.anchor, always_dump=True)
                return
            elif isinstance(node, SequenceNode):
                data3 = CommentedSeq()
                data3._yaml_set_line_col(node.start_mark.line, node.start_mark.column)
                if node.flow_style is True:
                    data3.fa.set_flow_style()
                elif node.flow_style is False:
                    data3.fa.set_block_style()
                data3.yaml_set_tag(node.tag)
                yield data3
                if node.anchor:
                    data3.yaml_set_anchor(node.anchor)
                data3.extend(self.construct_sequence(node))
                return
        except:  # NOQA
            pass
        raise ConstructorError(
            None,
            None,
            'could not determine a constructor for the tag %r' % utf8(node.tag),
            node.start_mark,
        )

    def construct_yaml_timestamp(self, node, values=None):
        # type: (Any, Any) -> Any
        try:
            match = self.timestamp_regexp.match(node.value)
        except TypeError:
            match = None
        if match is None:
            raise ConstructorError(
                None,
                None,
                'failed to construct timestamp from "{}"'.format(node.value),
                node.start_mark,
            )
        values = match.groupdict()
        if not values['hour']:
            return SafeConstructor.construct_yaml_timestamp(self, node, values)
        for part in ['t', 'tz_sign', 'tz_hour', 'tz_minute']:
            if values[part]:
                break
        else:
            return SafeConstructor.construct_yaml_timestamp(self, node, values)
        year = int(values['year'])
        month = int(values['month'])
        day = int(values['day'])
        hour = int(values['hour'])
        minute = int(values['minute'])
        second = int(values['second'])
        fraction = 0
        if values['fraction']:
            fraction_s = values['fraction'][:6]
            while len(fraction_s) < 6:
                fraction_s += '0'
            fraction = int(fraction_s)
            if len(values['fraction']) > 6 and int(values['fraction'][6]) > 4:
                fraction += 1
        delta = None
        if values['tz_sign']:
            tz_hour = int(values['tz_hour'])
            minutes = values['tz_minute']
            tz_minute = int(minutes) if minutes else 0
            delta = datetime.timedelta(hours=tz_hour, minutes=tz_minute)
            if values['tz_sign'] == '-':
                delta = -delta
        if delta:
            dt = datetime.datetime(year, month, day, hour, minute)
            dt -= delta
            data = TimeStamp(dt.year, dt.month, dt.day, dt.hour, dt.minute, second, fraction)
            data._yaml['delta'] = delta
            tz = values['tz_sign'] + values['tz_hour']
            if values['tz_minute']:
                tz += ':' + values['tz_minute']
            data._yaml['tz'] = tz
        else:
            data = TimeStamp(year, month, day, hour, minute, second, fraction)
            if values['tz']:  # no delta
                data._yaml['tz'] = values['tz']

        if values['t']:
            data._yaml['t'] = True
        return data

    def construct_yaml_bool(self, node):
        # type: (Any) -> Any
        b = SafeConstructor.construct_yaml_bool(self, node)
        if node.anchor:
            return ScalarBoolean(b, anchor=node.anchor)
        return b


RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:null', RoundTripConstructor.construct_yaml_null
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:bool', RoundTripConstructor.construct_yaml_bool
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:int', RoundTripConstructor.construct_yaml_int
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:float', RoundTripConstructor.construct_yaml_float
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:binary', RoundTripConstructor.construct_yaml_binary
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:timestamp', RoundTripConstructor.construct_yaml_timestamp
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:omap', RoundTripConstructor.construct_yaml_omap
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:pairs', RoundTripConstructor.construct_yaml_pairs
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:set', RoundTripConstructor.construct_yaml_set
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:str', RoundTripConstructor.construct_yaml_str
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:seq', RoundTripConstructor.construct_yaml_seq
)

RoundTripConstructor.add_constructor(
    u'tag:yaml.org,2002:map', RoundTripConstructor.construct_yaml_map
)

RoundTripConstructor.add_constructor(None, RoundTripConstructor.construct_undefined)
