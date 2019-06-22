# coding: utf-8

from __future__ import print_function, absolute_import, division


from ruamel.yaml.error import *  # NOQA
from ruamel.yaml.nodes import *  # NOQA
from ruamel.yaml.compat import text_type, binary_type, to_unicode, PY2, PY3, ordereddict
from ruamel.yaml.compat import nprint, nprintf  # NOQA
from ruamel.yaml.scalarstring import (
    LiteralScalarString,
    FoldedScalarString,
    SingleQuotedScalarString,
    DoubleQuotedScalarString,
    PlainScalarString,
)
from ruamel.yaml.scalarint import ScalarInt, BinaryInt, OctalInt, HexInt, HexCapsInt
from ruamel.yaml.scalarfloat import ScalarFloat
from ruamel.yaml.scalarbool import ScalarBoolean
from ruamel.yaml.timestamp import TimeStamp

import datetime
import sys
import types

if PY3:
    import copyreg
    import base64
else:
    import copy_reg as copyreg  # type: ignore

if False:  # MYPY
    from typing import Dict, List, Any, Union, Text, Optional  # NOQA

# fmt: off
__all__ = ['BaseRepresenter', 'SafeRepresenter', 'Representer',
           'RepresenterError', 'RoundTripRepresenter']
# fmt: on


class RepresenterError(YAMLError):
    pass


if PY2:

    def get_classobj_bases(cls):
        # type: (Any) -> Any
        bases = [cls]
        for base in cls.__bases__:
            bases.extend(get_classobj_bases(base))
        return bases


class BaseRepresenter(object):

    yaml_representers = {}  # type: Dict[Any, Any]
    yaml_multi_representers = {}  # type: Dict[Any, Any]

    def __init__(self, default_style=None, default_flow_style=None, dumper=None):
        # type: (Any, Any, Any, Any) -> None
        self.dumper = dumper
        if self.dumper is not None:
            self.dumper._representer = self
        self.default_style = default_style
        self.default_flow_style = default_flow_style
        self.represented_objects = {}  # type: Dict[Any, Any]
        self.object_keeper = []  # type: List[Any]
        self.alias_key = None  # type: Optional[int]
        self.sort_base_mapping_type_on_output = True

    @property
    def serializer(self):
        # type: () -> Any
        try:
            if hasattr(self.dumper, 'typ'):
                return self.dumper.serializer
            return self.dumper._serializer
        except AttributeError:
            return self  # cyaml

    def represent(self, data):
        # type: (Any) -> None
        node = self.represent_data(data)
        self.serializer.serialize(node)
        self.represented_objects = {}  # type: Dict[Any, Any]
        self.object_keeper = []  # type: List[Any]
        self.alias_key = None

    def represent_data(self, data):
        # type: (Any) -> Any
        if self.ignore_aliases(data):
            self.alias_key = None
        else:
            self.alias_key = id(data)
        if self.alias_key is not None:
            if self.alias_key in self.represented_objects:
                node = self.represented_objects[self.alias_key]
                # if node is None:
                #     raise RepresenterError(
                #          "recursive objects are not allowed: %r" % data)
                return node
            # self.represented_objects[alias_key] = None
            self.object_keeper.append(data)
        data_types = type(data).__mro__
        if PY2:
            # if type(data) is types.InstanceType:
            if isinstance(data, types.InstanceType):
                data_types = get_classobj_bases(data.__class__) + list(data_types)
        if data_types[0] in self.yaml_representers:
            node = self.yaml_representers[data_types[0]](self, data)
        else:
            for data_type in data_types:
                if data_type in self.yaml_multi_representers:
                    node = self.yaml_multi_representers[data_type](self, data)
                    break
            else:
                if None in self.yaml_multi_representers:
                    node = self.yaml_multi_representers[None](self, data)
                elif None in self.yaml_representers:
                    node = self.yaml_representers[None](self, data)
                else:
                    node = ScalarNode(None, text_type(data))
        # if alias_key is not None:
        #     self.represented_objects[alias_key] = node
        return node

    def represent_key(self, data):
        # type: (Any) -> Any
        """
        David Fraser: Extract a method to represent keys in mappings, so that
        a subclass can choose not to quote them (for example)
        used in represent_mapping
        https://bitbucket.org/davidfraser/pyyaml/commits/d81df6eb95f20cac4a79eed95ae553b5c6f77b8c
        """
        return self.represent_data(data)

    @classmethod
    def add_representer(cls, data_type, representer):
        # type: (Any, Any) -> None
        if 'yaml_representers' not in cls.__dict__:
            cls.yaml_representers = cls.yaml_representers.copy()
        cls.yaml_representers[data_type] = representer

    @classmethod
    def add_multi_representer(cls, data_type, representer):
        # type: (Any, Any) -> None
        if 'yaml_multi_representers' not in cls.__dict__:
            cls.yaml_multi_representers = cls.yaml_multi_representers.copy()
        cls.yaml_multi_representers[data_type] = representer

    def represent_scalar(self, tag, value, style=None, anchor=None):
        # type: (Any, Any, Any, Any) -> Any
        if style is None:
            style = self.default_style
        comment = None
        if style and style[0] in '|>':
            comment = getattr(value, 'comment', None)
            if comment:
                comment = [None, [comment]]
        node = ScalarNode(tag, value, style=style, comment=comment, anchor=anchor)
        if self.alias_key is not None:
            self.represented_objects[self.alias_key] = node
        return node

    def represent_sequence(self, tag, sequence, flow_style=None):
        # type: (Any, Any, Any) -> Any
        value = []  # type: List[Any]
        node = SequenceNode(tag, value, flow_style=flow_style)
        if self.alias_key is not None:
            self.represented_objects[self.alias_key] = node
        best_style = True
        for item in sequence:
            node_item = self.represent_data(item)
            if not (isinstance(node_item, ScalarNode) and not node_item.style):
                best_style = False
            value.append(node_item)
        if flow_style is None:
            if self.default_flow_style is not None:
                node.flow_style = self.default_flow_style
            else:
                node.flow_style = best_style
        return node

    def represent_omap(self, tag, omap, flow_style=None):
        # type: (Any, Any, Any) -> Any
        value = []  # type: List[Any]
        node = SequenceNode(tag, value, flow_style=flow_style)
        if self.alias_key is not None:
            self.represented_objects[self.alias_key] = node
        best_style = True
        for item_key in omap:
            item_val = omap[item_key]
            node_item = self.represent_data({item_key: item_val})
            # if not (isinstance(node_item, ScalarNode) \
            #    and not node_item.style):
            #     best_style = False
            value.append(node_item)
        if flow_style is None:
            if self.default_flow_style is not None:
                node.flow_style = self.default_flow_style
            else:
                node.flow_style = best_style
        return node

    def represent_mapping(self, tag, mapping, flow_style=None):
        # type: (Any, Any, Any) -> Any
        value = []  # type: List[Any]
        node = MappingNode(tag, value, flow_style=flow_style)
        if self.alias_key is not None:
            self.represented_objects[self.alias_key] = node
        best_style = True
        if hasattr(mapping, 'items'):
            mapping = list(mapping.items())
            if self.sort_base_mapping_type_on_output:
                try:
                    mapping = sorted(mapping)
                except TypeError:
                    pass
        for item_key, item_value in mapping:
            node_key = self.represent_key(item_key)
            node_value = self.represent_data(item_value)
            if not (isinstance(node_key, ScalarNode) and not node_key.style):
                best_style = False
            if not (isinstance(node_value, ScalarNode) and not node_value.style):
                best_style = False
            value.append((node_key, node_value))
        if flow_style is None:
            if self.default_flow_style is not None:
                node.flow_style = self.default_flow_style
            else:
                node.flow_style = best_style
        return node

    def ignore_aliases(self, data):
        # type: (Any) -> bool
        return False


class SafeRepresenter(BaseRepresenter):
    def ignore_aliases(self, data):
        # type: (Any) -> bool
        # https://docs.python.org/3/reference/expressions.html#parenthesized-forms :
        # "i.e. two occurrences of the empty tuple may or may not yield the same object"
        # so "data is ()" should not be used
        if data is None or (isinstance(data, tuple) and data == ()):
            return True
        if isinstance(data, (binary_type, text_type, bool, int, float)):
            return True
        return False

    def represent_none(self, data):
        # type: (Any) -> Any
        return self.represent_scalar(u'tag:yaml.org,2002:null', u'null')

    if PY3:

        def represent_str(self, data):
            # type: (Any) -> Any
            return self.represent_scalar(u'tag:yaml.org,2002:str', data)

        def represent_binary(self, data):
            # type: (Any) -> Any
            if hasattr(base64, 'encodebytes'):
                data = base64.encodebytes(data).decode('ascii')
            else:
                data = base64.encodestring(data).decode('ascii')
            return self.represent_scalar(u'tag:yaml.org,2002:binary', data, style='|')

    else:

        def represent_str(self, data):
            # type: (Any) -> Any
            tag = None
            style = None
            try:
                data = unicode(data, 'ascii')
                tag = u'tag:yaml.org,2002:str'
            except UnicodeDecodeError:
                try:
                    data = unicode(data, 'utf-8')
                    tag = u'tag:yaml.org,2002:str'
                except UnicodeDecodeError:
                    data = data.encode('base64')
                    tag = u'tag:yaml.org,2002:binary'
                    style = '|'
            return self.represent_scalar(tag, data, style=style)

        def represent_unicode(self, data):
            # type: (Any) -> Any
            return self.represent_scalar(u'tag:yaml.org,2002:str', data)

    def represent_bool(self, data, anchor=None):
        # type: (Any, Optional[Any]) -> Any
        try:
            value = self.dumper.boolean_representation[bool(data)]
        except AttributeError:
            if data:
                value = u'true'
            else:
                value = u'false'
        return self.represent_scalar(u'tag:yaml.org,2002:bool', value, anchor=anchor)

    def represent_int(self, data):
        # type: (Any) -> Any
        return self.represent_scalar(u'tag:yaml.org,2002:int', text_type(data))

    if PY2:

        def represent_long(self, data):
            # type: (Any) -> Any
            return self.represent_scalar(u'tag:yaml.org,2002:int', text_type(data))

    inf_value = 1e300
    while repr(inf_value) != repr(inf_value * inf_value):
        inf_value *= inf_value

    def represent_float(self, data):
        # type: (Any) -> Any
        if data != data or (data == 0.0 and data == 1.0):
            value = u'.nan'
        elif data == self.inf_value:
            value = u'.inf'
        elif data == -self.inf_value:
            value = u'-.inf'
        else:
            value = to_unicode(repr(data)).lower()
            if getattr(self.serializer, 'use_version', None) == (1, 1):
                if u'.' not in value and u'e' in value:
                    # Note that in some cases `repr(data)` represents a float number
                    # without the decimal parts.  For instance:
                    #   >>> repr(1e17)
                    #   '1e17'
                    # Unfortunately, this is not a valid float representation according
                    # to the definition of the `!!float` tag in YAML 1.1.  We fix
                    # this by adding '.0' before the 'e' symbol.
                    value = value.replace(u'e', u'.0e', 1)
        return self.represent_scalar(u'tag:yaml.org,2002:float', value)

    def represent_list(self, data):
        # type: (Any) -> Any
        # pairs = (len(data) > 0 and isinstance(data, list))
        # if pairs:
        #     for item in data:
        #         if not isinstance(item, tuple) or len(item) != 2:
        #             pairs = False
        #             break
        # if not pairs:
        return self.represent_sequence(u'tag:yaml.org,2002:seq', data)

    # value = []
    # for item_key, item_value in data:
    #     value.append(self.represent_mapping(u'tag:yaml.org,2002:map',
    #         [(item_key, item_value)]))
    # return SequenceNode(u'tag:yaml.org,2002:pairs', value)

    def represent_dict(self, data):
        # type: (Any) -> Any
        return self.represent_mapping(u'tag:yaml.org,2002:map', data)

    def represent_ordereddict(self, data):
        # type: (Any) -> Any
        return self.represent_omap(u'tag:yaml.org,2002:omap', data)

    def represent_set(self, data):
        # type: (Any) -> Any
        value = {}  # type: Dict[Any, None]
        for key in data:
            value[key] = None
        return self.represent_mapping(u'tag:yaml.org,2002:set', value)

    def represent_date(self, data):
        # type: (Any) -> Any
        value = to_unicode(data.isoformat())
        return self.represent_scalar(u'tag:yaml.org,2002:timestamp', value)

    def represent_datetime(self, data):
        # type: (Any) -> Any
        value = to_unicode(data.isoformat(' '))
        return self.represent_scalar(u'tag:yaml.org,2002:timestamp', value)

    def represent_yaml_object(self, tag, data, cls, flow_style=None):
        # type: (Any, Any, Any, Any) -> Any
        if hasattr(data, '__getstate__'):
            state = data.__getstate__()
        else:
            state = data.__dict__.copy()
        return self.represent_mapping(tag, state, flow_style=flow_style)

    def represent_undefined(self, data):
        # type: (Any) -> None
        raise RepresenterError('cannot represent an object: %s' % (data,))


SafeRepresenter.add_representer(type(None), SafeRepresenter.represent_none)

SafeRepresenter.add_representer(str, SafeRepresenter.represent_str)

if PY2:
    SafeRepresenter.add_representer(unicode, SafeRepresenter.represent_unicode)
else:
    SafeRepresenter.add_representer(bytes, SafeRepresenter.represent_binary)

SafeRepresenter.add_representer(bool, SafeRepresenter.represent_bool)

SafeRepresenter.add_representer(int, SafeRepresenter.represent_int)

if PY2:
    SafeRepresenter.add_representer(long, SafeRepresenter.represent_long)

SafeRepresenter.add_representer(float, SafeRepresenter.represent_float)

SafeRepresenter.add_representer(list, SafeRepresenter.represent_list)

SafeRepresenter.add_representer(tuple, SafeRepresenter.represent_list)

SafeRepresenter.add_representer(dict, SafeRepresenter.represent_dict)

SafeRepresenter.add_representer(set, SafeRepresenter.represent_set)

SafeRepresenter.add_representer(ordereddict, SafeRepresenter.represent_ordereddict)

if sys.version_info >= (2, 7):
    import collections

    SafeRepresenter.add_representer(
        collections.OrderedDict, SafeRepresenter.represent_ordereddict
    )

SafeRepresenter.add_representer(datetime.date, SafeRepresenter.represent_date)

SafeRepresenter.add_representer(datetime.datetime, SafeRepresenter.represent_datetime)

SafeRepresenter.add_representer(None, SafeRepresenter.represent_undefined)


class Representer(SafeRepresenter):
    if PY2:

        def represent_str(self, data):
            # type: (Any) -> Any
            tag = None
            style = None
            try:
                data = unicode(data, 'ascii')
                tag = u'tag:yaml.org,2002:str'
            except UnicodeDecodeError:
                try:
                    data = unicode(data, 'utf-8')
                    tag = u'tag:yaml.org,2002:python/str'
                except UnicodeDecodeError:
                    data = data.encode('base64')
                    tag = u'tag:yaml.org,2002:binary'
                    style = '|'
            return self.represent_scalar(tag, data, style=style)

        def represent_unicode(self, data):
            # type: (Any) -> Any
            tag = None
            try:
                data.encode('ascii')
                tag = u'tag:yaml.org,2002:python/unicode'
            except UnicodeEncodeError:
                tag = u'tag:yaml.org,2002:str'
            return self.represent_scalar(tag, data)

        def represent_long(self, data):
            # type: (Any) -> Any
            tag = u'tag:yaml.org,2002:int'
            if int(data) is not data:
                tag = u'tag:yaml.org,2002:python/long'
            return self.represent_scalar(tag, to_unicode(data))

    def represent_complex(self, data):
        # type: (Any) -> Any
        if data.imag == 0.0:
            data = u'%r' % data.real
        elif data.real == 0.0:
            data = u'%rj' % data.imag
        elif data.imag > 0:
            data = u'%r+%rj' % (data.real, data.imag)
        else:
            data = u'%r%rj' % (data.real, data.imag)
        return self.represent_scalar(u'tag:yaml.org,2002:python/complex', data)

    def represent_tuple(self, data):
        # type: (Any) -> Any
        return self.represent_sequence(u'tag:yaml.org,2002:python/tuple', data)

    def represent_name(self, data):
        # type: (Any) -> Any
        try:
            name = u'%s.%s' % (data.__module__, data.__qualname__)
        except AttributeError:
            # probably PY2
            name = u'%s.%s' % (data.__module__, data.__name__)
        return self.represent_scalar(u'tag:yaml.org,2002:python/name:' + name, "")

    def represent_module(self, data):
        # type: (Any) -> Any
        return self.represent_scalar(u'tag:yaml.org,2002:python/module:' + data.__name__, "")

    if PY2:

        def represent_instance(self, data):
            # type: (Any) -> Any
            # For instances of classic classes, we use __getinitargs__ and
            # __getstate__ to serialize the data.

            # If data.__getinitargs__ exists, the object must be reconstructed
            # by calling cls(**args), where args is a tuple returned by
            # __getinitargs__. Otherwise, the cls.__init__ method should never
            # be called and the class instance is created by instantiating a
            # trivial class and assigning to the instance's __class__ variable.

            # If data.__getstate__ exists, it returns the state of the object.
            # Otherwise, the state of the object is data.__dict__.

            # We produce either a !!python/object or !!python/object/new node.
            # If data.__getinitargs__ does not exist and state is a dictionary,
            # we produce a !!python/object node . Otherwise we produce a
            # !!python/object/new node.

            cls = data.__class__
            class_name = u'%s.%s' % (cls.__module__, cls.__name__)
            args = None
            state = None
            if hasattr(data, '__getinitargs__'):
                args = list(data.__getinitargs__())
            if hasattr(data, '__getstate__'):
                state = data.__getstate__()
            else:
                state = data.__dict__
            if args is None and isinstance(state, dict):
                return self.represent_mapping(
                    u'tag:yaml.org,2002:python/object:' + class_name, state
                )
            if isinstance(state, dict) and not state:
                return self.represent_sequence(
                    u'tag:yaml.org,2002:python/object/new:' + class_name, args
                )
            value = {}
            if bool(args):
                value['args'] = args
            value['state'] = state  # type: ignore
            return self.represent_mapping(
                u'tag:yaml.org,2002:python/object/new:' + class_name, value
            )

    def represent_object(self, data):
        # type: (Any) -> Any
        # We use __reduce__ API to save the data. data.__reduce__ returns
        # a tuple of length 2-5:
        #   (function, args, state, listitems, dictitems)

        # For reconstructing, we calls function(*args), then set its state,
        # listitems, and dictitems if they are not None.

        # A special case is when function.__name__ == '__newobj__'. In this
        # case we create the object with args[0].__new__(*args).

        # Another special case is when __reduce__ returns a string - we don't
        # support it.

        # We produce a !!python/object, !!python/object/new or
        # !!python/object/apply node.

        cls = type(data)
        if cls in copyreg.dispatch_table:
            reduce = copyreg.dispatch_table[cls](data)
        elif hasattr(data, '__reduce_ex__'):
            reduce = data.__reduce_ex__(2)
        elif hasattr(data, '__reduce__'):
            reduce = data.__reduce__()
        else:
            raise RepresenterError('cannot represent object: %r' % (data,))
        reduce = (list(reduce) + [None] * 5)[:5]
        function, args, state, listitems, dictitems = reduce
        args = list(args)
        if state is None:
            state = {}
        if listitems is not None:
            listitems = list(listitems)
        if dictitems is not None:
            dictitems = dict(dictitems)
        if function.__name__ == '__newobj__':
            function = args[0]
            args = args[1:]
            tag = u'tag:yaml.org,2002:python/object/new:'
            newobj = True
        else:
            tag = u'tag:yaml.org,2002:python/object/apply:'
            newobj = False
        try:
            function_name = u'%s.%s' % (function.__module__, function.__qualname__)
        except AttributeError:
            # probably PY2
            function_name = u'%s.%s' % (function.__module__, function.__name__)
        if not args and not listitems and not dictitems and isinstance(state, dict) and newobj:
            return self.represent_mapping(
                u'tag:yaml.org,2002:python/object:' + function_name, state
            )
        if not listitems and not dictitems and isinstance(state, dict) and not state:
            return self.represent_sequence(tag + function_name, args)
        value = {}
        if args:
            value['args'] = args
        if state or not isinstance(state, dict):
            value['state'] = state
        if listitems:
            value['listitems'] = listitems
        if dictitems:
            value['dictitems'] = dictitems
        return self.represent_mapping(tag + function_name, value)


if PY2:
    Representer.add_representer(str, Representer.represent_str)

    Representer.add_representer(unicode, Representer.represent_unicode)

    Representer.add_representer(long, Representer.represent_long)

Representer.add_representer(complex, Representer.represent_complex)

Representer.add_representer(tuple, Representer.represent_tuple)

Representer.add_representer(type, Representer.represent_name)

if PY2:
    Representer.add_representer(types.ClassType, Representer.represent_name)

Representer.add_representer(types.FunctionType, Representer.represent_name)

Representer.add_representer(types.BuiltinFunctionType, Representer.represent_name)

Representer.add_representer(types.ModuleType, Representer.represent_module)

if PY2:
    Representer.add_multi_representer(types.InstanceType, Representer.represent_instance)

Representer.add_multi_representer(object, Representer.represent_object)

Representer.add_multi_representer(type, Representer.represent_name)

from ruamel.yaml.comments import (
    CommentedMap,
    CommentedOrderedMap,
    CommentedSeq,
    CommentedKeySeq,
    CommentedKeyMap,
    CommentedSet,
    comment_attrib,
    merge_attrib,
    TaggedScalar,
)  # NOQA


class RoundTripRepresenter(SafeRepresenter):
    # need to add type here and write out the .comment
    # in serializer and emitter

    def __init__(self, default_style=None, default_flow_style=None, dumper=None):
        # type: (Any, Any, Any) -> None
        if not hasattr(dumper, 'typ') and default_flow_style is None:
            default_flow_style = False
        SafeRepresenter.__init__(
            self,
            default_style=default_style,
            default_flow_style=default_flow_style,
            dumper=dumper,
        )

    def ignore_aliases(self, data):
        # type: (Any) -> bool
        try:
            if data.anchor is not None and data.anchor.value is not None:
                return False
        except AttributeError:
            pass
        return SafeRepresenter.ignore_aliases(self, data)

    def represent_none(self, data):
        # type: (Any) -> Any
        if len(self.represented_objects) == 0 and not self.serializer.use_explicit_start:
            # this will be open ended (although it is not yet)
            return self.represent_scalar(u'tag:yaml.org,2002:null', u'null')
        return self.represent_scalar(u'tag:yaml.org,2002:null', "")

    def represent_literal_scalarstring(self, data):
        # type: (Any) -> Any
        tag = None
        style = '|'
        anchor = data.yaml_anchor(any=True)
        if PY2 and not isinstance(data, unicode):
            data = unicode(data, 'ascii')
        tag = u'tag:yaml.org,2002:str'
        return self.represent_scalar(tag, data, style=style, anchor=anchor)

    represent_preserved_scalarstring = represent_literal_scalarstring

    def represent_folded_scalarstring(self, data):
        # type: (Any) -> Any
        tag = None
        style = '>'
        anchor = data.yaml_anchor(any=True)
        for fold_pos in reversed(getattr(data, 'fold_pos', [])):
            if (
                data[fold_pos] == ' '
                and (fold_pos > 0 and not data[fold_pos - 1].isspace())
                and (fold_pos < len(data) and not data[fold_pos + 1].isspace())
            ):
                data = data[:fold_pos] + '\a' + data[fold_pos:]
        if PY2 and not isinstance(data, unicode):
            data = unicode(data, 'ascii')
        tag = u'tag:yaml.org,2002:str'
        return self.represent_scalar(tag, data, style=style, anchor=anchor)

    def represent_single_quoted_scalarstring(self, data):
        # type: (Any) -> Any
        tag = None
        style = "'"
        anchor = data.yaml_anchor(any=True)
        if PY2 and not isinstance(data, unicode):
            data = unicode(data, 'ascii')
        tag = u'tag:yaml.org,2002:str'
        return self.represent_scalar(tag, data, style=style, anchor=anchor)

    def represent_double_quoted_scalarstring(self, data):
        # type: (Any) -> Any
        tag = None
        style = '"'
        anchor = data.yaml_anchor(any=True)
        if PY2 and not isinstance(data, unicode):
            data = unicode(data, 'ascii')
        tag = u'tag:yaml.org,2002:str'
        return self.represent_scalar(tag, data, style=style, anchor=anchor)

    def represent_plain_scalarstring(self, data):
        # type: (Any) -> Any
        tag = None
        style = ''
        anchor = data.yaml_anchor(any=True)
        if PY2 and not isinstance(data, unicode):
            data = unicode(data, 'ascii')
        tag = u'tag:yaml.org,2002:str'
        return self.represent_scalar(tag, data, style=style, anchor=anchor)

    def insert_underscore(self, prefix, s, underscore, anchor=None):
        # type: (Any, Any, Any, Any) -> Any
        if underscore is None:
            return self.represent_scalar(u'tag:yaml.org,2002:int', prefix + s, anchor=anchor)
        if underscore[0]:
            sl = list(s)
            pos = len(s) - underscore[0]
            while pos > 0:
                sl.insert(pos, '_')
                pos -= underscore[0]
            s = "".join(sl)
        if underscore[1]:
            s = '_' + s
        if underscore[2]:
            s += '_'
        return self.represent_scalar(u'tag:yaml.org,2002:int', prefix + s, anchor=anchor)

    def represent_scalar_int(self, data):
        # type: (Any) -> Any
        if data._width is not None:
            s = '{:0{}d}'.format(data, data._width)
        else:
            s = format(data, 'd')
        anchor = data.yaml_anchor(any=True)
        return self.insert_underscore("", s, data._underscore, anchor=anchor)

    def represent_binary_int(self, data):
        # type: (Any) -> Any
        if data._width is not None:
            # cannot use '{:#0{}b}', that strips the zeros
            s = '{:0{}b}'.format(data, data._width)
        else:
            s = format(data, 'b')
        anchor = data.yaml_anchor(any=True)
        return self.insert_underscore('0b', s, data._underscore, anchor=anchor)

    def represent_octal_int(self, data):
        # type: (Any) -> Any
        if data._width is not None:
            # cannot use '{:#0{}o}', that strips the zeros
            s = '{:0{}o}'.format(data, data._width)
        else:
            s = format(data, 'o')
        anchor = data.yaml_anchor(any=True)
        return self.insert_underscore('0o', s, data._underscore, anchor=anchor)

    def represent_hex_int(self, data):
        # type: (Any) -> Any
        if data._width is not None:
            # cannot use '{:#0{}x}', that strips the zeros
            s = '{:0{}x}'.format(data, data._width)
        else:
            s = format(data, 'x')
        anchor = data.yaml_anchor(any=True)
        return self.insert_underscore('0x', s, data._underscore, anchor=anchor)

    def represent_hex_caps_int(self, data):
        # type: (Any) -> Any
        if data._width is not None:
            # cannot use '{:#0{}X}', that strips the zeros
            s = '{:0{}X}'.format(data, data._width)
        else:
            s = format(data, 'X')
        anchor = data.yaml_anchor(any=True)
        return self.insert_underscore('0x', s, data._underscore, anchor=anchor)

    def represent_scalar_float(self, data):
        # type: (Any) -> Any
        """ this is way more complicated """
        value = None
        anchor = data.yaml_anchor(any=True)
        if data != data or (data == 0.0 and data == 1.0):
            value = u'.nan'
        elif data == self.inf_value:
            value = u'.inf'
        elif data == -self.inf_value:
            value = u'-.inf'
        if value:
            return self.represent_scalar(u'tag:yaml.org,2002:float', value, anchor=anchor)
        if data._exp is None and data._prec > 0 and data._prec == data._width - 1:
            # no exponent, but trailing dot
            value = u'{}{:d}.'.format(data._m_sign if data._m_sign else "", abs(int(data)))
        elif data._exp is None:
            # no exponent, "normal" dot
            prec = data._prec
            ms = data._m_sign if data._m_sign else ""
            # -1 for the dot
            value = u'{}{:0{}.{}f}'.format(
                ms, abs(data), data._width - len(ms), data._width - prec - 1
            )
            if prec == 0 or (prec == 1 and ms != ""):
                value = value.replace(u'0.', u'.')
            while len(value) < data._width:
                value += u'0'
        else:
            # exponent
            m, es = u'{:{}.{}e}'.format(
                # data, data._width, data._width - data._prec + (1 if data._m_sign else 0)
                data,
                data._width,
                data._width + (1 if data._m_sign else 0),
            ).split('e')
            w = data._width if data._prec > 0 else (data._width + 1)
            if data < 0:
                w += 1
            m = m[:w]
            e = int(es)
            m1, m2 = m.split('.')  # always second?
            while len(m1) + len(m2) < data._width - (1 if data._prec >= 0 else 0):
                m2 += u'0'
            if data._m_sign and data > 0:
                m1 = '+' + m1
            esgn = u'+' if data._e_sign else ""
            if data._prec < 0:  # mantissa without dot
                if m2 != u'0':
                    e -= len(m2)
                else:
                    m2 = ""
                while (len(m1) + len(m2) - (1 if data._m_sign else 0)) < data._width:
                    m2 += u'0'
                    e -= 1
                value = m1 + m2 + data._exp + u'{:{}0{}d}'.format(e, esgn, data._e_width)
            elif data._prec == 0:  # mantissa with trailing dot
                e -= len(m2)
                value = (
                    m1 + m2 + u'.' + data._exp + u'{:{}0{}d}'.format(e, esgn, data._e_width)
                )
            else:
                if data._m_lead0 > 0:
                    m2 = u'0' * (data._m_lead0 - 1) + m1 + m2
                    m1 = u'0'
                    m2 = m2[: -data._m_lead0]  # these should be zeros
                    e += data._m_lead0
                while len(m1) < data._prec:
                    m1 += m2[0]
                    m2 = m2[1:]
                    e -= 1
                value = (
                    m1 + u'.' + m2 + data._exp + u'{:{}0{}d}'.format(e, esgn, data._e_width)
                )

        if value is None:
            value = to_unicode(repr(data)).lower()
        return self.represent_scalar(u'tag:yaml.org,2002:float', value, anchor=anchor)

    def represent_sequence(self, tag, sequence, flow_style=None):
        # type: (Any, Any, Any) -> Any
        value = []  # type: List[Any]
        # if the flow_style is None, the flow style tacked on to the object
        # explicitly will be taken. If that is None as well the default flow
        # style rules
        try:
            flow_style = sequence.fa.flow_style(flow_style)
        except AttributeError:
            flow_style = flow_style
        try:
            anchor = sequence.yaml_anchor()
        except AttributeError:
            anchor = None
        node = SequenceNode(tag, value, flow_style=flow_style, anchor=anchor)
        if self.alias_key is not None:
            self.represented_objects[self.alias_key] = node
        best_style = True
        try:
            comment = getattr(sequence, comment_attrib)
            node.comment = comment.comment
            # reset any comment already printed information
            if node.comment and node.comment[1]:
                for ct in node.comment[1]:
                    ct.reset()
            item_comments = comment.items
            for v in item_comments.values():
                if v and v[1]:
                    for ct in v[1]:
                        ct.reset()
            item_comments = comment.items
            node.comment = comment.comment
            try:
                node.comment.append(comment.end)
            except AttributeError:
                pass
        except AttributeError:
            item_comments = {}
        for idx, item in enumerate(sequence):
            node_item = self.represent_data(item)
            self.merge_comments(node_item, item_comments.get(idx))
            if not (isinstance(node_item, ScalarNode) and not node_item.style):
                best_style = False
            value.append(node_item)
        if flow_style is None:
            if len(sequence) != 0 and self.default_flow_style is not None:
                node.flow_style = self.default_flow_style
            else:
                node.flow_style = best_style
        return node

    def merge_comments(self, node, comments):
        # type: (Any, Any) -> Any
        if comments is None:
            assert hasattr(node, 'comment')
            return node
        if getattr(node, 'comment', None) is not None:
            for idx, val in enumerate(comments):
                if idx >= len(node.comment):
                    continue
                nc = node.comment[idx]
                if nc is not None:
                    assert val is None or val == nc
                    comments[idx] = nc
        node.comment = comments
        return node

    def represent_key(self, data):
        # type: (Any) -> Any
        if isinstance(data, CommentedKeySeq):
            self.alias_key = None
            return self.represent_sequence(u'tag:yaml.org,2002:seq', data, flow_style=True)
        if isinstance(data, CommentedKeyMap):
            self.alias_key = None
            return self.represent_mapping(u'tag:yaml.org,2002:map', data, flow_style=True)
        return SafeRepresenter.represent_key(self, data)

    def represent_mapping(self, tag, mapping, flow_style=None):
        # type: (Any, Any, Any) -> Any
        value = []  # type: List[Any]
        try:
            flow_style = mapping.fa.flow_style(flow_style)
        except AttributeError:
            flow_style = flow_style
        try:
            anchor = mapping.yaml_anchor()
        except AttributeError:
            anchor = None
        node = MappingNode(tag, value, flow_style=flow_style, anchor=anchor)
        if self.alias_key is not None:
            self.represented_objects[self.alias_key] = node
        best_style = True
        # no sorting! !!
        try:
            comment = getattr(mapping, comment_attrib)
            node.comment = comment.comment
            if node.comment and node.comment[1]:
                for ct in node.comment[1]:
                    ct.reset()
            item_comments = comment.items
            for v in item_comments.values():
                if v and v[1]:
                    for ct in v[1]:
                        ct.reset()
            try:
                node.comment.append(comment.end)
            except AttributeError:
                pass
        except AttributeError:
            item_comments = {}
        merge_list = [m[1] for m in getattr(mapping, merge_attrib, [])]
        try:
            merge_pos = getattr(mapping, merge_attrib, [[0]])[0][0]
        except IndexError:
            merge_pos = 0
        item_count = 0
        if bool(merge_list):
            items = mapping.non_merged_items()
        else:
            items = mapping.items()
        for item_key, item_value in items:
            item_count += 1
            node_key = self.represent_key(item_key)
            node_value = self.represent_data(item_value)
            item_comment = item_comments.get(item_key)
            if item_comment:
                assert getattr(node_key, 'comment', None) is None
                node_key.comment = item_comment[:2]
                nvc = getattr(node_value, 'comment', None)
                if nvc is not None:  # end comment already there
                    nvc[0] = item_comment[2]
                    nvc[1] = item_comment[3]
                else:
                    node_value.comment = item_comment[2:]
            if not (isinstance(node_key, ScalarNode) and not node_key.style):
                best_style = False
            if not (isinstance(node_value, ScalarNode) and not node_value.style):
                best_style = False
            value.append((node_key, node_value))
        if flow_style is None:
            if ((item_count != 0) or bool(merge_list)) and self.default_flow_style is not None:
                node.flow_style = self.default_flow_style
            else:
                node.flow_style = best_style
        if bool(merge_list):
            # because of the call to represent_data here, the anchors
            # are marked as being used and thereby created
            if len(merge_list) == 1:
                arg = self.represent_data(merge_list[0])
            else:
                arg = self.represent_data(merge_list)
                arg.flow_style = True
            value.insert(merge_pos, (ScalarNode(u'tag:yaml.org,2002:merge', '<<'), arg))
        return node

    def represent_omap(self, tag, omap, flow_style=None):
        # type: (Any, Any, Any) -> Any
        value = []  # type: List[Any]
        try:
            flow_style = omap.fa.flow_style(flow_style)
        except AttributeError:
            flow_style = flow_style
        try:
            anchor = omap.yaml_anchor()
        except AttributeError:
            anchor = None
        node = SequenceNode(tag, value, flow_style=flow_style, anchor=anchor)
        if self.alias_key is not None:
            self.represented_objects[self.alias_key] = node
        best_style = True
        try:
            comment = getattr(omap, comment_attrib)
            node.comment = comment.comment
            if node.comment and node.comment[1]:
                for ct in node.comment[1]:
                    ct.reset()
            item_comments = comment.items
            for v in item_comments.values():
                if v and v[1]:
                    for ct in v[1]:
                        ct.reset()
            try:
                node.comment.append(comment.end)
            except AttributeError:
                pass
        except AttributeError:
            item_comments = {}
        for item_key in omap:
            item_val = omap[item_key]
            node_item = self.represent_data({item_key: item_val})
            # node_item.flow_style = False
            # node item has two scalars in value: node_key and node_value
            item_comment = item_comments.get(item_key)
            if item_comment:
                if item_comment[1]:
                    node_item.comment = [None, item_comment[1]]
                assert getattr(node_item.value[0][0], 'comment', None) is None
                node_item.value[0][0].comment = [item_comment[0], None]
                nvc = getattr(node_item.value[0][1], 'comment', None)
                if nvc is not None:  # end comment already there
                    nvc[0] = item_comment[2]
                    nvc[1] = item_comment[3]
                else:
                    node_item.value[0][1].comment = item_comment[2:]
            # if not (isinstance(node_item, ScalarNode) \
            #    and not node_item.style):
            #     best_style = False
            value.append(node_item)
        if flow_style is None:
            if self.default_flow_style is not None:
                node.flow_style = self.default_flow_style
            else:
                node.flow_style = best_style
        return node

    def represent_set(self, setting):
        # type: (Any) -> Any
        flow_style = False
        tag = u'tag:yaml.org,2002:set'
        # return self.represent_mapping(tag, value)
        value = []  # type: List[Any]
        flow_style = setting.fa.flow_style(flow_style)
        try:
            anchor = setting.yaml_anchor()
        except AttributeError:
            anchor = None
        node = MappingNode(tag, value, flow_style=flow_style, anchor=anchor)
        if self.alias_key is not None:
            self.represented_objects[self.alias_key] = node
        best_style = True
        # no sorting! !!
        try:
            comment = getattr(setting, comment_attrib)
            node.comment = comment.comment
            if node.comment and node.comment[1]:
                for ct in node.comment[1]:
                    ct.reset()
            item_comments = comment.items
            for v in item_comments.values():
                if v and v[1]:
                    for ct in v[1]:
                        ct.reset()
            try:
                node.comment.append(comment.end)
            except AttributeError:
                pass
        except AttributeError:
            item_comments = {}
        for item_key in setting.odict:
            node_key = self.represent_key(item_key)
            node_value = self.represent_data(None)
            item_comment = item_comments.get(item_key)
            if item_comment:
                assert getattr(node_key, 'comment', None) is None
                node_key.comment = item_comment[:2]
            node_key.style = node_value.style = '?'
            if not (isinstance(node_key, ScalarNode) and not node_key.style):
                best_style = False
            if not (isinstance(node_value, ScalarNode) and not node_value.style):
                best_style = False
            value.append((node_key, node_value))
        best_style = best_style
        return node

    def represent_dict(self, data):
        # type: (Any) -> Any
        """write out tag if saved on loading"""
        try:
            t = data.tag.value
        except AttributeError:
            t = None
        if t:
            if t.startswith('!!'):
                tag = 'tag:yaml.org,2002:' + t[2:]
            else:
                tag = t
        else:
            tag = u'tag:yaml.org,2002:map'
        return self.represent_mapping(tag, data)

    def represent_list(self, data):
        # type: (Any) -> Any
        try:
            t = data.tag.value
        except AttributeError:
            t = None
        if t:
            if t.startswith('!!'):
                tag = 'tag:yaml.org,2002:' + t[2:]
            else:
                tag = t
        else:
            tag = u'tag:yaml.org,2002:seq'
        return self.represent_sequence(tag, data)

    def represent_datetime(self, data):
        # type: (Any) -> Any
        inter = 'T' if data._yaml['t'] else ' '
        _yaml = data._yaml
        if _yaml['delta']:
            data += _yaml['delta']
            value = data.isoformat(inter)
        else:
            value = data.isoformat(inter)
        if _yaml['tz']:
            value += _yaml['tz']
        return self.represent_scalar(u'tag:yaml.org,2002:timestamp', to_unicode(value))

    def represent_tagged_scalar(self, data):
        # type: (Any) -> Any
        try:
            tag = data.tag.value
        except AttributeError:
            tag = None
        try:
            anchor = data.yaml_anchor()
        except AttributeError:
            anchor = None
        return self.represent_scalar(tag, data.value, style=data.style, anchor=anchor)

    def represent_scalar_bool(self, data):
        # type: (Any) -> Any
        try:
            anchor = data.yaml_anchor()
        except AttributeError:
            anchor = None
        return SafeRepresenter.represent_bool(self, data, anchor=anchor)


RoundTripRepresenter.add_representer(type(None), RoundTripRepresenter.represent_none)

RoundTripRepresenter.add_representer(
    LiteralScalarString, RoundTripRepresenter.represent_literal_scalarstring
)

RoundTripRepresenter.add_representer(
    FoldedScalarString, RoundTripRepresenter.represent_folded_scalarstring
)

RoundTripRepresenter.add_representer(
    SingleQuotedScalarString, RoundTripRepresenter.represent_single_quoted_scalarstring
)

RoundTripRepresenter.add_representer(
    DoubleQuotedScalarString, RoundTripRepresenter.represent_double_quoted_scalarstring
)

RoundTripRepresenter.add_representer(
    PlainScalarString, RoundTripRepresenter.represent_plain_scalarstring
)

RoundTripRepresenter.add_representer(ScalarInt, RoundTripRepresenter.represent_scalar_int)

RoundTripRepresenter.add_representer(BinaryInt, RoundTripRepresenter.represent_binary_int)

RoundTripRepresenter.add_representer(OctalInt, RoundTripRepresenter.represent_octal_int)

RoundTripRepresenter.add_representer(HexInt, RoundTripRepresenter.represent_hex_int)

RoundTripRepresenter.add_representer(HexCapsInt, RoundTripRepresenter.represent_hex_caps_int)

RoundTripRepresenter.add_representer(ScalarFloat, RoundTripRepresenter.represent_scalar_float)

RoundTripRepresenter.add_representer(ScalarBoolean, RoundTripRepresenter.represent_scalar_bool)

RoundTripRepresenter.add_representer(CommentedSeq, RoundTripRepresenter.represent_list)

RoundTripRepresenter.add_representer(CommentedMap, RoundTripRepresenter.represent_dict)

RoundTripRepresenter.add_representer(
    CommentedOrderedMap, RoundTripRepresenter.represent_ordereddict
)

if sys.version_info >= (2, 7):
    import collections

    RoundTripRepresenter.add_representer(
        collections.OrderedDict, RoundTripRepresenter.represent_ordereddict
    )

RoundTripRepresenter.add_representer(CommentedSet, RoundTripRepresenter.represent_set)

RoundTripRepresenter.add_representer(
    TaggedScalar, RoundTripRepresenter.represent_tagged_scalar
)

RoundTripRepresenter.add_representer(TimeStamp, RoundTripRepresenter.represent_datetime)
