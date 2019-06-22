# coding: utf-8

from __future__ import absolute_import
from __future__ import print_function

# Emitter expects events obeying the following grammar:
# stream ::= STREAM-START document* STREAM-END
# document ::= DOCUMENT-START node DOCUMENT-END
# node ::= SCALAR | sequence | mapping
# sequence ::= SEQUENCE-START node* SEQUENCE-END
# mapping ::= MAPPING-START (node node)* MAPPING-END

import sys
from ruamel.yaml.error import YAMLError, YAMLStreamError
from ruamel.yaml.events import *  # NOQA

# fmt: off
from ruamel.yaml.compat import utf8, text_type, PY2, nprint, dbg, DBG_EVENT, \
    check_anchorname_char
# fmt: on

if False:  # MYPY
    from typing import Any, Dict, List, Union, Text, Tuple, Optional  # NOQA
    from ruamel.yaml.compat import StreamType  # NOQA

__all__ = ['Emitter', 'EmitterError']


class EmitterError(YAMLError):
    pass


class ScalarAnalysis(object):
    def __init__(
        self,
        scalar,
        empty,
        multiline,
        allow_flow_plain,
        allow_block_plain,
        allow_single_quoted,
        allow_double_quoted,
        allow_block,
    ):
        # type: (Any, Any, Any, bool, bool, bool, bool, bool) -> None
        self.scalar = scalar
        self.empty = empty
        self.multiline = multiline
        self.allow_flow_plain = allow_flow_plain
        self.allow_block_plain = allow_block_plain
        self.allow_single_quoted = allow_single_quoted
        self.allow_double_quoted = allow_double_quoted
        self.allow_block = allow_block


class Indents(object):
    # replacement for the list based stack of None/int
    def __init__(self):
        # type: () -> None
        self.values = []  # type: List[Tuple[int, bool]]

    def append(self, val, seq):
        # type: (Any, Any) -> None
        self.values.append((val, seq))

    def pop(self):
        # type: () -> Any
        return self.values.pop()[0]

    def last_seq(self):
        # type: () -> bool
        # return the seq(uence) value for the element added before the last one
        # in increase_indent()
        try:
            return self.values[-2][1]
        except IndexError:
            return False

    def seq_flow_align(self, seq_indent, column):
        # type: (int, int) -> int
        # extra spaces because of dash
        if len(self.values) < 2 or not self.values[-1][1]:
            return 0
        # -1 for the dash
        base = self.values[-1][0] if self.values[-1][0] is not None else 0
        return base + seq_indent - column - 1

    def __len__(self):
        # type: () -> int
        return len(self.values)


class Emitter(object):
    # fmt: off
    DEFAULT_TAG_PREFIXES = {
        u'!': u'!',
        u'tag:yaml.org,2002:': u'!!',
    }
    # fmt: on

    MAX_SIMPLE_KEY_LENGTH = 128

    def __init__(
        self,
        stream,
        canonical=None,
        indent=None,
        width=None,
        allow_unicode=None,
        line_break=None,
        block_seq_indent=None,
        top_level_colon_align=None,
        prefix_colon=None,
        brace_single_entry_mapping_in_flow_sequence=None,
        dumper=None,
    ):
        # type: (StreamType, Any, Optional[int], Optional[int], Optional[bool], Any, Optional[int], Optional[bool], Any, Optional[bool], Any) -> None  # NOQA
        self.dumper = dumper
        if self.dumper is not None and getattr(self.dumper, '_emitter', None) is None:
            self.dumper._emitter = self
        self.stream = stream

        # Encoding can be overriden by STREAM-START.
        self.encoding = None  # type: Optional[Text]
        self.allow_space_break = None

        # Emitter is a state machine with a stack of states to handle nested
        # structures.
        self.states = []  # type: List[Any]
        self.state = self.expect_stream_start  # type: Any

        # Current event and the event queue.
        self.events = []  # type: List[Any]
        self.event = None  # type: Any

        # The current indentation level and the stack of previous indents.
        self.indents = Indents()
        self.indent = None  # type: Optional[int]

        # flow_context is an expanding/shrinking list consisting of '{' and '['
        # for each unclosed flow context. If empty list that means block context
        self.flow_context = []  # type: List[Text]

        # Contexts.
        self.root_context = False
        self.sequence_context = False
        self.mapping_context = False
        self.simple_key_context = False

        # Characteristics of the last emitted character:
        #  - current position.
        #  - is it a whitespace?
        #  - is it an indention character
        #    (indentation space, '-', '?', or ':')?
        self.line = 0
        self.column = 0
        self.whitespace = True
        self.indention = True
        self.compact_seq_seq = True  # dash after dash
        self.compact_seq_map = True  # key after dash
        # self.compact_ms = False   # dash after key, only when excplicit key with ?
        self.no_newline = None  # type: Optional[bool]  # set if directly after `- `

        # Whether the document requires an explicit document indicator
        self.open_ended = False

        # colon handling
        self.colon = u':'
        self.prefixed_colon = self.colon if prefix_colon is None else prefix_colon + self.colon
        # single entry mappings in flow sequence
        self.brace_single_entry_mapping_in_flow_sequence = (
            brace_single_entry_mapping_in_flow_sequence
        )  # NOQA

        # Formatting details.
        self.canonical = canonical
        self.allow_unicode = allow_unicode
        # set to False to get "\Uxxxxxxxx" for non-basic unicode like emojis
        self.unicode_supplementary = sys.maxunicode > 0xffff
        self.sequence_dash_offset = block_seq_indent if block_seq_indent else 0
        self.top_level_colon_align = top_level_colon_align
        self.best_sequence_indent = 2
        self.requested_indent = indent  # specific for literal zero indent
        if indent and 1 < indent < 10:
            self.best_sequence_indent = indent
        self.best_map_indent = self.best_sequence_indent
        # if self.best_sequence_indent < self.sequence_dash_offset + 1:
        #     self.best_sequence_indent = self.sequence_dash_offset + 1
        self.best_width = 80
        if width and width > self.best_sequence_indent * 2:
            self.best_width = width
        self.best_line_break = u'\n'  # type: Any
        if line_break in [u'\r', u'\n', u'\r\n']:
            self.best_line_break = line_break

        # Tag prefixes.
        self.tag_prefixes = None  # type: Any

        # Prepared anchor and tag.
        self.prepared_anchor = None  # type: Any
        self.prepared_tag = None  # type: Any

        # Scalar analysis and style.
        self.analysis = None  # type: Any
        self.style = None  # type: Any

    @property
    def stream(self):
        # type: () -> Any
        try:
            return self._stream
        except AttributeError:
            raise YAMLStreamError('output stream needs to specified')

    @stream.setter
    def stream(self, val):
        # type: (Any) -> None
        if val is None:
            return
        if not hasattr(val, 'write'):
            raise YAMLStreamError('stream argument needs to have a write() method')
        self._stream = val

    @property
    def serializer(self):
        # type: () -> Any
        try:
            if hasattr(self.dumper, 'typ'):
                return self.dumper.serializer
            return self.dumper._serializer
        except AttributeError:
            return self  # cyaml

    @property
    def flow_level(self):
        # type: () -> int
        return len(self.flow_context)

    def dispose(self):
        # type: () -> None
        # Reset the state attributes (to clear self-references)
        self.states = []  # type: List[Any]
        self.state = None

    def emit(self, event):
        # type: (Any) -> None
        if dbg(DBG_EVENT):
            nprint(event)
        self.events.append(event)
        while not self.need_more_events():
            self.event = self.events.pop(0)
            self.state()
            self.event = None

    # In some cases, we wait for a few next events before emitting.

    def need_more_events(self):
        # type: () -> bool
        if not self.events:
            return True
        event = self.events[0]
        if isinstance(event, DocumentStartEvent):
            return self.need_events(1)
        elif isinstance(event, SequenceStartEvent):
            return self.need_events(2)
        elif isinstance(event, MappingStartEvent):
            return self.need_events(3)
        else:
            return False

    def need_events(self, count):
        # type: (int) -> bool
        level = 0
        for event in self.events[1:]:
            if isinstance(event, (DocumentStartEvent, CollectionStartEvent)):
                level += 1
            elif isinstance(event, (DocumentEndEvent, CollectionEndEvent)):
                level -= 1
            elif isinstance(event, StreamEndEvent):
                level = -1
            if level < 0:
                return False
        return len(self.events) < count + 1

    def increase_indent(self, flow=False, sequence=None, indentless=False):
        # type: (bool, Optional[bool], bool) -> None
        self.indents.append(self.indent, sequence)
        if self.indent is None:  # top level
            if flow:
                # self.indent = self.best_sequence_indent if self.indents.last_seq() else \
                #              self.best_map_indent
                # self.indent = self.best_sequence_indent
                self.indent = self.requested_indent
            else:
                self.indent = 0
        elif not indentless:
            self.indent += (
                self.best_sequence_indent if self.indents.last_seq() else self.best_map_indent
            )
            # if self.indents.last_seq():
            #     if self.indent == 0: # top level block sequence
            #         self.indent = self.best_sequence_indent - self.sequence_dash_offset
            #     else:
            #         self.indent += self.best_sequence_indent
            # else:
            #     self.indent += self.best_map_indent

    # States.

    # Stream handlers.

    def expect_stream_start(self):
        # type: () -> None
        if isinstance(self.event, StreamStartEvent):
            if PY2:
                if self.event.encoding and not getattr(self.stream, 'encoding', None):
                    self.encoding = self.event.encoding
            else:
                if self.event.encoding and not hasattr(self.stream, 'encoding'):
                    self.encoding = self.event.encoding
            self.write_stream_start()
            self.state = self.expect_first_document_start
        else:
            raise EmitterError('expected StreamStartEvent, but got %s' % (self.event,))

    def expect_nothing(self):
        # type: () -> None
        raise EmitterError('expected nothing, but got %s' % (self.event,))

    # Document handlers.

    def expect_first_document_start(self):
        # type: () -> Any
        return self.expect_document_start(first=True)

    def expect_document_start(self, first=False):
        # type: (bool) -> None
        if isinstance(self.event, DocumentStartEvent):
            if (self.event.version or self.event.tags) and self.open_ended:
                self.write_indicator(u'...', True)
                self.write_indent()
            if self.event.version:
                version_text = self.prepare_version(self.event.version)
                self.write_version_directive(version_text)
            self.tag_prefixes = self.DEFAULT_TAG_PREFIXES.copy()
            if self.event.tags:
                handles = sorted(self.event.tags.keys())
                for handle in handles:
                    prefix = self.event.tags[handle]
                    self.tag_prefixes[prefix] = handle
                    handle_text = self.prepare_tag_handle(handle)
                    prefix_text = self.prepare_tag_prefix(prefix)
                    self.write_tag_directive(handle_text, prefix_text)
            implicit = (
                first
                and not self.event.explicit
                and not self.canonical
                and not self.event.version
                and not self.event.tags
                and not self.check_empty_document()
            )
            if not implicit:
                self.write_indent()
                self.write_indicator(u'---', True)
                if self.canonical:
                    self.write_indent()
            self.state = self.expect_document_root
        elif isinstance(self.event, StreamEndEvent):
            if self.open_ended:
                self.write_indicator(u'...', True)
                self.write_indent()
            self.write_stream_end()
            self.state = self.expect_nothing
        else:
            raise EmitterError('expected DocumentStartEvent, but got %s' % (self.event,))

    def expect_document_end(self):
        # type: () -> None
        if isinstance(self.event, DocumentEndEvent):
            self.write_indent()
            if self.event.explicit:
                self.write_indicator(u'...', True)
                self.write_indent()
            self.flush_stream()
            self.state = self.expect_document_start
        else:
            raise EmitterError('expected DocumentEndEvent, but got %s' % (self.event,))

    def expect_document_root(self):
        # type: () -> None
        self.states.append(self.expect_document_end)
        self.expect_node(root=True)

    # Node handlers.

    def expect_node(self, root=False, sequence=False, mapping=False, simple_key=False):
        # type: (bool, bool, bool, bool) -> None
        self.root_context = root
        self.sequence_context = sequence  # not used in PyYAML
        self.mapping_context = mapping
        self.simple_key_context = simple_key
        if isinstance(self.event, AliasEvent):
            self.expect_alias()
        elif isinstance(self.event, (ScalarEvent, CollectionStartEvent)):
            if (
                self.process_anchor(u'&')
                and isinstance(self.event, ScalarEvent)
                and self.sequence_context
            ):
                self.sequence_context = False
            self.process_tag()
            if isinstance(self.event, ScalarEvent):
                # nprint('@', self.indention, self.no_newline, self.column)
                self.expect_scalar()
            elif isinstance(self.event, SequenceStartEvent):
                # nprint('@', self.indention, self.no_newline, self.column)
                i2, n2 = self.indention, self.no_newline  # NOQA
                if self.event.comment:
                    if self.event.flow_style is False and self.event.comment:
                        if self.write_post_comment(self.event):
                            self.indention = False
                            self.no_newline = True
                    if self.write_pre_comment(self.event):
                        self.indention = i2
                        self.no_newline = not self.indention
                if (
                    self.flow_level
                    or self.canonical
                    or self.event.flow_style
                    or self.check_empty_sequence()
                ):
                    self.expect_flow_sequence()
                else:
                    self.expect_block_sequence()
            elif isinstance(self.event, MappingStartEvent):
                if self.event.flow_style is False and self.event.comment:
                    self.write_post_comment(self.event)
                if self.event.comment and self.event.comment[1]:
                    self.write_pre_comment(self.event)
                if (
                    self.flow_level
                    or self.canonical
                    or self.event.flow_style
                    or self.check_empty_mapping()
                ):
                    self.expect_flow_mapping(single=self.event.nr_items == 1)
                else:
                    self.expect_block_mapping()
        else:
            raise EmitterError('expected NodeEvent, but got %s' % (self.event,))

    def expect_alias(self):
        # type: () -> None
        if self.event.anchor is None:
            raise EmitterError('anchor is not specified for alias')
        self.process_anchor(u'*')
        self.state = self.states.pop()

    def expect_scalar(self):
        # type: () -> None
        self.increase_indent(flow=True)
        self.process_scalar()
        self.indent = self.indents.pop()
        self.state = self.states.pop()

    # Flow sequence handlers.

    def expect_flow_sequence(self):
        # type: () -> None
        ind = self.indents.seq_flow_align(self.best_sequence_indent, self.column)
        self.write_indicator(u' ' * ind + u'[', True, whitespace=True)
        self.increase_indent(flow=True, sequence=True)
        self.flow_context.append('[')
        self.state = self.expect_first_flow_sequence_item

    def expect_first_flow_sequence_item(self):
        # type: () -> None
        if isinstance(self.event, SequenceEndEvent):
            self.indent = self.indents.pop()
            popped = self.flow_context.pop()
            assert popped == '['
            self.write_indicator(u']', False)
            if self.event.comment and self.event.comment[0]:
                # eol comment on empty flow sequence
                self.write_post_comment(self.event)
            elif self.flow_level == 0:
                self.write_line_break()
            self.state = self.states.pop()
        else:
            if self.canonical or self.column > self.best_width:
                self.write_indent()
            self.states.append(self.expect_flow_sequence_item)
            self.expect_node(sequence=True)

    def expect_flow_sequence_item(self):
        # type: () -> None
        if isinstance(self.event, SequenceEndEvent):
            self.indent = self.indents.pop()
            popped = self.flow_context.pop()
            assert popped == '['
            if self.canonical:
                self.write_indicator(u',', False)
                self.write_indent()
            self.write_indicator(u']', False)
            if self.event.comment and self.event.comment[0]:
                # eol comment on flow sequence
                self.write_post_comment(self.event)
            else:
                self.no_newline = False
            self.state = self.states.pop()
        else:
            self.write_indicator(u',', False)
            if self.canonical or self.column > self.best_width:
                self.write_indent()
            self.states.append(self.expect_flow_sequence_item)
            self.expect_node(sequence=True)

    # Flow mapping handlers.

    def expect_flow_mapping(self, single=False):
        # type: (Optional[bool]) -> None
        ind = self.indents.seq_flow_align(self.best_sequence_indent, self.column)
        map_init = u'{'
        if (
            single
            and self.flow_level
            and self.flow_context[-1] == '['
            and not self.canonical
            and not self.brace_single_entry_mapping_in_flow_sequence
        ):
            # single map item with flow context, no curly braces necessary
            map_init = u''
        self.write_indicator(u' ' * ind + map_init, True, whitespace=True)
        self.flow_context.append(map_init)
        self.increase_indent(flow=True, sequence=False)
        self.state = self.expect_first_flow_mapping_key

    def expect_first_flow_mapping_key(self):
        # type: () -> None
        if isinstance(self.event, MappingEndEvent):
            self.indent = self.indents.pop()
            popped = self.flow_context.pop()
            assert popped == '{'  # empty flow mapping
            self.write_indicator(u'}', False)
            if self.event.comment and self.event.comment[0]:
                # eol comment on empty mapping
                self.write_post_comment(self.event)
            elif self.flow_level == 0:
                self.write_line_break()
            self.state = self.states.pop()
        else:
            if self.canonical or self.column > self.best_width:
                self.write_indent()
            if not self.canonical and self.check_simple_key():
                self.states.append(self.expect_flow_mapping_simple_value)
                self.expect_node(mapping=True, simple_key=True)
            else:
                self.write_indicator(u'?', True)
                self.states.append(self.expect_flow_mapping_value)
                self.expect_node(mapping=True)

    def expect_flow_mapping_key(self):
        # type: () -> None
        if isinstance(self.event, MappingEndEvent):
            # if self.event.comment and self.event.comment[1]:
            #     self.write_pre_comment(self.event)
            self.indent = self.indents.pop()
            popped = self.flow_context.pop()
            assert popped in [u'{', u'']
            if self.canonical:
                self.write_indicator(u',', False)
                self.write_indent()
            if popped != u'':
                self.write_indicator(u'}', False)
            if self.event.comment and self.event.comment[0]:
                # eol comment on flow mapping, never reached on empty mappings
                self.write_post_comment(self.event)
            else:
                self.no_newline = False
            self.state = self.states.pop()
        else:
            self.write_indicator(u',', False)
            if self.canonical or self.column > self.best_width:
                self.write_indent()
            if not self.canonical and self.check_simple_key():
                self.states.append(self.expect_flow_mapping_simple_value)
                self.expect_node(mapping=True, simple_key=True)
            else:
                self.write_indicator(u'?', True)
                self.states.append(self.expect_flow_mapping_value)
                self.expect_node(mapping=True)

    def expect_flow_mapping_simple_value(self):
        # type: () -> None
        self.write_indicator(self.prefixed_colon, False)
        self.states.append(self.expect_flow_mapping_key)
        self.expect_node(mapping=True)

    def expect_flow_mapping_value(self):
        # type: () -> None
        if self.canonical or self.column > self.best_width:
            self.write_indent()
        self.write_indicator(self.prefixed_colon, True)
        self.states.append(self.expect_flow_mapping_key)
        self.expect_node(mapping=True)

    # Block sequence handlers.

    def expect_block_sequence(self):
        # type: () -> None
        if self.mapping_context:
            indentless = not self.indention
        else:
            indentless = False
            if not self.compact_seq_seq and self.column != 0:
                self.write_line_break()
        self.increase_indent(flow=False, sequence=True, indentless=indentless)
        self.state = self.expect_first_block_sequence_item

    def expect_first_block_sequence_item(self):
        # type: () -> Any
        return self.expect_block_sequence_item(first=True)

    def expect_block_sequence_item(self, first=False):
        # type: (bool) -> None
        if not first and isinstance(self.event, SequenceEndEvent):
            if self.event.comment and self.event.comment[1]:
                # final comments on a block list e.g. empty line
                self.write_pre_comment(self.event)
            self.indent = self.indents.pop()
            self.state = self.states.pop()
            self.no_newline = False
        else:
            if self.event.comment and self.event.comment[1]:
                self.write_pre_comment(self.event)
            nonl = self.no_newline if self.column == 0 else False
            self.write_indent()
            ind = self.sequence_dash_offset  # if  len(self.indents) > 1 else 0
            self.write_indicator(u' ' * ind + u'-', True, indention=True)
            if nonl or self.sequence_dash_offset + 2 > self.best_sequence_indent:
                self.no_newline = True
            self.states.append(self.expect_block_sequence_item)
            self.expect_node(sequence=True)

    # Block mapping handlers.

    def expect_block_mapping(self):
        # type: () -> None
        if not self.mapping_context and not (self.compact_seq_map or self.column == 0):
            self.write_line_break()
        self.increase_indent(flow=False, sequence=False)
        self.state = self.expect_first_block_mapping_key

    def expect_first_block_mapping_key(self):
        # type: () -> None
        return self.expect_block_mapping_key(first=True)

    def expect_block_mapping_key(self, first=False):
        # type: (Any) -> None
        if not first and isinstance(self.event, MappingEndEvent):
            if self.event.comment and self.event.comment[1]:
                # final comments from a doc
                self.write_pre_comment(self.event)
            self.indent = self.indents.pop()
            self.state = self.states.pop()
        else:
            if self.event.comment and self.event.comment[1]:
                # final comments from a doc
                self.write_pre_comment(self.event)
            self.write_indent()
            if self.check_simple_key():
                if not isinstance(
                    self.event, (SequenceStartEvent, MappingStartEvent)
                ):  # sequence keys
                    try:
                        if self.event.style == '?':
                            self.write_indicator(u'?', True, indention=True)
                    except AttributeError:  # aliases have no style
                        pass
                self.states.append(self.expect_block_mapping_simple_value)
                self.expect_node(mapping=True, simple_key=True)
                if isinstance(self.event, AliasEvent):
                    self.stream.write(u' ')
            else:
                self.write_indicator(u'?', True, indention=True)
                self.states.append(self.expect_block_mapping_value)
                self.expect_node(mapping=True)

    def expect_block_mapping_simple_value(self):
        # type: () -> None
        if getattr(self.event, 'style', None) != '?':
            # prefix = u''
            if self.indent == 0 and self.top_level_colon_align is not None:
                # write non-prefixed colon
                c = u' ' * (self.top_level_colon_align - self.column) + self.colon
            else:
                c = self.prefixed_colon
            self.write_indicator(c, False)
        self.states.append(self.expect_block_mapping_key)
        self.expect_node(mapping=True)

    def expect_block_mapping_value(self):
        # type: () -> None
        self.write_indent()
        self.write_indicator(self.prefixed_colon, True, indention=True)
        self.states.append(self.expect_block_mapping_key)
        self.expect_node(mapping=True)

    # Checkers.

    def check_empty_sequence(self):
        # type: () -> bool
        return (
            isinstance(self.event, SequenceStartEvent)
            and bool(self.events)
            and isinstance(self.events[0], SequenceEndEvent)
        )

    def check_empty_mapping(self):
        # type: () -> bool
        return (
            isinstance(self.event, MappingStartEvent)
            and bool(self.events)
            and isinstance(self.events[0], MappingEndEvent)
        )

    def check_empty_document(self):
        # type: () -> bool
        if not isinstance(self.event, DocumentStartEvent) or not self.events:
            return False
        event = self.events[0]
        return (
            isinstance(event, ScalarEvent)
            and event.anchor is None
            and event.tag is None
            and event.implicit
            and event.value == ""
        )

    def check_simple_key(self):
        # type: () -> bool
        length = 0
        if isinstance(self.event, NodeEvent) and self.event.anchor is not None:
            if self.prepared_anchor is None:
                self.prepared_anchor = self.prepare_anchor(self.event.anchor)
            length += len(self.prepared_anchor)
        if (
            isinstance(self.event, (ScalarEvent, CollectionStartEvent))
            and self.event.tag is not None
        ):
            if self.prepared_tag is None:
                self.prepared_tag = self.prepare_tag(self.event.tag)
            length += len(self.prepared_tag)
        if isinstance(self.event, ScalarEvent):
            if self.analysis is None:
                self.analysis = self.analyze_scalar(self.event.value)
            length += len(self.analysis.scalar)
        return length < self.MAX_SIMPLE_KEY_LENGTH and (
            isinstance(self.event, AliasEvent)
            or (isinstance(self.event, SequenceStartEvent) and self.event.flow_style is True)
            or (isinstance(self.event, MappingStartEvent) and self.event.flow_style is True)
            or (
                isinstance(self.event, ScalarEvent)
                and not self.analysis.empty
                and not self.analysis.multiline
            )
            or self.check_empty_sequence()
            or self.check_empty_mapping()
        )

    # Anchor, Tag, and Scalar processors.

    def process_anchor(self, indicator):
        # type: (Any) -> bool
        if self.event.anchor is None:
            self.prepared_anchor = None
            return False
        if self.prepared_anchor is None:
            self.prepared_anchor = self.prepare_anchor(self.event.anchor)
        if self.prepared_anchor:
            self.write_indicator(indicator + self.prepared_anchor, True)
            # issue 288
            self.no_newline = False
        self.prepared_anchor = None
        return True

    def process_tag(self):
        # type: () -> None
        tag = self.event.tag
        if isinstance(self.event, ScalarEvent):
            if self.style is None:
                self.style = self.choose_scalar_style()
            if (not self.canonical or tag is None) and (
                (self.style == "" and self.event.implicit[0])
                or (self.style != "" and self.event.implicit[1])
            ):
                self.prepared_tag = None
                return
            if self.event.implicit[0] and tag is None:
                tag = u'!'
                self.prepared_tag = None
        else:
            if (not self.canonical or tag is None) and self.event.implicit:
                self.prepared_tag = None
                return
        if tag is None:
            raise EmitterError('tag is not specified')
        if self.prepared_tag is None:
            self.prepared_tag = self.prepare_tag(tag)
        if self.prepared_tag:
            self.write_indicator(self.prepared_tag, True)
            if (
                self.sequence_context
                and not self.flow_level
                and isinstance(self.event, ScalarEvent)
            ):
                self.no_newline = True
        self.prepared_tag = None

    def choose_scalar_style(self):
        # type: () -> Any
        if self.analysis is None:
            self.analysis = self.analyze_scalar(self.event.value)
        if self.event.style == '"' or self.canonical:
            return '"'
        if (not self.event.style or self.event.style == '?') and (
            self.event.implicit[0] or not self.event.implicit[2]
        ):
            if not (
                self.simple_key_context and (self.analysis.empty or self.analysis.multiline)
            ) and (
                self.flow_level
                and self.analysis.allow_flow_plain
                or (not self.flow_level and self.analysis.allow_block_plain)
            ):
                return ""
        self.analysis.allow_block = True
        if self.event.style and self.event.style in '|>':
            if (
                not self.flow_level
                and not self.simple_key_context
                and self.analysis.allow_block
            ):
                return self.event.style
        if not self.event.style and self.analysis.allow_double_quoted:
            if "'" in self.event.value or '\n' in self.event.value:
                return '"'
        if not self.event.style or self.event.style == "'":
            if self.analysis.allow_single_quoted and not (
                self.simple_key_context and self.analysis.multiline
            ):
                return "'"
        return '"'

    def process_scalar(self):
        # type: () -> None
        if self.analysis is None:
            self.analysis = self.analyze_scalar(self.event.value)
        if self.style is None:
            self.style = self.choose_scalar_style()
        split = not self.simple_key_context
        # if self.analysis.multiline and split    \
        #         and (not self.style or self.style in '\'\"'):
        #     self.write_indent()
        # nprint('xx', self.sequence_context, self.flow_level)
        if self.sequence_context and not self.flow_level:
            self.write_indent()
        if self.style == '"':
            self.write_double_quoted(self.analysis.scalar, split)
        elif self.style == "'":
            self.write_single_quoted(self.analysis.scalar, split)
        elif self.style == '>':
            self.write_folded(self.analysis.scalar)
        elif self.style == '|':
            self.write_literal(self.analysis.scalar, self.event.comment)
        else:
            self.write_plain(self.analysis.scalar, split)
        self.analysis = None
        self.style = None
        if self.event.comment:
            self.write_post_comment(self.event)

    # Analyzers.

    def prepare_version(self, version):
        # type: (Any) -> Any
        major, minor = version
        if major != 1:
            raise EmitterError('unsupported YAML version: %d.%d' % (major, minor))
        return u'%d.%d' % (major, minor)

    def prepare_tag_handle(self, handle):
        # type: (Any) -> Any
        if not handle:
            raise EmitterError('tag handle must not be empty')
        if handle[0] != u'!' or handle[-1] != u'!':
            raise EmitterError("tag handle must start and end with '!': %r" % (utf8(handle)))
        for ch in handle[1:-1]:
            if not (
                u'0' <= ch <= u'9' or u'A' <= ch <= u'Z' or u'a' <= ch <= u'z' or ch in u'-_'
            ):
                raise EmitterError(
                    'invalid character %r in the tag handle: %r' % (utf8(ch), utf8(handle))
                )
        return handle

    def prepare_tag_prefix(self, prefix):
        # type: (Any) -> Any
        if not prefix:
            raise EmitterError('tag prefix must not be empty')
        chunks = []  # type: List[Any]
        start = end = 0
        if prefix[0] == u'!':
            end = 1
        while end < len(prefix):
            ch = prefix[end]
            if (
                u'0' <= ch <= u'9'
                or u'A' <= ch <= u'Z'
                or u'a' <= ch <= u'z'
                or ch in u"-;/?!:@&=+$,_.~*'()[]"
            ):
                end += 1
            else:
                if start < end:
                    chunks.append(prefix[start:end])
                start = end = end + 1
                data = utf8(ch)
                for ch in data:
                    chunks.append(u'%%%02X' % ord(ch))
        if start < end:
            chunks.append(prefix[start:end])
        return "".join(chunks)

    def prepare_tag(self, tag):
        # type: (Any) -> Any
        if not tag:
            raise EmitterError('tag must not be empty')
        if tag == u'!':
            return tag
        handle = None
        suffix = tag
        prefixes = sorted(self.tag_prefixes.keys())
        for prefix in prefixes:
            if tag.startswith(prefix) and (prefix == u'!' or len(prefix) < len(tag)):
                handle = self.tag_prefixes[prefix]
                suffix = tag[len(prefix) :]
        chunks = []  # type: List[Any]
        start = end = 0
        while end < len(suffix):
            ch = suffix[end]
            if (
                u'0' <= ch <= u'9'
                or u'A' <= ch <= u'Z'
                or u'a' <= ch <= u'z'
                or ch in u"-;/?:@&=+$,_.~*'()[]"
                or (ch == u'!' and handle != u'!')
            ):
                end += 1
            else:
                if start < end:
                    chunks.append(suffix[start:end])
                start = end = end + 1
                data = utf8(ch)
                for ch in data:
                    chunks.append(u'%%%02X' % ord(ch))
        if start < end:
            chunks.append(suffix[start:end])
        suffix_text = "".join(chunks)
        if handle:
            return u'%s%s' % (handle, suffix_text)
        else:
            return u'!<%s>' % suffix_text

    def prepare_anchor(self, anchor):
        # type: (Any) -> Any
        if not anchor:
            raise EmitterError('anchor must not be empty')
        for ch in anchor:
            if not check_anchorname_char(ch):
                raise EmitterError(
                    'invalid character %r in the anchor: %r' % (utf8(ch), utf8(anchor))
                )
        return anchor

    def analyze_scalar(self, scalar):
        # type: (Any) -> Any
        # Empty scalar is a special case.
        if not scalar:
            return ScalarAnalysis(
                scalar=scalar,
                empty=True,
                multiline=False,
                allow_flow_plain=False,
                allow_block_plain=True,
                allow_single_quoted=True,
                allow_double_quoted=True,
                allow_block=False,
            )

        # Indicators and special characters.
        block_indicators = False
        flow_indicators = False
        line_breaks = False
        special_characters = False

        # Important whitespace combinations.
        leading_space = False
        leading_break = False
        trailing_space = False
        trailing_break = False
        break_space = False
        space_break = False

        # Check document indicators.
        if scalar.startswith(u'---') or scalar.startswith(u'...'):
            block_indicators = True
            flow_indicators = True

        # First character or preceded by a whitespace.
        preceeded_by_whitespace = True

        # Last character or followed by a whitespace.
        followed_by_whitespace = len(scalar) == 1 or scalar[1] in u'\0 \t\r\n\x85\u2028\u2029'

        # The previous character is a space.
        previous_space = False

        # The previous character is a break.
        previous_break = False

        index = 0
        while index < len(scalar):
            ch = scalar[index]

            # Check for indicators.
            if index == 0:
                # Leading indicators are special characters.
                if ch in u'#,[]{}&*!|>\'"%@`':
                    flow_indicators = True
                    block_indicators = True
                if ch in u'?:':  # ToDo
                    if self.serializer.use_version == (1, 1):
                        flow_indicators = True
                    elif len(scalar) == 1:  # single character
                        flow_indicators = True
                    if followed_by_whitespace:
                        block_indicators = True
                if ch == u'-' and followed_by_whitespace:
                    flow_indicators = True
                    block_indicators = True
            else:
                # Some indicators cannot appear within a scalar as well.
                if ch in u',[]{}':  # http://yaml.org/spec/1.2/spec.html#id2788859
                    flow_indicators = True
                if ch == u'?' and self.serializer.use_version == (1, 1):
                    flow_indicators = True
                if ch == u':':
                    if followed_by_whitespace:
                        flow_indicators = True
                        block_indicators = True
                if ch == u'#' and preceeded_by_whitespace:
                    flow_indicators = True
                    block_indicators = True

            # Check for line breaks, special, and unicode characters.
            if ch in u'\n\x85\u2028\u2029':
                line_breaks = True
            if not (ch == u'\n' or u'\x20' <= ch <= u'\x7E'):
                if (
                    ch == u'\x85'
                    or u'\xA0' <= ch <= u'\uD7FF'
                    or u'\uE000' <= ch <= u'\uFFFD'
                    or (self.unicode_supplementary and (u'\U00010000' <= ch <= u'\U0010FFFF'))
                ) and ch != u'\uFEFF':
                    # unicode_characters = True
                    if not self.allow_unicode:
                        special_characters = True
                else:
                    special_characters = True

            # Detect important whitespace combinations.
            if ch == u' ':
                if index == 0:
                    leading_space = True
                if index == len(scalar) - 1:
                    trailing_space = True
                if previous_break:
                    break_space = True
                previous_space = True
                previous_break = False
            elif ch in u'\n\x85\u2028\u2029':
                if index == 0:
                    leading_break = True
                if index == len(scalar) - 1:
                    trailing_break = True
                if previous_space:
                    space_break = True
                previous_space = False
                previous_break = True
            else:
                previous_space = False
                previous_break = False

            # Prepare for the next character.
            index += 1
            preceeded_by_whitespace = ch in u'\0 \t\r\n\x85\u2028\u2029'
            followed_by_whitespace = (
                index + 1 >= len(scalar) or scalar[index + 1] in u'\0 \t\r\n\x85\u2028\u2029'
            )

        # Let's decide what styles are allowed.
        allow_flow_plain = True
        allow_block_plain = True
        allow_single_quoted = True
        allow_double_quoted = True
        allow_block = True

        # Leading and trailing whitespaces are bad for plain scalars.
        if leading_space or leading_break or trailing_space or trailing_break:
            allow_flow_plain = allow_block_plain = False

        # We do not permit trailing spaces for block scalars.
        if trailing_space:
            allow_block = False

        # Spaces at the beginning of a new line are only acceptable for block
        # scalars.
        if break_space:
            allow_flow_plain = allow_block_plain = allow_single_quoted = False

        # Spaces followed by breaks, as well as special character are only
        # allowed for double quoted scalars.
        if special_characters:
            allow_flow_plain = allow_block_plain = allow_single_quoted = allow_block = False
        elif space_break:
            allow_flow_plain = allow_block_plain = allow_single_quoted = False
            if not self.allow_space_break:
                allow_block = False

        # Although the plain scalar writer supports breaks, we never emit
        # multiline plain scalars.
        if line_breaks:
            allow_flow_plain = allow_block_plain = False

        # Flow indicators are forbidden for flow plain scalars.
        if flow_indicators:
            allow_flow_plain = False

        # Block indicators are forbidden for block plain scalars.
        if block_indicators:
            allow_block_plain = False

        return ScalarAnalysis(
            scalar=scalar,
            empty=False,
            multiline=line_breaks,
            allow_flow_plain=allow_flow_plain,
            allow_block_plain=allow_block_plain,
            allow_single_quoted=allow_single_quoted,
            allow_double_quoted=allow_double_quoted,
            allow_block=allow_block,
        )

    # Writers.

    def flush_stream(self):
        # type: () -> None
        if hasattr(self.stream, 'flush'):
            self.stream.flush()

    def write_stream_start(self):
        # type: () -> None
        # Write BOM if needed.
        if self.encoding and self.encoding.startswith('utf-16'):
            self.stream.write(u'\uFEFF'.encode(self.encoding))

    def write_stream_end(self):
        # type: () -> None
        self.flush_stream()

    def write_indicator(self, indicator, need_whitespace, whitespace=False, indention=False):
        # type: (Any, Any, bool, bool) -> None
        if self.whitespace or not need_whitespace:
            data = indicator
        else:
            data = u' ' + indicator
        self.whitespace = whitespace
        self.indention = self.indention and indention
        self.column += len(data)
        self.open_ended = False
        if bool(self.encoding):
            data = data.encode(self.encoding)
        self.stream.write(data)

    def write_indent(self):
        # type: () -> None
        indent = self.indent or 0
        if (
            not self.indention
            or self.column > indent
            or (self.column == indent and not self.whitespace)
        ):
            if bool(self.no_newline):
                self.no_newline = False
            else:
                self.write_line_break()
        if self.column < indent:
            self.whitespace = True
            data = u' ' * (indent - self.column)
            self.column = indent
            if self.encoding:
                data = data.encode(self.encoding)
            self.stream.write(data)

    def write_line_break(self, data=None):
        # type: (Any) -> None
        if data is None:
            data = self.best_line_break
        self.whitespace = True
        self.indention = True
        self.line += 1
        self.column = 0
        if bool(self.encoding):
            data = data.encode(self.encoding)
        self.stream.write(data)

    def write_version_directive(self, version_text):
        # type: (Any) -> None
        data = u'%%YAML %s' % version_text
        if self.encoding:
            data = data.encode(self.encoding)
        self.stream.write(data)
        self.write_line_break()

    def write_tag_directive(self, handle_text, prefix_text):
        # type: (Any, Any) -> None
        data = u'%%TAG %s %s' % (handle_text, prefix_text)
        if self.encoding:
            data = data.encode(self.encoding)
        self.stream.write(data)
        self.write_line_break()

    # Scalar streams.

    def write_single_quoted(self, text, split=True):
        # type: (Any, Any) -> None
        if self.root_context:
            if self.requested_indent is not None:
                self.write_line_break()
                if self.requested_indent != 0:
                    self.write_indent()
        self.write_indicator(u"'", True)
        spaces = False
        breaks = False
        start = end = 0
        while end <= len(text):
            ch = None
            if end < len(text):
                ch = text[end]
            if spaces:
                if ch is None or ch != u' ':
                    if (
                        start + 1 == end
                        and self.column > self.best_width
                        and split
                        and start != 0
                        and end != len(text)
                    ):
                        self.write_indent()
                    else:
                        data = text[start:end]
                        self.column += len(data)
                        if bool(self.encoding):
                            data = data.encode(self.encoding)
                        self.stream.write(data)
                    start = end
            elif breaks:
                if ch is None or ch not in u'\n\x85\u2028\u2029':
                    if text[start] == u'\n':
                        self.write_line_break()
                    for br in text[start:end]:
                        if br == u'\n':
                            self.write_line_break()
                        else:
                            self.write_line_break(br)
                    self.write_indent()
                    start = end
            else:
                if ch is None or ch in u' \n\x85\u2028\u2029' or ch == u"'":
                    if start < end:
                        data = text[start:end]
                        self.column += len(data)
                        if bool(self.encoding):
                            data = data.encode(self.encoding)
                        self.stream.write(data)
                        start = end
            if ch == u"'":
                data = u"''"
                self.column += 2
                if bool(self.encoding):
                    data = data.encode(self.encoding)
                self.stream.write(data)
                start = end + 1
            if ch is not None:
                spaces = ch == u' '
                breaks = ch in u'\n\x85\u2028\u2029'
            end += 1
        self.write_indicator(u"'", False)

    ESCAPE_REPLACEMENTS = {
        u'\0': u'0',
        u'\x07': u'a',
        u'\x08': u'b',
        u'\x09': u't',
        u'\x0A': u'n',
        u'\x0B': u'v',
        u'\x0C': u'f',
        u'\x0D': u'r',
        u'\x1B': u'e',
        u'"': u'"',
        u'\\': u'\\',
        u'\x85': u'N',
        u'\xA0': u'_',
        u'\u2028': u'L',
        u'\u2029': u'P',
    }

    def write_double_quoted(self, text, split=True):
        # type: (Any, Any) -> None
        if self.root_context:
            if self.requested_indent is not None:
                self.write_line_break()
                if self.requested_indent != 0:
                    self.write_indent()
        self.write_indicator(u'"', True)
        start = end = 0
        while end <= len(text):
            ch = None
            if end < len(text):
                ch = text[end]
            if (
                ch is None
                or ch in u'"\\\x85\u2028\u2029\uFEFF'
                or not (
                    u'\x20' <= ch <= u'\x7E'
                    or (
                        self.allow_unicode
                        and (u'\xA0' <= ch <= u'\uD7FF' or u'\uE000' <= ch <= u'\uFFFD')
                    )
                )
            ):
                if start < end:
                    data = text[start:end]
                    self.column += len(data)
                    if bool(self.encoding):
                        data = data.encode(self.encoding)
                    self.stream.write(data)
                    start = end
                if ch is not None:
                    if ch in self.ESCAPE_REPLACEMENTS:
                        data = u'\\' + self.ESCAPE_REPLACEMENTS[ch]
                    elif ch <= u'\xFF':
                        data = u'\\x%02X' % ord(ch)
                    elif ch <= u'\uFFFF':
                        data = u'\\u%04X' % ord(ch)
                    else:
                        data = u'\\U%08X' % ord(ch)
                    self.column += len(data)
                    if bool(self.encoding):
                        data = data.encode(self.encoding)
                    self.stream.write(data)
                    start = end + 1
            if (
                0 < end < len(text) - 1
                and (ch == u' ' or start >= end)
                and self.column + (end - start) > self.best_width
                and split
            ):
                data = text[start:end] + u'\\'
                if start < end:
                    start = end
                self.column += len(data)
                if bool(self.encoding):
                    data = data.encode(self.encoding)
                self.stream.write(data)
                self.write_indent()
                self.whitespace = False
                self.indention = False
                if text[start] == u' ':
                    data = u'\\'
                    self.column += len(data)
                    if bool(self.encoding):
                        data = data.encode(self.encoding)
                    self.stream.write(data)
            end += 1
        self.write_indicator(u'"', False)

    def determine_block_hints(self, text):
        # type: (Any) -> Any
        indent = 0
        indicator = u''
        hints = u''
        if text:
            if text[0] in u' \n\x85\u2028\u2029':
                indent = self.best_sequence_indent
                hints += text_type(indent)
            elif self.root_context:
                for end in ['\n---', '\n...']:
                    pos = 0
                    while True:
                        pos = text.find(end, pos)
                        if pos == -1:
                            break
                        try:
                            if text[pos + 4] in ' \r\n':
                                break
                        except IndexError:
                            pass
                        pos += 1
                    if pos > -1:
                        break
                if pos > 0:
                    indent = self.best_sequence_indent
            if text[-1] not in u'\n\x85\u2028\u2029':
                indicator = u'-'
            elif len(text) == 1 or text[-2] in u'\n\x85\u2028\u2029':
                indicator = u'+'
        hints += indicator
        return hints, indent, indicator

    def write_folded(self, text):
        # type: (Any) -> None
        hints, _indent, _indicator = self.determine_block_hints(text)
        self.write_indicator(u'>' + hints, True)
        if _indicator == u'+':
            self.open_ended = True
        self.write_line_break()
        leading_space = True
        spaces = False
        breaks = True
        start = end = 0
        while end <= len(text):
            ch = None
            if end < len(text):
                ch = text[end]
            if breaks:
                if ch is None or ch not in u'\n\x85\u2028\u2029\a':
                    if (
                        not leading_space
                        and ch is not None
                        and ch != u' '
                        and text[start] == u'\n'
                    ):
                        self.write_line_break()
                    leading_space = ch == u' '
                    for br in text[start:end]:
                        if br == u'\n':
                            self.write_line_break()
                        else:
                            self.write_line_break(br)
                    if ch is not None:
                        self.write_indent()
                    start = end
            elif spaces:
                if ch != u' ':
                    if start + 1 == end and self.column > self.best_width:
                        self.write_indent()
                    else:
                        data = text[start:end]
                        self.column += len(data)
                        if bool(self.encoding):
                            data = data.encode(self.encoding)
                        self.stream.write(data)
                    start = end
            else:
                if ch is None or ch in u' \n\x85\u2028\u2029\a':
                    data = text[start:end]
                    self.column += len(data)
                    if bool(self.encoding):
                        data = data.encode(self.encoding)
                    self.stream.write(data)
                    if ch == u'\a':
                        if end < (len(text) - 1) and not text[end + 2].isspace():
                            self.write_line_break()
                            self.write_indent()
                            end += 2  # \a and the space that is inserted on the fold
                        else:
                            raise EmitterError('unexcpected fold indicator \\a before space')
                    if ch is None:
                        self.write_line_break()
                    start = end
            if ch is not None:
                breaks = ch in u'\n\x85\u2028\u2029'
                spaces = ch == u' '
            end += 1

    def write_literal(self, text, comment=None):
        # type: (Any, Any) -> None
        hints, _indent, _indicator = self.determine_block_hints(text)
        self.write_indicator(u'|' + hints, True)
        try:
            comment = comment[1][0]
            if comment:
                self.stream.write(comment)
        except (TypeError, IndexError):
            pass
        if _indicator == u'+':
            self.open_ended = True
        self.write_line_break()
        breaks = True
        start = end = 0
        while end <= len(text):
            ch = None
            if end < len(text):
                ch = text[end]
            if breaks:
                if ch is None or ch not in u'\n\x85\u2028\u2029':
                    for br in text[start:end]:
                        if br == u'\n':
                            self.write_line_break()
                        else:
                            self.write_line_break(br)
                    if ch is not None:
                        if self.root_context:
                            idnx = self.indent if self.indent is not None else 0
                            self.stream.write(u' ' * (_indent + idnx))
                        else:
                            self.write_indent()
                    start = end
            else:
                if ch is None or ch in u'\n\x85\u2028\u2029':
                    data = text[start:end]
                    if bool(self.encoding):
                        data = data.encode(self.encoding)
                    self.stream.write(data)
                    if ch is None:
                        self.write_line_break()
                    start = end
            if ch is not None:
                breaks = ch in u'\n\x85\u2028\u2029'
            end += 1

    def write_plain(self, text, split=True):
        # type: (Any, Any) -> None
        if self.root_context:
            if self.requested_indent is not None:
                self.write_line_break()
                if self.requested_indent != 0:
                    self.write_indent()
            else:
                self.open_ended = True
        if not text:
            return
        if not self.whitespace:
            data = u' '
            self.column += len(data)
            if self.encoding:
                data = data.encode(self.encoding)
            self.stream.write(data)
        self.whitespace = False
        self.indention = False
        spaces = False
        breaks = False
        start = end = 0
        while end <= len(text):
            ch = None
            if end < len(text):
                ch = text[end]
            if spaces:
                if ch != u' ':
                    if start + 1 == end and self.column > self.best_width and split:
                        self.write_indent()
                        self.whitespace = False
                        self.indention = False
                    else:
                        data = text[start:end]
                        self.column += len(data)
                        if self.encoding:
                            data = data.encode(self.encoding)
                        self.stream.write(data)
                    start = end
            elif breaks:
                if ch not in u'\n\x85\u2028\u2029':
                    if text[start] == u'\n':
                        self.write_line_break()
                    for br in text[start:end]:
                        if br == u'\n':
                            self.write_line_break()
                        else:
                            self.write_line_break(br)
                    self.write_indent()
                    self.whitespace = False
                    self.indention = False
                    start = end
            else:
                if ch is None or ch in u' \n\x85\u2028\u2029':
                    data = text[start:end]
                    self.column += len(data)
                    if self.encoding:
                        data = data.encode(self.encoding)
                    try:
                        self.stream.write(data)
                    except:  # NOQA
                        sys.stdout.write(repr(data) + '\n')
                        raise
                    start = end
            if ch is not None:
                spaces = ch == u' '
                breaks = ch in u'\n\x85\u2028\u2029'
            end += 1

    def write_comment(self, comment, pre=False):
        # type: (Any, bool) -> None
        value = comment.value
        # nprintf('{:02d} {:02d} {!r}'.format(self.column, comment.start_mark.column, value))
        if not pre and value[-1] == '\n':
            value = value[:-1]
        try:
            # get original column position
            col = comment.start_mark.column
            if comment.value and comment.value.startswith('\n'):
                # never inject extra spaces if the comment starts with a newline
                # and not a real comment (e.g. if you have an empty line following a key-value
                col = self.column
            elif col < self.column + 1:
                ValueError
        except ValueError:
            col = self.column + 1
        # nprint('post_comment', self.line, self.column, value)
        try:
            # at least one space if the current column >= the start column of the comment
            # but not at the start of a line
            nr_spaces = col - self.column
            if self.column and value.strip() and nr_spaces < 1 and value[0] != '\n':
                nr_spaces = 1
            value = ' ' * nr_spaces + value
            try:
                if bool(self.encoding):
                    value = value.encode(self.encoding)
            except UnicodeDecodeError:
                pass
            self.stream.write(value)
        except TypeError:
            raise
        if not pre:
            self.write_line_break()

    def write_pre_comment(self, event):
        # type: (Any) -> bool
        comments = event.comment[1]
        if comments is None:
            return False
        try:
            start_events = (MappingStartEvent, SequenceStartEvent)
            for comment in comments:
                if isinstance(event, start_events) and getattr(comment, 'pre_done', None):
                    continue
                if self.column != 0:
                    self.write_line_break()
                self.write_comment(comment, pre=True)
                if isinstance(event, start_events):
                    comment.pre_done = True
        except TypeError:
            sys.stdout.write('eventtt {} {}'.format(type(event), event))
            raise
        return True

    def write_post_comment(self, event):
        # type: (Any) -> bool
        if self.event.comment[0] is None:
            return False
        comment = event.comment[0]
        self.write_comment(comment)
        return True
