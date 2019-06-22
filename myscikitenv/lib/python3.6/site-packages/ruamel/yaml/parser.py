# coding: utf-8

from __future__ import absolute_import

# The following YAML grammar is LL(1) and is parsed by a recursive descent
# parser.
#
# stream            ::= STREAM-START implicit_document? explicit_document*
#                                                                   STREAM-END
# implicit_document ::= block_node DOCUMENT-END*
# explicit_document ::= DIRECTIVE* DOCUMENT-START block_node? DOCUMENT-END*
# block_node_or_indentless_sequence ::=
#                       ALIAS
#                       | properties (block_content |
#                                                   indentless_block_sequence)?
#                       | block_content
#                       | indentless_block_sequence
# block_node        ::= ALIAS
#                       | properties block_content?
#                       | block_content
# flow_node         ::= ALIAS
#                       | properties flow_content?
#                       | flow_content
# properties        ::= TAG ANCHOR? | ANCHOR TAG?
# block_content     ::= block_collection | flow_collection | SCALAR
# flow_content      ::= flow_collection | SCALAR
# block_collection  ::= block_sequence | block_mapping
# flow_collection   ::= flow_sequence | flow_mapping
# block_sequence    ::= BLOCK-SEQUENCE-START (BLOCK-ENTRY block_node?)*
#                                                                   BLOCK-END
# indentless_sequence   ::= (BLOCK-ENTRY block_node?)+
# block_mapping     ::= BLOCK-MAPPING_START
#                       ((KEY block_node_or_indentless_sequence?)?
#                       (VALUE block_node_or_indentless_sequence?)?)*
#                       BLOCK-END
# flow_sequence     ::= FLOW-SEQUENCE-START
#                       (flow_sequence_entry FLOW-ENTRY)*
#                       flow_sequence_entry?
#                       FLOW-SEQUENCE-END
# flow_sequence_entry   ::= flow_node | KEY flow_node? (VALUE flow_node?)?
# flow_mapping      ::= FLOW-MAPPING-START
#                       (flow_mapping_entry FLOW-ENTRY)*
#                       flow_mapping_entry?
#                       FLOW-MAPPING-END
# flow_mapping_entry    ::= flow_node | KEY flow_node? (VALUE flow_node?)?
#
# FIRST sets:
#
# stream: { STREAM-START }
# explicit_document: { DIRECTIVE DOCUMENT-START }
# implicit_document: FIRST(block_node)
# block_node: { ALIAS TAG ANCHOR SCALAR BLOCK-SEQUENCE-START
#                  BLOCK-MAPPING-START FLOW-SEQUENCE-START FLOW-MAPPING-START }
# flow_node: { ALIAS ANCHOR TAG SCALAR FLOW-SEQUENCE-START FLOW-MAPPING-START }
# block_content: { BLOCK-SEQUENCE-START BLOCK-MAPPING-START
#                               FLOW-SEQUENCE-START FLOW-MAPPING-START SCALAR }
# flow_content: { FLOW-SEQUENCE-START FLOW-MAPPING-START SCALAR }
# block_collection: { BLOCK-SEQUENCE-START BLOCK-MAPPING-START }
# flow_collection: { FLOW-SEQUENCE-START FLOW-MAPPING-START }
# block_sequence: { BLOCK-SEQUENCE-START }
# block_mapping: { BLOCK-MAPPING-START }
# block_node_or_indentless_sequence: { ALIAS ANCHOR TAG SCALAR
#               BLOCK-SEQUENCE-START BLOCK-MAPPING-START FLOW-SEQUENCE-START
#               FLOW-MAPPING-START BLOCK-ENTRY }
# indentless_sequence: { ENTRY }
# flow_collection: { FLOW-SEQUENCE-START FLOW-MAPPING-START }
# flow_sequence: { FLOW-SEQUENCE-START }
# flow_mapping: { FLOW-MAPPING-START }
# flow_sequence_entry: { ALIAS ANCHOR TAG SCALAR FLOW-SEQUENCE-START
#                                                    FLOW-MAPPING-START KEY }
# flow_mapping_entry: { ALIAS ANCHOR TAG SCALAR FLOW-SEQUENCE-START
#                                                    FLOW-MAPPING-START KEY }

# need to have full path with import, as pkg_resources tries to load parser.py in __init__.py
# only to not do anything with the package afterwards
# and for Jython too


from ruamel.yaml.error import MarkedYAMLError
from ruamel.yaml.tokens import *  # NOQA
from ruamel.yaml.events import *  # NOQA
from ruamel.yaml.scanner import Scanner, RoundTripScanner, ScannerError  # NOQA
from ruamel.yaml.compat import utf8, nprint, nprintf  # NOQA

if False:  # MYPY
    from typing import Any, Dict, Optional, List  # NOQA

__all__ = ['Parser', 'RoundTripParser', 'ParserError']


class ParserError(MarkedYAMLError):
    pass


class Parser(object):
    # Since writing a recursive-descendant parser is a straightforward task, we
    # do not give many comments here.

    DEFAULT_TAGS = {u'!': u'!', u'!!': u'tag:yaml.org,2002:'}

    def __init__(self, loader):
        # type: (Any) -> None
        self.loader = loader
        if self.loader is not None and getattr(self.loader, '_parser', None) is None:
            self.loader._parser = self
        self.reset_parser()

    def reset_parser(self):
        # type: () -> None
        # Reset the state attributes (to clear self-references)
        self.current_event = None
        self.yaml_version = None
        self.tag_handles = {}  # type: Dict[Any, Any]
        self.states = []  # type: List[Any]
        self.marks = []  # type: List[Any]
        self.state = self.parse_stream_start  # type: Any

    def dispose(self):
        # type: () -> None
        self.reset_parser()

    @property
    def scanner(self):
        # type: () -> Any
        if hasattr(self.loader, 'typ'):
            return self.loader.scanner
        return self.loader._scanner

    @property
    def resolver(self):
        # type: () -> Any
        if hasattr(self.loader, 'typ'):
            return self.loader.resolver
        return self.loader._resolver

    def check_event(self, *choices):
        # type: (Any) -> bool
        # Check the type of the next event.
        if self.current_event is None:
            if self.state:
                self.current_event = self.state()
        if self.current_event is not None:
            if not choices:
                return True
            for choice in choices:
                if isinstance(self.current_event, choice):
                    return True
        return False

    def peek_event(self):
        # type: () -> Any
        # Get the next event.
        if self.current_event is None:
            if self.state:
                self.current_event = self.state()
        return self.current_event

    def get_event(self):
        # type: () -> Any
        # Get the next event and proceed further.
        if self.current_event is None:
            if self.state:
                self.current_event = self.state()
        value = self.current_event
        self.current_event = None
        return value

    # stream    ::= STREAM-START implicit_document? explicit_document*
    #                                                               STREAM-END
    # implicit_document ::= block_node DOCUMENT-END*
    # explicit_document ::= DIRECTIVE* DOCUMENT-START block_node? DOCUMENT-END*

    def parse_stream_start(self):
        # type: () -> Any
        # Parse the stream start.
        token = self.scanner.get_token()
        token.move_comment(self.scanner.peek_token())
        event = StreamStartEvent(token.start_mark, token.end_mark, encoding=token.encoding)

        # Prepare the next state.
        self.state = self.parse_implicit_document_start

        return event

    def parse_implicit_document_start(self):
        # type: () -> Any
        # Parse an implicit document.
        if not self.scanner.check_token(DirectiveToken, DocumentStartToken, StreamEndToken):
            self.tag_handles = self.DEFAULT_TAGS
            token = self.scanner.peek_token()
            start_mark = end_mark = token.start_mark
            event = DocumentStartEvent(start_mark, end_mark, explicit=False)

            # Prepare the next state.
            self.states.append(self.parse_document_end)
            self.state = self.parse_block_node

            return event

        else:
            return self.parse_document_start()

    def parse_document_start(self):
        # type: () -> Any
        # Parse any extra document end indicators.
        while self.scanner.check_token(DocumentEndToken):
            self.scanner.get_token()
        # Parse an explicit document.
        if not self.scanner.check_token(StreamEndToken):
            token = self.scanner.peek_token()
            start_mark = token.start_mark
            version, tags = self.process_directives()
            if not self.scanner.check_token(DocumentStartToken):
                raise ParserError(
                    None,
                    None,
                    "expected '<document start>', but found %r" % self.scanner.peek_token().id,
                    self.scanner.peek_token().start_mark,
                )
            token = self.scanner.get_token()
            end_mark = token.end_mark
            event = DocumentStartEvent(
                start_mark, end_mark, explicit=True, version=version, tags=tags
            )  # type: Any
            self.states.append(self.parse_document_end)
            self.state = self.parse_document_content
        else:
            # Parse the end of the stream.
            token = self.scanner.get_token()
            event = StreamEndEvent(token.start_mark, token.end_mark, comment=token.comment)
            assert not self.states
            assert not self.marks
            self.state = None
        return event

    def parse_document_end(self):
        # type: () -> Any
        # Parse the document end.
        token = self.scanner.peek_token()
        start_mark = end_mark = token.start_mark
        explicit = False
        if self.scanner.check_token(DocumentEndToken):
            token = self.scanner.get_token()
            end_mark = token.end_mark
            explicit = True
        event = DocumentEndEvent(start_mark, end_mark, explicit=explicit)

        # Prepare the next state.
        if self.resolver.processing_version == (1, 1):
            self.state = self.parse_document_start
        else:
            self.state = self.parse_implicit_document_start

        return event

    def parse_document_content(self):
        # type: () -> Any
        if self.scanner.check_token(
            DirectiveToken, DocumentStartToken, DocumentEndToken, StreamEndToken
        ):
            event = self.process_empty_scalar(self.scanner.peek_token().start_mark)
            self.state = self.states.pop()
            return event
        else:
            return self.parse_block_node()

    def process_directives(self):
        # type: () -> Any
        self.yaml_version = None
        self.tag_handles = {}
        while self.scanner.check_token(DirectiveToken):
            token = self.scanner.get_token()
            if token.name == u'YAML':
                if self.yaml_version is not None:
                    raise ParserError(
                        None, None, 'found duplicate YAML directive', token.start_mark
                    )
                major, minor = token.value
                if major != 1:
                    raise ParserError(
                        None,
                        None,
                        'found incompatible YAML document (version 1.* is ' 'required)',
                        token.start_mark,
                    )
                self.yaml_version = token.value
            elif token.name == u'TAG':
                handle, prefix = token.value
                if handle in self.tag_handles:
                    raise ParserError(
                        None, None, 'duplicate tag handle %r' % utf8(handle), token.start_mark
                    )
                self.tag_handles[handle] = prefix
        if bool(self.tag_handles):
            value = self.yaml_version, self.tag_handles.copy()  # type: Any
        else:
            value = self.yaml_version, None
        for key in self.DEFAULT_TAGS:
            if key not in self.tag_handles:
                self.tag_handles[key] = self.DEFAULT_TAGS[key]
        return value

    # block_node_or_indentless_sequence ::= ALIAS
    #               | properties (block_content | indentless_block_sequence)?
    #               | block_content
    #               | indentless_block_sequence
    # block_node    ::= ALIAS
    #                   | properties block_content?
    #                   | block_content
    # flow_node     ::= ALIAS
    #                   | properties flow_content?
    #                   | flow_content
    # properties    ::= TAG ANCHOR? | ANCHOR TAG?
    # block_content     ::= block_collection | flow_collection | SCALAR
    # flow_content      ::= flow_collection | SCALAR
    # block_collection  ::= block_sequence | block_mapping
    # flow_collection   ::= flow_sequence | flow_mapping

    def parse_block_node(self):
        # type: () -> Any
        return self.parse_node(block=True)

    def parse_flow_node(self):
        # type: () -> Any
        return self.parse_node()

    def parse_block_node_or_indentless_sequence(self):
        # type: () -> Any
        return self.parse_node(block=True, indentless_sequence=True)

    def transform_tag(self, handle, suffix):
        # type: (Any, Any) -> Any
        return self.tag_handles[handle] + suffix

    def parse_node(self, block=False, indentless_sequence=False):
        # type: (bool, bool) -> Any
        if self.scanner.check_token(AliasToken):
            token = self.scanner.get_token()
            event = AliasEvent(token.value, token.start_mark, token.end_mark)  # type: Any
            self.state = self.states.pop()
            return event

        anchor = None
        tag = None
        start_mark = end_mark = tag_mark = None
        if self.scanner.check_token(AnchorToken):
            token = self.scanner.get_token()
            start_mark = token.start_mark
            end_mark = token.end_mark
            anchor = token.value
            if self.scanner.check_token(TagToken):
                token = self.scanner.get_token()
                tag_mark = token.start_mark
                end_mark = token.end_mark
                tag = token.value
        elif self.scanner.check_token(TagToken):
            token = self.scanner.get_token()
            start_mark = tag_mark = token.start_mark
            end_mark = token.end_mark
            tag = token.value
            if self.scanner.check_token(AnchorToken):
                token = self.scanner.get_token()
                start_mark = tag_mark = token.start_mark
                end_mark = token.end_mark
                anchor = token.value
        if tag is not None:
            handle, suffix = tag
            if handle is not None:
                if handle not in self.tag_handles:
                    raise ParserError(
                        'while parsing a node',
                        start_mark,
                        'found undefined tag handle %r' % utf8(handle),
                        tag_mark,
                    )
                tag = self.transform_tag(handle, suffix)
            else:
                tag = suffix
        # if tag == u'!':
        #     raise ParserError("while parsing a node", start_mark,
        #             "found non-specific tag '!'", tag_mark,
        #      "Please check 'http://pyyaml.org/wiki/YAMLNonSpecificTag'
        #     and share your opinion.")
        if start_mark is None:
            start_mark = end_mark = self.scanner.peek_token().start_mark
        event = None
        implicit = tag is None or tag == u'!'
        if indentless_sequence and self.scanner.check_token(BlockEntryToken):
            comment = None
            pt = self.scanner.peek_token()
            if pt.comment and pt.comment[0]:
                comment = [pt.comment[0], []]
                pt.comment[0] = None
            end_mark = self.scanner.peek_token().end_mark
            event = SequenceStartEvent(
                anchor, tag, implicit, start_mark, end_mark, flow_style=False, comment=comment
            )
            self.state = self.parse_indentless_sequence_entry
            return event

        if self.scanner.check_token(ScalarToken):
            token = self.scanner.get_token()
            # self.scanner.peek_token_same_line_comment(token)
            end_mark = token.end_mark
            if (token.plain and tag is None) or tag == u'!':
                implicit = (True, False)
            elif tag is None:
                implicit = (False, True)
            else:
                implicit = (False, False)
            # nprint('se', token.value, token.comment)
            event = ScalarEvent(
                anchor,
                tag,
                implicit,
                token.value,
                start_mark,
                end_mark,
                style=token.style,
                comment=token.comment,
            )
            self.state = self.states.pop()
        elif self.scanner.check_token(FlowSequenceStartToken):
            pt = self.scanner.peek_token()
            end_mark = pt.end_mark
            event = SequenceStartEvent(
                anchor,
                tag,
                implicit,
                start_mark,
                end_mark,
                flow_style=True,
                comment=pt.comment,
            )
            self.state = self.parse_flow_sequence_first_entry
        elif self.scanner.check_token(FlowMappingStartToken):
            pt = self.scanner.peek_token()
            end_mark = pt.end_mark
            event = MappingStartEvent(
                anchor,
                tag,
                implicit,
                start_mark,
                end_mark,
                flow_style=True,
                comment=pt.comment,
            )
            self.state = self.parse_flow_mapping_first_key
        elif block and self.scanner.check_token(BlockSequenceStartToken):
            end_mark = self.scanner.peek_token().start_mark
            # should inserting the comment be dependent on the
            # indentation?
            pt = self.scanner.peek_token()
            comment = pt.comment
            # nprint('pt0', type(pt))
            if comment is None or comment[1] is None:
                comment = pt.split_comment()
            # nprint('pt1', comment)
            event = SequenceStartEvent(
                anchor, tag, implicit, start_mark, end_mark, flow_style=False, comment=comment
            )
            self.state = self.parse_block_sequence_first_entry
        elif block and self.scanner.check_token(BlockMappingStartToken):
            end_mark = self.scanner.peek_token().start_mark
            comment = self.scanner.peek_token().comment
            event = MappingStartEvent(
                anchor, tag, implicit, start_mark, end_mark, flow_style=False, comment=comment
            )
            self.state = self.parse_block_mapping_first_key
        elif anchor is not None or tag is not None:
            # Empty scalars are allowed even if a tag or an anchor is
            # specified.
            event = ScalarEvent(anchor, tag, (implicit, False), "", start_mark, end_mark)
            self.state = self.states.pop()
        else:
            if block:
                node = 'block'
            else:
                node = 'flow'
            token = self.scanner.peek_token()
            raise ParserError(
                'while parsing a %s node' % node,
                start_mark,
                'expected the node content, but found %r' % token.id,
                token.start_mark,
            )
        return event

    # block_sequence ::= BLOCK-SEQUENCE-START (BLOCK-ENTRY block_node?)*
    #                                                               BLOCK-END

    def parse_block_sequence_first_entry(self):
        # type: () -> Any
        token = self.scanner.get_token()
        # move any comment from start token
        # token.move_comment(self.scanner.peek_token())
        self.marks.append(token.start_mark)
        return self.parse_block_sequence_entry()

    def parse_block_sequence_entry(self):
        # type: () -> Any
        if self.scanner.check_token(BlockEntryToken):
            token = self.scanner.get_token()
            token.move_comment(self.scanner.peek_token())
            if not self.scanner.check_token(BlockEntryToken, BlockEndToken):
                self.states.append(self.parse_block_sequence_entry)
                return self.parse_block_node()
            else:
                self.state = self.parse_block_sequence_entry
                return self.process_empty_scalar(token.end_mark)
        if not self.scanner.check_token(BlockEndToken):
            token = self.scanner.peek_token()
            raise ParserError(
                'while parsing a block collection',
                self.marks[-1],
                'expected <block end>, but found %r' % token.id,
                token.start_mark,
            )
        token = self.scanner.get_token()  # BlockEndToken
        event = SequenceEndEvent(token.start_mark, token.end_mark, comment=token.comment)
        self.state = self.states.pop()
        self.marks.pop()
        return event

    # indentless_sequence ::= (BLOCK-ENTRY block_node?)+

    # indentless_sequence?
    # sequence:
    # - entry
    #  - nested

    def parse_indentless_sequence_entry(self):
        # type: () -> Any
        if self.scanner.check_token(BlockEntryToken):
            token = self.scanner.get_token()
            token.move_comment(self.scanner.peek_token())
            if not self.scanner.check_token(
                BlockEntryToken, KeyToken, ValueToken, BlockEndToken
            ):
                self.states.append(self.parse_indentless_sequence_entry)
                return self.parse_block_node()
            else:
                self.state = self.parse_indentless_sequence_entry
                return self.process_empty_scalar(token.end_mark)
        token = self.scanner.peek_token()
        event = SequenceEndEvent(token.start_mark, token.start_mark, comment=token.comment)
        self.state = self.states.pop()
        return event

    # block_mapping     ::= BLOCK-MAPPING_START
    #                       ((KEY block_node_or_indentless_sequence?)?
    #                       (VALUE block_node_or_indentless_sequence?)?)*
    #                       BLOCK-END

    def parse_block_mapping_first_key(self):
        # type: () -> Any
        token = self.scanner.get_token()
        self.marks.append(token.start_mark)
        return self.parse_block_mapping_key()

    def parse_block_mapping_key(self):
        # type: () -> Any
        if self.scanner.check_token(KeyToken):
            token = self.scanner.get_token()
            token.move_comment(self.scanner.peek_token())
            if not self.scanner.check_token(KeyToken, ValueToken, BlockEndToken):
                self.states.append(self.parse_block_mapping_value)
                return self.parse_block_node_or_indentless_sequence()
            else:
                self.state = self.parse_block_mapping_value
                return self.process_empty_scalar(token.end_mark)
        if self.resolver.processing_version > (1, 1) and self.scanner.check_token(ValueToken):
            self.state = self.parse_block_mapping_value
            return self.process_empty_scalar(self.scanner.peek_token().start_mark)
        if not self.scanner.check_token(BlockEndToken):
            token = self.scanner.peek_token()
            raise ParserError(
                'while parsing a block mapping',
                self.marks[-1],
                'expected <block end>, but found %r' % token.id,
                token.start_mark,
            )
        token = self.scanner.get_token()
        token.move_comment(self.scanner.peek_token())
        event = MappingEndEvent(token.start_mark, token.end_mark, comment=token.comment)
        self.state = self.states.pop()
        self.marks.pop()
        return event

    def parse_block_mapping_value(self):
        # type: () -> Any
        if self.scanner.check_token(ValueToken):
            token = self.scanner.get_token()
            # value token might have post comment move it to e.g. block
            if self.scanner.check_token(ValueToken):
                token.move_comment(self.scanner.peek_token())
            else:
                if not self.scanner.check_token(KeyToken):
                    token.move_comment(self.scanner.peek_token(), empty=True)
                # else: empty value for this key cannot move token.comment
            if not self.scanner.check_token(KeyToken, ValueToken, BlockEndToken):
                self.states.append(self.parse_block_mapping_key)
                return self.parse_block_node_or_indentless_sequence()
            else:
                self.state = self.parse_block_mapping_key
                comment = token.comment
                if comment is None:
                    token = self.scanner.peek_token()
                    comment = token.comment
                    if comment:
                        token._comment = [None, comment[1]]
                        comment = [comment[0], None]
                return self.process_empty_scalar(token.end_mark, comment=comment)
        else:
            self.state = self.parse_block_mapping_key
            token = self.scanner.peek_token()
            return self.process_empty_scalar(token.start_mark)

    # flow_sequence     ::= FLOW-SEQUENCE-START
    #                       (flow_sequence_entry FLOW-ENTRY)*
    #                       flow_sequence_entry?
    #                       FLOW-SEQUENCE-END
    # flow_sequence_entry   ::= flow_node | KEY flow_node? (VALUE flow_node?)?
    #
    # Note that while production rules for both flow_sequence_entry and
    # flow_mapping_entry are equal, their interpretations are different.
    # For `flow_sequence_entry`, the part `KEY flow_node? (VALUE flow_node?)?`
    # generate an inline mapping (set syntax).

    def parse_flow_sequence_first_entry(self):
        # type: () -> Any
        token = self.scanner.get_token()
        self.marks.append(token.start_mark)
        return self.parse_flow_sequence_entry(first=True)

    def parse_flow_sequence_entry(self, first=False):
        # type: (bool) -> Any
        if not self.scanner.check_token(FlowSequenceEndToken):
            if not first:
                if self.scanner.check_token(FlowEntryToken):
                    self.scanner.get_token()
                else:
                    token = self.scanner.peek_token()
                    raise ParserError(
                        'while parsing a flow sequence',
                        self.marks[-1],
                        "expected ',' or ']', but got %r" % token.id,
                        token.start_mark,
                    )

            if self.scanner.check_token(KeyToken):
                token = self.scanner.peek_token()
                event = MappingStartEvent(
                    None, None, True, token.start_mark, token.end_mark, flow_style=True
                )  # type: Any
                self.state = self.parse_flow_sequence_entry_mapping_key
                return event
            elif not self.scanner.check_token(FlowSequenceEndToken):
                self.states.append(self.parse_flow_sequence_entry)
                return self.parse_flow_node()
        token = self.scanner.get_token()
        event = SequenceEndEvent(token.start_mark, token.end_mark, comment=token.comment)
        self.state = self.states.pop()
        self.marks.pop()
        return event

    def parse_flow_sequence_entry_mapping_key(self):
        # type: () -> Any
        token = self.scanner.get_token()
        if not self.scanner.check_token(ValueToken, FlowEntryToken, FlowSequenceEndToken):
            self.states.append(self.parse_flow_sequence_entry_mapping_value)
            return self.parse_flow_node()
        else:
            self.state = self.parse_flow_sequence_entry_mapping_value
            return self.process_empty_scalar(token.end_mark)

    def parse_flow_sequence_entry_mapping_value(self):
        # type: () -> Any
        if self.scanner.check_token(ValueToken):
            token = self.scanner.get_token()
            if not self.scanner.check_token(FlowEntryToken, FlowSequenceEndToken):
                self.states.append(self.parse_flow_sequence_entry_mapping_end)
                return self.parse_flow_node()
            else:
                self.state = self.parse_flow_sequence_entry_mapping_end
                return self.process_empty_scalar(token.end_mark)
        else:
            self.state = self.parse_flow_sequence_entry_mapping_end
            token = self.scanner.peek_token()
            return self.process_empty_scalar(token.start_mark)

    def parse_flow_sequence_entry_mapping_end(self):
        # type: () -> Any
        self.state = self.parse_flow_sequence_entry
        token = self.scanner.peek_token()
        return MappingEndEvent(token.start_mark, token.start_mark)

    # flow_mapping  ::= FLOW-MAPPING-START
    #                   (flow_mapping_entry FLOW-ENTRY)*
    #                   flow_mapping_entry?
    #                   FLOW-MAPPING-END
    # flow_mapping_entry    ::= flow_node | KEY flow_node? (VALUE flow_node?)?

    def parse_flow_mapping_first_key(self):
        # type: () -> Any
        token = self.scanner.get_token()
        self.marks.append(token.start_mark)
        return self.parse_flow_mapping_key(first=True)

    def parse_flow_mapping_key(self, first=False):
        # type: (Any) -> Any
        if not self.scanner.check_token(FlowMappingEndToken):
            if not first:
                if self.scanner.check_token(FlowEntryToken):
                    self.scanner.get_token()
                else:
                    token = self.scanner.peek_token()
                    raise ParserError(
                        'while parsing a flow mapping',
                        self.marks[-1],
                        "expected ',' or '}', but got %r" % token.id,
                        token.start_mark,
                    )
            if self.scanner.check_token(KeyToken):
                token = self.scanner.get_token()
                if not self.scanner.check_token(
                    ValueToken, FlowEntryToken, FlowMappingEndToken
                ):
                    self.states.append(self.parse_flow_mapping_value)
                    return self.parse_flow_node()
                else:
                    self.state = self.parse_flow_mapping_value
                    return self.process_empty_scalar(token.end_mark)
            elif self.resolver.processing_version > (1, 1) and self.scanner.check_token(
                ValueToken
            ):
                self.state = self.parse_flow_mapping_value
                return self.process_empty_scalar(self.scanner.peek_token().end_mark)
            elif not self.scanner.check_token(FlowMappingEndToken):
                self.states.append(self.parse_flow_mapping_empty_value)
                return self.parse_flow_node()
        token = self.scanner.get_token()
        event = MappingEndEvent(token.start_mark, token.end_mark, comment=token.comment)
        self.state = self.states.pop()
        self.marks.pop()
        return event

    def parse_flow_mapping_value(self):
        # type: () -> Any
        if self.scanner.check_token(ValueToken):
            token = self.scanner.get_token()
            if not self.scanner.check_token(FlowEntryToken, FlowMappingEndToken):
                self.states.append(self.parse_flow_mapping_key)
                return self.parse_flow_node()
            else:
                self.state = self.parse_flow_mapping_key
                return self.process_empty_scalar(token.end_mark)
        else:
            self.state = self.parse_flow_mapping_key
            token = self.scanner.peek_token()
            return self.process_empty_scalar(token.start_mark)

    def parse_flow_mapping_empty_value(self):
        # type: () -> Any
        self.state = self.parse_flow_mapping_key
        return self.process_empty_scalar(self.scanner.peek_token().start_mark)

    def process_empty_scalar(self, mark, comment=None):
        # type: (Any, Any) -> Any
        return ScalarEvent(None, None, (True, False), "", mark, mark, comment=comment)


class RoundTripParser(Parser):
    """roundtrip is a safe loader, that wants to see the unmangled tag"""

    def transform_tag(self, handle, suffix):
        # type: (Any, Any) -> Any
        # return self.tag_handles[handle]+suffix
        if handle == '!!' and suffix in (
            u'null',
            u'bool',
            u'int',
            u'float',
            u'binary',
            u'timestamp',
            u'omap',
            u'pairs',
            u'set',
            u'str',
            u'seq',
            u'map',
        ):
            return Parser.transform_tag(self, handle, suffix)
        return handle + suffix
