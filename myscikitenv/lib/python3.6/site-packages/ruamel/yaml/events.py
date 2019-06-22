# coding: utf-8

# Abstract classes.

if False:  # MYPY
    from typing import Any, Dict, Optional, List  # NOQA


def CommentCheck():
    # type: () -> None
    pass


class Event(object):
    __slots__ = 'start_mark', 'end_mark', 'comment'

    def __init__(self, start_mark=None, end_mark=None, comment=CommentCheck):
        # type: (Any, Any, Any) -> None
        self.start_mark = start_mark
        self.end_mark = end_mark
        # assert comment is not CommentCheck
        if comment is CommentCheck:
            comment = None
        self.comment = comment

    def __repr__(self):
        # type: () -> Any
        attributes = [
            key
            for key in ['anchor', 'tag', 'implicit', 'value', 'flow_style', 'style']
            if hasattr(self, key)
        ]
        arguments = ', '.join(['%s=%r' % (key, getattr(self, key)) for key in attributes])
        if self.comment not in [None, CommentCheck]:
            arguments += ', comment={!r}'.format(self.comment)
        return '%s(%s)' % (self.__class__.__name__, arguments)


class NodeEvent(Event):
    __slots__ = ('anchor',)

    def __init__(self, anchor, start_mark=None, end_mark=None, comment=None):
        # type: (Any, Any, Any, Any) -> None
        Event.__init__(self, start_mark, end_mark, comment)
        self.anchor = anchor


class CollectionStartEvent(NodeEvent):
    __slots__ = 'tag', 'implicit', 'flow_style', 'nr_items'

    def __init__(
        self,
        anchor,
        tag,
        implicit,
        start_mark=None,
        end_mark=None,
        flow_style=None,
        comment=None,
        nr_items=None,
    ):
        # type: (Any, Any, Any, Any, Any, Any, Any, Optional[int]) -> None
        NodeEvent.__init__(self, anchor, start_mark, end_mark, comment)
        self.tag = tag
        self.implicit = implicit
        self.flow_style = flow_style
        self.nr_items = nr_items


class CollectionEndEvent(Event):
    __slots__ = ()


# Implementations.


class StreamStartEvent(Event):
    __slots__ = ('encoding',)

    def __init__(self, start_mark=None, end_mark=None, encoding=None, comment=None):
        # type: (Any, Any, Any, Any) -> None
        Event.__init__(self, start_mark, end_mark, comment)
        self.encoding = encoding


class StreamEndEvent(Event):
    __slots__ = ()


class DocumentStartEvent(Event):
    __slots__ = 'explicit', 'version', 'tags'

    def __init__(
        self,
        start_mark=None,
        end_mark=None,
        explicit=None,
        version=None,
        tags=None,
        comment=None,
    ):
        # type: (Any, Any, Any, Any, Any, Any) -> None
        Event.__init__(self, start_mark, end_mark, comment)
        self.explicit = explicit
        self.version = version
        self.tags = tags


class DocumentEndEvent(Event):
    __slots__ = ('explicit',)

    def __init__(self, start_mark=None, end_mark=None, explicit=None, comment=None):
        # type: (Any, Any, Any, Any) -> None
        Event.__init__(self, start_mark, end_mark, comment)
        self.explicit = explicit


class AliasEvent(NodeEvent):
    __slots__ = ()


class ScalarEvent(NodeEvent):
    __slots__ = 'tag', 'implicit', 'value', 'style'

    def __init__(
        self,
        anchor,
        tag,
        implicit,
        value,
        start_mark=None,
        end_mark=None,
        style=None,
        comment=None,
    ):
        # type: (Any, Any, Any, Any, Any, Any, Any, Any) -> None
        NodeEvent.__init__(self, anchor, start_mark, end_mark, comment)
        self.tag = tag
        self.implicit = implicit
        self.value = value
        self.style = style


class SequenceStartEvent(CollectionStartEvent):
    __slots__ = ()


class SequenceEndEvent(CollectionEndEvent):
    __slots__ = ()


class MappingStartEvent(CollectionStartEvent):
    __slots__ = ()


class MappingEndEvent(CollectionEndEvent):
    __slots__ = ()
