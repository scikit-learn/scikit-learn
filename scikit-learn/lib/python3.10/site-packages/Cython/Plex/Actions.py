"""
Python Lexical Analyser

Actions for use in token specifications
"""

class Action:
    def perform(self, token_stream, text):
        pass  # abstract

    def __copy__(self):
        return self  # immutable, no need to copy

    def __deepcopy__(self, memo):
        return self  # immutable, no need to copy


class Return(Action):
    """
    Internal Plex action which causes |value| to
    be returned as the value of the associated token
    """

    def __init__(self, value):
        self.value = value

    def perform(self, token_stream, text):
        return self.value

    def __repr__(self):
        return "Return(%r)" % self.value


class Call(Action):
    """
    Internal Plex action which causes a function to be called.
    """

    def __init__(self, function):
        self.function = function

    def perform(self, token_stream, text):
        return self.function(token_stream, text)

    def __repr__(self):
        return "Call(%s)" % self.function.__name__


class Method(Action):
    """
    Plex action that calls a specific method on the token stream,
    passing the matched text and any provided constant keyword arguments.
    """

    def __init__(self, name, **kwargs):
        self.name = name
        self.kwargs = kwargs or None

    def perform(self, token_stream, text):
        method = getattr(token_stream, self.name)
        # self.kwargs is almost always unused => avoid call overhead
        return method(text, **self.kwargs) if self.kwargs is not None else method(text)

    def __repr__(self):
        kwargs = (
            ', '.join(sorted(['%s=%r' % item for item in self.kwargs.items()]))
            if self.kwargs is not None else '')
        return "Method(%s%s%s)" % (self.name, ', ' if kwargs else '', kwargs)


class Begin(Action):
    """
    Begin(state_name) is a Plex action which causes the Scanner to
    enter the state |state_name|. See the docstring of Plex.Lexicon
    for more information.
    """

    def __init__(self, state_name):
        self.state_name = state_name

    def perform(self, token_stream, text):
        token_stream.begin(self.state_name)

    def __repr__(self):
        return "Begin(%s)" % self.state_name


class Ignore(Action):
    """
    IGNORE is a Plex action which causes its associated token
    to be ignored. See the docstring of Plex.Lexicon  for more
    information.
    """

    def perform(self, token_stream, text):
        return None

    def __repr__(self):
        return "IGNORE"


IGNORE = Ignore()


class Text(Action):
    """
    TEXT is a Plex action which causes the text of a token to
    be returned as the value of the token. See the docstring of
    Plex.Lexicon  for more information.
    """

    def perform(self, token_stream, text):
        return text

    def __repr__(self):
        return "TEXT"


TEXT = Text()
