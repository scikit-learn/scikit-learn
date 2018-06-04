from __future__ import unicode_literals

from six import text_type

from prompt_toolkit.enums import IncrementalSearchDirection, SEARCH_BUFFER
from prompt_toolkit.token import Token

from .utils import token_list_len
from .processors import Processor, Transformation

__all__ = (
    'DefaultPrompt',
)


class DefaultPrompt(Processor):
    """
    Default prompt. This one shows the 'arg' and reverse search like
    Bash/readline normally do.

    There are two ways to instantiate a ``DefaultPrompt``. For a prompt
    with a static message, do for instance::

        prompt = DefaultPrompt.from_message('prompt> ')

    For a dynamic prompt, generated from a token list function::

        def get_tokens(cli):
            return [(Token.A, 'text'), (Token.B, 'text2')]

        prompt = DefaultPrompt(get_tokens)
    """
    def __init__(self, get_tokens):
        assert callable(get_tokens)
        self.get_tokens = get_tokens

    @classmethod
    def from_message(cls, message='> '):
        """
        Create a default prompt with a static message text.
        """
        assert isinstance(message, text_type)

        def get_message_tokens(cli):
            return [(Token.Prompt, message)]
        return cls(get_message_tokens)

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        # Get text before cursor.
        if cli.is_searching:
            before = _get_isearch_tokens(cli)

        elif cli.input_processor.arg is not None:
            before = _get_arg_tokens(cli)

        else:
            before = self.get_tokens(cli)

        # Insert before buffer text.
        shift_position = token_list_len(before)

        # Only show the prompt before the first line. For the following lines,
        # only indent using spaces.
        if lineno != 0:
            before = [(Token.Prompt, ' ' * shift_position)]

        return Transformation(
                tokens=before + tokens,
                source_to_display=lambda i: i + shift_position,
                display_to_source=lambda i: i - shift_position)

    def has_focus(self, cli):
        # Obtain focus when the CLI is searching.

        # Usually, when using this `DefaultPrompt`, we don't have a
        # `BufferControl` instance that displays the content of the search
        # buffer. Instead the search text is displayed before the current text.
        # So, we can still show the cursor here, while it's actually not this
        # buffer that's focussed.
        return cli.is_searching


def _get_isearch_tokens(cli):
    def before():
        if cli.search_state.direction == IncrementalSearchDirection.BACKWARD:
            text = 'reverse-i-search'
        else:
            text = 'i-search'

        return [(Token.Prompt.Search, '(%s)`' % text)]

    def text():
        return [(Token.Prompt.Search.Text, cli.buffers[SEARCH_BUFFER].text)]

    def after():
        return [(Token.Prompt.Search, '`: ')]

    return before() + text() + after()


def _get_arg_tokens(cli):
    """
    Tokens for the arg-prompt.
    """
    arg = cli.input_processor.arg

    return [
        (Token.Prompt.Arg, '(arg: '),
        (Token.Prompt.Arg.Text, str(arg)),
        (Token.Prompt.Arg, ') '),
    ]
