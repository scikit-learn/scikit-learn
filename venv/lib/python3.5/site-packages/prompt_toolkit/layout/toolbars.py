from __future__ import unicode_literals

from ..enums import IncrementalSearchDirection

from .processors import BeforeInput

from .lexers import SimpleLexer
from .dimension import LayoutDimension
from .controls import BufferControl, TokenListControl, UIControl, UIContent
from .containers import Window, ConditionalContainer
from .screen import Char
from .utils import token_list_len
from prompt_toolkit.enums import SEARCH_BUFFER, SYSTEM_BUFFER
from prompt_toolkit.filters import HasFocus, HasArg, HasCompletions, HasValidationError, HasSearch, Always, IsDone
from prompt_toolkit.token import Token

__all__ = (
    'TokenListToolbar',
    'ArgToolbar',
    'CompletionsToolbar',
    'SearchToolbar',
    'SystemToolbar',
    'ValidationToolbar',
)


class TokenListToolbar(ConditionalContainer):
    def __init__(self, get_tokens, filter=Always(), **kw):
        super(TokenListToolbar, self).__init__(
            content=Window(
                TokenListControl(get_tokens, **kw),
                height=LayoutDimension.exact(1)),
            filter=filter)


class SystemToolbarControl(BufferControl):
    def __init__(self):
        token = Token.Toolbar.System

        super(SystemToolbarControl, self).__init__(
            buffer_name=SYSTEM_BUFFER,
            default_char=Char(token=token),
            lexer=SimpleLexer(token=token.Text),
            input_processors=[BeforeInput.static('Shell command: ', token)],)


class SystemToolbar(ConditionalContainer):
    def __init__(self):
        super(SystemToolbar, self).__init__(
            content=Window(
                SystemToolbarControl(),
                height=LayoutDimension.exact(1)),
            filter=HasFocus(SYSTEM_BUFFER) & ~IsDone())


class ArgToolbarControl(TokenListControl):
    def __init__(self):
        def get_tokens(cli):
            arg = cli.input_processor.arg
            if arg == '-':
                arg = '-1'

            return [
                (Token.Toolbar.Arg, 'Repeat: '),
                (Token.Toolbar.Arg.Text, arg),
            ]

        super(ArgToolbarControl, self).__init__(get_tokens)


class ArgToolbar(ConditionalContainer):
    def __init__(self):
        super(ArgToolbar, self).__init__(
            content=Window(
                ArgToolbarControl(),
                height=LayoutDimension.exact(1)),
            filter=HasArg())


class SearchToolbarControl(BufferControl):
    """
    :param vi_mode: Display '/' and '?' instead of I-search.
    """
    def __init__(self, vi_mode=False):
        token = Token.Toolbar.Search

        def get_before_input(cli):
            if not cli.is_searching:
                text = ''
            elif cli.search_state.direction == IncrementalSearchDirection.BACKWARD:
                text = ('?' if vi_mode else 'I-search backward: ')
            else:
                text = ('/' if vi_mode else 'I-search: ')

            return [(token, text)]

        super(SearchToolbarControl, self).__init__(
            buffer_name=SEARCH_BUFFER,
            input_processors=[BeforeInput(get_before_input)],
            default_char=Char(token=token),
            lexer=SimpleLexer(token=token.Text))


class SearchToolbar(ConditionalContainer):
    def __init__(self, vi_mode=False):
        super(SearchToolbar, self).__init__(
            content=Window(
                SearchToolbarControl(vi_mode=vi_mode),
                height=LayoutDimension.exact(1)),
            filter=HasSearch() & ~IsDone())


class CompletionsToolbarControl(UIControl):
    token = Token.Toolbar.Completions

    def create_content(self, cli, width, height):
        complete_state = cli.current_buffer.complete_state
        if complete_state:
            completions = complete_state.current_completions
            index = complete_state.complete_index  # Can be None!

            # Width of the completions without the left/right arrows in the margins.
            content_width = width - 6

            # Booleans indicating whether we stripped from the left/right
            cut_left = False
            cut_right = False

            # Create Menu content.
            tokens = []

            for i, c in enumerate(completions):
                # When there is no more place for the next completion
                if token_list_len(tokens) + len(c.display) >= content_width:
                    # If the current one was not yet displayed, page to the next sequence.
                    if i <= (index or 0):
                        tokens = []
                        cut_left = True
                    # If the current one is visible, stop here.
                    else:
                        cut_right = True
                        break

                tokens.append((self.token.Completion.Current if i == index else self.token.Completion, c.display))
                tokens.append((self.token, ' '))

            # Extend/strip until the content width.
            tokens.append((self.token, ' ' * (content_width - token_list_len(tokens))))
            tokens = tokens[:content_width]

            # Return tokens
            all_tokens = [
                (self.token, ' '),
                (self.token.Arrow, '<' if cut_left else ' '),
                (self.token, ' '),
            ] + tokens + [
                (self.token, ' '),
                (self.token.Arrow, '>' if cut_right else ' '),
                (self.token, ' '),
            ]
        else:
            all_tokens = []

        def get_line(i):
            return all_tokens

        return UIContent(get_line=get_line, line_count=1)


class CompletionsToolbar(ConditionalContainer):
    def __init__(self, extra_filter=Always()):
        super(CompletionsToolbar, self).__init__(
            content=Window(
                CompletionsToolbarControl(),
                height=LayoutDimension.exact(1)),
            filter=HasCompletions() & ~IsDone() & extra_filter)


class ValidationToolbarControl(TokenListControl):
    def __init__(self, show_position=False):
        token = Token.Toolbar.Validation

        def get_tokens(cli):
            buffer = cli.current_buffer

            if buffer.validation_error:
                row, column = buffer.document.translate_index_to_position(
                    buffer.validation_error.cursor_position)

                if show_position:
                    text = '%s (line=%s column=%s)' % (
                        buffer.validation_error.message, row + 1, column + 1)
                else:
                    text = buffer.validation_error.message

                return [(token, text)]
            else:
                return []

        super(ValidationToolbarControl, self).__init__(get_tokens)


class ValidationToolbar(ConditionalContainer):
    def __init__(self, show_position=False):
        super(ValidationToolbar, self).__init__(
            content=Window(
                ValidationToolbarControl(show_position=show_position),
                height=LayoutDimension.exact(1)),
            filter=HasValidationError() & ~IsDone())
