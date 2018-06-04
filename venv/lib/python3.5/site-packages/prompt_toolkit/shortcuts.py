"""
Shortcuts for retrieving input from the user.

If you are using this library for retrieving some input from the user (as a
pure Python replacement for GNU readline), probably for 90% of the use cases,
the :func:`.prompt` function is all you need. It's the easiest shortcut which
does a lot of the underlying work like creating a
:class:`~prompt_toolkit.interface.CommandLineInterface` instance for you.

When is this not sufficient:
    - When you want to have more complicated layouts (maybe with sidebars or
      multiple toolbars. Or visibility of certain user interface controls
      according to some conditions.)
    - When you wish to have multiple input buffers. (If you would create an
      editor like a Vi clone.)
    - Something else that requires more customization than what is possible
      with the parameters of `prompt`.

In that case, study the code in this file and build your own
`CommandLineInterface` instance. It's not too complicated.
"""
from __future__ import unicode_literals

from .buffer import Buffer, AcceptAction
from .document import Document
from .enums import DEFAULT_BUFFER, SEARCH_BUFFER, EditingMode
from .filters import IsDone, HasFocus, RendererHeightIsKnown, to_simple_filter, to_cli_filter, Condition
from .history import InMemoryHistory
from .interface import CommandLineInterface, Application, AbortAction
from .key_binding.defaults import load_key_bindings_for_prompt
from .key_binding.registry import Registry
from .keys import Keys
from .layout import Window, HSplit, FloatContainer, Float
from .layout.containers import ConditionalContainer
from .layout.controls import BufferControl, TokenListControl
from .layout.dimension import LayoutDimension
from .layout.lexers import PygmentsLexer
from .layout.margins import PromptMargin, ConditionalMargin
from .layout.menus import CompletionsMenu, MultiColumnCompletionsMenu
from .layout.processors import PasswordProcessor, ConditionalProcessor, AppendAutoSuggestion, HighlightSearchProcessor, HighlightSelectionProcessor, DisplayMultipleCursors
from .layout.prompt import DefaultPrompt
from .layout.screen import Char
from .layout.toolbars import ValidationToolbar, SystemToolbar, ArgToolbar, SearchToolbar
from .layout.utils import explode_tokens
from .renderer import print_tokens as renderer_print_tokens
from .styles import DEFAULT_STYLE, Style, style_from_dict
from .token import Token
from .utils import is_conemu_ansi, is_windows, DummyContext

from six import text_type, exec_, PY2

import os
import sys
import textwrap
import threading
import time

try:
    from pygments.lexer import Lexer as pygments_Lexer
    from pygments.style import Style as pygments_Style
except ImportError:
    pygments_Lexer = None
    pygments_Style = None

if is_windows():
    from .terminal.win32_output import Win32Output
    from .terminal.conemu_output import ConEmuOutput
else:
    from .terminal.vt100_output import Vt100_Output


__all__ = (
    'create_eventloop',
    'create_output',
    'create_prompt_layout',
    'create_prompt_application',
    'prompt',
    'prompt_async',
    'create_confirm_application',
    'run_application',
    'confirm',
    'print_tokens',
    'clear',
)


def create_eventloop(inputhook=None, recognize_win32_paste=True):
    """
    Create and return an
    :class:`~prompt_toolkit.eventloop.base.EventLoop` instance for a
    :class:`~prompt_toolkit.interface.CommandLineInterface`.
    """
    if is_windows():
        from prompt_toolkit.eventloop.win32 import Win32EventLoop as Loop
        return Loop(inputhook=inputhook, recognize_paste=recognize_win32_paste)
    else:
        from prompt_toolkit.eventloop.posix import PosixEventLoop as Loop
        return Loop(inputhook=inputhook)


def create_output(stdout=None, true_color=False, ansi_colors_only=None):
    """
    Return an :class:`~prompt_toolkit.output.Output` instance for the command
    line.

    :param true_color: When True, use 24bit colors instead of 256 colors.
        (`bool` or :class:`~prompt_toolkit.filters.SimpleFilter`.)
    :param ansi_colors_only: When True, restrict to 16 ANSI colors only.
        (`bool` or :class:`~prompt_toolkit.filters.SimpleFilter`.)
    """
    stdout = stdout or sys.__stdout__
    true_color = to_simple_filter(true_color)

    if is_windows():
        if is_conemu_ansi():
            return ConEmuOutput(stdout)
        else:
            return Win32Output(stdout)
    else:
        term = os.environ.get('TERM', '')
        if PY2:
            term = term.decode('utf-8')

        return Vt100_Output.from_pty(
            stdout, true_color=true_color,
            ansi_colors_only=ansi_colors_only, term=term)


def create_asyncio_eventloop(loop=None):
    """
    Returns an asyncio :class:`~prompt_toolkit.eventloop.EventLoop` instance
    for usage in a :class:`~prompt_toolkit.interface.CommandLineInterface`. It
    is a wrapper around an asyncio loop.

    :param loop: The asyncio eventloop (or `None` if the default asyncioloop
                 should be used.)
    """
    # Inline import, to make sure the rest doesn't break on Python 2. (Where
    # asyncio is not available.)
    if is_windows():
        from prompt_toolkit.eventloop.asyncio_win32 import Win32AsyncioEventLoop as AsyncioEventLoop
    else:
        from prompt_toolkit.eventloop.asyncio_posix import PosixAsyncioEventLoop as AsyncioEventLoop

    return AsyncioEventLoop(loop)


def _split_multiline_prompt(get_prompt_tokens):
    """
    Take a `get_prompt_tokens` function and return three new functions instead.
    One that tells whether this prompt consists of multiple lines; one that
    returns the tokens to be shown on the lines above the input; and another
    one with the tokens to be shown at the first line of the input.
    """
    def has_before_tokens(cli):
        for token, char in get_prompt_tokens(cli):
            if '\n' in char:
                return True
        return False

    def before(cli):
        result = []
        found_nl = False
        for token, char in reversed(explode_tokens(get_prompt_tokens(cli))):
            if found_nl:
                result.insert(0, (token, char))
            elif char == '\n':
                found_nl = True
        return result

    def first_input_line(cli):
        result = []
        for token, char in reversed(explode_tokens(get_prompt_tokens(cli))):
            if char == '\n':
                break
            else:
                result.insert(0, (token, char))
        return result

    return has_before_tokens, before, first_input_line


class _RPrompt(Window):
    " The prompt that is displayed on the right side of the Window. "
    def __init__(self, get_tokens=None):
        get_tokens = get_tokens or (lambda cli: [])

        super(_RPrompt, self).__init__(
            TokenListControl(get_tokens, align_right=True))


def create_prompt_layout(message='', lexer=None, is_password=False,
                         reserve_space_for_menu=8,
                         get_prompt_tokens=None, get_continuation_tokens=None,
                         get_rprompt_tokens=None,
                         get_bottom_toolbar_tokens=None,
                         display_completions_in_columns=False,
                         extra_input_processors=None, multiline=False,
                         wrap_lines=True):
    """
    Create a :class:`.Container` instance for a prompt.

    :param message: Text to be used as prompt.
    :param lexer: :class:`~prompt_toolkit.layout.lexers.Lexer` to be used for
        the highlighting.
    :param is_password: `bool` or :class:`~prompt_toolkit.filters.CLIFilter`.
        When True, display input as '*'.
    :param reserve_space_for_menu: Space to be reserved for the menu. When >0,
        make sure that a minimal height is allocated in the terminal, in order
        to display the completion menu.
    :param get_prompt_tokens: An optional callable that returns the tokens to be
        shown in the menu. (To be used instead of a `message`.)
    :param get_continuation_tokens: An optional callable that takes a
        CommandLineInterface and width as input and returns a list of (Token,
        text) tuples to be used for the continuation.
    :param get_bottom_toolbar_tokens: An optional callable that returns the
        tokens for a toolbar at the bottom.
    :param display_completions_in_columns: `bool` or
        :class:`~prompt_toolkit.filters.CLIFilter`. Display the completions in
        multiple columns.
    :param multiline: `bool` or :class:`~prompt_toolkit.filters.CLIFilter`.
        When True, prefer a layout that is more adapted for multiline input.
        Text after newlines is automatically indented, and search/arg input is
        shown below the input, instead of replacing the prompt.
    :param wrap_lines: `bool` or :class:`~prompt_toolkit.filters.CLIFilter`.
        When True (the default), automatically wrap long lines instead of
        scrolling horizontally.
    """
    assert isinstance(message, text_type), 'Please provide a unicode string.'
    assert get_bottom_toolbar_tokens is None or callable(get_bottom_toolbar_tokens)
    assert get_prompt_tokens is None or callable(get_prompt_tokens)
    assert get_rprompt_tokens is None or callable(get_rprompt_tokens)
    assert not (message and get_prompt_tokens)

    display_completions_in_columns = to_cli_filter(display_completions_in_columns)
    multiline = to_cli_filter(multiline)

    if get_prompt_tokens is None:
        get_prompt_tokens = lambda _: [(Token.Prompt, message)]

    has_before_tokens, get_prompt_tokens_1, get_prompt_tokens_2 = \
        _split_multiline_prompt(get_prompt_tokens)

    # `lexer` is supposed to be a `Lexer` instance. But if a Pygments lexer
    # class is given, turn it into a PygmentsLexer. (Important for
    # backwards-compatibility.)
    try:
        if pygments_Lexer and issubclass(lexer, pygments_Lexer):
            lexer = PygmentsLexer(lexer, sync_from_start=True)
    except TypeError: # Happens when lexer is `None` or an instance of something else.
        pass

    # Create processors list.
    input_processors = [
        ConditionalProcessor(
            # By default, only highlight search when the search
            # input has the focus. (Note that this doesn't mean
            # there is no search: the Vi 'n' binding for instance
            # still allows to jump to the next match in
            # navigation mode.)
            HighlightSearchProcessor(preview_search=True),
            HasFocus(SEARCH_BUFFER)),
        HighlightSelectionProcessor(),
        ConditionalProcessor(AppendAutoSuggestion(), HasFocus(DEFAULT_BUFFER) & ~IsDone()),
        ConditionalProcessor(PasswordProcessor(), is_password),
        DisplayMultipleCursors(DEFAULT_BUFFER),
    ]

    if extra_input_processors:
        input_processors.extend(extra_input_processors)

    # Show the prompt before the input (using the DefaultPrompt processor.
    # This also replaces it with reverse-i-search and 'arg' when required.
    # (Only for single line mode.)
    # (DefaultPrompt should always be at the end of the processors.)
    input_processors.append(ConditionalProcessor(
        DefaultPrompt(get_prompt_tokens_2), ~multiline))

    # Create bottom toolbar.
    if get_bottom_toolbar_tokens:
        toolbars = [ConditionalContainer(
            Window(TokenListControl(get_bottom_toolbar_tokens,
                                    default_char=Char(' ', Token.Toolbar)),
                                    height=LayoutDimension.exact(1)),
            filter=~IsDone() & RendererHeightIsKnown())]
    else:
        toolbars = []

    def get_height(cli):
        # If there is an autocompletion menu to be shown, make sure that our
        # layout has at least a minimal height in order to display it.
        if reserve_space_for_menu and not cli.is_done:
            buff = cli.current_buffer

            # Reserve the space, either when there are completions, or when
            # `complete_while_typing` is true and we expect completions very
            # soon.
            if buff.complete_while_typing() or buff.complete_state is not None:
                return LayoutDimension(min=reserve_space_for_menu)

        return LayoutDimension()

    # Create and return Container instance.
    return HSplit([
        # The main input, with completion menus floating on top of it.
        FloatContainer(
            HSplit([
                ConditionalContainer(
                    Window(
                        TokenListControl(get_prompt_tokens_1),
                        dont_extend_height=True),
                    Condition(has_before_tokens)
                ),
                Window(
                    BufferControl(
                        input_processors=input_processors,
                        lexer=lexer,
                        # Enable preview_search, we want to have immediate feedback
                        # in reverse-i-search mode.
                        preview_search=True),
                    get_height=get_height,
                    left_margins=[
                        # In multiline mode, use the window margin to display
                        # the prompt and continuation tokens.
                        ConditionalMargin(
                            PromptMargin(get_prompt_tokens_2, get_continuation_tokens),
                            filter=multiline
                        )
                    ],
                    wrap_lines=wrap_lines,
                ),
            ]),
            [
                # Completion menus.
                Float(xcursor=True,
                      ycursor=True,
                      content=CompletionsMenu(
                          max_height=16,
                          scroll_offset=1,
                          extra_filter=HasFocus(DEFAULT_BUFFER) &
                                       ~display_completions_in_columns)),
                Float(xcursor=True,
                      ycursor=True,
                      content=MultiColumnCompletionsMenu(
                          extra_filter=HasFocus(DEFAULT_BUFFER) &
                                       display_completions_in_columns,
                          show_meta=True)),

                # The right prompt.
                Float(right=0, top=0, hide_when_covering_content=True,
                      content=_RPrompt(get_rprompt_tokens)),
            ]
        ),
        ValidationToolbar(),
        SystemToolbar(),

        # In multiline mode, we use two toolbars for 'arg' and 'search'.
        ConditionalContainer(ArgToolbar(), multiline),
        ConditionalContainer(SearchToolbar(), multiline),
    ] + toolbars)


def create_prompt_application(
        message='',
        multiline=False,
        wrap_lines=True,
        is_password=False,
        vi_mode=False,
        editing_mode=EditingMode.EMACS,
        complete_while_typing=True,
        enable_history_search=False,
        lexer=None,
        enable_system_bindings=False,
        enable_open_in_editor=False,
        validator=None,
        completer=None,
        reserve_space_for_menu=8,
        auto_suggest=None,
        style=None,
        history=None,
        clipboard=None,
        get_prompt_tokens=None,
        get_continuation_tokens=None,
        get_rprompt_tokens=None,
        get_bottom_toolbar_tokens=None,
        display_completions_in_columns=False,
        get_title=None,
        mouse_support=False,
        extra_input_processors=None,
        key_bindings_registry=None,
        on_abort=AbortAction.RAISE_EXCEPTION,
        on_exit=AbortAction.RAISE_EXCEPTION,
        accept_action=AcceptAction.RETURN_DOCUMENT,
        erase_when_done=False,
        default=''):
    """
    Create an :class:`~Application` instance for a prompt.

    (It is meant to cover 90% of the prompt use cases, where no extreme
    customization is required. For more complex input, it is required to create
    a custom :class:`~Application` instance.)

    :param message: Text to be shown before the prompt.
    :param mulitiline: Allow multiline input. Pressing enter will insert a
                       newline. (This requires Meta+Enter to accept the input.)
    :param wrap_lines: `bool` or :class:`~prompt_toolkit.filters.CLIFilter`.
        When True (the default), automatically wrap long lines instead of
        scrolling horizontally.
    :param is_password: Show asterisks instead of the actual typed characters.
    :param editing_mode: ``EditingMode.VI`` or ``EditingMode.EMACS``.
    :param vi_mode: `bool`, if True, Identical to ``editing_mode=EditingMode.VI``.
    :param complete_while_typing: `bool` or
        :class:`~prompt_toolkit.filters.SimpleFilter`. Enable autocompletion
        while typing.
    :param enable_history_search: `bool` or
        :class:`~prompt_toolkit.filters.SimpleFilter`. Enable up-arrow parting
        string matching.
    :param lexer: :class:`~prompt_toolkit.layout.lexers.Lexer` to be used for
        the syntax highlighting.
    :param validator: :class:`~prompt_toolkit.validation.Validator` instance
        for input validation.
    :param completer: :class:`~prompt_toolkit.completion.Completer` instance
        for input completion.
    :param reserve_space_for_menu: Space to be reserved for displaying the menu.
        (0 means that no space needs to be reserved.)
    :param auto_suggest: :class:`~prompt_toolkit.auto_suggest.AutoSuggest`
        instance for input suggestions.
    :param style: :class:`.Style` instance for the color scheme.
    :param enable_system_bindings: `bool` or
        :class:`~prompt_toolkit.filters.CLIFilter`. Pressing Meta+'!' will show
        a system prompt.
    :param enable_open_in_editor: `bool` or
        :class:`~prompt_toolkit.filters.CLIFilter`. Pressing 'v' in Vi mode or
        C-X C-E in emacs mode will open an external editor.
    :param history: :class:`~prompt_toolkit.history.History` instance.
    :param clipboard: :class:`~prompt_toolkit.clipboard.base.Clipboard` instance.
        (e.g. :class:`~prompt_toolkit.clipboard.in_memory.InMemoryClipboard`)
    :param get_bottom_toolbar_tokens: Optional callable which takes a
        :class:`~prompt_toolkit.interface.CommandLineInterface` and returns a
        list of tokens for the bottom toolbar.
    :param display_completions_in_columns: `bool` or
        :class:`~prompt_toolkit.filters.CLIFilter`. Display the completions in
        multiple columns.
    :param get_title: Callable that returns the title to be displayed in the
        terminal.
    :param mouse_support: `bool` or :class:`~prompt_toolkit.filters.CLIFilter`
        to enable mouse support.
    :param default: The default text to be shown in the input buffer. (This can
        be edited by the user.)
    """
    if key_bindings_registry is None:
        key_bindings_registry = load_key_bindings_for_prompt(
            enable_system_bindings=enable_system_bindings,
            enable_open_in_editor=enable_open_in_editor)

    # Ensure backwards-compatibility, when `vi_mode` is passed.
    if vi_mode:
        editing_mode = EditingMode.VI

    # Make sure that complete_while_typing is disabled when enable_history_search
    # is enabled. (First convert to SimpleFilter, to avoid doing bitwise operations
    # on bool objects.)
    complete_while_typing = to_simple_filter(complete_while_typing)
    enable_history_search = to_simple_filter(enable_history_search)
    multiline = to_simple_filter(multiline)

    complete_while_typing = complete_while_typing & ~enable_history_search

    # Accept Pygments styles as well for backwards compatibility.
    try:
        if pygments_Style and issubclass(style, pygments_Style):
            style = style_from_dict(style.styles)
    except TypeError:  # Happens when style is `None` or an instance of something else.
        pass

    # Create application
    return Application(
        layout=create_prompt_layout(
            message=message,
            lexer=lexer,
            is_password=is_password,
            reserve_space_for_menu=(reserve_space_for_menu if completer is not None else 0),
            multiline=Condition(lambda cli: multiline()),
            get_prompt_tokens=get_prompt_tokens,
            get_continuation_tokens=get_continuation_tokens,
            get_rprompt_tokens=get_rprompt_tokens,
            get_bottom_toolbar_tokens=get_bottom_toolbar_tokens,
            display_completions_in_columns=display_completions_in_columns,
            extra_input_processors=extra_input_processors,
            wrap_lines=wrap_lines),
        buffer=Buffer(
            enable_history_search=enable_history_search,
            complete_while_typing=complete_while_typing,
            is_multiline=multiline,
            history=(history or InMemoryHistory()),
            validator=validator,
            completer=completer,
            auto_suggest=auto_suggest,
            accept_action=accept_action,
            initial_document=Document(default),
        ),
        style=style or DEFAULT_STYLE,
        clipboard=clipboard,
        key_bindings_registry=key_bindings_registry,
        get_title=get_title,
        mouse_support=mouse_support,
        editing_mode=editing_mode,
        erase_when_done=erase_when_done,
        reverse_vi_search_direction=True,
        on_abort=on_abort,
        on_exit=on_exit)


def prompt(message='', **kwargs):
    """
    Get input from the user and return it.

    This is a wrapper around a lot of ``prompt_toolkit`` functionality and can
    be a replacement for `raw_input`. (or GNU readline.)

    If you want to keep your history across several calls, create one
    :class:`~prompt_toolkit.history.History` instance and pass it every time.

    This function accepts many keyword arguments. Except for the following,
    they are a proxy to the arguments of :func:`.create_prompt_application`.

    :param patch_stdout: Replace ``sys.stdout`` by a proxy that ensures that
            print statements from other threads won't destroy the prompt. (They
            will be printed above the prompt instead.)
    :param return_asyncio_coroutine: When True, return a asyncio coroutine. (Python >3.3)
    :param true_color: When True, use 24bit colors instead of 256 colors.
    :param refresh_interval: (number; in seconds) When given, refresh the UI
        every so many seconds.
    """
    patch_stdout = kwargs.pop('patch_stdout', False)
    return_asyncio_coroutine = kwargs.pop('return_asyncio_coroutine', False)
    true_color = kwargs.pop('true_color', False)
    refresh_interval = kwargs.pop('refresh_interval', 0)
    eventloop = kwargs.pop('eventloop', None)

    application = create_prompt_application(message, **kwargs)

    return run_application(application,
        patch_stdout=patch_stdout,
        return_asyncio_coroutine=return_asyncio_coroutine,
        true_color=true_color,
        refresh_interval=refresh_interval,
        eventloop=eventloop)


def run_application(
        application, patch_stdout=False, return_asyncio_coroutine=False,
        true_color=False, refresh_interval=0, eventloop=None):
    """
    Run a prompt toolkit application.

    :param patch_stdout: Replace ``sys.stdout`` by a proxy that ensures that
            print statements from other threads won't destroy the prompt. (They
            will be printed above the prompt instead.)
    :param return_asyncio_coroutine: When True, return a asyncio coroutine. (Python >3.3)
    :param true_color: When True, use 24bit colors instead of 256 colors.
    :param refresh_interval: (number; in seconds) When given, refresh the UI
        every so many seconds.
    """
    assert isinstance(application, Application)

    if return_asyncio_coroutine:
        eventloop = create_asyncio_eventloop()
    else:
        eventloop = eventloop or create_eventloop()

    # Create CommandLineInterface.
    cli = CommandLineInterface(
        application=application,
        eventloop=eventloop,
        output=create_output(true_color=true_color))

    # Set up refresh interval.
    if refresh_interval:
        done = [False]
        def start_refresh_loop(cli):
            def run():
                while not done[0]:
                    time.sleep(refresh_interval)
                    cli.request_redraw()
            t = threading.Thread(target=run)
            t.daemon = True
            t.start()

        def stop_refresh_loop(cli):
            done[0] = True

        cli.on_start += start_refresh_loop
        cli.on_stop += stop_refresh_loop

    # Replace stdout.
    patch_context = cli.patch_stdout_context(raw=True) if patch_stdout else DummyContext()

    # Read input and return it.
    if return_asyncio_coroutine:
        # Create an asyncio coroutine and call it.
        exec_context = {'patch_context': patch_context, 'cli': cli,
                        'Document': Document}
        exec_(textwrap.dedent('''
        def prompt_coro():
            # Inline import, because it slows down startup when asyncio is not
            # needed.
            import asyncio

            @asyncio.coroutine
            def run():
                with patch_context:
                    result = yield from cli.run_async()

                if isinstance(result, Document):  # Backwards-compatibility.
                    return result.text
                return result
            return run()
        '''), exec_context)

        return exec_context['prompt_coro']()
    else:
        try:
            with patch_context:
                result = cli.run()

            if isinstance(result, Document):  # Backwards-compatibility.
                return result.text
            return result
        finally:
            eventloop.close()


def prompt_async(message='', **kwargs):
    """
    Similar to :func:`.prompt`, but return an asyncio coroutine instead.
    """
    kwargs['return_asyncio_coroutine'] = True
    return prompt(message, **kwargs)


def create_confirm_application(message):
    """
    Create a confirmation `Application` that returns True/False.
    """
    registry = Registry()

    @registry.add_binding('y')
    @registry.add_binding('Y')
    def _(event):
        event.cli.buffers[DEFAULT_BUFFER].text = 'y'
        event.cli.set_return_value(True)

    @registry.add_binding('n')
    @registry.add_binding('N')
    @registry.add_binding(Keys.ControlC)
    def _(event):
        event.cli.buffers[DEFAULT_BUFFER].text = 'n'
        event.cli.set_return_value(False)

    return create_prompt_application(message, key_bindings_registry=registry)


def confirm(message='Confirm (y or n) '):
    """
    Display a confirmation prompt.
    """
    assert isinstance(message, text_type)

    app = create_confirm_application(message)
    return run_application(app)


def print_tokens(tokens, style=None, true_color=False, file=None):
    """
    Print a list of (Token, text) tuples in the given style to the output.
    E.g.::

        style = style_from_dict({
            Token.Hello: '#ff0066',
            Token.World: '#884444 italic',
        })
        tokens = [
            (Token.Hello, 'Hello'),
            (Token.World, 'World'),
        ]
        print_tokens(tokens, style=style)

    :param tokens: List of ``(Token, text)`` tuples.
    :param style: :class:`.Style` instance for the color scheme.
    :param true_color: When True, use 24bit colors instead of 256 colors.
    :param file: The output file. This can be `sys.stdout` or `sys.stderr`.
    """
    if style is None:
        style = DEFAULT_STYLE
    assert isinstance(style, Style)

    output = create_output(true_color=true_color, stdout=file)
    renderer_print_tokens(output, tokens, style)


def clear():
    """
    Clear the screen.
    """
    out = create_output()
    out.erase_screen()
    out.cursor_goto(0, 0)
    out.flush()


# Deprecated alias for `prompt`.
get_input = prompt
# Deprecated alias for create_prompt_layout
create_default_layout = create_prompt_layout
# Deprecated alias for create_prompt_application
create_default_application = create_prompt_application
