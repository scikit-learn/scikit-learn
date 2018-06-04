from __future__ import unicode_literals

from .buffer import Buffer, AcceptAction
from .buffer_mapping import BufferMapping
from .clipboard import Clipboard, InMemoryClipboard
from .enums import DEFAULT_BUFFER, EditingMode
from .filters import CLIFilter, to_cli_filter
from .key_binding.bindings.basic import load_basic_bindings
from .key_binding.bindings.emacs import load_emacs_bindings
from .key_binding.bindings.vi import load_vi_bindings
from .key_binding.registry import BaseRegistry
from .key_binding.defaults import load_key_bindings
from .layout import Window
from .layout.containers import Container
from .layout.controls import BufferControl
from .styles import DEFAULT_STYLE, Style
import six

__all__ = (
    'AbortAction',
    'Application',
)


class AbortAction(object):
    """
    Actions to take on an Exit or Abort exception.
    """
    RETRY = 'retry'
    RAISE_EXCEPTION = 'raise-exception'
    RETURN_NONE = 'return-none'

    _all = (RETRY, RAISE_EXCEPTION, RETURN_NONE)


class Application(object):
    """
    Application class to be passed to a
    :class:`~prompt_toolkit.interface.CommandLineInterface`.

    This contains all customizable logic that is not I/O dependent.
    (So, what is independent of event loops, input and output.)

    This way, such an :class:`.Application` can run easily on several
    :class:`~prompt_toolkit.interface.CommandLineInterface` instances, each
    with a different I/O backends. that runs for instance over telnet, SSH or
    any other I/O backend.

    :param layout: A :class:`~prompt_toolkit.layout.containers.Container` instance.
    :param buffer: A :class:`~prompt_toolkit.buffer.Buffer` instance for the default buffer.
    :param initial_focussed_buffer: Name of the buffer that is focussed during start-up.
    :param key_bindings_registry:
        :class:`~prompt_toolkit.key_binding.registry.BaseRegistry` instance for
        the key bindings.
    :param clipboard: :class:`~prompt_toolkit.clipboard.base.Clipboard` to use.
    :param on_abort: What to do when Control-C is pressed.
    :param on_exit: What to do when Control-D is pressed.
    :param use_alternate_screen: When True, run the application on the alternate screen buffer.
    :param get_title: Callable that returns the current title to be displayed in the terminal.
    :param erase_when_done: (bool) Clear the application output when it finishes.
    :param reverse_vi_search_direction: Normally, in Vi mode, a '/' searches
        forward and a '?' searches backward. In readline mode, this is usually
        reversed.

    Filters:

    :param mouse_support: (:class:`~prompt_toolkit.filters.CLIFilter` or
        boolean). When True, enable mouse support.
    :param paste_mode: :class:`~prompt_toolkit.filters.CLIFilter` or boolean.
    :param ignore_case: :class:`~prompt_toolkit.filters.CLIFilter` or boolean.
    :param editing_mode: :class:`~prompt_toolkit.enums.EditingMode`.

    Callbacks (all of these should accept a
    :class:`~prompt_toolkit.interface.CommandLineInterface` object as input.)

    :param on_input_timeout: Called when there is no input for x seconds.
                    (Fired when any eventloop.onInputTimeout is fired.)
    :param on_start: Called when reading input starts.
    :param on_stop: Called when reading input ends.
    :param on_reset: Called during reset.
    :param on_buffer_changed: Called when the content of a buffer has been changed.
    :param on_initialize: Called after the
        :class:`~prompt_toolkit.interface.CommandLineInterface` initializes.
    :param on_render: Called right after rendering.
    :param on_invalidate: Called when the UI has been invalidated.
    """
    def __init__(self, layout=None, buffer=None, buffers=None,
                 initial_focussed_buffer=DEFAULT_BUFFER,
                 style=None,
                 key_bindings_registry=None, clipboard=None,
                 on_abort=AbortAction.RAISE_EXCEPTION, on_exit=AbortAction.RAISE_EXCEPTION,
                 use_alternate_screen=False, mouse_support=False,
                 get_title=None,

                 paste_mode=False, ignore_case=False, editing_mode=EditingMode.EMACS,
                 erase_when_done=False,
                 reverse_vi_search_direction=False,

                 on_input_timeout=None, on_start=None, on_stop=None,
                 on_reset=None, on_initialize=None, on_buffer_changed=None,
                 on_render=None, on_invalidate=None):

        paste_mode = to_cli_filter(paste_mode)
        ignore_case = to_cli_filter(ignore_case)
        mouse_support = to_cli_filter(mouse_support)
        reverse_vi_search_direction = to_cli_filter(reverse_vi_search_direction)

        assert layout is None or isinstance(layout, Container)
        assert buffer is None or isinstance(buffer, Buffer)
        assert buffers is None or isinstance(buffers, (dict, BufferMapping))
        assert key_bindings_registry is None or isinstance(key_bindings_registry, BaseRegistry)
        assert clipboard is None or isinstance(clipboard, Clipboard)
        assert on_abort in AbortAction._all
        assert on_exit in AbortAction._all
        assert isinstance(use_alternate_screen, bool)
        assert get_title is None or callable(get_title)
        assert isinstance(paste_mode, CLIFilter)
        assert isinstance(ignore_case, CLIFilter)
        assert isinstance(editing_mode, six.string_types)
        assert on_input_timeout is None or callable(on_input_timeout)
        assert style is None or isinstance(style, Style)
        assert isinstance(erase_when_done, bool)

        assert on_start is None or callable(on_start)
        assert on_stop is None or callable(on_stop)
        assert on_reset is None or callable(on_reset)
        assert on_buffer_changed is None or callable(on_buffer_changed)
        assert on_initialize is None or callable(on_initialize)
        assert on_render is None or callable(on_render)
        assert on_invalidate is None or callable(on_invalidate)

        self.layout = layout or Window(BufferControl())

        # Make sure that the 'buffers' dictionary is a BufferMapping.
        # NOTE: If no buffer is given, we create a default Buffer, with IGNORE as
        #       default accept_action. This is what makes sense for most users
        #       creating full screen applications. Doing nothing is the obvious
        #       default. Those creating a REPL would use the shortcuts module that
        #       passes in RETURN_DOCUMENT.
        self.buffer = buffer or Buffer(accept_action=AcceptAction.IGNORE)
        if not buffers or not isinstance(buffers, BufferMapping):
            self.buffers = BufferMapping(buffers, initial=initial_focussed_buffer)
        else:
            self.buffers = buffers

        if buffer:
            self.buffers[DEFAULT_BUFFER] = buffer

        self.initial_focussed_buffer = initial_focussed_buffer

        self.style = style or DEFAULT_STYLE

        if key_bindings_registry is None:
            key_bindings_registry = load_key_bindings()

        if get_title is None:
            get_title = lambda: None

        self.key_bindings_registry = key_bindings_registry
        self.clipboard = clipboard or InMemoryClipboard()
        self.on_abort = on_abort
        self.on_exit = on_exit
        self.use_alternate_screen = use_alternate_screen
        self.mouse_support = mouse_support
        self.get_title = get_title

        self.paste_mode = paste_mode
        self.ignore_case = ignore_case
        self.editing_mode = editing_mode
        self.erase_when_done = erase_when_done
        self.reverse_vi_search_direction = reverse_vi_search_direction

        def dummy_handler(cli):
            " Dummy event handler. "

        self.on_input_timeout = on_input_timeout or dummy_handler
        self.on_start = on_start or dummy_handler
        self.on_stop = on_stop or dummy_handler
        self.on_reset = on_reset or dummy_handler
        self.on_initialize = on_initialize or dummy_handler
        self.on_buffer_changed = on_buffer_changed or dummy_handler
        self.on_render = on_render or dummy_handler
        self.on_invalidate = on_invalidate or dummy_handler

        # List of 'extra' functions to execute before a CommandLineInterface.run.
        # Note: It's important to keep this here, and not in the
        #       CommandLineInterface itself. shortcuts.run_application creates
        #       a new Application instance everytime. (Which is correct, it
        #       could be that we want to detach from one IO backend and attach
        #       the UI on a different backend.) But important is to keep as
        #       much state as possible between runs.
        self.pre_run_callables = []
