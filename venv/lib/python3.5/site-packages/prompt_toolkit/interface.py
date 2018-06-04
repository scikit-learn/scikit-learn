"""
The main `CommandLineInterface` class and logic.
"""
from __future__ import unicode_literals

import functools
import os
import signal
import six
import sys
import textwrap
import threading
import time
import types
import weakref

from subprocess import Popen

from .application import Application, AbortAction
from .buffer import Buffer
from .buffer_mapping import BufferMapping
from .completion import CompleteEvent, get_common_complete_suffix
from .enums import SEARCH_BUFFER
from .eventloop.base import EventLoop
from .eventloop.callbacks import EventLoopCallbacks
from .filters import Condition
from .input import StdinInput, Input
from .key_binding.input_processor import InputProcessor
from .key_binding.input_processor import KeyPress
from .key_binding.registry import Registry
from .key_binding.vi_state import ViState
from .keys import Keys
from .output import Output
from .renderer import Renderer, print_tokens
from .search_state import SearchState
from .utils import Event

# Following import is required for backwards compatibility.
from .buffer import AcceptAction

__all__ = (
    'AbortAction',
    'CommandLineInterface',
)


class CommandLineInterface(object):
    """
    Wrapper around all the other classes, tying everything together.

    Typical usage::

        application = Application(...)
        cli = CommandLineInterface(application, eventloop)
        result = cli.run()
        print(result)

    :param application: :class:`~prompt_toolkit.application.Application` instance.
    :param eventloop: The :class:`~prompt_toolkit.eventloop.base.EventLoop` to
                      be used when `run` is called. The easiest way to create
                      an eventloop is by calling
                      :meth:`~prompt_toolkit.shortcuts.create_eventloop`.
    :param input: :class:`~prompt_toolkit.input.Input` instance.
    :param output: :class:`~prompt_toolkit.output.Output` instance. (Probably
                   Vt100_Output or Win32Output.)
    """
    def __init__(self, application, eventloop=None, input=None, output=None):
        assert isinstance(application, Application)
        assert isinstance(eventloop, EventLoop), 'Passing an eventloop is required.'
        assert output is None or isinstance(output, Output)
        assert input is None or isinstance(input, Input)

        from .shortcuts import create_output

        self.application = application
        self.eventloop = eventloop
        self._is_running = False

        # Inputs and outputs.
        self.output = output or create_output()
        self.input = input or StdinInput(sys.stdin)

        #: The input buffers.
        assert isinstance(application.buffers, BufferMapping)
        self.buffers = application.buffers

        #: EditingMode.VI or EditingMode.EMACS
        self.editing_mode = application.editing_mode

        #: Quoted insert. This flag is set if we go into quoted insert mode.
        self.quoted_insert = False

        #: Vi state. (For Vi key bindings.)
        self.vi_state = ViState()

        #: The `Renderer` instance.
        # Make sure that the same stdout is used, when a custom renderer has been passed.
        self.renderer = Renderer(
            self.application.style,
            self.output,
            use_alternate_screen=application.use_alternate_screen,
            mouse_support=application.mouse_support)

        #: Render counter. This one is increased every time the UI is rendered.
        #: It can be used as a key for caching certain information during one
        #: rendering.
        self.render_counter = 0

        #: When there is high CPU, postpone the renderering max x seconds.
        #: '0' means: don't postpone. '.5' means: try to draw at least twice a second.
        self.max_render_postpone_time = 0  # E.g. .5

        # Invalidate flag. When 'True', a repaint has been scheduled.
        self._invalidated = False

        #: The `InputProcessor` instance.
        self.input_processor = InputProcessor(application.key_bindings_registry, weakref.ref(self))

        self._async_completers = {}  # Map buffer name to completer function.

        # Pointer to sub CLI. (In chain of CLI instances.)
        self._sub_cli = None  # None or other CommandLineInterface instance.

        # Call `add_buffer` for each buffer.
        for name, b in self.buffers.items():
            self.add_buffer(name, b)

        # Events.
        self.on_buffer_changed = Event(self, application.on_buffer_changed)
        self.on_initialize = Event(self, application.on_initialize)
        self.on_input_timeout = Event(self, application.on_input_timeout)
        self.on_invalidate = Event(self, application.on_invalidate)
        self.on_render = Event(self, application.on_render)
        self.on_reset = Event(self, application.on_reset)
        self.on_start = Event(self, application.on_start)
        self.on_stop = Event(self, application.on_stop)

        # Trigger initialize callback.
        self.reset()
        self.on_initialize += self.application.on_initialize
        self.on_initialize.fire()

    @property
    def layout(self):
        return self.application.layout

    @property
    def clipboard(self):
        return self.application.clipboard

    @property
    def pre_run_callables(self):
        return self.application.pre_run_callables

    def add_buffer(self, name, buffer, focus=False):
        """
        Insert a new buffer.
        """
        assert isinstance(buffer, Buffer)
        self.buffers[name] = buffer

        if focus:
            self.buffers.focus(name)

        # Create asynchronous completer / auto suggestion.
        auto_suggest_function = self._create_auto_suggest_function(buffer)
        completer_function = self._create_async_completer(buffer)
        self._async_completers[name] = completer_function

        # Complete/suggest on text insert.
        def create_on_insert_handler():
            """
            Wrapper around the asynchronous completer and auto suggestion, that
            ensures that it's only called while typing if the
            `complete_while_typing` filter is enabled.
            """
            def on_text_insert(_):
                # Only complete when "complete_while_typing" is enabled.
                if buffer.completer and buffer.complete_while_typing():
                    completer_function()

                # Call auto_suggest.
                if buffer.auto_suggest:
                    auto_suggest_function()

            return on_text_insert

        buffer.on_text_insert += create_on_insert_handler()

        def buffer_changed(_):
            """
            When the text in a buffer changes.
            (A paste event is also a change, but not an insert. So we don't
            want to do autocompletions in this case, but we want to propagate
            the on_buffer_changed event.)
            """
            # Trigger on_buffer_changed.
            self.on_buffer_changed.fire()

        buffer.on_text_changed += buffer_changed

    def start_completion(self, buffer_name=None, select_first=False,
                         select_last=False, insert_common_part=False,
                         complete_event=None):
        """
        Start asynchronous autocompletion of this buffer.
        (This will do nothing if a previous completion was still in progress.)
        """
        buffer_name = buffer_name or self.current_buffer_name
        completer = self._async_completers.get(buffer_name)

        if completer:
            completer(select_first=select_first,
                      select_last=select_last,
                      insert_common_part=insert_common_part,
                      complete_event=CompleteEvent(completion_requested=True))

    @property
    def current_buffer_name(self):
        """
        The name of the current  :class:`.Buffer`. (Or `None`.)
        """
        return self.buffers.current_name(self)

    @property
    def current_buffer(self):
        """
        The currently focussed :class:`~.Buffer`.

        (This returns a dummy :class:`.Buffer` when none of the actual buffers
        has the focus. In this case, it's really not practical to check for
        `None` values or catch exceptions every time.)
        """
        return self.buffers.current(self)

    def focus(self, buffer_name):
        """
        Focus the buffer with the given name on the focus stack.
        """
        self.buffers.focus(self, buffer_name)

    def push_focus(self, buffer_name):
        """
        Push to the focus stack.
        """
        self.buffers.push_focus(self, buffer_name)

    def pop_focus(self):
        """
        Pop from the focus stack.
        """
        self.buffers.pop_focus(self)

    @property
    def terminal_title(self):
        """
        Return the current title to be displayed in the terminal.
        When this in `None`, the terminal title remains the original.
        """
        result = self.application.get_title()

        # Make sure that this function returns a unicode object,
        # and not a byte string.
        assert result is None or isinstance(result, six.text_type)
        return result

    @property
    def is_searching(self):
        """
        True when we are searching.
        """
        return self.current_buffer_name == SEARCH_BUFFER

    def reset(self, reset_current_buffer=False):
        """
        Reset everything, for reading the next input.

        :param reset_current_buffer: XXX: not used anymore. The reason for
            having this option in the past was when this CommandLineInterface
            is run multiple times, that we could reset the buffer content from
            the previous run. This is now handled in the AcceptAction.
        """
        # Notice that we don't reset the buffers. (This happens just before
        # returning, and when we have multiple buffers, we clearly want the
        # content in the other buffers to remain unchanged between several
        # calls of `run`. (And the same is true for the focus stack.)

        self._exit_flag = False
        self._abort_flag = False

        self._return_value = None

        self.renderer.reset()
        self.input_processor.reset()
        self.layout.reset()
        self.vi_state.reset()

        # Search new search state. (Does also remember what has to be
        # highlighted.)
        self.search_state = SearchState(ignore_case=Condition(lambda: self.is_ignoring_case))

        # Trigger reset event.
        self.on_reset.fire()

    @property
    def in_paste_mode(self):
        """ True when we are in paste mode. """
        return self.application.paste_mode(self)

    @property
    def is_ignoring_case(self):
        """ True when we currently ignore casing. """
        return self.application.ignore_case(self)

    def invalidate(self):
        """
        Thread safe way of sending a repaint trigger to the input event loop.
        """
        # Never schedule a second redraw, when a previous one has not yet been
        # executed. (This should protect against other threads calling
        # 'invalidate' many times, resulting in 100% CPU.)
        if self._invalidated:
            return
        else:
            self._invalidated = True

        # Trigger event.
        self.on_invalidate.fire()

        if self.eventloop is not None:
            def redraw():
                self._invalidated = False
                self._redraw()

            # Call redraw in the eventloop (thread safe).
            # Usually with the high priority, in order to make the application
            # feel responsive, but this can be tuned by changing the value of
            # `max_render_postpone_time`.
            if self.max_render_postpone_time:
                _max_postpone_until = time.time() + self.max_render_postpone_time
            else:
                _max_postpone_until = None

            self.eventloop.call_from_executor(
                redraw, _max_postpone_until=_max_postpone_until)

    # Depracated alias for 'invalidate'.
    request_redraw = invalidate

    def _redraw(self):
        """
        Render the command line again. (Not thread safe!) (From other threads,
        or if unsure, use :meth:`.CommandLineInterface.invalidate`.)
        """
        # Only draw when no sub application was started.
        if self._is_running and self._sub_cli is None:
            self.render_counter += 1
            self.renderer.render(self, self.layout, is_done=self.is_done)

            # Fire render event.
            self.on_render.fire()

    def _on_resize(self):
        """
        When the window size changes, we erase the current output and request
        again the cursor position. When the CPR answer arrives, the output is
        drawn again.
        """
        # Erase, request position (when cursor is at the start position)
        # and redraw again. -- The order is important.
        self.renderer.erase(leave_alternate_screen=False, erase_title=False)
        self.renderer.request_absolute_cursor_position()
        self._redraw()

    def _load_next_buffer_indexes(self):
        for buff, index in self._next_buffer_indexes.items():
            if buff in self.buffers:
                self.buffers[buff].working_index = index

    def _pre_run(self, pre_run=None):
        " Called during `run`. "
        if pre_run:
            pre_run()

        # Process registered "pre_run_callables" and clear list.
        for c in self.pre_run_callables:
            c()
        del self.pre_run_callables[:]

    def run(self, reset_current_buffer=False, pre_run=None):
        """
        Read input from the command line.
        This runs the eventloop until a return value has been set.

        :param reset_current_buffer: XXX: Not used anymore.
        :param pre_run: Callable that is called right after the reset has taken
            place. This allows custom initialisation.
        """
        assert pre_run is None or callable(pre_run)

        try:
            self._is_running = True

            self.on_start.fire()
            self.reset()

            # Call pre_run.
            self._pre_run(pre_run)

            # Run eventloop in raw mode.
            with self.input.raw_mode():
                self.renderer.request_absolute_cursor_position()
                self._redraw()

                self.eventloop.run(self.input, self.create_eventloop_callbacks())
        finally:
            # Clean up renderer. (This will leave the alternate screen, if we use
            # that.)

            # If exit/abort haven't been called set, but another exception was
            # thrown instead for some reason, make sure that we redraw in exit
            # mode.
            if not self.is_done:
                self._exit_flag = True
                self._redraw()

            self.renderer.reset()
            self.on_stop.fire()
            self._is_running = False

        # Return result.
        return self.return_value()

    try:
        # The following `run_async` function is compiled at runtime
        # because it contains syntax which is not supported on older Python
        # versions. (A 'return' inside a generator.)
        six.exec_(textwrap.dedent('''
        def run_async(self, reset_current_buffer=True, pre_run=None):
            """
            Same as `run`, but this returns a coroutine.

            This is only available on Python >3.3, with asyncio.
            """
            # Inline import, because it slows down startup when asyncio is not
            # needed.
            import asyncio

            @asyncio.coroutine
            def run():
                assert pre_run is None or callable(pre_run)

                try:
                    self._is_running = True

                    self.on_start.fire()
                    self.reset()

                    # Call pre_run.
                    self._pre_run(pre_run)

                    with self.input.raw_mode():
                        self.renderer.request_absolute_cursor_position()
                        self._redraw()

                        yield from self.eventloop.run_as_coroutine(
                                self.input, self.create_eventloop_callbacks())

                    return self.return_value()
                finally:
                    if not self.is_done:
                        self._exit_flag = True
                        self._redraw()

                    self.renderer.reset()
                    self.on_stop.fire()
                    self._is_running = False

            return run()
        '''))
    except SyntaxError:
        # Python2, or early versions of Python 3.
        def run_async(self, reset_current_buffer=True, pre_run=None):
            """
            Same as `run`, but this returns a coroutine.

            This is only available on Python >3.3, with asyncio.
            """
            raise NotImplementedError

    def run_sub_application(self, application, done_callback=None, erase_when_done=False,
                            _from_application_generator=False):
        # `erase_when_done` is deprecated, set Application.erase_when_done instead.
        """
        Run a sub :class:`~prompt_toolkit.application.Application`.

        This will suspend the main application and display the sub application
        until that one returns a value. The value is returned by calling
        `done_callback` with the result.

        The sub application will share the same I/O of the main application.
        That means, it uses the same input and output channels and it shares
        the same event loop.

        .. note:: Technically, it gets another Eventloop instance, but that is
            only a proxy to our main event loop. The reason is that calling
            'stop' --which returns the result of an application when it's
            done-- is handled differently.
        """
        assert isinstance(application, Application)
        assert done_callback is None or callable(done_callback)

        if self._sub_cli is not None:
            raise RuntimeError('Another sub application started already.')

        # Erase current application.
        if not _from_application_generator:
            self.renderer.erase()

        # Callback when the sub app is done.
        def done():
            # Redraw sub app in done state.
            # and reset the renderer. (This reset will also quit the alternate
            # screen, if the sub application used that.)
            sub_cli._redraw()
            if erase_when_done or application.erase_when_done:
                sub_cli.renderer.erase()
            sub_cli.renderer.reset()
            sub_cli._is_running = False  # Don't render anymore.

            self._sub_cli = None

            # Restore main application.
            if not _from_application_generator:
                self.renderer.request_absolute_cursor_position()
                self._redraw()

            # Deliver result.
            if done_callback:
                done_callback(sub_cli.return_value())

        # Create sub CommandLineInterface.
        sub_cli = CommandLineInterface(
            application=application,
            eventloop=_SubApplicationEventLoop(self, done),
            input=self.input,
            output=self.output)
        sub_cli._is_running = True  # Allow rendering of sub app.

        sub_cli._redraw()
        self._sub_cli = sub_cli

    def exit(self):
        """
        Set exit. When Control-D has been pressed.
        """
        on_exit = self.application.on_exit
        self._exit_flag = True
        self._redraw()

        if on_exit == AbortAction.RAISE_EXCEPTION:
            def eof_error():
                raise EOFError()
            self._set_return_callable(eof_error)

        elif on_exit == AbortAction.RETRY:
            self.reset()
            self.renderer.request_absolute_cursor_position()
            self.current_buffer.reset()

        elif on_exit == AbortAction.RETURN_NONE:
            self.set_return_value(None)

    def abort(self):
        """
        Set abort. When Control-C has been pressed.
        """
        on_abort = self.application.on_abort
        self._abort_flag = True
        self._redraw()

        if on_abort == AbortAction.RAISE_EXCEPTION:
            def keyboard_interrupt():
                raise KeyboardInterrupt()
            self._set_return_callable(keyboard_interrupt)

        elif on_abort == AbortAction.RETRY:
            self.reset()
            self.renderer.request_absolute_cursor_position()
            self.current_buffer.reset()

        elif on_abort == AbortAction.RETURN_NONE:
            self.set_return_value(None)

    # Deprecated aliase for exit/abort.
    set_exit = exit
    set_abort = abort

    def set_return_value(self, document):
        """
        Set a return value. The eventloop can retrieve the result it by calling
        `return_value`.
        """
        self._set_return_callable(lambda: document)
        self._redraw()  # Redraw in "done" state, after the return value has been set.

    def _set_return_callable(self, value):
        assert callable(value)
        self._return_value = value

        if self.eventloop:
            self.eventloop.stop()

    def run_in_terminal(self, func, render_cli_done=False, cooked_mode=True):
        """
        Run function on the terminal above the prompt.

        What this does is first hiding the prompt, then running this callable
        (which can safely output to the terminal), and then again rendering the
        prompt which causes the output of this function to scroll above the
        prompt.

        :param func: The callable to execute.
        :param render_cli_done: When True, render the interface in the
                'Done' state first, then execute the function. If False,
                erase the interface first.
        :param cooked_mode: When True (the default), switch the input to
                cooked mode while executing the function.

        :returns: the result of `func`.
        """
        # Draw interface in 'done' state, or erase.
        if render_cli_done:
            self._return_value = True
            self._redraw()
            self.renderer.reset()  # Make sure to disable mouse mode, etc...
        else:
            self.renderer.erase()
        self._return_value = None

        # Run system command.
        if cooked_mode:
            with self.input.cooked_mode():
                result = func()
        else:
            result = func()

        # Redraw interface again.
        self.renderer.reset()
        self.renderer.request_absolute_cursor_position()
        self._redraw()

        return result

    def run_application_generator(self, coroutine, render_cli_done=False):
        """
        EXPERIMENTAL
        Like `run_in_terminal`, but takes a generator that can yield Application instances.

        Example:

            def f():
                yield Application1(...)
                print('...')
                yield Application2(...)
            cli.run_in_terminal_async(f)

        The values which are yielded by the given coroutine are supposed to be
        `Application` instances that run in the current CLI, all other code is
        supposed to be CPU bound, so except for yielding the applications,
        there should not be any user interaction or I/O in the given function.
        """
        # Draw interface in 'done' state, or erase.
        if render_cli_done:
            self._return_value = True
            self._redraw()
            self.renderer.reset()  # Make sure to disable mouse mode, etc...
        else:
            self.renderer.erase()
        self._return_value = None

        # Loop through the generator.
        g = coroutine()
        assert isinstance(g, types.GeneratorType)

        def step_next(send_value=None):
            " Execute next step of the coroutine."
            try:
                # Run until next yield, in cooked mode.
                with self.input.cooked_mode():
                    result = g.send(send_value)
            except StopIteration:
                done()
            except:
                done()
                raise
            else:
                # Process yielded value from coroutine.
                assert isinstance(result, Application)
                self.run_sub_application(result, done_callback=step_next,
                                         _from_application_generator=True)

        def done():
            # Redraw interface again.
            self.renderer.reset()
            self.renderer.request_absolute_cursor_position()
            self._redraw()

        # Start processing coroutine.
        step_next()

    def run_system_command(self, command):
        """
        Run system command (While hiding the prompt. When finished, all the
        output will scroll above the prompt.)

        :param command: Shell command to be executed.
        """
        def wait_for_enter():
            """
            Create a sub application to wait for the enter key press.
            This has two advantages over using 'input'/'raw_input':
            - This will share the same input/output I/O.
            - This doesn't block the event loop.
            """
            from .shortcuts import create_prompt_application

            registry = Registry()

            @registry.add_binding(Keys.ControlJ)
            @registry.add_binding(Keys.ControlM)
            def _(event):
                event.cli.set_return_value(None)

            application = create_prompt_application(
                message='Press ENTER to continue...',
                key_bindings_registry=registry)
            self.run_sub_application(application)

        def run():
            # Try to use the same input/output file descriptors as the one,
            # used to run this application.
            try:
                input_fd = self.input.fileno()
            except AttributeError:
                input_fd = sys.stdin.fileno()
            try:
                output_fd = self.output.fileno()
            except AttributeError:
                output_fd = sys.stdout.fileno()

            # Run sub process.
            # XXX: This will still block the event loop.
            p = Popen(command, shell=True,
                      stdin=input_fd, stdout=output_fd)
            p.wait()

            # Wait for the user to press enter.
            wait_for_enter()

        self.run_in_terminal(run)

    def suspend_to_background(self, suspend_group=True):
        """
        (Not thread safe -- to be called from inside the key bindings.)
        Suspend process.

        :param suspend_group: When true, suspend the whole process group.
            (This is the default, and probably what you want.)
        """
        # Only suspend when the opperating system supports it.
        # (Not on Windows.)
        if hasattr(signal, 'SIGTSTP'):
            def run():
                # Send `SIGSTP` to own process.
                # This will cause it to suspend.

                # Usually we want the whole process group to be suspended. This
                # handles the case when input is piped from another process.
                if suspend_group:
                    os.kill(0, signal.SIGTSTP)
                else:
                    os.kill(os.getpid(), signal.SIGTSTP)

            self.run_in_terminal(run)

    def print_tokens(self, tokens, style=None):
        """
        Print a list of (Token, text) tuples to the output.
        (When the UI is running, this method has to be called through
        `run_in_terminal`, otherwise it will destroy the UI.)

        :param style: Style class to use. Defaults to the active style in the CLI.
        """
        print_tokens(self.output, tokens, style or self.application.style)

    @property
    def is_exiting(self):
        """
        ``True`` when the exit flag as been set.
        """
        return self._exit_flag

    @property
    def is_aborting(self):
        """
        ``True`` when the abort flag as been set.
        """
        return self._abort_flag

    @property
    def is_returning(self):
        """
        ``True`` when a return value has been set.
        """
        return self._return_value is not None

    def return_value(self):
        """
        Get the return value. Not that this method can throw an exception.
        """
        # Note that it's a method, not a property, because it can throw
        # exceptions.
        if self._return_value:
            return self._return_value()

    @property
    def is_done(self):
        return self.is_exiting or self.is_aborting or self.is_returning

    def _create_async_completer(self, buffer):
        """
        Create function for asynchronous autocompletion.
        (Autocomplete in other thread.)
        """
        complete_thread_running = [False]  # By ref.

        def completion_does_nothing(document, completion):
            """
            Return `True` if applying this completion doesn't have any effect.
            (When it doesn't insert any new text.
            """
            text_before_cursor = document.text_before_cursor
            replaced_text = text_before_cursor[
                len(text_before_cursor) + completion.start_position:]
            return replaced_text == completion.text

        def async_completer(select_first=False, select_last=False,
                            insert_common_part=False, complete_event=None):
            document = buffer.document
            complete_event = complete_event or CompleteEvent(text_inserted=True)

            # Don't start two threads at the same time.
            if complete_thread_running[0]:
                return

            # Don't complete when we already have completions.
            if buffer.complete_state or not buffer.completer:
                return

            # Otherwise, get completions in other thread.
            complete_thread_running[0] = True

            def run():
                completions = list(buffer.completer.get_completions(document, complete_event))

                def callback():
                    """
                    Set the new complete_state in a safe way. Don't replace an
                    existing complete_state if we had one. (The user could have
                    pressed 'Tab' in the meantime. Also don't set it if the text
                    was changed in the meantime.
                    """
                    complete_thread_running[0] = False

                    # When there is only one completion, which has nothing to add, ignore it.
                    if (len(completions) == 1 and
                            completion_does_nothing(document, completions[0])):
                        del completions[:]

                    # Set completions if the text was not yet changed.
                    if buffer.text == document.text and \
                            buffer.cursor_position == document.cursor_position and \
                            not buffer.complete_state:

                        set_completions = True
                        select_first_anyway = False

                        # When the common part has to be inserted, and there
                        # is a common part.
                        if insert_common_part:
                            common_part = get_common_complete_suffix(document, completions)
                            if common_part:
                                # Insert the common part, update completions.
                                buffer.insert_text(common_part)
                                if len(completions) > 1:
                                    # (Don't call `async_completer` again, but
                                    # recalculate completions. See:
                                    # https://github.com/ipython/ipython/issues/9658)
                                    completions[:] = [
                                        c.new_completion_from_position(len(common_part))
                                        for c in completions]
                                else:
                                    set_completions = False
                            else:
                                # When we were asked to insert the "common"
                                # prefix, but there was no common suffix but
                                # still exactly one match, then select the
                                # first. (It could be that we have a completion
                                # which does * expansion, like '*.py', with
                                # exactly one match.)
                                if len(completions) == 1:
                                    select_first_anyway = True

                        if set_completions:
                            buffer.set_completions(
                                completions=completions,
                                go_to_first=select_first or select_first_anyway,
                                go_to_last=select_last)
                        self.invalidate()
                    elif not buffer.complete_state:
                        # Otherwise, restart thread.
                        async_completer()

                if self.eventloop:
                    self.eventloop.call_from_executor(callback)

            self.eventloop.run_in_executor(run)
        return async_completer

    def _create_auto_suggest_function(self, buffer):
        """
        Create function for asynchronous auto suggestion.
        (AutoSuggest in other thread.)
        """
        suggest_thread_running = [False]  # By ref.

        def async_suggestor():
            document = buffer.document

            # Don't start two threads at the same time.
            if suggest_thread_running[0]:
                return

            # Don't suggest when we already have a suggestion.
            if buffer.suggestion or not buffer.auto_suggest:
                return

            # Otherwise, get completions in other thread.
            suggest_thread_running[0] = True

            def run():
                suggestion = buffer.auto_suggest.get_suggestion(self, buffer, document)

                def callback():
                    suggest_thread_running[0] = False

                    # Set suggestion only if the text was not yet changed.
                    if buffer.text == document.text and \
                            buffer.cursor_position == document.cursor_position:

                        # Set suggestion and redraw interface.
                        buffer.suggestion = suggestion
                        self.invalidate()
                    else:
                        # Otherwise, restart thread.
                        async_suggestor()

                if self.eventloop:
                    self.eventloop.call_from_executor(callback)

            self.eventloop.run_in_executor(run)
        return async_suggestor

    def stdout_proxy(self, raw=False):
        """
        Create an :class:`_StdoutProxy` class which can be used as a patch for
        `sys.stdout`. Writing to this proxy will make sure that the text
        appears above the prompt, and that it doesn't destroy the output from
        the renderer.

        :param raw: (`bool`) When True, vt100 terminal escape sequences are not
                    removed/escaped.
        """
        return _StdoutProxy(self, raw=raw)

    def patch_stdout_context(self, raw=False, patch_stdout=True, patch_stderr=True):
        """
        Return a context manager that will replace ``sys.stdout`` with a proxy
        that makes sure that all printed text will appear above the prompt, and
        that it doesn't destroy the output from the renderer.

        :param patch_stdout: Replace `sys.stdout`.
        :param patch_stderr: Replace `sys.stderr`.
        """
        return _PatchStdoutContext(
            self.stdout_proxy(raw=raw),
            patch_stdout=patch_stdout, patch_stderr=patch_stderr)

    def create_eventloop_callbacks(self):
        return _InterfaceEventLoopCallbacks(self)


class _InterfaceEventLoopCallbacks(EventLoopCallbacks):
    """
    Callbacks on the :class:`.CommandLineInterface` object, to which an
    eventloop can talk.
    """
    def __init__(self, cli):
        assert isinstance(cli, CommandLineInterface)
        self.cli = cli

    @property
    def _active_cli(self):
        """
        Return the active `CommandLineInterface`.
        """
        cli = self.cli

        # If there is a sub CLI. That one is always active.
        while cli._sub_cli:
            cli = cli._sub_cli

        return cli

    def terminal_size_changed(self):
        """
        Report terminal size change. This will trigger a redraw.
        """
        self._active_cli._on_resize()

    def input_timeout(self):
        cli = self._active_cli
        cli.on_input_timeout.fire()

    def feed_key(self, key_press):
        """
        Feed a key press to the CommandLineInterface.
        """
        assert isinstance(key_press, KeyPress)
        cli = self._active_cli

        # Feed the key and redraw.
        # (When the CLI is in 'done' state, it should return to the event loop
        # as soon as possible. Ignore all key presses beyond this point.)
        if not cli.is_done:
            cli.input_processor.feed(key_press)
            cli.input_processor.process_keys()


class _PatchStdoutContext(object):
    def __init__(self, new_stdout, patch_stdout=True, patch_stderr=True):
        self.new_stdout = new_stdout
        self.patch_stdout = patch_stdout
        self.patch_stderr = patch_stderr

    def __enter__(self):
        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr

        if self.patch_stdout:
            sys.stdout = self.new_stdout
        if self.patch_stderr:
            sys.stderr = self.new_stdout

    def __exit__(self, *a, **kw):
        if self.patch_stdout:
            sys.stdout = self.original_stdout

        if self.patch_stderr:
            sys.stderr = self.original_stderr


class _StdoutProxy(object):
    """
    Proxy for stdout, as returned by
    :class:`CommandLineInterface.stdout_proxy`.
    """
    def __init__(self, cli, raw=False):
        assert isinstance(cli, CommandLineInterface)
        assert isinstance(raw, bool)

        self._lock = threading.RLock()
        self._cli = cli
        self._raw = raw
        self._buffer = []

        self.errors = sys.__stdout__.errors
        self.encoding = sys.__stdout__.encoding

    def _do(self, func):
        if self._cli._is_running:
            run_in_terminal = functools.partial(self._cli.run_in_terminal, func)
            self._cli.eventloop.call_from_executor(run_in_terminal)
        else:
            func()

    def _write(self, data):
        """
        Note: print()-statements cause to multiple write calls.
              (write('line') and write('\n')). Of course we don't want to call
              `run_in_terminal` for every individual call, because that's too
              expensive, and as long as the newline hasn't been written, the
              text itself is again overwritter by the rendering of the input
              command line. Therefor, we have a little buffer which holds the
              text until a newline is written to stdout.
        """
        if '\n' in data:
            # When there is a newline in the data, write everything before the
            # newline, including the newline itself.
            before, after = data.rsplit('\n', 1)
            to_write = self._buffer + [before, '\n']
            self._buffer = [after]

            def run():
                for s in to_write:
                    if self._raw:
                        self._cli.output.write_raw(s)
                    else:
                        self._cli.output.write(s)
            self._do(run)
        else:
            # Otherwise, cache in buffer.
            self._buffer.append(data)

    def write(self, data):
        with self._lock:
            self._write(data)

    def _flush(self):
        def run():
            for s in self._buffer:
                if self._raw:
                    self._cli.output.write_raw(s)
                else:
                    self._cli.output.write(s)
            self._buffer = []
            self._cli.output.flush()
        self._do(run)

    def flush(self):
        """
        Flush buffered output.
        """
        with self._lock:
            self._flush()


class _SubApplicationEventLoop(EventLoop):
    """
    Eventloop used by sub applications.

    A sub application is an `Application` that is "spawned" by a parent
    application. The parent application is suspended temporarily and the sub
    application is displayed instead.

    It doesn't need it's own event loop. The `EventLoopCallbacks` from the
    parent application are redirected to the sub application. So if the event
    loop that is run by the parent application detects input, the callbacks
    will make sure that it's forwarded to the sub application.

    When the sub application has a return value set, it will terminate
    by calling the `stop` method of this event loop. This is used to
    transfer control back to the parent application.
    """
    def __init__(self, cli, stop_callback):
        assert isinstance(cli, CommandLineInterface)
        assert callable(stop_callback)

        self.cli = cli
        self.stop_callback = stop_callback

    def stop(self):
        self.stop_callback()

    def close(self):
        pass

    def run_in_executor(self, callback):
        self.cli.eventloop.run_in_executor(callback)

    def call_from_executor(self, callback, _max_postpone_until=None):
        self.cli.eventloop.call_from_executor(
            callback, _max_postpone_until=_max_postpone_until)

    def add_reader(self, fd, callback):
        self.cli.eventloop.add_reader(fd, callback)

    def remove_reader(self, fd):
        self.cli.eventloop.remove_reader(fd)
