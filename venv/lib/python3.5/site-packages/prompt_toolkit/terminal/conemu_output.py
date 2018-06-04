from __future__ import unicode_literals

from prompt_toolkit.renderer import Output

from .win32_output import Win32Output
from .vt100_output import Vt100_Output

__all__ = (
    'ConEmuOutput',
)


class ConEmuOutput(object):
    """
    ConEmu (Windows) output abstraction.

    ConEmu is a Windows console application, but it also supports ANSI escape
    sequences. This output class is actually a proxy to both `Win32Output` and
    `Vt100_Output`. It uses `Win32Output` for console sizing and scrolling, but
    all cursor movements and scrolling happens through the `Vt100_Output`.

    This way, we can have 256 colors in ConEmu and Cmder. Rendering will be
    even a little faster as well.

    http://conemu.github.io/
    http://gooseberrycreative.com/cmder/
    """
    def __init__(self, stdout):
        self.win32_output = Win32Output(stdout)
        self.vt100_output = Vt100_Output(stdout, lambda: None)

    def __getattr__(self, name):
        if name in ('get_size', 'get_rows_below_cursor_position',
                    'enable_mouse_support', 'disable_mouse_support',
                    'scroll_buffer_to_prompt', 'get_win32_screen_buffer_info',
                    'enable_bracketed_paste', 'disable_bracketed_paste'):
            return getattr(self.win32_output, name)
        else:
            return getattr(self.vt100_output, name)


Output.register(ConEmuOutput)
