"""
Command line layout definitions
-------------------------------

The layout of a command line interface is defined by a Container instance.
There are two main groups of classes here. Containers and controls:

- A container can contain other containers or controls, it can have multiple
  children and it decides about the dimensions.
- A control is responsible for rendering the actual content to a screen.
  A control can propose some dimensions, but it's the container who decides
  about the dimensions -- or when the control consumes more space -- which part
  of the control will be visible.


Container classes::

    - Container (Abstract base class)
       |- HSplit (Horizontal split)
       |- VSplit (Vertical split)
       |- FloatContainer (Container which can also contain menus and other floats)
       `- Window (Container which contains one actual control

Control classes::

    - UIControl (Abstract base class)
       |- TokenListControl (Renders a simple list of tokens)
       |- FillControl (Fills control with one token/character.)
       `- BufferControl (Renders an input buffer.)


Usually, you end up wrapping every control inside a `Window` object, because
that's the only way to render it in a layout.

There are some prepared toolbars which are ready to use::

- SystemToolbar (Shows the 'system' input buffer, for entering system commands.)
- ArgToolbar (Shows the input 'arg', for repetition of input commands.)
- SearchToolbar (Shows the 'search' input buffer, for incremental search.)
- CompletionsToolbar (Shows the completions of the current buffer.)
- ValidationToolbar (Shows validation errors of the current buffer.)

And one prepared menu:

- CompletionsMenu

"""
from __future__ import unicode_literals

from .containers import Float, FloatContainer, HSplit, VSplit, Window, ConditionalContainer
from .controls import TokenListControl, FillControl, BufferControl
