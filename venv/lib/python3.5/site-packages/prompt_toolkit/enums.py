from __future__ import unicode_literals


class IncrementalSearchDirection(object):
    FORWARD = 'FORWARD'
    BACKWARD = 'BACKWARD'


class EditingMode(object):
    # The set of key bindings that is active.
    VI = 'VI'
    EMACS = 'EMACS'


#: Name of the search buffer.
SEARCH_BUFFER = 'SEARCH_BUFFER'

#: Name of the default buffer.
DEFAULT_BUFFER = 'DEFAULT_BUFFER'

#: Name of the system buffer.
SYSTEM_BUFFER = 'SYSTEM_BUFFER'

# Dummy buffer. This is the buffer returned by
# `CommandLineInterface.current_buffer` when the top of the `FocusStack` is
# `None`. This could be the case when there is some widget has the focus and no
# actual text editing is possible. This buffer should also never be displayed.
# (It will never contain any actual text.)
DUMMY_BUFFER = 'DUMMY_BUFFER'
