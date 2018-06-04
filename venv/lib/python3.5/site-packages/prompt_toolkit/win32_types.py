from ctypes import Union, Structure, c_char, c_short, c_long, c_ulong
from ctypes.wintypes import DWORD, BOOL, LPVOID, WORD, WCHAR


# Input/Output standard device numbers. Note that these are not handle objects.
# It's the `windll.kernel32.GetStdHandle` system call that turns them into a
# real handle object.
STD_INPUT_HANDLE = c_ulong(-10)
STD_OUTPUT_HANDLE = c_ulong(-11)
STD_ERROR_HANDLE = c_ulong(-12)


class COORD(Structure):
    """
    Struct in wincon.h
    http://msdn.microsoft.com/en-us/library/windows/desktop/ms682119(v=vs.85).aspx
    """
    _fields_ = [
        ('X', c_short),  # Short
        ('Y', c_short),  # Short
    ]

    def __repr__(self):
        return '%s(X=%r, Y=%r, type_x=%r, type_y=%r)' % (
            self.__class__.__name__, self.X, self.Y, type(self.X), type(self.Y))


class UNICODE_OR_ASCII(Union):
    _fields_ = [
        ('AsciiChar', c_char),
        ('UnicodeChar', WCHAR),
    ]


class KEY_EVENT_RECORD(Structure):
    """
    http://msdn.microsoft.com/en-us/library/windows/desktop/ms684166(v=vs.85).aspx
    """
    _fields_ = [
        ('KeyDown', c_long),  # bool
        ('RepeatCount', c_short),  # word
        ('VirtualKeyCode', c_short),  # word
        ('VirtualScanCode', c_short),  # word
        ('uChar', UNICODE_OR_ASCII),  # Unicode or ASCII.
        ('ControlKeyState', c_long)  # double word
    ]


class MOUSE_EVENT_RECORD(Structure):
    """
    http://msdn.microsoft.com/en-us/library/windows/desktop/ms684239(v=vs.85).aspx
    """
    _fields_ = [
        ('MousePosition', COORD),
        ('ButtonState', c_long),  # dword
        ('ControlKeyState', c_long),  # dword
        ('EventFlags', c_long)  # dword
    ]


class WINDOW_BUFFER_SIZE_RECORD(Structure):
    """
    http://msdn.microsoft.com/en-us/library/windows/desktop/ms687093(v=vs.85).aspx
    """
    _fields_ = [
        ('Size', COORD)
    ]


class MENU_EVENT_RECORD(Structure):
    """
    http://msdn.microsoft.com/en-us/library/windows/desktop/ms684213(v=vs.85).aspx
    """
    _fields_ = [
        ('CommandId', c_long)  # uint
    ]


class FOCUS_EVENT_RECORD(Structure):
    """
    http://msdn.microsoft.com/en-us/library/windows/desktop/ms683149(v=vs.85).aspx
    """
    _fields_ = [
        ('SetFocus', c_long)   # bool
    ]


class EVENT_RECORD(Union):
    _fields_ = [
        ('KeyEvent', KEY_EVENT_RECORD),
        ('MouseEvent', MOUSE_EVENT_RECORD),
        ('WindowBufferSizeEvent', WINDOW_BUFFER_SIZE_RECORD),
        ('MenuEvent', MENU_EVENT_RECORD),
        ('FocusEvent', FOCUS_EVENT_RECORD)
    ]


class INPUT_RECORD(Structure):
    """
    http://msdn.microsoft.com/en-us/library/windows/desktop/ms683499(v=vs.85).aspx
    """
    _fields_ = [
        ('EventType', c_short),  # word
        ('Event', EVENT_RECORD)  # Union.
    ]


EventTypes = {
    1: 'KeyEvent',
    2: 'MouseEvent',
    4: 'WindowBufferSizeEvent',
    8: 'MenuEvent',
    16: 'FocusEvent'
}


class SMALL_RECT(Structure):
    """struct in wincon.h."""
    _fields_ = [
        ("Left", c_short),
        ("Top", c_short),
        ("Right", c_short),
        ("Bottom", c_short),
    ]


class CONSOLE_SCREEN_BUFFER_INFO(Structure):
    """struct in wincon.h."""
    _fields_ = [
        ("dwSize", COORD),
        ("dwCursorPosition", COORD),
        ("wAttributes", WORD),
        ("srWindow", SMALL_RECT),
        ("dwMaximumWindowSize", COORD),
    ]

    def __str__(self):
        return '(%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d)' % (
            self.dwSize.Y, self.dwSize.X,
            self.dwCursorPosition.Y, self.dwCursorPosition.X,
            self.wAttributes,
            self.srWindow.Top, self.srWindow.Left, self.srWindow.Bottom, self.srWindow.Right,
            self.dwMaximumWindowSize.Y, self.dwMaximumWindowSize.X,
        )


class SECURITY_ATTRIBUTES(Structure):
    """
    http://msdn.microsoft.com/en-us/library/windows/desktop/aa379560(v=vs.85).aspx
    """
    _fields_ = [
        ('nLength', DWORD),
        ('lpSecurityDescriptor', LPVOID),
        ('bInheritHandle', BOOL),
    ]
