from __future__ import unicode_literals

__all__ = (
    'Key',
    'Keys',
)


class Key(object):
    def __init__(self, name):

        #: Descriptive way of writing keys in configuration files. e.g. <C-A>
        #: for ``Control-A``.
        self.name = name

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.name)


class Keys(object):
    Escape = Key('<Escape>')

    ControlA = Key('<C-A>')
    ControlB = Key('<C-B>')
    ControlC = Key('<C-C>')
    ControlD = Key('<C-D>')
    ControlE = Key('<C-E>')
    ControlF = Key('<C-F>')
    ControlG = Key('<C-G>')
    ControlH = Key('<C-H>')
    ControlI = Key('<C-I>')  # Tab
    ControlJ = Key('<C-J>')  # Enter
    ControlK = Key('<C-K>')
    ControlL = Key('<C-L>')
    ControlM = Key('<C-M>')  # Enter
    ControlN = Key('<C-N>')
    ControlO = Key('<C-O>')
    ControlP = Key('<C-P>')
    ControlQ = Key('<C-Q>')
    ControlR = Key('<C-R>')
    ControlS = Key('<C-S>')
    ControlT = Key('<C-T>')
    ControlU = Key('<C-U>')
    ControlV = Key('<C-V>')
    ControlW = Key('<C-W>')
    ControlX = Key('<C-X>')
    ControlY = Key('<C-Y>')
    ControlZ = Key('<C-Z>')

    ControlSpace       = Key('<C-Space>')
    ControlBackslash   = Key('<C-Backslash>')
    ControlSquareClose = Key('<C-SquareClose>')
    ControlCircumflex  = Key('<C-Circumflex>')
    ControlUnderscore  = Key('<C-Underscore>')
    ControlLeft        = Key('<C-Left>')
    ControlRight       = Key('<C-Right>')
    ControlUp          = Key('<C-Up>')
    ControlDown        = Key('<C-Down>')

    Up          = Key('<Up>')
    Down        = Key('<Down>')
    Right       = Key('<Right>')
    Left        = Key('<Left>')

    ShiftLeft   = Key('<ShiftLeft>')
    ShiftUp     = Key('<ShiftUp>')
    ShiftDown   = Key('<ShiftDown>')
    ShiftRight  = Key('<ShiftRight>')

    Home        = Key('<Home>')
    End         = Key('<End>')
    Delete      = Key('<Delete>')
    ShiftDelete = Key('<ShiftDelete>')
    ControlDelete = Key('<C-Delete>')
    PageUp      = Key('<PageUp>')
    PageDown    = Key('<PageDown>')
    BackTab     = Key('<BackTab>')  # shift + tab
    Insert      = Key('<Insert>')
    Backspace   = Key('<Backspace>')

    # Aliases.
    Tab         = ControlI
    Enter       = ControlJ
        # XXX: Actually Enter equals ControlM, not ControlJ,
        #      However, in prompt_toolkit, we made the mistake of translating
        #      \r into \n during the input, so everyone is now handling the
        #      enter key by binding ControlJ.

        #      From now on, it's better to bind `Keys.Enter` everywhere,
        #      because that's future compatible, and will still work when we
        #      stop replacing \r by \n.

    F1 = Key('<F1>')
    F2 = Key('<F2>')
    F3 = Key('<F3>')
    F4 = Key('<F4>')
    F5 = Key('<F5>')
    F6 = Key('<F6>')
    F7 = Key('<F7>')
    F8 = Key('<F8>')
    F9 = Key('<F9>')
    F10 = Key('<F10>')
    F11 = Key('<F11>')
    F12 = Key('<F12>')
    F13 = Key('<F13>')
    F14 = Key('<F14>')
    F15 = Key('<F15>')
    F16 = Key('<F16>')
    F17 = Key('<F17>')
    F18 = Key('<F18>')
    F19 = Key('<F19>')
    F20 = Key('<F20>')
    F21 = Key('<F21>')
    F22 = Key('<F22>')
    F23 = Key('<F23>')
    F24 = Key('<F24>')

    # Matches any key.
    Any = Key('<Any>')

    # Special
    CPRResponse = Key('<Cursor-Position-Response>')
    Vt100MouseEvent = Key('<Vt100-Mouse-Event>')
    WindowsMouseEvent = Key('<Windows-Mouse-Event>')
    BracketedPaste = Key('<Bracketed-Paste>')

    # Key which is ignored. (The key binding for this key should not do
    # anything.)
    Ignore = Key('<Ignore>')
