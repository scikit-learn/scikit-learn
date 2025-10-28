from pygments.style import Style
from pygments.token import (
    Comment,
    Error,
    Generic,
    Keyword,
    Literal,
    Name,
    Number,
    Operator,
    Other,
    Punctuation,
    String,
    Text,
)


class Colors:
    comment = "#d9dee3"  # grey[2]
    red = "#ff9492"  # red[3]
    orange = "#ffb757"  # orange[2]
    green = "#72f088"  # green[1]
    blue = "#91cbff"  # blue[2]
    purple = "#dbb7ff"  # purple[2]
    black = "#C9D1D9"  # fg.default


class Theme(Style):
    """
    This style mimics the github dark high contrast theme from vs code themes.
    """

    default_style = ""

    background_color = "#0d1117"  # canvas.default
    highlight_color = "#58a6ff70"  # accent.fg

    styles = {
        Text: Colors.black,  # class:  ''
        Error: Colors.red,  # class: 'err'
        Other: "",  # class 'x'
        Comment: Colors.comment,  # class: 'c'
        Keyword: Colors.red,  # class: 'k'
        Keyword.Constant: Colors.blue,  # class: 'kc'
        # Keyword.Declaration:       "",            # class: 'kd'
        # Keyword.Namespace:         "",            # class: 'kn'
        # Keyword.Pseudo:            "",            # class: 'kp'
        # Keyword.Reserved:          "",            # class: 'kr'
        Keyword.Type: Colors.red,  # class: 'kt'
        Operator: Colors.green,  # class: 'o'
        Operator.Word: Colors.purple,  # class: 'ow'
        Punctuation: Colors.black,  # class: 'p'
        Name: Colors.purple,  # class: 'n'
        Name.Attribute: Colors.orange,  # class: 'na'
        Name.Builtin: Colors.orange,  # class: 'nb'
        Name.Builtin.Pseudo: Colors.orange,  # class: 'bp'
        Name.Class: Colors.blue,  # class: 'nc'
        Name.Constant: Colors.blue,  # class: 'no'
        Name.Decorator: Colors.orange,  # class: 'nd'
        Name.Entity: Colors.green,  # class: 'ni'
        Name.Exception: Colors.purple,  # class: 'ne'
        Name.Function: Colors.blue,  # class: 'nf'
        Name.Property: Colors.blue,  # class: 'py'
        Name.Label: Colors.orange,  # class: 'nl'
        Name.Namespace: Colors.black,  # class: 'nn'
        # Name.Other:                "",            # class: 'nx'
        Name.Tag: Colors.green,  # class: 'nt'
        Name.Variable: Colors.orange,  # class: 'nv'
        Name.Variable.Magic: Colors.orange,
        # Name.Variable.Class:       "",            # class: 'vc'
        # Name.Variable.Global:      "",            # class: 'vg'
        # Name.Variable.Instance:    "",            # class: 'vi'
        Number: Colors.orange,  # class: 'm'
        # Number.Float:              "",            # class: 'mf'
        # Number.Hex:                "",            # class: 'mh'
        # Number.Integer:            "",            # class: 'mi'
        # Number.Integer.Long:       "",            # class: 'il'
        # Number.Oct:                "",            # class: 'mo'
        Literal: Colors.orange,  # class: 'l'
        # Literal.Date:              "",            # class: 'ld'
        String: Colors.blue,  # class: 's'
        String.Backtick: Colors.blue,  # class: 'sb'
        # String.Char:               "",            # class: 'sc'
        # String.Doc:                "",            # class: 'sd'
        # String.Double:             "",            # class: 's2'
        # String.Escape:             "",            # class: 'se'
        # String.Heredoc:            "",            # class: 'sh'
        # String.Interpol:           "",            # class: 'si'
        # String.Other:              "",            # class: 'sx'
        String.Regex: Colors.blue,  # class: 'sr'
        # String.Single:             "",            # class: 's1'
        String.Symbol: Colors.blue,  # class: 'ss'
        # Generic:                   "",            # class: 'g'
        Generic.Deleted: Colors.blue,  # class: 'gd',
        Generic.Emph: "italic",  # class: 'ge'
        Generic.Error: Colors.red,  # class: 'gr'
        Generic.Heading: Colors.blue,  # class: 'gh'
        Generic.Subheading: Colors.blue,  # class: 'gu'
        # Generic.Inserted:          "",            # class: 'gi'
        # Generic.Output:            "",            # class: 'go'
        # Generic.Prompt:            "",            # class: 'gp'
        Generic.Strong: "bold",  # class: 'gs'
        # Generic.Traceback:         "",            # class: 'gt'
    }
