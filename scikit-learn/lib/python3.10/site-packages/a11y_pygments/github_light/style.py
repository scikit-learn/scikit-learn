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
    comment = "#6e7781"  # grey[5]
    red = "#cf222e"  # red[5]
    orange = "#953800"  # orange[6]
    green = "#116329"  # green[6]
    blue = "#0550ae"  # blue[6]
    purple = "#8250df"  # purple[5]
    black = "#24292f"  # fg.default


class Theme(Style):
    """
    This style mimics the github light theme from vscode themes.
    """

    default_style = ""

    background_color = "#ffffff"  # canvas.default
    highlight_color = "#0969da4a"  # accent.fg

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
