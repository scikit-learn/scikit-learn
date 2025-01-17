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
    comment = "#141414"
    red = "#9F4E55"
    orange = "#A25E53"
    yellow = "#98661B"
    green = "#437A6B"
    blue = "#3D73A9"
    purple = "#974EB7"
    black = "#141414"


class Theme(Style):
    """
    This style mimics the gotthard light theme from vscode.
    """

    default_style = ""

    background_color = "#F5F5F5"
    highlight_color = "#E1E1E1"

    styles = {
        Text: Colors.black,  # class:  ''
        Error: Colors.red,  # class: 'err'
        Other: "",  # class 'x'
        Comment: Colors.purple,  # class: 'c'
        Keyword: Colors.purple,  # class: 'k'
        Keyword.Constant: Colors.red,  # class: 'kc'
        # Keyword.Declaration:       "",            # class: 'kd'
        # Keyword.Namespace:         "",            # class: 'kn'
        # Keyword.Pseudo:            "",            # class: 'kp'
        # Keyword.Reserved:          "",            # class: 'kr'
        Keyword.Type: Colors.green,  # class: 'kt'
        Operator: Colors.blue,  # class: 'o'
        Operator.Word: Colors.purple,  # class: 'ow'
        Punctuation: Colors.black,  # class: 'p'
        Name: Colors.black,  # class: 'n'
        Name.Attribute: Colors.purple,  # class: 'na'
        Name.Builtin: Colors.green,  # class: 'nb'
        Name.Builtin.Pseudo: Colors.green,  # class: 'bp'
        Name.Class: Colors.yellow,  # class: 'nc'
        Name.Constant: Colors.red,  # class: 'no'
        Name.Decorator: Colors.green,  # class: 'nd'
        Name.Entity: Colors.green,  # class: 'ni'
        Name.Exception: Colors.red,  # class: 'ne'
        Name.Function: Colors.purple,  # class: 'nf'
        Name.Property: Colors.purple,  # class: 'py'
        Name.Label: Colors.green,  # class: 'nl'
        Name.Namespace: Colors.yellow,  # class: 'nn'
        # Name.Other:                "",            # class: 'nx'
        Name.Tag: Colors.red,  # class: 'nt'
        Name.Variable: Colors.comment,  # class: 'nv'
        Name.Variable.Magic: Colors.comment,
        # Name.Variable.Class:       "",            # class: 'vc'
        # Name.Variable.Global:      "",            # class: 'vg'
        # Name.Variable.Instance:    "",            # class: 'vi'
        Number: Colors.red,  # class: 'm'
        # Number.Float:              "",            # class: 'mf'
        # Number.Hex:                "",            # class: 'mh'
        # Number.Integer:            "",            # class: 'mi'
        # Number.Integer.Long:       "",            # class: 'il'
        # Number.Oct:                "",            # class: 'mo'
        Literal: Colors.purple,  # class: 'l'
        # Literal.Date:              "",            # class: 'ld'
        String: Colors.green,  # class: 's'
        String.Backtick: Colors.yellow,  # class: 'sb'
        # String.Char:               "",            # class: 'sc'
        # String.Doc:                "",            # class: 'sd'
        # String.Double:             "",            # class: 's2'
        String.Escape: Colors.blue,  # class: 'se'
        # String.Heredoc:            "",            # class: 'sh'
        # String.Interpol:           "",            # class: 'si'
        # String.Other:              "",            # class: 'sx'
        String.Regex: Colors.blue,  # class: 'sr'
        # String.Single:             "",            # class: 's1'
        String.Symbol: Colors.green,  # class: 'ss'
        # Generic:                   "",            # class: 'g'
        Generic.Deleted: Colors.red,  # class: 'gd',
        # Generic.Emph:              "italic",      # class: 'ge'
        # Generic.Error:             "",            # class: 'gr'
        Generic.Heading: Colors.green,  # class: 'gh'
        Generic.Subheading: Colors.green,  # class: 'gu'
        # Generic.Inserted:          "",            # class: 'gi'
        # Generic.Output:            "",            # class: 'go'
        # Generic.Prompt:            "",            # class: 'gp'
        Generic.Strong: "bold",  # class: 'gs'
        # Generic.Traceback:         "",            # class: 'gt'
    }
