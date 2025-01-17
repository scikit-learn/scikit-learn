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
    comment = "#797979"
    orange = "#F78C6C"
    yellow = "#9e8741"
    green = "#c5e478"
    blue = "#a2bffc"
    cyan = "#5ca7e4"
    purple = "#9e86c8"
    black = "#a2bffc"


class Theme(Style):
    """
    This style mimics greative theme from vscode themes.
    """

    default_style = ""

    background_color = "#010726"
    highlight_color = "#473d18"

    styles = {
        Text: Colors.black,  # class:  ''
        Error: Colors.blue,  # class: 'err'
        Other: "",  # class 'x'
        Comment: Colors.comment,  # class: 'c'
        Keyword: Colors.purple,  # class: 'k'
        Keyword.Constant: Colors.purple,  # class: 'kc'
        # Keyword.Declaration:       "",            # class: 'kd'
        # Keyword.Namespace:         "",            # class: 'kn'
        # Keyword.Pseudo:            "",            # class: 'kp'
        # Keyword.Reserved:          "",            # class: 'kr'
        Keyword.Type: Colors.green,  # class: 'kt'
        Operator: Colors.orange,  # class: 'o'
        Operator.Word: Colors.purple,  # class: 'ow'
        Punctuation: Colors.black,  # class: 'p'
        Name: Colors.blue,  # class: 'n'
        Name.Attribute: Colors.purple,  # class: 'na'
        Name.Builtin: Colors.green,  # class: 'nb'
        Name.Builtin.Pseudo: Colors.green,  # class: 'bp'
        Name.Class: Colors.orange,  # class: 'nc'
        Name.Constant: Colors.blue,  # class: 'no'
        Name.Decorator: Colors.yellow,  # class: 'nd'
        Name.Entity: Colors.yellow,  # class: 'ni'
        Name.Exception: Colors.blue,  # class: 'ne'
        Name.Function: Colors.purple,  # class: 'nf'
        Name.Property: Colors.blue,  # class: 'py'
        Name.Label: Colors.orange,  # class: 'nl'
        Name.Namespace: Colors.green,  # class: 'nn'
        # Name.Other:                "",            # class: 'nx'
        Name.Tag: Colors.cyan,  # class: 'nt'
        Name.Variable: Colors.green,  # class: 'nv'
        Name.Variable.Magic: Colors.orange,
        # Name.Variable.Class:       "",            # class: 'vc'
        # Name.Variable.Global:      "",            # class: 'vg'
        # Name.Variable.Instance:    "",            # class: 'vi'
        Number: Colors.purple,  # class: 'm'
        # Number.Float:              "",            # class: 'mf'
        # Number.Hex:                "",            # class: 'mh'
        # Number.Integer:            "",            # class: 'mi'
        # Number.Integer.Long:       "",            # class: 'il'
        # Number.Oct:                "",            # class: 'mo'
        Literal: Colors.blue,  # class: 'l'
        # Literal.Date:              "",            # class: 'ld'
        String: Colors.yellow,  # class: 's'
        # String.Backtick:           Colors.green,   # class: 'sb'
        # String.Char:               "",            # class: 'sc'
        # String.Doc:                "",            # class: 'sd'
        # String.Double:             "",            # class: 's2'
        String.Escape: Colors.orange,  # class: 'se'
        # String.Heredoc:            "",            # class: 'sh'
        # String.Interpol:           "",            # class: 'si'
        # String.Other:              "",            # class: 'sx'
        String.Regex: Colors.cyan,  # class: 'sr'
        # String.Single:             "",            # class: 's1'
        String.Symbol: Colors.orange,  # class: 'ss'
        # Generic:                   "",            # class: 'g'
        Generic.Deleted: Colors.blue,  # class: 'gd',
        Generic.Emph: "italic",  # class: 'ge'
        # Generic.Error:             "",            # class: 'gr'
        Generic.Heading: Colors.blue,  # class: 'gh'
        Generic.Subheading: Colors.blue,  # class: 'gu'
        # Generic.Inserted:          "",            # class: 'gi'
        # Generic.Output:            "",            # class: 'go'
        # Generic.Prompt:            "",            # class: 'gp'
        Generic.Strong: "bold",  # class: 'gs'
        # Generic.Traceback:         "",            # class: 'gt'
    }
