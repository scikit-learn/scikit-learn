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
    comment = "#515151"
    red = "#d71835"
    orange = "#7f4707"
    yellow = "#7f4707"
    green = "#116633"
    blue = "#00749c"
    purple = "#8045e5"
    black = "#1e1e1e"


class Theme(Style):
    """
    This style inspired by the a11y-light theme from eric bailey's accessible themes.
    """

    default_style = ""

    background_color = "#f2f2f2"
    highlight_color = "#fdf2e2"

    styles = {
        Text: Colors.black,  # class:  ''
        Error: Colors.red,  # class: 'err'
        Other: "",  # class 'x'
        Comment: Colors.comment,  # class: 'c'
        Keyword: Colors.purple,  # class: 'k'
        Keyword.Constant: Colors.purple,  # class: 'kc'
        # Keyword.Declaration:       "",            # class: 'kd'
        # Keyword.Namespace:         "",            # class: 'kn'
        # Keyword.Pseudo:            "",            # class: 'kp'
        # Keyword.Reserved:          "",            # class: 'kr'
        Keyword.Type: Colors.orange,  # class: 'kt'
        Operator: Colors.green,  # class: 'o'
        Operator.Word: Colors.purple,  # class: 'ow'
        Punctuation: Colors.black,  # class: 'p'
        Name: Colors.black,  # class: 'n'
        Name.Attribute: Colors.yellow,  # class: 'na'
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
        Name.Tag: Colors.blue,  # class: 'nt'
        Name.Variable: Colors.red,  # class: 'nv'
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
        String: Colors.green,  # class: 's'
        String.Backtick: Colors.green,  # class: 'sb'
        # String.Char:               "",            # class: 'sc'
        # String.Doc:                "",            # class: 'sd'
        # String.Double:             "",            # class: 's2'
        # String.Escape:             "",            # class: 'se'
        # String.Heredoc:            "",            # class: 'sh'
        # String.Interpol:           "",            # class: 'si'
        # String.Other:              "",            # class: 'sx'
        String.Regex: Colors.red,  # class: 'sr'
        # String.Single:             "",            # class: 's1'
        String.Symbol: Colors.blue,  # class: 'ss'
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
