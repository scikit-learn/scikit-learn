import sys

__all__ = [
    "AMPER",
    "AMPEREQUAL",
    "ASYNC",
    "AT",
    "ATEQUAL",
    "AWAIT",
    "CIRCUMFLEX",
    "CIRCUMFLEXEQUAL",
    "COLON",
    "COLONEQUAL",
    "COMMA",
    "DEDENT",
    "DOT",
    "DOUBLESLASH",
    "DOUBLESLASHEQUAL",
    "DOUBLESTAR",
    "DOUBLESTAREQUAL",
    "ELLIPSIS",
    "ENDMARKER",
    "EQEQUAL",
    "EQUAL",
    "ERRORTOKEN",
    "GREATER",
    "GREATEREQUAL",
    "INDENT",
    "ISEOF",
    "ISNONTERMINAL",
    "ISTERMINAL",
    "LBRACE",
    "LEFTSHIFT",
    "LEFTSHIFTEQUAL",
    "LESS",
    "LESSEQUAL",
    "LPAR",
    "LSQB",
    "MINEQUAL",
    "MINUS",
    "NAME",
    "NEWLINE",
    "NOTEQUAL",
    "NT_OFFSET",
    "NUMBER",
    "N_TOKENS",
    "OP",
    "PERCENT",
    "PERCENTEQUAL",
    "PLUS",
    "PLUSEQUAL",
    "RARROW",
    "RBRACE",
    "RIGHTSHIFT",
    "RIGHTSHIFTEQUAL",
    "RPAR",
    "RSQB",
    "SEMI",
    "SLASH",
    "SLASHEQUAL",
    "STAR",
    "STAREQUAL",
    "STRING",
    "TILDE",
    "TYPE_COMMENT",
    "TYPE_IGNORE",
    "VBAR",
    "VBAREQUAL",
    "tok_name",
    "ENCODING",
    "NL",
    "COMMENT",
]

if sys.version_info >= (3, 10):
    __all__ += ["SOFT_KEYWORD"]

if sys.version_info >= (3, 12):
    __all__ += ["EXCLAMATION", "FSTRING_END", "FSTRING_MIDDLE", "FSTRING_START"]

ENDMARKER: int
NAME: int
NUMBER: int
STRING: int
NEWLINE: int
INDENT: int
DEDENT: int
LPAR: int
RPAR: int
LSQB: int
RSQB: int
COLON: int
COMMA: int
SEMI: int
PLUS: int
MINUS: int
STAR: int
SLASH: int
VBAR: int
AMPER: int
LESS: int
GREATER: int
EQUAL: int
DOT: int
PERCENT: int
LBRACE: int
RBRACE: int
EQEQUAL: int
NOTEQUAL: int
LESSEQUAL: int
GREATEREQUAL: int
TILDE: int
CIRCUMFLEX: int
LEFTSHIFT: int
RIGHTSHIFT: int
DOUBLESTAR: int
PLUSEQUAL: int
MINEQUAL: int
STAREQUAL: int
SLASHEQUAL: int
PERCENTEQUAL: int
AMPEREQUAL: int
VBAREQUAL: int
CIRCUMFLEXEQUAL: int
LEFTSHIFTEQUAL: int
RIGHTSHIFTEQUAL: int
DOUBLESTAREQUAL: int
DOUBLESLASH: int
DOUBLESLASHEQUAL: int
AT: int
RARROW: int
ELLIPSIS: int
ATEQUAL: int
AWAIT: int
ASYNC: int
OP: int
ERRORTOKEN: int
N_TOKENS: int
NT_OFFSET: int
tok_name: dict[int, str]
COMMENT: int
NL: int
ENCODING: int
TYPE_COMMENT: int
TYPE_IGNORE: int
COLONEQUAL: int
EXACT_TOKEN_TYPES: dict[str, int]
if sys.version_info >= (3, 10):
    SOFT_KEYWORD: int

if sys.version_info >= (3, 12):
    EXCLAMATION: int
    FSTRING_END: int
    FSTRING_MIDDLE: int
    FSTRING_START: int

def ISTERMINAL(x: int) -> bool: ...
def ISNONTERMINAL(x: int) -> bool: ...
def ISEOF(x: int) -> bool: ...
