"""
Important note on ids
----------------------------------------------------------------------------

Multiple id generation schemes are used due to backwards compatibility.
- v1: 1.2.3 <= version < 1.3
      The style used before the rewrite.
      It is not the actual old code, but a replication of the behaviour.
- v2: 1.3 <= version < now
      Standardised mangling scheme from
      https://itanium-cxx-abi.github.io/cxx-abi/abi.html#mangling
      though not completely implemented.
All versions are generated and attached to elements. The newest is used for
the index. All of the versions should work as permalinks.


Signature Nodes and Tagnames
----------------------------------------------------------------------------

Each signature is in a desc_signature node, where all children are
desc_signature_line nodes. Each of these lines will have the attribute
'sphinx_line_type' set to one of the following (prioritized):
- 'declarator', if the line contains the name of the declared object.
- 'templateParams', if the line starts a template parameter list,
- 'templateParams', if the line has template parameters
  Note: such lines might get a new tag in the future.
- 'templateIntroduction, if the line is on the form 'conceptName{...}'
No other desc_signature nodes should exist (so far).


Grammar
----------------------------------------------------------------------------

See https://www.nongnu.org/hcb/ for the grammar,
and https://github.com/cplusplus/draft/blob/master/source/grammar.tex,
and https://github.com/cplusplus/concepts-ts
for the newest grammar.

common grammar things:
    template-declaration ->
        "template" "<" template-parameter-list ">" declaration
    template-parameter-list ->
          template-parameter
        | template-parameter-list "," template-parameter
    template-parameter ->
          type-parameter
        | parameter-declaration # i.e., same as a function argument

    type-parameter ->
          "class"    "..."[opt] identifier[opt]
        | "class"               identifier[opt] "=" type-id
        | "typename" "..."[opt] identifier[opt]
        | "typename"            identifier[opt] "=" type-id
        | "template" "<" template-parameter-list ">"
            "class"  "..."[opt] identifier[opt]
        | "template" "<" template-parameter-list ">"
            "class"             identifier[opt] "=" id-expression
        # also, from C++17 we can have "typename" in template templates
    templateDeclPrefix ->
        "template" "<" template-parameter-list ">"

    simple-declaration ->
        attribute-specifier-seq[opt] decl-specifier-seq[opt]
            init-declarator-list[opt] ;
    # Make the semicolon optional.
    # For now: drop the attributes (TODO).
    # Use at most 1 init-declarator.
    -> decl-specifier-seq init-declarator
    -> decl-specifier-seq declarator initializer

    decl-specifier ->
          storage-class-specifier ->
             (  "static" (only for member_object and function_object)
              | "extern" (only for member_object and function_object)
              | "register"
             )
             thread_local[opt] (only for member_object)
                               (it can also appear before the others)

        | type-specifier -> trailing-type-specifier
        | function-specifier -> "inline" | "virtual" | "explicit" (only
          for function_object)
        | "friend" (only for function_object)
        | "constexpr" (only for member_object and function_object)
    trailing-type-specifier ->
          simple-type-specifier
        | elaborated-type-specifier
        | typename-specifier
        | cv-qualifier -> "const" | "volatile"
    stricter grammar for decl-specifier-seq (with everything, each object
    uses a subset):
        visibility storage-class-specifier function-specifier "friend"
        "constexpr" "volatile" "const" trailing-type-specifier
        # where trailing-type-specifier can no be cv-qualifier
    # Inside e.g., template parameters a strict subset is used
    # (see type-specifier-seq)
    trailing-type-specifier ->
          simple-type-specifier ->
            ::[opt] nested-name-specifier[opt] type-name
          | ::[opt] nested-name-specifier "template" simple-template-id
          | "char" | "bool" | etc.
          | decltype-specifier
        | elaborated-type-specifier ->
            class-key attribute-specifier-seq[opt] ::[opt]
            nested-name-specifier[opt] identifier
          | class-key ::[opt] nested-name-specifier[opt] template[opt]
            simple-template-id
          | "enum" ::[opt] nested-name-specifier[opt] identifier
        | typename-specifier ->
            "typename" ::[opt] nested-name-specifier identifier
          | "typename" ::[opt] nested-name-specifier template[opt]
            simple-template-id
    class-key -> "class" | "struct" | "union"
    type-name ->* identifier | simple-template-id
    # ignoring attributes and decltype, and then some left-factoring
    trailing-type-specifier ->
        rest-of-trailing
        ("class" | "struct" | "union" | "typename") rest-of-trailing
        built-in -> "char" | "bool" | etc.
        decltype-specifier
    rest-of-trailing -> (with some simplification)
        "::"[opt] list-of-elements-separated-by-::
    element ->
        "template"[opt] identifier ("<" template-argument-list ">")[opt]
    template-argument-list ->
          template-argument "..."[opt]
        | template-argument-list "," template-argument "..."[opt]
    template-argument ->
          constant-expression
        | type-specifier-seq abstract-declarator
        | id-expression


    declarator ->
          ptr-declarator
        | noptr-declarator parameters-and-qualifiers trailing-return-type
    ptr-declarator ->
          noptr-declarator
        | ptr-operator ptr-declarator
    noptr-declarator ->
          declarator-id attribute-specifier-seq[opt] ->
                "..."[opt] id-expression
              | rest-of-trailing
        | noptr-declarator parameters-and-qualifiers
        | noptr-declarator "[" constant-expression[opt] "]"
          attribute-specifier-seq[opt]
        | "(" ptr-declarator ")"
    ptr-operator ->
          "*"  attribute-specifier-seq[opt] cv-qualifier-seq[opt]
        | "&   attribute-specifier-seq[opt]
        | "&&" attribute-specifier-seq[opt]
        | "::"[opt] nested-name-specifier "*" attribute-specifier-seq[opt]
            cv-qualifier-seq[opt]
    # function_object must use a parameters-and-qualifiers, the others may
    # use it (e.g., function pointers)
    parameters-and-qualifiers ->
        "(" parameter-clause ")" attribute-specifier-seq[opt]
        cv-qualifier-seq[opt] ref-qualifier[opt]
        exception-specification[opt]
    ref-qualifier -> "&" | "&&"
    exception-specification ->
        "noexcept" ("(" constant-expression ")")[opt]
        "throw" ("(" type-id-list ")")[opt]
    # TODO: we don't implement attributes
    # member functions can have initializers, but we fold them into here
    memberFunctionInit -> "=" "0"
    # (note: only "0" is allowed as the value, according to the standard,
    # right?)

    enum-head ->
        enum-key attribute-specifier-seq[opt] nested-name-specifier[opt]
            identifier enum-base[opt]
    enum-key -> "enum" | "enum struct" | "enum class"
    enum-base ->
        ":" type
    enumerator-definition ->
          identifier
        | identifier "=" constant-expression

We additionally add the possibility for specifying the visibility as the
first thing.

concept_object:
    goal:
        just a declaration of the name (for now)

    grammar: only a single template parameter list, and the nested name
        may not have any template argument lists

        "template" "<" template-parameter-list ">"
        nested-name-specifier

type_object:
    goal:
        either a single type (e.g., "MyClass:Something_T" or a typedef-like
        thing (e.g. "Something Something_T" or "int I_arr[]"
    grammar, single type: based on a type in a function parameter, but
    without a name:
           parameter-declaration
        -> attribute-specifier-seq[opt] decl-specifier-seq
           abstract-declarator[opt]
        # Drop the attributes
        -> decl-specifier-seq abstract-declarator[opt]
    grammar, typedef-like: no initializer
        decl-specifier-seq declarator
    Can start with a templateDeclPrefix.

member_object:
    goal: as a type_object which must have a declarator, and optionally
    with a initializer
    grammar:
        decl-specifier-seq declarator initializer
    Can start with a templateDeclPrefix.

function_object:
    goal: a function declaration, TODO: what about templates? for now: skip
    grammar: no initializer
       decl-specifier-seq declarator
    Can start with a templateDeclPrefix.

class_object:
    goal: a class declaration, but with specification of a base class
    grammar:
          attribute-specifier-seq[opt]
              nested-name "final"[opt] (":" base-specifier-list)[opt]
        base-specifier-list ->
          base-specifier "..."[opt]
        | base-specifier-list, base-specifier "..."[opt]
        base-specifier ->
          base-type-specifier
        | "virtual" access-spe"cifier[opt]    base-type-specifier
        | access-specifier[opt] "virtual"[opt] base-type-specifier
    Can start with a templateDeclPrefix.

enum_object:
    goal: an unscoped enum or a scoped enum, optionally with the underlying
          type specified
    grammar:
        ("class" | "struct")[opt] visibility[opt]
            attribute-specifier-seq[opt] nested-name (":" type)[opt]
enumerator_object:
    goal: an element in a scoped or unscoped enum. The name should be
          injected according to the scopedness.
    grammar:
        nested-name ("=" constant-expression)

namespace_object:
    goal: a directive to put all following declarations in a specific scope
    grammar:
        nested-name
"""

from __future__ import annotations

import re

udl_identifier_re = re.compile(r'''
    [a-zA-Z_][a-zA-Z0-9_]*\b   # note, no word boundary in the beginning
''', re.VERBOSE)
_string_re = re.compile(r"[LuU8]?('([^'\\]*(?:\\.[^'\\]*)*)'"
                        r'|"([^"\\]*(?:\\.[^"\\]*)*)")', re.DOTALL)
_visibility_re = re.compile(r'\b(public|private|protected)\b')
_operator_re = re.compile(r'''
        \[\s*\]
    |   \(\s*\)
    |   \+\+ | --
    |   ->\*? | \,
    |   (<<|>>)=? | && | \|\|
    |   <=>
    |   [!<>=/*%+|&^~-]=?
    |   (\b(and|and_eq|bitand|bitor|compl|not|not_eq|or|or_eq|xor|xor_eq)\b)
''', re.VERBOSE)
_fold_operator_re = re.compile(r'''
        ->\*    |    \.\*    |    \,
    |   (<<|>>)=?    |    &&    |    \|\|
    |   !=
    |   [<>=/*%+|&^~-]=?
''', re.VERBOSE)
# see https://en.cppreference.com/w/cpp/keyword
_keywords = [
    'alignas', 'alignof', 'and', 'and_eq', 'asm', 'auto', 'bitand', 'bitor',
    'bool', 'break', 'case', 'catch', 'char', 'char8_t', 'char16_t', 'char32_t',
    'class', 'compl', 'concept', 'const', 'consteval', 'constexpr', 'constinit',
    'const_cast', 'continue',
    'decltype', 'default', 'delete', 'do', 'double', 'dynamic_cast', 'else',
    'enum', 'explicit', 'export', 'extern', 'false', 'float', 'for', 'friend',
    'goto', 'if', 'inline', 'int', 'long', 'mutable', 'namespace', 'new',
    'noexcept', 'not', 'not_eq', 'nullptr', 'operator', 'or', 'or_eq',
    'private', 'protected', 'public', 'register', 'reinterpret_cast',
    'requires', 'return', 'short', 'signed', 'sizeof', 'static',
    'static_assert', 'static_cast', 'struct', 'switch', 'template', 'this',
    'thread_local', 'throw', 'true', 'try', 'typedef', 'typeid', 'typename',
    'union', 'unsigned', 'using', 'virtual', 'void', 'volatile', 'wchar_t',
    'while', 'xor', 'xor_eq',
]


_simple_type_specifiers_re = re.compile(r"""
    \b(
    auto|void|bool
    |signed|unsigned
    |short|long
    |char|wchar_t|char(8|16|32)_t
    |int
    |__int(64|128)  # extension
    |float|double
    |__float80|_Float64x|__float128|_Float128  # extension
    |_Complex|_Imaginary  # extension
    )\b
""", re.VERBOSE)

_max_id = 4
_id_prefix = [None, '', '_CPPv2', '_CPPv3', '_CPPv4']
# Ids are used in lookup keys which are used across pickled files,
# so when _max_id changes, make sure to update the ENV_VERSION.

# ------------------------------------------------------------------------------
# Id v1 constants
# ------------------------------------------------------------------------------

_id_fundamental_v1 = {
    'char': 'c',
    'signed char': 'c',
    'unsigned char': 'C',
    'int': 'i',
    'signed int': 'i',
    'unsigned int': 'U',
    'long': 'l',
    'signed long': 'l',
    'unsigned long': 'L',
    'bool': 'b',
}
_id_shorthands_v1 = {
    'std::string': 'ss',
    'std::ostream': 'os',
    'std::istream': 'is',
    'std::iostream': 'ios',
    'std::vector': 'v',
    'std::map': 'm',
}
_id_operator_v1 = {
    'new': 'new-operator',
    'new[]': 'new-array-operator',
    'delete': 'delete-operator',
    'delete[]': 'delete-array-operator',
    # the arguments will make the difference between unary and binary
    # '+(unary)' : 'ps',
    # '-(unary)' : 'ng',
    # '&(unary)' : 'ad',
    # '*(unary)' : 'de',
    '~': 'inv-operator',
    '+': 'add-operator',
    '-': 'sub-operator',
    '*': 'mul-operator',
    '/': 'div-operator',
    '%': 'mod-operator',
    '&': 'and-operator',
    '|': 'or-operator',
    '^': 'xor-operator',
    '=': 'assign-operator',
    '+=': 'add-assign-operator',
    '-=': 'sub-assign-operator',
    '*=': 'mul-assign-operator',
    '/=': 'div-assign-operator',
    '%=': 'mod-assign-operator',
    '&=': 'and-assign-operator',
    '|=': 'or-assign-operator',
    '^=': 'xor-assign-operator',
    '<<': 'lshift-operator',
    '>>': 'rshift-operator',
    '<<=': 'lshift-assign-operator',
    '>>=': 'rshift-assign-operator',
    '==': 'eq-operator',
    '!=': 'neq-operator',
    '<': 'lt-operator',
    '>': 'gt-operator',
    '<=': 'lte-operator',
    '>=': 'gte-operator',
    '!': 'not-operator',
    '&&': 'sand-operator',
    '||': 'sor-operator',
    '++': 'inc-operator',
    '--': 'dec-operator',
    ',': 'comma-operator',
    '->*': 'pointer-by-pointer-operator',
    '->': 'pointer-operator',
    '()': 'call-operator',
    '[]': 'subscript-operator',
}

# ------------------------------------------------------------------------------
# Id v > 1 constants
# ------------------------------------------------------------------------------

_id_fundamental_v2 = {
    # not all of these are actually parsed as fundamental types, TODO: do that
    'void': 'v',
    'bool': 'b',
    'char': 'c',
    'signed char': 'a',
    'unsigned char': 'h',
    'wchar_t': 'w',
    'char32_t': 'Di',
    'char16_t': 'Ds',
    'char8_t': 'Du',
    'short': 's',
    'short int': 's',
    'signed short': 's',
    'signed short int': 's',
    'unsigned short': 't',
    'unsigned short int': 't',
    'int': 'i',
    'signed': 'i',
    'signed int': 'i',
    'unsigned': 'j',
    'unsigned int': 'j',
    'long': 'l',
    'long int': 'l',
    'signed long': 'l',
    'signed long int': 'l',
    'unsigned long': 'm',
    'unsigned long int': 'm',
    'long long': 'x',
    'long long int': 'x',
    'signed long long': 'x',
    'signed long long int': 'x',
    '__int64': 'x',
    'unsigned long long': 'y',
    'unsigned long long int': 'y',
    '__int128': 'n',
    'signed __int128': 'n',
    'unsigned __int128': 'o',
    'float': 'f',
    'double': 'd',
    'long double': 'e',
    '__float80': 'e', '_Float64x': 'e',
    '__float128': 'g', '_Float128': 'g',
    '_Complex float': 'Cf',
    '_Complex double': 'Cd',
    '_Complex long double': 'Ce',
    '_Imaginary float': 'f',
    '_Imaginary double': 'd',
    '_Imaginary long double': 'e',
    'auto': 'Da',
    'decltype(auto)': 'Dc',
    'std::nullptr_t': 'Dn',
}
_id_operator_v2 = {
    'new': 'nw',
    'new[]': 'na',
    'delete': 'dl',
    'delete[]': 'da',
    # the arguments will make the difference between unary and binary
    # in operator definitions
    # '+(unary)' : 'ps',
    # '-(unary)' : 'ng',
    # '&(unary)' : 'ad',
    # '*(unary)' : 'de',
    '~': 'co', 'compl': 'co',
    '+': 'pl',
    '-': 'mi',
    '*': 'ml',
    '/': 'dv',
    '%': 'rm',
    '&': 'an', 'bitand': 'an',
    '|': 'or', 'bitor': 'or',
    '^': 'eo', 'xor': 'eo',
    '=': 'aS',
    '+=': 'pL',
    '-=': 'mI',
    '*=': 'mL',
    '/=': 'dV',
    '%=': 'rM',
    '&=': 'aN', 'and_eq': 'aN',
    '|=': 'oR', 'or_eq': 'oR',
    '^=': 'eO', 'xor_eq': 'eO',
    '<<': 'ls',
    '>>': 'rs',
    '<<=': 'lS',
    '>>=': 'rS',
    '==': 'eq',
    '!=': 'ne', 'not_eq': 'ne',
    '<': 'lt',
    '>': 'gt',
    '<=': 'le',
    '>=': 'ge',
    '<=>': 'ss',
    '!': 'nt', 'not': 'nt',
    '&&': 'aa', 'and': 'aa',
    '||': 'oo', 'or': 'oo',
    '++': 'pp',
    '--': 'mm',
    ',': 'cm',
    '->*': 'pm',
    '->': 'pt',
    '()': 'cl',
    '[]': 'ix',
    '.*': 'ds',  # this one is not overloadable, but we need it for expressions
    '?': 'qu',
}
_id_operator_unary_v2 = {
    '++': 'pp_',
    '--': 'mm_',
    '*': 'de',
    '&': 'ad',
    '+': 'ps',
    '-': 'ng',
    '!': 'nt', 'not': 'nt',
    '~': 'co', 'compl': 'co',
}
_id_char_from_prefix: dict[str | None, str] = {
    None: 'c', 'u8': 'c',
    'u': 'Ds', 'U': 'Di', 'L': 'w',
}
# these are ordered by preceedence
_expression_bin_ops = [
    ['||', 'or'],
    ['&&', 'and'],
    ['|', 'bitor'],
    ['^', 'xor'],
    ['&', 'bitand'],
    ['==', '!=', 'not_eq'],
    ['<=>', '<=', '>=', '<', '>'],
    ['<<', '>>'],
    ['+', '-'],
    ['*', '/', '%'],
    ['.*', '->*'],
]
_expression_unary_ops = ["++", "--", "*", "&", "+", "-", "!", "not", "~", "compl"]
_expression_assignment_ops = ["=", "*=", "/=", "%=", "+=", "-=",
                              ">>=", "<<=", "&=", "and_eq", "^=", "|=", "xor_eq", "or_eq"]
_id_explicit_cast = {
    'dynamic_cast': 'dc',
    'static_cast': 'sc',
    'const_cast': 'cc',
    'reinterpret_cast': 'rc',
}
