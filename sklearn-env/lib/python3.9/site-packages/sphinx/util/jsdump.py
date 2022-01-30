"""
    sphinx.util.jsdump
    ~~~~~~~~~~~~~~~~~~

    This module implements a simple JavaScript serializer.
    Uses the basestring encode function from simplejson by Bob Ippolito.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import re
from typing import IO, Any, Dict, List, Match, Union

_str_re = re.compile(r'"(\\\\|\\"|[^"])*"')
_int_re = re.compile(r'\d+')
_name_re = re.compile(r'[a-zA-Z_]\w*')
_nameonly_re = re.compile(r'[a-zA-Z_][a-zA-Z0-9_]*$')

# escape \, ", control characters and everything outside ASCII
ESCAPE_ASCII = re.compile(r'([\\"]|[^\ -~])')
ESCAPE_DICT = {
    '\\': '\\\\',
    '"': '\\"',
    '\b': '\\b',
    '\f': '\\f',
    '\n': '\\n',
    '\r': '\\r',
    '\t': '\\t',
}

ESCAPED = re.compile(r'\\u.{4}|\\.')


def encode_string(s: str) -> str:
    def replace(match: Match) -> str:
        s = match.group(0)
        try:
            return ESCAPE_DICT[s]
        except KeyError:
            n = ord(s)
            if n < 0x10000:
                return '\\u%04x' % (n,)
            else:
                # surrogate pair
                n -= 0x10000
                s1 = 0xd800 | ((n >> 10) & 0x3ff)
                s2 = 0xdc00 | (n & 0x3ff)
                return '\\u%04x\\u%04x' % (s1, s2)
    return '"' + str(ESCAPE_ASCII.sub(replace, s)) + '"'


def decode_string(s: str) -> str:
    return ESCAPED.sub(lambda m: eval('"' + m.group() + '"'), s)


reswords = set("""\
abstract   else   instanceof   switch
boolean   enum   int   synchronized
break   export   interface   this
byte   extends   long   throw
case   false   native   throws
catch   final   new   transient
char   finally   null   true
class   float   package   try
const   for   private   typeof
continue   function   protected   var
debugger   goto   public   void
default   if   return   volatile
delete   implements   short   while
do   import   static   with
double   in   super""".split())


def dumps(obj: Any, key: bool = False) -> str:
    if key:
        if not isinstance(obj, str):
            obj = str(obj)
        if _nameonly_re.match(obj) and obj not in reswords:
            return obj  # return it as a bare word
        else:
            return encode_string(obj)
    if obj is None:
        return 'null'
    elif obj is True or obj is False:
        return 'true' if obj else 'false'
    elif isinstance(obj, (int, float)):
        return str(obj)
    elif isinstance(obj, dict):
        return '{%s}' % ','.join(sorted('%s:%s' % (
            dumps(key, True),
            dumps(value)
        ) for key, value in obj.items()))
    elif isinstance(obj, set):
        return '[%s]' % ','.join(sorted(dumps(x) for x in obj))
    elif isinstance(obj, (tuple, list)):
        return '[%s]' % ','.join(dumps(x) for x in obj)
    elif isinstance(obj, str):
        return encode_string(obj)
    raise TypeError(type(obj))


def dump(obj: Any, f: IO) -> None:
    f.write(dumps(obj))


def loads(x: str) -> Any:
    """Loader that can read the JS subset the indexer produces."""
    nothing = object()
    i = 0
    n = len(x)
    stack: List[Union[List, Dict]] = []
    obj: Any = nothing
    key = False
    keys = []
    while i < n:
        c = x[i]
        if c == '{':
            obj = {}
            stack.append(obj)
            key = True
            keys.append(nothing)
            i += 1
        elif c == '[':
            obj = []
            stack.append(obj)
            key = False
            keys.append(nothing)
            i += 1
        elif c in '}]':
            if key:
                if keys[-1] is not nothing:
                    raise ValueError("unfinished dict")
                # empty dict
                key = False
            oldobj = stack.pop()
            keys.pop()
            if stack:
                obj = stack[-1]
                if isinstance(obj, dict):
                    if keys[-1] is nothing:
                        raise ValueError("invalid key object", oldobj)
                    obj[keys[-1]] = oldobj
                else:
                    obj.append(oldobj)
            else:
                break
            i += 1
        elif c == ',':
            if key:
                raise ValueError("multiple keys")
            if isinstance(obj, dict):
                key = True
            i += 1
        elif c == ':':
            if not isinstance(obj, dict):
                raise ValueError("colon in list")
            i += 1
            if not key:
                raise ValueError("multiple values")
            key = False
        else:
            y: Any = None
            m = _str_re.match(x, i)
            if m:
                y = decode_string(m.group()[1:-1])
            else:
                m = _int_re.match(x, i)
                if m:
                    y = int(m.group())
                else:
                    m = _name_re.match(x, i)
                    if m:
                        y = m.group()
                        if y == 'true':
                            y = True
                        elif y == 'false':
                            y = False
                        elif y == 'null':
                            y = None
                        elif not key:
                            raise ValueError("bareword as value")
                    else:
                        raise ValueError("read error at pos %d" % i)
            i = m.end()
            if isinstance(obj, dict):
                if key:
                    keys[-1] = y
                else:
                    obj[keys[-1]] = y
                    key = False
            else:
                obj.append(y)
    if obj is nothing:
        raise ValueError("nothing loaded from string")
    return obj


def load(f: IO) -> Any:
    return loads(f.read())
