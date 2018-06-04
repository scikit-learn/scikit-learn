"""
Tools to open .py files as Unicode, using the encoding specified within the file,
as per PEP 263.

Much of the code is taken from the tokenize module in Python 3.2.
"""

import io
from io import TextIOWrapper, BytesIO
import re
from tokenize import open, detect_encoding

cookie_re = re.compile(r"coding[:=]\s*([-\w.]+)", re.UNICODE)
cookie_comment_re = re.compile(r"^\s*#.*coding[:=]\s*([-\w.]+)", re.UNICODE)

def source_to_unicode(txt, errors='replace', skip_encoding_cookie=True):
    """Converts a bytes string with python source code to unicode.

    Unicode strings are passed through unchanged. Byte strings are checked
    for the python source file encoding cookie to determine encoding.
    txt can be either a bytes buffer or a string containing the source
    code.
    """
    if isinstance(txt, str):
        return txt
    if isinstance(txt, bytes):
        buffer = BytesIO(txt)
    else:
        buffer = txt
    try:
        encoding, _ = detect_encoding(buffer.readline)
    except SyntaxError:
        encoding = "ascii"
    buffer.seek(0)
    text = TextIOWrapper(buffer, encoding, errors=errors, line_buffering=True)
    text.mode = 'r'
    if skip_encoding_cookie:
        return u"".join(strip_encoding_cookie(text))
    else:
        return text.read()

def strip_encoding_cookie(filelike):
    """Generator to pull lines from a text-mode file, skipping the encoding
    cookie if it is found in the first two lines.
    """
    it = iter(filelike)
    try:
        first = next(it)
        if not cookie_comment_re.match(first):
            yield first
        second = next(it)
        if not cookie_comment_re.match(second):
            yield second
    except StopIteration:
        return
    
    for line in it:
        yield line

def read_py_file(filename, skip_encoding_cookie=True):
    """Read a Python file, using the encoding declared inside the file.
    
    Parameters
    ----------
    filename : str
      The path to the file to read.
    skip_encoding_cookie : bool
      If True (the default), and the encoding declaration is found in the first
      two lines, that line will be excluded from the output - compiling a
      unicode string with an encoding declaration is a SyntaxError in Python 2.
    
    Returns
    -------
    A unicode string containing the contents of the file.
    """
    with open(filename) as f:   # the open function defined in this module.
        if skip_encoding_cookie:
            return "".join(strip_encoding_cookie(f))
        else:
            return f.read()

def read_py_url(url, errors='replace', skip_encoding_cookie=True):
    """Read a Python file from a URL, using the encoding declared inside the file.
    
    Parameters
    ----------
    url : str
      The URL from which to fetch the file.
    errors : str
      How to handle decoding errors in the file. Options are the same as for
      bytes.decode(), but here 'replace' is the default.
    skip_encoding_cookie : bool
      If True (the default), and the encoding declaration is found in the first
      two lines, that line will be excluded from the output - compiling a
      unicode string with an encoding declaration is a SyntaxError in Python 2.
    
    Returns
    -------
    A unicode string containing the contents of the file.
    """
    # Deferred import for faster start
    from urllib.request import urlopen 
    response = urlopen(url)
    buffer = io.BytesIO(response.read())
    return source_to_unicode(buffer, errors, skip_encoding_cookie)

def _list_readline(x):
    """Given a list, returns a readline() function that returns the next element
    with each call.
    """
    x = iter(x)
    def readline():
        return next(x)
    return readline

# Code for going between .py files and cached .pyc files ----------------------
try: 
    from importlib.util import source_from_cache, cache_from_source
except ImportError :
    ## deprecated since 3.4
    from imp import source_from_cache, cache_from_source
