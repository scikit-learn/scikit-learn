from collections import namedtuple
import re
import sys
from ast import literal_eval

from parso._compatibility import unicode, total_ordering


Version = namedtuple('Version', 'major, minor, micro')


def split_lines(string, keepends=False):
    r"""
    Intended for Python code. In contrast to Python's :py:meth:`str.splitlines`,
    looks at form feeds and other special characters as normal text. Just
    splits ``\n`` and ``\r\n``.
    Also different: Returns ``[""]`` for an empty string input.

    In Python 2.7 form feeds are used as normal characters when using
    str.splitlines. However in Python 3 somewhere there was a decision to split
    also on form feeds.
    """
    if keepends:
        lst = string.splitlines(True)

        # We have to merge lines that were broken by form feed characters.
        merge = []
        for i, line in enumerate(lst):
            if line.endswith('\f'):
                merge.append(i)

        for index in reversed(merge):
            try:
                lst[index] = lst[index] + lst[index + 1]
                del lst[index + 1]
            except IndexError:
                # index + 1 can be empty and therefore there's no need to
                # merge.
                pass

        # The stdlib's implementation of the end is inconsistent when calling
        # it with/without keepends. One time there's an empty string in the
        # end, one time there's none.
        if string.endswith('\n') or string == '':
            lst.append('')
        return lst
    else:
        return re.split('\n|\r\n', string)


def python_bytes_to_unicode(source, encoding='utf-8', errors='strict'):
    """
    Checks for unicode BOMs and PEP 263 encoding declarations. Then returns a
    unicode object like in :py:meth:`bytes.decode`.

    :param encoding: See :py:meth:`bytes.decode` documentation.
    :param errors: See :py:meth:`bytes.decode` documentation. ``errors`` can be
        ``'strict'``, ``'replace'`` or ``'ignore'``.
    """
    def detect_encoding():
        """
        For the implementation of encoding definitions in Python, look at:
        - http://www.python.org/dev/peps/pep-0263/
        - http://docs.python.org/2/reference/lexical_analysis.html#encoding-declarations
        """
        byte_mark = literal_eval(r"b'\xef\xbb\xbf'")
        if source.startswith(byte_mark):
            # UTF-8 byte-order mark
            return 'utf-8'

        first_two_lines = re.match(br'(?:[^\n]*\n){0,2}', source).group(0)
        possible_encoding = re.search(br"coding[=:]\s*([-\w.]+)",
                                      first_two_lines)
        if possible_encoding:
            return possible_encoding.group(1)
        else:
            # the default if nothing else has been set -> PEP 263
            return encoding

    if isinstance(source, unicode):
        # only cast str/bytes
        return source

    encoding = detect_encoding()
    if not isinstance(encoding, unicode):
        encoding = unicode(encoding, 'utf-8', 'replace')

    # Cast to unicode
    return unicode(source, encoding, errors)


def version_info():
    """
    Returns a namedtuple of parso's version, similar to Python's
    ``sys.version_info``.
    """
    from parso import __version__
    tupl = re.findall(r'[a-z]+|\d+', __version__)
    return Version(*[x if i == 3 else int(x) for i, x in enumerate(tupl)])


def _parse_version(version):
    match = re.match(r'(\d+)(?:\.(\d)(?:\.\d+)?)?$', version)
    if match is None:
        raise ValueError('The given version is not in the right format. '
                         'Use something like "3.2" or "3".')

    major = int(match.group(1))
    minor = match.group(2)
    if minor is None:
        # Use the latest Python in case it's not exactly defined, because the
        # grammars are typically backwards compatible?
        if major == 2:
            minor = "7"
        elif major == 3:
            minor = "6"
        else:
            raise NotImplementedError("Sorry, no support yet for those fancy new/old versions.")
    minor = int(minor)
    return PythonVersionInfo(major, minor)


@total_ordering
class PythonVersionInfo(namedtuple('Version', 'major, minor')):
    def __gt__(self, other):
        if isinstance(other, tuple):
            if len(other) != 2:
                raise ValueError("Can only compare to tuples of length 2.")
            return (self.major, self.minor) > other
        super(PythonVersionInfo, self).__gt__(other)

        return (self.major, self.minor)

    def __eq__(self, other):
        if isinstance(other, tuple):
            if len(other) != 2:
                raise ValueError("Can only compare to tuples of length 2.")
            return (self.major, self.minor) == other
        super(PythonVersionInfo, self).__eq__(other)

    def __ne__(self, other):
        return not self.__eq__(other)


def parse_version_string(version=None):
    """
    Checks for a valid version number (e.g. `3.2` or `2.7.1` or `3`) and
    returns a corresponding version info that is always two characters long in
    decimal.
    """
    if version is None:
        version = '%s.%s' % sys.version_info[:2]
    if not isinstance(version, (unicode, str)):
        raise TypeError("version must be a string like 3.2.")

    return _parse_version(version)
