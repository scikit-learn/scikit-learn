import base64
import decimal
import json as _json
import sys
import re
from functools import reduce

from _plotly_utils.optional_imports import get_module
from _plotly_utils.basevalidators import (
    ImageUriValidator,
    copy_to_readonly_numpy_array,
    is_homogeneous_array,
)


int8min = -128
int8max = 127
int16min = -32768
int16max = 32767
int32min = -2147483648
int32max = 2147483647

uint8max = 255
uint16max = 65535
uint32max = 4294967295

plotlyjsShortTypes = {
    "int8": "i1",
    "uint8": "u1",
    "int16": "i2",
    "uint16": "u2",
    "int32": "i4",
    "uint32": "u4",
    "float32": "f4",
    "float64": "f8",
}


def to_typed_array_spec(v):
    """
    Convert numpy array to plotly.js typed array spec
    If not possible return the original value
    """
    v = copy_to_readonly_numpy_array(v)

    # Skip b64 encoding if numpy is not installed,
    # or if v is not a numpy array, or if v is empty
    np = get_module("numpy", should_load=False)
    if not np or not isinstance(v, np.ndarray) or v.size == 0:
        return v

    dtype = str(v.dtype)

    # convert default Big Ints until we could support them in plotly.js
    if dtype == "int64":
        max = v.max()
        min = v.min()
        if max <= int8max and min >= int8min:
            v = v.astype("int8")
        elif max <= int16max and min >= int16min:
            v = v.astype("int16")
        elif max <= int32max and min >= int32min:
            v = v.astype("int32")
        else:
            return v

    elif dtype == "uint64":
        max = v.max()
        min = v.min()
        if max <= uint8max and min >= 0:
            v = v.astype("uint8")
        elif max <= uint16max and min >= 0:
            v = v.astype("uint16")
        elif max <= uint32max and min >= 0:
            v = v.astype("uint32")
        else:
            return v

    dtype = str(v.dtype)

    if dtype in plotlyjsShortTypes:
        arrObj = {
            "dtype": plotlyjsShortTypes[dtype],
            "bdata": base64.b64encode(v).decode("ascii"),
        }

        if v.ndim > 1:
            arrObj["shape"] = str(v.shape)[1:-1]

        return arrObj

    return v


def is_skipped_key(key):
    """
    Return whether the key is skipped for conversion to the typed array spec
    """
    skipped_keys = ["geojson", "layer", "layers", "range"]
    return any(skipped_key == key for skipped_key in skipped_keys)


def convert_to_base64(obj):
    if isinstance(obj, dict):
        for key, value in obj.items():
            if is_skipped_key(key):
                continue
            elif is_homogeneous_array(value):
                obj[key] = to_typed_array_spec(value)
            else:
                convert_to_base64(value)
    elif isinstance(obj, list) or isinstance(obj, tuple):
        for value in obj:
            convert_to_base64(value)


def cumsum(x):
    """
    Custom cumsum to avoid a numpy import.
    """

    def _reducer(a, x):
        if len(a) == 0:
            return [x]
        return a + [a[-1] + x]

    ret = reduce(_reducer, x, [])
    return ret


class PlotlyJSONEncoder(_json.JSONEncoder):
    """
    Meant to be passed as the `cls` kwarg to json.dumps(obj, cls=..)

    See PlotlyJSONEncoder.default for more implementation information.

    Additionally, this encoder overrides nan functionality so that 'Inf',
    'NaN' and '-Inf' encode to 'null'. Which is stricter JSON than the Python
    version.

    """

    def coerce_to_strict(self, const):
        """
        This is used to ultimately *encode* into strict JSON, see `encode`

        """
        # before python 2.7, 'true', 'false', 'null', were include here.
        if const in ("Infinity", "-Infinity", "NaN"):
            return None
        else:
            return const

    def encode(self, o):
        """
        Load and then dump the result using parse_constant kwarg

        Note that setting invalid separators will cause a failure at this step.

        """
        # this will raise errors in a normal-expected way
        encoded_o = super(PlotlyJSONEncoder, self).encode(o)
        # Brute force guessing whether NaN or Infinity values are in the string
        # We catch false positive cases (e.g. strings such as titles, labels etc.)
        # but this is ok since the intention is to skip the decoding / reencoding
        # step when it's completely safe

        if not ("NaN" in encoded_o or "Infinity" in encoded_o):
            return encoded_o

        # now:
        #    1. `loads` to switch Infinity, -Infinity, NaN to None
        #    2. `dumps` again so you get 'null' instead of extended JSON
        try:
            new_o = _json.loads(encoded_o, parse_constant=self.coerce_to_strict)
        except ValueError:
            # invalid separators will fail here. raise a helpful exception
            raise ValueError(
                "Encoding into strict JSON failed. Did you set the separators "
                "valid JSON separators?"
            )
        else:
            return _json.dumps(
                new_o,
                sort_keys=self.sort_keys,
                indent=self.indent,
                separators=(self.item_separator, self.key_separator),
            )

    def default(self, obj):
        """
        Accept an object (of unknown type) and try to encode with priority:
        1. builtin:     user-defined objects
        2. sage:        sage math cloud
        3. pandas:      dataframes/series
        4. numpy:       ndarrays
        5. datetime:    time/datetime objects

        Each method throws a NotEncoded exception if it fails.

        The default method will only get hit if the object is not a type that
        is naturally encoded by json:

            Normal objects:
                dict                object
                list, tuple         array
                str, unicode        string
                int, long, float    number
                True                true
                False               false
                None                null

            Extended objects:
                float('nan')        'NaN'
                float('infinity')   'Infinity'
                float('-infinity')  '-Infinity'

        Therefore, we only anticipate either unknown iterables or values here.

        """
        # TODO: The ordering if these methods is *very* important. Is this OK?
        encoding_methods = (
            self.encode_as_plotly,
            self.encode_as_sage,
            self.encode_as_numpy,
            self.encode_as_pandas,
            self.encode_as_datetime,
            self.encode_as_date,
            self.encode_as_list,  # because some values have `tolist` do last.
            self.encode_as_decimal,
            self.encode_as_pil,
        )
        for encoding_method in encoding_methods:
            try:
                return encoding_method(obj)
            except NotEncodable:
                pass
        return _json.JSONEncoder.default(self, obj)

    @staticmethod
    def encode_as_plotly(obj):
        """Attempt to use a builtin `to_plotly_json` method."""
        try:
            return obj.to_plotly_json()
        except AttributeError:
            raise NotEncodable

    @staticmethod
    def encode_as_list(obj):
        """Attempt to use `tolist` method to convert to normal Python list."""
        if hasattr(obj, "tolist"):
            return obj.tolist()
        else:
            raise NotEncodable

    @staticmethod
    def encode_as_sage(obj):
        """Attempt to convert sage.all.RR to floats and sage.all.ZZ to ints"""
        sage_all = get_module("sage.all")
        if not sage_all:
            raise NotEncodable

        if obj in sage_all.RR:
            return float(obj)
        elif obj in sage_all.ZZ:
            return int(obj)
        else:
            raise NotEncodable

    @staticmethod
    def encode_as_pandas(obj):
        """Attempt to convert pandas.NaT / pandas.NA"""
        pandas = get_module("pandas", should_load=False)
        if not pandas:
            raise NotEncodable

        if obj is pandas.NaT:
            return None

        # pandas.NA was introduced in pandas 1.0
        if hasattr(pandas, "NA") and obj is pandas.NA:
            return None

        raise NotEncodable

    @staticmethod
    def encode_as_numpy(obj):
        """Attempt to convert numpy.ma.core.masked"""
        numpy = get_module("numpy", should_load=False)
        if not numpy:
            raise NotEncodable

        if obj is numpy.ma.core.masked:
            return float("nan")
        elif isinstance(obj, numpy.ndarray) and obj.dtype.kind == "M":
            try:
                return numpy.datetime_as_string(obj).tolist()
            except TypeError:
                pass

        raise NotEncodable

    @staticmethod
    def encode_as_datetime(obj):
        """Convert datetime objects to iso-format strings"""
        try:
            return obj.isoformat()
        except AttributeError:
            raise NotEncodable

    @staticmethod
    def encode_as_date(obj):
        """Attempt to convert to utc-iso time string using date methods."""
        try:
            time_string = obj.isoformat()
        except AttributeError:
            raise NotEncodable
        else:
            return iso_to_plotly_time_string(time_string)

    @staticmethod
    def encode_as_decimal(obj):
        """Attempt to encode decimal by converting it to float"""
        if isinstance(obj, decimal.Decimal):
            return float(obj)
        else:
            raise NotEncodable

    @staticmethod
    def encode_as_pil(obj):
        """Attempt to convert PIL.Image.Image to base64 data uri"""
        image = get_module("PIL.Image")
        if image is not None and isinstance(obj, image.Image):
            return ImageUriValidator.pil_image_to_uri(obj)
        else:
            raise NotEncodable


class NotEncodable(Exception):
    pass


def iso_to_plotly_time_string(iso_string):
    """Remove timezone info and replace 'T' delimeter with ' ' (ws)."""
    # make sure we don't send timezone info to plotly
    if (iso_string.split("-")[:3] == "00:00") or (iso_string.split("+")[0] == "00:00"):
        raise Exception(
            "Plotly won't accept timestrings with timezone info.\n"
            "All timestrings are assumed to be in UTC."
        )

    iso_string = iso_string.replace("-00:00", "").replace("+00:00", "")

    if iso_string.endswith("T00:00:00"):
        return iso_string.replace("T00:00:00", "")
    else:
        return iso_string.replace("T", " ")


def template_doc(**names):
    def _decorator(func):
        if not sys.version_info[:2] == (3, 2):
            if func.__doc__ is not None:
                func.__doc__ = func.__doc__.format(**names)
        return func

    return _decorator


def _natural_sort_strings(vals, reverse=False):
    def key(v):
        v_parts = re.split(r"(\d+)", v)
        for i in range(len(v_parts)):
            try:
                v_parts[i] = int(v_parts[i])
            except ValueError:
                # not an int
                pass
        return tuple(v_parts)

    return sorted(vals, key=key, reverse=reverse)


def _get_int_type():
    np = get_module("numpy", should_load=False)
    if np:
        int_type = (int, np.integer)
    else:
        int_type = (int,)
    return int_type


def split_multichar(ss, chars):
    """
    Split all the strings in ss at any of the characters in chars.
    Example:

        >>> ss = ["a.string[0].with_separators"]
        >>> chars = list(".[]_")
        >>> split_multichar(ss, chars)
        ['a', 'string', '0', '', 'with', 'separators']

    :param (list) ss: A list of strings.
    :param (list) chars: Is a list of chars (note: not a string).
    """
    if len(chars) == 0:
        return ss
    c = chars.pop()
    ss = reduce(lambda x, y: x + y, map(lambda x: x.split(c), ss))
    return split_multichar(ss, chars)


def split_string_positions(ss):
    """
    Given a list of strings split using split_multichar, return a list of
    integers representing the indices of the first character of every string in
    the original string.
    Example:

        >>> ss = ["a.string[0].with_separators"]
        >>> chars = list(".[]_")
        >>> ss_split = split_multichar(ss, chars)
        >>> ss_split
        ['a', 'string', '0', '', 'with', 'separators']
        >>> split_string_positions(ss_split)
        [0, 2, 9, 11, 12, 17]

    :param (list) ss: A list of strings.
    """
    return list(
        map(
            lambda t: t[0] + t[1],
            zip(range(len(ss)), cumsum([0] + list(map(len, ss[:-1])))),
        )
    )


def display_string_positions(p, i=None, offset=0, length=1, char="^", trim=True):
    """
    Return a string that is whitespace except at p[i] which is replaced with char.
    If i is None then all the indices of the string in p are replaced with char.

    Example:

        >>> ss = ["a.string[0].with_separators"]
        >>> chars = list(".[]_")
        >>> ss_split = split_multichar(ss, chars)
        >>> ss_split
        ['a', 'string', '0', '', 'with', 'separators']
        >>> ss_pos = split_string_positions(ss_split)
        >>> ss[0]
        'a.string[0].with_separators'
        >>> display_string_positions(ss_pos,4)
        '            ^'
        >>> display_string_positions(ss_pos,4,offset=1,length=3,char="~",trim=False)
        '             ~~~      '
        >>> display_string_positions(ss_pos)
        '^ ^      ^ ^^    ^'
    :param (list) p: A list of integers.
    :param (integer|None) i: Optional index of p to display.
    :param (integer) offset: Allows adding a number of spaces to the replacement.
    :param (integer) length: Allows adding a replacement that is the char
                             repeated length times.
    :param (str) char: allows customizing the replacement character.
    :param (boolean) trim: trims the remaining whitespace if True.
    """
    s = [" " for _ in range(max(p) + 1 + offset + length)]
    maxaddr = 0
    if i is None:
        for p_ in p:
            for temp in range(length):
                maxaddr = p_ + offset + temp
                s[maxaddr] = char
    else:
        for temp in range(length):
            maxaddr = p[i] + offset + temp
            s[maxaddr] = char
    ret = "".join(s)
    if trim:
        ret = ret[: maxaddr + 1]
    return ret


def chomp_empty_strings(strings, c, reverse=False):
    """
    Given a list of strings, some of which are the empty string "", replace the
    empty strings with c and combine them with the closest non-empty string on
    the left or "" if it is the first string.
    Examples:
    for c="_"
    ['hey', '', 'why', '', '', 'whoa', '', ''] -> ['hey_', 'why__', 'whoa__']
    ['', 'hi', '', "I'm", 'bob', '', ''] -> ['_', 'hi_', "I'm", 'bob__']
    ['hi', "i'm", 'a', 'good', 'string'] -> ['hi', "i'm", 'a', 'good', 'string']
    Some special cases are:
    [] -> []
    [''] -> ['']
    ['', ''] -> ['_']
    ['', '', '', ''] -> ['___']
    If reverse is true, empty strings are combined with closest non-empty string
    on the right or "" if it is the last string.
    """

    def _rev(vals):
        return [s[::-1] for s in vals][::-1]

    if reverse:
        return _rev(chomp_empty_strings(_rev(strings), c))
    if not len(strings):
        return strings
    if sum(map(len, strings)) == 0:
        return [c * (len(strings) - 1)]

    class _Chomper:
        def __init__(self, c):
            self.c = c

        def __call__(self, x, y):
            # x is list up to now
            # y is next item in list
            # x should be [""] initially, and then empty strings filtered out at the
            # end
            if len(y) == 0:
                return x[:-1] + [x[-1] + self.c]
            else:
                return x + [y]

    return list(filter(len, reduce(_Chomper(c), strings, [""])))


# taken from
# https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
def levenshtein(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)  # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            # j+1 instead of j since previous_row and current_row are one character longer
            # than s2
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def find_closest_string(string, strings):
    def _key(s):
        # sort by levenshtein distance and lexographically to maintain a stable
        # sort for different keys with the same levenshtein distance
        return (levenshtein(s, string), s)

    return sorted(strings, key=_key)[0]
