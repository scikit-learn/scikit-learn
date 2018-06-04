# Last Change: Mon Aug 20 08:00 PM 2007 J
from __future__ import division, print_function, absolute_import

import re
import itertools
import datetime
from functools import partial

import numpy as np

from scipy._lib.six import next

"""A module to read arff files."""

__all__ = ['MetaData', 'loadarff', 'ArffError', 'ParseArffError']

# An Arff file is basically two parts:
#   - header
#   - data
#
# A header has each of its components starting by @META where META is one of
# the keyword (attribute of relation, for now).

# TODO:
#   - both integer and reals are treated as numeric -> the integer info
#    is lost!
#   - Replace ValueError by ParseError or something

# We know can handle the following:
#   - numeric and nominal attributes
#   - missing values for numeric attributes

r_meta = re.compile(r'^\s*@')
# Match a comment
r_comment = re.compile(r'^%')
# Match an empty line
r_empty = re.compile(r'^\s+$')
# Match a header line, that is a line which starts by @ + a word
r_headerline = re.compile(r'^@\S*')
r_datameta = re.compile(r'^@[Dd][Aa][Tt][Aa]')
r_relation = re.compile(r'^@[Rr][Ee][Ll][Aa][Tt][Ii][Oo][Nn]\s*(\S*)')
r_attribute = re.compile(r'^@[Aa][Tt][Tt][Rr][Ii][Bb][Uu][Tt][Ee]\s*(..*$)')

# To get attributes name enclosed with ''
r_comattrval = re.compile(r"'(..+)'\s+(..+$)")
# To get normal attributes
r_wcomattrval = re.compile(r"(\S+)\s+(..+$)")

#-------------------------
# Module defined exception
#-------------------------


class ArffError(IOError):
    pass


class ParseArffError(ArffError):
    pass

#------------------
# Various utilities
#------------------

# An attribute  is defined as @attribute name value


def parse_type(attrtype):
    """Given an arff attribute value (meta data), returns its type.

    Expect the value to be a name."""
    uattribute = attrtype.lower().strip()
    if uattribute[0] == '{':
        return 'nominal'
    elif uattribute[:len('real')] == 'real':
        return 'numeric'
    elif uattribute[:len('integer')] == 'integer':
        return 'numeric'
    elif uattribute[:len('numeric')] == 'numeric':
        return 'numeric'
    elif uattribute[:len('string')] == 'string':
        return 'string'
    elif uattribute[:len('relational')] == 'relational':
        return 'relational'
    elif uattribute[:len('date')] == 'date':
        return 'date'
    else:
        raise ParseArffError("unknown attribute %s" % uattribute)


def get_nominal(attribute):
    """If attribute is nominal, returns a list of the values"""
    return attribute.split(',')


def read_data_list(ofile):
    """Read each line of the iterable and put it in a list."""
    data = [next(ofile)]
    if data[0].strip()[0] == '{':
        raise ValueError("This looks like a sparse ARFF: not supported yet")
    data.extend([i for i in ofile])
    return data


def get_ndata(ofile):
    """Read the whole file to get number of data attributes."""
    data = [next(ofile)]
    loc = 1
    if data[0].strip()[0] == '{':
        raise ValueError("This looks like a sparse ARFF: not supported yet")
    for i in ofile:
        loc += 1
    return loc


def maxnomlen(atrv):
    """Given a string containing a nominal type definition, returns the
    string len of the biggest component.

    A nominal type is defined as seomthing framed between brace ({}).

    Parameters
    ----------
    atrv : str
       Nominal type definition

    Returns
    -------
    slen : int
       length of longest component

    Examples
    --------
    maxnomlen("{floup, bouga, fl, ratata}") returns 6 (the size of
    ratata, the longest nominal value).

    >>> maxnomlen("{floup, bouga, fl, ratata}")
    6
    """
    nomtp = get_nom_val(atrv)
    return max(len(i) for i in nomtp)


def get_nom_val(atrv):
    """Given a string containing a nominal type, returns a tuple of the
    possible values.

    A nominal type is defined as something framed between braces ({}).

    Parameters
    ----------
    atrv : str
       Nominal type definition

    Returns
    -------
    poss_vals : tuple
       possible values

    Examples
    --------
    >>> get_nom_val("{floup, bouga, fl, ratata}")
    ('floup', 'bouga', 'fl', 'ratata')
    """
    r_nominal = re.compile('{(.+)}')
    m = r_nominal.match(atrv)
    if m:
        return tuple(i.strip() for i in m.group(1).split(','))
    else:
        raise ValueError("This does not look like a nominal string")


def get_date_format(atrv):
    r_date = re.compile(r"[Dd][Aa][Tt][Ee]\s+[\"']?(.+?)[\"']?$")
    m = r_date.match(atrv)
    if m:
        pattern = m.group(1).strip()
        # convert time pattern from Java's SimpleDateFormat to C's format
        datetime_unit = None
        if "yyyy" in pattern:
            pattern = pattern.replace("yyyy", "%Y")
            datetime_unit = "Y"
        elif "yy":
            pattern = pattern.replace("yy", "%y")
            datetime_unit = "Y"
        if "MM" in pattern:
            pattern = pattern.replace("MM", "%m")
            datetime_unit = "M"
        if "dd" in pattern:
            pattern = pattern.replace("dd", "%d")
            datetime_unit = "D"
        if "HH" in pattern:
            pattern = pattern.replace("HH", "%H")
            datetime_unit = "h"
        if "mm" in pattern:
            pattern = pattern.replace("mm", "%M")
            datetime_unit = "m"
        if "ss" in pattern:
            pattern = pattern.replace("ss", "%S")
            datetime_unit = "s"
        if "z" in pattern or "Z" in pattern:
            raise ValueError("Date type attributes with time zone not "
                             "supported, yet")

        if datetime_unit is None:
            raise ValueError("Invalid or unsupported date format")

        return pattern, datetime_unit
    else:
        raise ValueError("Invalid or no date format")


def go_data(ofile):
    """Skip header.

    the first next() call of the returned iterator will be the @data line"""
    return itertools.dropwhile(lambda x: not r_datameta.match(x), ofile)


#----------------
# Parsing header
#----------------
def tokenize_attribute(iterable, attribute):
    """Parse a raw string in header (eg starts by @attribute).

    Given a raw string attribute, try to get the name and type of the
    attribute. Constraints:

    * The first line must start with @attribute (case insensitive, and
      space like characters before @attribute are allowed)
    * Works also if the attribute is spread on multilines.
    * Works if empty lines or comments are in between

    Parameters
    ----------
    attribute : str
       the attribute string.

    Returns
    -------
    name : str
       name of the attribute
    value : str
       value of the attribute
    next : str
       next line to be parsed

    Examples
    --------
    If attribute is a string defined in python as r"floupi real", will
    return floupi as name, and real as value.

    >>> iterable = iter([0] * 10) # dummy iterator
    >>> tokenize_attribute(iterable, r"@attribute floupi real")
    ('floupi', 'real', 0)

    If attribute is r"'floupi 2' real", will return 'floupi 2' as name,
    and real as value.

    >>> tokenize_attribute(iterable, r"  @attribute 'floupi 2' real   ")
    ('floupi 2', 'real', 0)

    """
    sattr = attribute.strip()
    mattr = r_attribute.match(sattr)
    if mattr:
        # atrv is everything after @attribute
        atrv = mattr.group(1)
        if r_comattrval.match(atrv):
            name, type = tokenize_single_comma(atrv)
            next_item = next(iterable)
        elif r_wcomattrval.match(atrv):
            name, type = tokenize_single_wcomma(atrv)
            next_item = next(iterable)
        else:
            # Not sure we should support this, as it does not seem supported by
            # weka.
            raise ValueError("multi line not supported yet")
            #name, type, next_item = tokenize_multilines(iterable, atrv)
    else:
        raise ValueError("First line unparsable: %s" % sattr)

    if type == 'relational':
        raise ValueError("relational attributes not supported yet")
    return name, type, next_item


def tokenize_single_comma(val):
    # XXX we match twice the same string (here and at the caller level). It is
    # stupid, but it is easier for now...
    m = r_comattrval.match(val)
    if m:
        try:
            name = m.group(1).strip()
            type = m.group(2).strip()
        except IndexError:
            raise ValueError("Error while tokenizing attribute")
    else:
        raise ValueError("Error while tokenizing single %s" % val)
    return name, type


def tokenize_single_wcomma(val):
    # XXX we match twice the same string (here and at the caller level). It is
    # stupid, but it is easier for now...
    m = r_wcomattrval.match(val)
    if m:
        try:
            name = m.group(1).strip()
            type = m.group(2).strip()
        except IndexError:
            raise ValueError("Error while tokenizing attribute")
    else:
        raise ValueError("Error while tokenizing single %s" % val)
    return name, type


def read_header(ofile):
    """Read the header of the iterable ofile."""
    i = next(ofile)

    # Pass first comments
    while r_comment.match(i):
        i = next(ofile)

    # Header is everything up to DATA attribute ?
    relation = None
    attributes = []
    while not r_datameta.match(i):
        m = r_headerline.match(i)
        if m:
            isattr = r_attribute.match(i)
            if isattr:
                name, type, i = tokenize_attribute(ofile, i)
                attributes.append((name, type))
            else:
                isrel = r_relation.match(i)
                if isrel:
                    relation = isrel.group(1)
                else:
                    raise ValueError("Error parsing line %s" % i)
                i = next(ofile)
        else:
            i = next(ofile)

    return relation, attributes


#--------------------
# Parsing actual data
#--------------------
def safe_float(x):
    """given a string x, convert it to a float. If the stripped string is a ?,
    return a Nan (missing value).

    Parameters
    ----------
    x : str
       string to convert

    Returns
    -------
    f : float
       where float can be nan

    Examples
    --------
    >>> safe_float('1')
    1.0
    >>> safe_float('1\\n')
    1.0
    >>> safe_float('?\\n')
    nan
    """
    if '?' in x:
        return np.nan
    else:
        return float(x)


def safe_nominal(value, pvalue):
    svalue = value.strip()
    if svalue in pvalue:
        return svalue
    elif svalue == '?':
        return svalue
    else:
        raise ValueError("%s value not in %s" % (str(svalue), str(pvalue)))


def safe_date(value, date_format, datetime_unit):
    date_str = value.strip().strip("'").strip('"')
    if date_str == '?':
        return np.datetime64('NaT', datetime_unit)
    else:
        dt = datetime.datetime.strptime(date_str, date_format)
        return np.datetime64(dt).astype("datetime64[%s]" % datetime_unit)


class MetaData(object):
    """Small container to keep useful information on a ARFF dataset.

    Knows about attributes names and types.

    Examples
    --------
    ::

        data, meta = loadarff('iris.arff')
        # This will print the attributes names of the iris.arff dataset
        for i in meta:
            print(i)
        # This works too
        meta.names()
        # Getting attribute type
        types = meta.types()

    Notes
    -----
    Also maintains the list of attributes in order, i.e. doing for i in
    meta, where meta is an instance of MetaData, will return the
    different attribute names in the order they were defined.
    """
    def __init__(self, rel, attr):
        self.name = rel
        # We need the dictionary to be ordered
        # XXX: may be better to implement an ordered dictionary
        self._attributes = {}
        self._attrnames = []
        for name, value in attr:
            tp = parse_type(value)
            self._attrnames.append(name)
            if tp == 'nominal':
                self._attributes[name] = (tp, get_nom_val(value))
            elif tp == 'date':
                self._attributes[name] = (tp, get_date_format(value)[0])
            else:
                self._attributes[name] = (tp, None)

    def __repr__(self):
        msg = ""
        msg += "Dataset: %s\n" % self.name
        for i in self._attrnames:
            msg += "\t%s's type is %s" % (i, self._attributes[i][0])
            if self._attributes[i][1]:
                msg += ", range is %s" % str(self._attributes[i][1])
            msg += '\n'
        return msg

    def __iter__(self):
        return iter(self._attrnames)

    def __getitem__(self, key):
        return self._attributes[key]

    def names(self):
        """Return the list of attribute names."""
        return self._attrnames

    def types(self):
        """Return the list of attribute types."""
        attr_types = [self._attributes[name][0] for name in self._attrnames]
        return attr_types


def loadarff(f):
    """
    Read an arff file.

    The data is returned as a record array, which can be accessed much like
    a dictionary of numpy arrays.  For example, if one of the attributes is
    called 'pressure', then its first 10 data points can be accessed from the
    ``data`` record array like so: ``data['pressure'][0:10]``


    Parameters
    ----------
    f : file-like or str
       File-like object to read from, or filename to open.

    Returns
    -------
    data : record array
       The data of the arff file, accessible by attribute names.
    meta : `MetaData`
       Contains information about the arff file such as name and
       type of attributes, the relation (name of the dataset), etc...

    Raises
    ------
    ParseArffError
        This is raised if the given file is not ARFF-formatted.
    NotImplementedError
        The ARFF file has an attribute which is not supported yet.

    Notes
    -----

    This function should be able to read most arff files. Not
    implemented functionality include:

    * date type attributes
    * string type attributes

    It can read files with numeric and nominal attributes.  It cannot read
    files with sparse data ({} in the file).  However, this function can
    read files with missing data (? in the file), representing the data
    points as NaNs.

    Examples
    --------
    >>> from scipy.io import arff
    >>> from io import StringIO
    >>> content = \"\"\"
    ... @relation foo
    ... @attribute width  numeric
    ... @attribute height numeric
    ... @attribute color  {red,green,blue,yellow,black}
    ... @data
    ... 5.0,3.25,blue
    ... 4.5,3.75,green
    ... 3.0,4.00,red
    ... \"\"\"
    >>> f = StringIO(content)
    >>> data, meta = arff.loadarff(f)
    >>> data
    array([(5.0, 3.25, 'blue'), (4.5, 3.75, 'green'), (3.0, 4.0, 'red')],
          dtype=[('width', '<f8'), ('height', '<f8'), ('color', '|S6')])
    >>> meta
    Dataset: foo
    \twidth's type is numeric
    \theight's type is numeric
    \tcolor's type is nominal, range is ('red', 'green', 'blue', 'yellow', 'black')

    """
    if hasattr(f, 'read'):
        ofile = f
    else:
        ofile = open(f, 'rt')
    try:
        return _loadarff(ofile)
    finally:
        if ofile is not f:  # only close what we opened
            ofile.close()


def _loadarff(ofile):
    # Parse the header file
    try:
        rel, attr = read_header(ofile)
    except ValueError as e:
        msg = "Error while parsing header, error was: " + str(e)
        raise ParseArffError(msg)

    # Check whether we have a string attribute (not supported yet)
    hasstr = False
    for name, value in attr:
        type = parse_type(value)
        if type == 'string':
            hasstr = True

    meta = MetaData(rel, attr)

    # XXX The following code is not great
    # Build the type descriptor descr and the list of convertors to convert
    # each attribute to the suitable type (which should match the one in
    # descr).

    # This can be used once we want to support integer as integer values and
    # not as numeric anymore (using masked arrays ?).
    acls2dtype = {'real': float, 'integer': float, 'numeric': float}
    acls2conv = {'real': safe_float,
                 'integer': safe_float,
                 'numeric': safe_float}
    descr = []
    convertors = []
    if not hasstr:
        for name, value in attr:
            type = parse_type(value)
            if type == 'date':
                date_format, datetime_unit = get_date_format(value)
                descr.append((name, "datetime64[%s]" % datetime_unit))
                convertors.append(partial(safe_date, date_format=date_format,
                                          datetime_unit=datetime_unit))
            elif type == 'nominal':
                n = maxnomlen(value)
                descr.append((name, 'S%d' % n))
                pvalue = get_nom_val(value)
                convertors.append(partial(safe_nominal, pvalue=pvalue))
            else:
                descr.append((name, acls2dtype[type]))
                convertors.append(safe_float)
                #dc.append(acls2conv[type])
                #sdescr.append((name, acls2sdtype[type]))
    else:
        # How to support string efficiently ? Ideally, we should know the max
        # size of the string before allocating the numpy array.
        raise NotImplementedError("String attributes not supported yet, sorry")

    ni = len(convertors)

    def generator(row_iter, delim=','):
        # TODO: this is where we are spending times (~80%). I think things
        # could be made more efficiently:
        #   - We could for example "compile" the function, because some values
        #   do not change here.
        #   - The function to convert a line to dtyped values could also be
        #   generated on the fly from a string and be executed instead of
        #   looping.
        #   - The regex are overkill: for comments, checking that a line starts
        #   by % should be enough and faster, and for empty lines, same thing
        #   --> this does not seem to change anything.

        # 'compiling' the range since it does not change
        # Note, I have already tried zipping the converters and
        # row elements and got slightly worse performance.
        elems = list(range(ni))

        for raw in row_iter:
            # We do not abstract skipping comments and empty lines for
            # performance reasons.
            if r_comment.match(raw) or r_empty.match(raw):
                continue
            row = raw.split(delim)
            yield tuple([convertors[i](row[i]) for i in elems])

    a = generator(ofile)
    # No error should happen here: it is a bug otherwise
    data = np.fromiter(a, descr)
    return data, meta


#-----
# Misc
#-----
def basic_stats(data):
    nbfac = data.size * 1. / (data.size - 1)
    return np.nanmin(data), np.nanmax(data), np.mean(data), np.std(data) * nbfac


def print_attribute(name, tp, data):
    type = tp[0]
    if type == 'numeric' or type == 'real' or type == 'integer':
        min, max, mean, std = basic_stats(data)
        print("%s,%s,%f,%f,%f,%f" % (name, type, min, max, mean, std))
    else:
        msg = name + ",{"
        for i in range(len(tp[1])-1):
            msg += tp[1][i] + ","
        msg += tp[1][-1]
        msg += "}"
        print(msg)


def test_weka(filename):
    data, meta = loadarff(filename)
    print(len(data.dtype))
    print(data.size)
    for i in meta:
        print_attribute(i, meta[i], data[i])


# make sure nose does not find this as a test
test_weka.__test__ = False


if __name__ == '__main__':
    import sys
    filename = sys.argv[1]
    test_weka(filename)
