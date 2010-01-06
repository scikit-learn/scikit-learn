#! /usr/bin/env python
# Last Change: Mon Aug 20 02:00 PM 2007 J
import re
import itertools
import sys

from scikits.learn.common import partial

import numpy as N

"""A module to read arff files."""

# An Arff file is basically two parts: 
#   - header
#   - data
#
# A header has each of its components starting by @META where META is one of
# the keyword (attribute of relation, for now).

# TODO:
#   - both integer and reals are treated as numeric -> the integer info is lost !
#   - Replace ValueError by ParseError or something

# We know can handle the following:
#   - numeric and nominal attributes
#   - missing values for numeric attributes

r_meta = re.compile('^\s*@')
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
# To get attributes name enclosed with '', possibly spread accross multilines
r_mcomattrval = re.compile(r"'([..\n]+)'\s+(..+$)")
# To get normal attributes 
r_wcomattrval = re.compile(r"(\S+)\s+(..+$)")

class ArffError(IOError):
    pass

class ParseArffError(ArffError):
    pass

# An attribute  is defined as @attribute name value
def parse_type(attrtype):
    """Given an arff attribute value, returns its type.
    
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
    else:
        raise ValueError("unknown attribute %s" % uattribute)


def get_nominal(attribute):
    """If attribute is nominal, returns a list of the values"""
    return attribute.split(',')
        

#-------------------------------------
# Functions to parse lines into tokens
#-------------------------------------
def tokenize_attribute(iterable, attribute):
    """Parse a raw string in header (eg starts by @attribute).
    
    Given a raw string attribute, try to get the name and type of the
    attribute. Constraints:
        - The first line must start with @attribute (case insensitive, and
          space like characters begore @attribute are allowed)
        - Works also if the attribute is spread on multilines. 
        - Works if empty lines or comments are in between
    
    :Parameters:
        attribute : str
            the attribute string. 
    
    :Returns:
        name : str
            name of the attribute
        value : str
            value of the attribute
        next : str
            next line to be parsed

    Example:
        - if attribute is a string defined in python as r"floupi real", will
          return floupi as name, and real as value.
        - if attribute is r"'floupi 2' real", will return 'floupi 2' as name,
          and real as value. """
    sattr = attribute.strip()
    mattr = r_attribute.match(sattr)
    if mattr:
        # atrv is everything after @attribute
        atrv = mattr.group(1)
        if r_comattrval.match(atrv):
            name, type = tokenize_single_comma(atrv)
            next = iterable.next()
        elif r_wcomattrval.match(atrv):
            name, type = tokenize_single_wcomma(atrv)
            next = iterable.next()
        else:
            # Not sure we should support this, as it does not seem supported by
            # weka.
            raise ValueError("multi line not supported yet")
            #name, type, next = tokenize_multilines(iterable, atrv)
    else:
        raise ValueError("First line unparsable: %s" % sattr)

    if type == 'relational':
        raise ValueError("relational attributes not supported yet")
    return name, type, next

def tokenize_multilines(iterable, val):
    """Can tokenize an attribute spread over several lines."""
    ## skip empty lines
    #while r_empty.match(val):
    #    val = iterable.next()

    # If one line does not match, read all the following lines up to next
    # line with meta character, and try to parse everything up to there.
    if not r_mcomattrval.match(val):
        all = [val]
        i = iterable.next()
        while not r_meta.match(i):
            all.append(i)
            i = iterable.next()
        if r_mend.search(i):
            raise ValueError("relational attribute not supported yet")
        print "".join(all[:-1])
        m = r_comattrval.match("".join(all[:-1]))
        return m.group(1), m.group(2), i
    else:
        raise ValueError("Cannot parse attribute names spread over multi "\
                        "lines yet")
    
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
    """Reader header of the iterable ofile."""
    i = ofile.next()

    # Pass first comments
    while r_comment.match(i):
        i = ofile.next()

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
                i = ofile.next()
        else:
            i = ofile.next()

    return relation, attributes

def attribute2dtype(typestr):
    """From a type string, returns the corresponding dtype.

    For example, if typestr is 'numeric', returns float dtype, etc..."""

def read_data_list(ofile):
    """Read each line of the iterable and put it in a list."""
    data = [ofile.next()]
    if data[0].strip()[0] == '{':
        raise ValueError("This looks like a sparse ARFF: not supported yet")
    data.extend([i for i in ofile])
    return data

def get_ndata(ofile):
    """Read the whole file to get number of data attributes."""
    #data = [ofile.next()]
    #if data[0].strip()[0] == '{':
    #    raise ValueError("This looks like a sparse ARFF: not supported yet")
    #data.extend([i for i in ofile])
    #return data
    data = [ofile.next()]
    loc = 1
    if data[0].strip()[0] == '{':
        raise ValueError("This looks like a sparse ARFF: not supported yet")
    for i in ofile:
        loc += 1
    return loc

def maxnomlen(atrv):
    """Given a string contening a nominal type, returns the string len of the
    biggest component.
    
    A nominal type is defined as seomthing framed between brace ({})."""
    r_nominal = re.compile('{(..+)}')
    m = r_nominal.match(atrv)
    if m:
        values = m.group(1).split(',')
        return max([len(i.strip()) for i in values])
    else:
        raise ValueError("This does not look like a nominal string")

def get_nom_val(atrv):
    """Given a string contening a nominal type, returns a tuple of the possible
    values.    
    
    A nominal type is defined as something framed between brace ({})."""
    r_nominal = re.compile('{(..+)}')
    m = r_nominal.match(atrv)
    if m:
        return tuple(i.strip() for i in m.group(1).split(','))
    else:
        raise ValueError("This does not look like a nominal string")

def go_data(ofile):
    """Skip header.
    
    the first next() call of the returned iterator will be the @data line"""
    return itertools.dropwhile(lambda x : not r_datameta.match(x), ofile)

#def get_header(ofile):
#    """Get the header as a list of lines."""
#    return itertools.takewhile(lambda x : not r_datameta.match(x), ofile)

def safe_float(x):
    """given a string x, convert it to a float. If the stripped string is a ?,
    return a Nan (missing value)."""
    if x.strip() == '?':
        return N.nan
    else:
        return N.float(x)

def safe_nominal(value, pvalue):
    svalue = value.strip()
    if svalue in pvalue:
        return svalue
    elif svalue == '?':
        return svalue
    else:
        raise ValueError("%s value not in %s" % (str(svalue), str(pvalue)))

def get_delim(line):
    """Given a string representing a line of data, check whether the
    delimiter is ',' or space."""
    l = line.split(',')
    if len(l) > 1:
        return ','
    else:
        l = line.split(' ')
        if len(l) > 1:
            return ' '
        else:
            raise ValueError("delimiter not understood: " + line)

def read_arff(filename):
    ofile = open(filename)

    # Parse the header file 
    try:
        rel, attr = read_header(ofile)
    except ValueError, e:
        msg = "Error while parsing header, error was: " + str(e)
        raise ParseArffError(msg)

    # Check whether we have a string attribute (not supported yet)
    hasstr = False
    for name, value in attr:
        type = parse_type(value)
        if type == 'string':
            hasstr = True

    # XXX The following code is ugly
    # Build the type descr from the attributes
    #acls2sdtype = {'real' : 'N.float', 'integer' : 'N.float', 'numeric' : 'N.float'}
    acls2dtype = {'real' : N.float, 'integer' : N.float, 'numeric' : N.float}
    acls2conv = {'real' : safe_float, 'integer' : safe_float, 'numeric' : safe_float}
    descr = []
    dc = []
    if not hasstr:
        for name, value in attr:
            type = parse_type(value)
            if type == 'date':
                raise ValueError("date type not supported yet, sorry")
            elif type == 'nominal':
                n = maxnomlen(value)
                descr.append((name, 'S%d' % n))
                pvalue = get_nom_val(value)
                dc.append(partial(safe_nominal, pvalue = pvalue))
            else:
                descr.append((name, acls2dtype[type]))
                dc.append(safe_float)
                #dc.append(acls2conv[type])
                #sdescr.append((name, acls2sdtype[type]))
    else:
        raise NotImplementedError("String attributes not supported yet, sorry")

    # dc[i] returns a callable which can convert the ith element of a row of
    # data
    ni = len(dc)

    # Get the delimiter from the first line of data:
    def next_data_line(row_iter):
        """Assumes we are already in the data part (eg after @data)."""
        raw = row_iter.next()
        while r_empty.match(raw):
            raw = row_iter.next()
        while r_comment.match(raw):
            raw = row_iter.next()
        return raw

    try:
        dtline = next_data_line(ofile)
        delim = get_delim(dtline)
    except ValueError, e:
        raise ParseArffError("Error while parsing delimiter: " + str(e))
    finally:
        ofile.seek(0, 0)
        ofile = go_data(ofile)
        # skip the @data line
        ofile.next()

    def generator(row_iter, delim = ','):
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

        # We do not abstract skipping comments and empty lines for performances
        # reason.
        raw = row_iter.next()
        while r_empty.match(raw):
            raw = row_iter.next()
        while r_comment.match(raw):
            raw = row_iter.next()

        row = raw.split(delim)
        yield tuple([dc[i](row[i]) for i in range(ni)])
        for raw in row_iter:
            while r_comment.match(raw):
                raw = row_iter.next()
            while r_empty.match(raw):
                raw = row_iter.next()
            row = raw.split(delim)
            yield tuple([dc[i](row[i]) for i in range(ni)])

    a = generator(ofile, delim = delim)
    # No error should happen here: it is a bug otherwise
    data = N.fromiter(a, descr)
    return data, rel, [parse_type(j) for i, j in attr]

def basic_stats(data):
    return N.min(data), N.max(data), N.mean(data), N.std(data)

def floupi(filename):
    data, rel, types = read_arff(filename)
    from attrselect import print_dataset_info
    print_dataset_info(data)
    #print "relation %s, has %d instances" % (rel, data.size)
    #itp = iter(types)
    #for i in data.dtype.names:
    #    tp = itp.next()
    #    if tp == 'numeric' or tp == 'real' or tp == 'integer':
    #        min, max, mean, std = basic_stats(data[i])
    #        print "\tinstance %s: min %f, max %f, mean %f, std %f" % \
    #                (i, min, max, mean, std)
    #    else:
    #        print "\tinstance %s is nominal" % i

if __name__ == '__main__':
    import glob
    for i in glob.glob('arff.bak/data/*'):
        relation, attributes = read_header(open(i))
        print "Parsing header of %s: relation %s, %d attributes" % (i,
                relation, len(attributes))

    #import sys
    #filename = sys.argv[1]
    #filename = 'arff.bak/data/pharynx.arff'
    #floupi(filename)

    gf = []
    wf = []
    for i in glob.glob('arff.bak/data/*'):
        try:
            print "=============== reading %s ======================" % i
            floupi(i)
            gf.append(i)
        except ValueError, e:
            print "!!!! Error parsing the file !!!!!"
            print e
            wf.append(i)
        except IndexError, e:
            print "!!!! Error parsing the file !!!!!"
            print e
            wf.append(i)
        except ArffError, e:
            print "!!!! Error parsing the file !!!!!"
            print e
            wf.append(i)

    print "%d good files" % len(gf)
    print "%d bad files" % len(wf)
