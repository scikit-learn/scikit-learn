#! /usr/bin/env python
# Last Change: Thu Aug 02 07:00 PM 2007 J
import re

import numpy as N

# Match a comment
r_comment = re.compile(r'^%')
# Match a header line, that is a line which starts by @ + a word
r_headerline = re.compile(r'^@\S*')
r_datameta = re.compile(r'^@[Dd][Aa][Tt][Aa]')
r_relation = re.compile(r'^@[Rr][Ee][Ll][Aa][Tt][Ii][Oo][Nn]\s*(\S*)')
r_attribute = re.compile(r'^@[Aa][Tt][Tt][Rr][Ii][Bb][Uu][Tt][Ee]\s*(..*$)')

# To get attributes name enclosed with ''
r_comattrval = re.compile(r"'(..*)'\s*(..*$)")
# To get normal attributes 
r_wcomattrval = re.compile(r"(\S*)\s*(..*$)")
_arff_aclass = {
    'numeric' : 0,
    'nominal' : 1,
    'string' : 2,
    'date' : 4,
    'relational' : 8,
}

acls2id = dict(_arff_aclass)
acls2id['real'] = _arff_aclass['numeric']
acls2id['integer'] = _arff_aclass['numeric']
id2acls = N.empty(len(acls2id), 'S%d' % N.max([len(i) for i in acls2id.keys()]))
id2dtype = {acls2id['numeric'] : N.float,
    acls2id['nominal'] : N.int}

def parse_type(attrtype):
    """Given an arff attribute type string, returns its type"""
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
        

def parse_attribute(attribute):
    """Parse a raw string attribute.
    
    Given a raw string attribute (eg everything after @attribute), try to get
    the name and type of the attribute.
    
    :Parameters:
        attribute : str
            the attribute string (defined as everything after @attribute)
    
    :Returns:
        name : str
            name of the attribute
        type : str
            type of the attribute

    Example:
        - if attribute is a string defined in python as r"floupi real", will
          return floupi as name, and real as type.
        - if attribute is r"'floupi 2' real", will return 'floupi 2' as name,
          and real as type.
        """
    sattr = attribute.strip()
    m = r_comattrval.search(sattr)
    if m:
        try:
            name = m.group(1).strip()
            type = m.group(2).strip()
        except IndexError:
            raise ValueError("Error while parsing attribute")
    else:
        m = r_wcomattrval.search(sattr)
        if not m:
            raise ValueError("Error while parsing attribute")
        try:
            name = m.group(1).strip()
            type = m.group(2).strip()
        except IndexError:
            raise ValueError("Error while parsing attribute")

    return name, type

class FLOUPI:
    attr1 = "floupi numeric"
    attr2 = " floupi numeric "
    attr3 = " 'floupi  2' numeric "

def read_header(ofile):
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
                name, type = parse_attribute(isattr.group(1))
                if type == 'relational':
                    raise ValueError("relational attribute not supported yet")
                attributes.append((name, type))
            else:
                isrel = r_relation.match(i)
                if isrel:
                    relation = isrel.group(1)
        i = ofile.next()

    return relation, attributes

def read_data(ofile, attributes, sb = None):
    """If sb is none, the whole file is read to get the size of returned
    data."""
    if sb:
        raise ValueError("Offline reading not supported yet")
    data = [ofile.next()]
    if data[0].strip()[0] == '{':
        raise ValueError("This looks like a sparse ARFF: not supported yet")
    data.extend([i for i in ofile])
    descr = _descr_from_attr(attributes)
    ndata = N.empty(len(data), descr)
    for i in range(10):
        ndata[i] = [descr[j][1](data[i][j]) for j in range(len(data[i]))]
    return ndata

def _descr_from_attr(attr):
    return [(name, id2dtype[acls2id[parse_type(type)]]) for name, type in attr]

def _read_data_raw(ofile, ar, sb = None):
    if sb is None:
        data = [i for i in ofile.readlines()]
    else:
        data = [i for i in ofile.readlines(sb)]
    ar[:len(data)] = data
    
if __name__ == '__main__':
    import sys
    filename = sys.argv[1]
    a = open(filename)
    rel, attr = read_header(a)
    for name, type in attr:
        if name.lower() == 'class':
            if parse_type(type) == 'nominal':
                print "Class has %d values" % len(get_nominal(type))
            else:
                print "Class is type %s" % parse_type(type)
        else:
            print "attribute %s is class %s" % (name, parse_type(type))
    print "%d attributes " % len(attr)
    data = read_data(a, attr)
    print data
