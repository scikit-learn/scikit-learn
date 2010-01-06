#! /usr/bin/env python
# Last Change: Sun Jul 22 03:00 PM 2007 J

"""This module implements function to extract attributes and/or classes from
datasets."""
import numpy as N

def print_dataset_info(data, label = None, cl = None):
    # Attributes info
    attr = get_attributes(data)
    msg = "data has: \n" + "\t%d attributes: " % len(attr)
    if len(attr) > 0:
        msg += ", ".join([i for i in attr[:-1]])
        msg += " and " + attr[-1]
    msg += '\n'

    if label is not None:
        # Classes info
        if cl is not None:
            if len(cl) > 0:
                msg += "\t%d classes: " % len(cl)
                msg += ", ".join([i for i in cl])
    else:
        msg += "\tNo classes"

    msg += '\n'

    # Number of samples
    ns = len(data)
    if label is not None:
        if cl is not None and len(cl) > 0:
            msg += "\t%d samples in the dataset:\n" % ns
            c2ind = get_c2ind(cl, label)
            msg += "".join(["\t\t%d samples for class %s\n" \
                              % (len(c2ind[cname]), cname) \
                              for cname in cl])
    else:
        msg += "\t%d samples in the dataset\n" % ns

    print msg

def get_attributes(data):
    """Given a numpy array, tries to find its fields as a list of names. If  no
    fields, returns a empty list."""
    if not isinstance(data, N.ndarray):
        raise ValueError("Expect a numpy array")

    f = data.dtype.fields
    try:
        attr = f.keys()
    except AttributeError, e:
        attr = []

    return attr

def get_data_from_class(data, label, cl, cname):
    """return all data samples of the class cname."""
    if cname not in cl:
        raise ValueError("class %s not found in this dataset" % cname)
    c2lab = get_c2lab(cl)
    return data[label==c2lab[cname]]

def get_subset_from_class(data, label, cl, cname, ind):
    """Function to index inside one class.
    
    In other workds, return all d[ind] where d are the data from the class
    cname."""
    c2ind = get_c2ind(cl, label)
    try:
       cind = c2ind[cname]
       find = cind[ind]
       return data[find]
    except KeyError:
        raise ValueError("class %s is not found in the dataset" % cname)
    except IndexError:
        raise IndexError("Given index is outside the class %s range" % cname)

def get_c2lab(cl):
    """Given an array of string, return a dictionary which keys are the
    strings, and the values the index."""
    return dict([(cl[i], i) for i in range(len(cl))])

def select_data_attr(data, attributes):
    """Return only the given attributes of the data."""
    return N.hstack([data[i][:, N.newaxis] for i in attributes])

def get_c2ind(cl, label):
    """Given a sequence of strings and an array of label, returns a dictionary
    which keys are the strings, and the values the index arrays corresponding
    to the classes. """
    c2lab = get_c2lab(cl)
    c2ind = dict([(cname, N.where(label == c2lab[cname])[0]) for cname in cl])
    return c2ind

if __name__ == '__main__':
    from scikits.learn.datasets import iris, german, pendigits, oldfaithful
    d = iris.load()
    data, lab, cl = d['data'], d['label'], d['class']
    print_dataset_info(data, lab, cl)

    d = german.load()
    data, lab, cl = d['data'], d['label'], d['class']
    print_dataset_info(data, lab, cl)

    d = oldfaithful.load()
    data = d['data']
    print_dataset_info(data)

    d = pendigits.testing.load()
    data, lab, cl = d['data'], d['label'], d['class']
    print_dataset_info(data, lab, cl)

    d = pendigits.training.load()
    data, lab, cl = d['data'], d['label'], d['class']
    print_dataset_info(data, lab, cl)
