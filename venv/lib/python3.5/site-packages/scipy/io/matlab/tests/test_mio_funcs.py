''' Jottings to work out format for __function_workspace__ matrix at end
of mat file.

'''
from __future__ import division, print_function, absolute_import

import os.path
import sys
import io

from numpy.compat import asstr

from scipy.io.matlab.mio5 import (MatlabObject, MatFile5Writer,
                                  MatFile5Reader, MatlabFunction)

test_data_path = os.path.join(os.path.dirname(__file__), 'data')


def read_minimat_vars(rdr):
    rdr.initialize_read()
    mdict = {'__globals__': []}
    i = 0
    while not rdr.end_of_stream():
        hdr, next_position = rdr.read_var_header()
        name = asstr(hdr.name)
        if name == '':
            name = 'var_%d' % i
            i += 1
        res = rdr.read_var_array(hdr, process=False)
        rdr.mat_stream.seek(next_position)
        mdict[name] = res
        if hdr.is_global:
            mdict['__globals__'].append(name)
    return mdict


def read_workspace_vars(fname):
    fp = open(fname, 'rb')
    rdr = MatFile5Reader(fp, struct_as_record=True)
    vars = rdr.get_variables()
    fws = vars['__function_workspace__']
    ws_bs = io.BytesIO(fws.tostring())
    ws_bs.seek(2)
    rdr.mat_stream = ws_bs
    # Guess byte order.
    mi = rdr.mat_stream.read(2)
    rdr.byte_order = mi == b'IM' and '<' or '>'
    rdr.mat_stream.read(4)  # presumably byte padding
    mdict = read_minimat_vars(rdr)
    fp.close()
    return mdict


def test_jottings():
    # example
    fname = os.path.join(test_data_path, 'parabola.mat')
    ws_vars = read_workspace_vars(fname)
