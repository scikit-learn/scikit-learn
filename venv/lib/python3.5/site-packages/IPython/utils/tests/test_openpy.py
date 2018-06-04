import io
import os.path
import nose.tools as nt

from IPython.utils import openpy

mydir = os.path.dirname(__file__)
nonascii_path = os.path.join(mydir, '../../core/tests/nonascii.py')

def test_detect_encoding():
    f = open(nonascii_path, 'rb')
    enc, lines = openpy.detect_encoding(f.readline)
    nt.assert_equal(enc, 'iso-8859-5')

def test_read_file():
    read_specified_enc = io.open(nonascii_path, encoding='iso-8859-5').read()
    read_detected_enc = openpy.read_py_file(nonascii_path, skip_encoding_cookie=False)
    nt.assert_equal(read_detected_enc, read_specified_enc)
    assert u'coding: iso-8859-5' in read_detected_enc
    
    read_strip_enc_cookie = openpy.read_py_file(nonascii_path, skip_encoding_cookie=True)
    assert u'coding: iso-8859-5' not in read_strip_enc_cookie

def test_source_to_unicode():
    with io.open(nonascii_path, 'rb') as f:
        source_bytes = f.read()
    nt.assert_equal(openpy.source_to_unicode(source_bytes, skip_encoding_cookie=False).splitlines(),
                    source_bytes.decode('iso-8859-5').splitlines())

    source_no_cookie = openpy.source_to_unicode(source_bytes, skip_encoding_cookie=True)
    nt.assert_not_in(u'coding: iso-8859-5', source_no_cookie)

def test_list_readline():
    l = ['a', 'b']
    readline = openpy._list_readline(l)
    nt.assert_equal(readline(), 'a')
    nt.assert_equal(readline(), 'b')
    with nt.assert_raises(StopIteration):
        readline()