# encoding: utf-8
"""Tests for IPython.utils.capture"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2013 The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


import sys

import nose.tools as nt

from IPython.utils import capture

#-----------------------------------------------------------------------------
# Globals
#-----------------------------------------------------------------------------

_mime_map = dict(
    _repr_png_="image/png",
    _repr_jpeg_="image/jpeg",
    _repr_svg_="image/svg+xml",
    _repr_html_="text/html",
    _repr_json_="application/json",
    _repr_javascript_="application/javascript",
)

basic_data = {
    'image/png' : b'binarydata',
    'text/html' : "<b>bold</b>",
}
basic_metadata = {
    'image/png' : {
        'width' : 10,
        'height' : 20,
    },
}

full_data = {
    'image/png' : b'binarydata',
    'image/jpeg' : b'binarydata',
    'image/svg+xml' : "<svg>",
    'text/html' : "<b>bold</b>",
    'application/javascript' : "alert();",
    'application/json' : "{}",
}
full_metadata = {
    'image/png' : {"png" : "exists"},
    'image/jpeg' : {"jpeg" : "exists"},
    'image/svg+xml' : {"svg" : "exists"},
    'text/html' : {"html" : "exists"},
    'application/javascript' : {"js" : "exists"},
    'application/json' : {"json" : "exists"},
}

hello_stdout = "hello, stdout"
hello_stderr = "hello, stderr"

#-----------------------------------------------------------------------------
# Test Functions
#-----------------------------------------------------------------------------

def test_rich_output_empty():
    """RichOutput with no args"""
    rich = capture.RichOutput()
    for method, mime in _mime_map.items():
        yield nt.assert_equal, getattr(rich, method)(), None
    
def test_rich_output():
    """test RichOutput basics"""
    data = basic_data
    metadata = basic_metadata
    rich = capture.RichOutput(data=data, metadata=metadata)
    yield nt.assert_equal, rich._repr_html_(), data['text/html']
    yield nt.assert_equal, rich._repr_png_(), (data['image/png'], metadata['image/png'])
    yield nt.assert_equal, rich._repr_latex_(), None
    yield nt.assert_equal, rich._repr_javascript_(), None
    yield nt.assert_equal, rich._repr_svg_(), None

def test_rich_output_no_metadata():
    """test RichOutput with no metadata"""
    data = full_data
    rich = capture.RichOutput(data=data)
    for method, mime in _mime_map.items():
        yield nt.assert_equal, getattr(rich, method)(), data[mime]

def test_rich_output_metadata():
    """test RichOutput with metadata"""
    data = full_data
    metadata = full_metadata
    rich = capture.RichOutput(data=data, metadata=metadata)
    for method, mime in _mime_map.items():
        yield nt.assert_equal, getattr(rich, method)(), (data[mime], metadata[mime])

def test_rich_output_display():
    """test RichOutput.display
    
    This is a bit circular, because we are actually using the capture code we are testing
    to test itself.
    """
    data = full_data
    rich = capture.RichOutput(data=data)
    with capture.capture_output() as cap:
        rich.display()
    yield nt.assert_equal, len(cap.outputs), 1
    rich2 = cap.outputs[0]
    yield nt.assert_equal, rich2.data, rich.data
    yield nt.assert_equal, rich2.metadata, rich.metadata

def test_capture_output():
    """capture_output works"""
    rich = capture.RichOutput(data=full_data)
    with capture.capture_output() as cap:
        print(hello_stdout, end="")
        print(hello_stderr, end="", file=sys.stderr)
        rich.display()
    yield nt.assert_equal, hello_stdout, cap.stdout
    yield nt.assert_equal, hello_stderr, cap.stderr

def test_capture_output_no_stdout():
    """test capture_output(stdout=False)"""
    rich = capture.RichOutput(data=full_data)
    with capture.capture_output(stdout=False) as cap:
        print(hello_stdout, end="")
        print(hello_stderr, end="", file=sys.stderr)
        rich.display()
    yield nt.assert_equal, "", cap.stdout
    yield nt.assert_equal, hello_stderr, cap.stderr
    yield nt.assert_equal, len(cap.outputs), 1

def test_capture_output_no_stderr():
    """test capture_output(stderr=False)"""
    rich = capture.RichOutput(data=full_data)
    # add nested capture_output so stderr doesn't make it to nose output
    with capture.capture_output(), capture.capture_output(stderr=False) as cap:
        print(hello_stdout, end="")
        print(hello_stderr, end="", file=sys.stderr)
        rich.display()
    yield nt.assert_equal, hello_stdout, cap.stdout
    yield nt.assert_equal, "", cap.stderr
    yield nt.assert_equal, len(cap.outputs), 1

def test_capture_output_no_display():
    """test capture_output(display=False)"""
    rich = capture.RichOutput(data=full_data)
    with capture.capture_output(display=False) as cap:
        print(hello_stdout, end="")
        print(hello_stderr, end="", file=sys.stderr)
        rich.display()
    yield nt.assert_equal, hello_stdout, cap.stdout
    yield nt.assert_equal, hello_stderr, cap.stderr
    yield nt.assert_equal, cap.outputs, []