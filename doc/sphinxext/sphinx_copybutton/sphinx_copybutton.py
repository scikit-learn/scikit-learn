# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2018 Chris Holdgraf

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
# This sphinx extension is to add a button to copy the code from examples.
# Details at https://github.com/choldgraf/sphinx-copybutton
import os


def scb_static_path(app):
    static_path = os.path.abspath(os.path.join(
                                  os.path.dirname(__file__), '_static'))
    app.config.html_static_path.append(static_path)


# Clipboard.js script for sphinx_copybutton:
clipboard_js_url = ('https://cdnjs.cloudflare.com/ajax/'
                    'libs/clipboard.js/2.0.0/clipboard.min.js')


def setup(app):
    # Add our static path
    app.connect('builder-inited', scb_static_path)

    # Add relevant code to headers
    app.add_stylesheet('sphinx_copybutton.css')
    app.add_javascript('sphinx_copybutton.js')
    app.add_javascript(clipboard_js_url)
