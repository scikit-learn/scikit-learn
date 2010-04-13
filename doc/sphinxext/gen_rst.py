"""
generate the rst files for the examples by iterating over the pylab examples.

Images must be created manually and dropped in directory auto_examples/images. Files that generate images should start with 'plot'

"""
import os, glob

import os
import re
import sys
fileList = []

                          
import token, tokenize      

                                                                                                                    

rst_template = """%(docstring)s

.. literalinclude:: %(short_fname)s
    :lines: %(end_row)s-
    """

plot_rst_template = """%(docstring)s

.. image:: images/%(image_name)s

.. literalinclude:: %(short_fname)s
    :lines: %(end_row)s-
    """


def extract_docstring(filename):
    # Extract a module-level docstring, if any
    lines = file(filename).readlines()
    start_row = 0
    if lines[0].startswith('#!'):
        lines.pop(0)
        start_row = 1

    docstring = ''
    first_par = ''
    tokens = tokenize.generate_tokens(lines.__iter__().next)
    for tok_type, tok_content, _, (erow, _), _ in tokens:
        tok_type = token.tok_name[tok_type]    
        if tok_type in ('NEWLINE', 'COMMENT', 'NL', 'INDENT', 'DEDENT'):
            continue
        elif tok_type == 'STRING':
            docstring = eval(tok_content)
            # If the docstring is formatted with several paragraphs, extract
            # the first one:
            paragraphs = '\n'.join(line.rstrip() 
                                for line in docstring.split('\n')).split('\n\n')
            if len(paragraphs) > 0:
                first_par = paragraphs[0]
        break
    return docstring, first_par, erow+1+start_row


def generate_example_rst(app):
    rootdir = os.path.join(app.builder.srcdir, 'auto_examples')
    exampledir = os.path.abspath(app.builder.srcdir +  '/../' + 'examples')
    if not os.path.exists(exampledir):
        os.makedirs(exampledir)

    datad = []

    for root, dirs, files in os.walk(exampledir):
        for fname in files:
            global rst_template, plot_rst_template
            if  not fname.endswith('py'): continue
            if fname.startswith('plot'): rst_template = plot_rst_template
            complete_fname = os.path.join(exampledir, fname)
            docstring, short_desc, end_row = extract_docstring(complete_fname)

            short_fname = '../../examples/' + fname

            image_name = fname[:-2] + 'png'

            f = open(os.path.join(rootdir, fname[:-2] + 'rst'),'w')
            print rst_template
            f.write( rst_template % locals())
            f.flush()
            datad.append(fname)

    # we create an index.rst with all examples
    fhindex = file(os.path.join(rootdir, 'index.rst'), 'w')
    fhindex.write("""\
.. _examples-index:

    :Release: |version|
    :Date: |today|

.. toctree::

""")
    
    for fname in datad:
        fhindex.write('    %s\n' % (fname[:-2]+'rst'))
    fhindex.flush()

def setup(app):
    app.connect('builder-inited', generate_example_rst)

