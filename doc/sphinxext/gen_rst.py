"""
Example generation for the scikit learn

Generate the rst files for the examples by iterating over the python
example files.

Files that generate images should start with 'plot'

"""
import os

fileList = []

import matplotlib
matplotlib.use('Agg')
import IPython.Shell
mplshell = IPython.Shell.MatplotlibShell('mpl')
                          
import token, tokenize      

rst_template = """

.. %(short_fname)s_example:

%(docstring)s

**Source code:** :download:`%(fname)s <%(short_fname)s>`

.. literalinclude:: %(short_fname)s
    :lines: %(end_row)s-
    """

plot_rst_template = """

.. %(short_fname)s_example:

%(docstring)s

.. image:: images/%(image_name)s
    :align: center

**Source code:** :download:`%(fname)s <%(short_fname)s>`

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
            image_name = fname[:-2] + 'png'
            global rst_template, plot_rst_template
            this_template = rst_template
            short_fname = '../../examples/' + fname
            if  not fname.endswith('py'): 
                continue
            example_file = os.path.join(exampledir, fname)
            if fname.startswith('plot'):
                # generate the plot as png image if file name
                # starts with plot and if it is more recent than an
                # existing image.
                if not os.path.exists(
                                    os.path.join(rootdir, 'images')):
                    os.makedirs(os.path.join(rootdir, 'images'))
                image_file = os.path.join(rootdir, 'images', image_name)
                if (not os.path.exists(image_file) or
                      os.stat(image_file).st_mtime <= 
                            os.stat(example_file).st_mtime):
                    print 'plotting %s' % fname
                    import matplotlib.pyplot as plt
                    plt.close('all')
                    mplshell.magic_run(example_file)
                    plt.savefig(image_file)
                this_template = plot_rst_template

            docstring, short_desc, end_row = extract_docstring(example_file)

            f = open(os.path.join(rootdir, fname[:-2] + 'rst'),'w')
            f.write( this_template % locals())
            f.flush()
            datad.append(fname)

    # we create an index.rst with all examples
    fhindex = file(os.path.join(rootdir, 'index.rst'), 'w')
    fhindex.write("""\
.. _examples-index:

Examples
==========

    :Release: |version|
    :Date: |today|

.. toctree::

""")
    
    for fname in datad:
        fhindex.write('    %s\n' % (fname[:-3]))
    fhindex.flush()

def setup(app):
    app.connect('builder-inited', generate_example_rst)

