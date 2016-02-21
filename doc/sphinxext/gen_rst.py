"""
Example generation for the scikit learn

Generate the rst files for the examples by iterating over the python
example files.

Files that generate images should start with 'plot'

"""
from __future__ import division, print_function
from time import time
import ast
import os
import re
import shutil
import traceback
import glob
import sys
import gzip
import posixpath
import subprocess
import warnings
from sklearn.externals import six


# Try Python 2 first, otherwise load from Python 3
try:
    from StringIO import StringIO
    import cPickle as pickle
    import urllib2 as urllib
    from urllib2 import HTTPError, URLError
except ImportError:
    from io import StringIO
    import pickle
    import urllib.request
    import urllib.error
    import urllib.parse
    from urllib.error import HTTPError, URLError


try:
    # Python 2 built-in
    execfile
except NameError:
    def execfile(filename, global_vars=None, local_vars=None):
        with open(filename, encoding='utf-8') as f:
            code = compile(f.read(), filename, 'exec')
            exec(code, global_vars, local_vars)

try:
    basestring
except NameError:
    basestring = str

import token
import tokenize
import numpy as np

try:
    # make sure that the Agg backend is set before importing any
    # matplotlib
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    # this script can be imported by nosetest to find tests to run: we should not
    # impose the matplotlib requirement in that case.
    pass


from sklearn.externals import joblib

###############################################################################
# A tee object to redict streams to multiple outputs

class Tee(object):

    def __init__(self, file1, file2):
        self.file1 = file1
        self.file2 = file2

    def write(self, data):
        self.file1.write(data)
        self.file2.write(data)

    def flush(self):
        self.file1.flush()
        self.file2.flush()

###############################################################################
# Documentation link resolver objects


def _get_data(url):
    """Helper function to get data over http or from a local file"""
    if url.startswith('http://'):
        # Try Python 2, use Python 3 on exception
        try:
            resp = urllib.urlopen(url)
            encoding = resp.headers.dict.get('content-encoding', 'plain')
        except AttributeError:
            resp = urllib.request.urlopen(url)
            encoding = resp.headers.get('content-encoding', 'plain')
        data = resp.read()
        if encoding == 'plain':
            pass
        elif encoding == 'gzip':
            data = StringIO(data)
            data = gzip.GzipFile(fileobj=data).read()
        else:
            raise RuntimeError('unknown encoding')
    else:
        with open(url, 'r') as fid:
            data = fid.read()
        fid.close()

    return data

mem = joblib.Memory(cachedir='_build')
get_data = mem.cache(_get_data)


def parse_sphinx_searchindex(searchindex):
    """Parse a Sphinx search index

    Parameters
    ----------
    searchindex : str
        The Sphinx search index (contents of searchindex.js)

    Returns
    -------
    filenames : list of str
        The file names parsed from the search index.
    objects : dict
        The objects parsed from the search index.
    """
    def _select_block(str_in, start_tag, end_tag):
        """Select first block delimited by start_tag and end_tag"""
        start_pos = str_in.find(start_tag)
        if start_pos < 0:
            raise ValueError('start_tag not found')
        depth = 0
        for pos in range(start_pos, len(str_in)):
            if str_in[pos] == start_tag:
                depth += 1
            elif str_in[pos] == end_tag:
                depth -= 1

            if depth == 0:
                break
        sel = str_in[start_pos + 1:pos]
        return sel

    def _parse_dict_recursive(dict_str):
        """Parse a dictionary from the search index"""
        dict_out = dict()
        pos_last = 0
        pos = dict_str.find(':')
        while pos >= 0:
            key = dict_str[pos_last:pos]
            if dict_str[pos + 1] == '[':
                # value is a list
                pos_tmp = dict_str.find(']', pos + 1)
                if pos_tmp < 0:
                    raise RuntimeError('error when parsing dict')
                value = dict_str[pos + 2: pos_tmp].split(',')
                # try to convert elements to int
                for i in range(len(value)):
                    try:
                        value[i] = int(value[i])
                    except ValueError:
                        pass
            elif dict_str[pos + 1] == '{':
                # value is another dictionary
                subdict_str = _select_block(dict_str[pos:], '{', '}')
                value = _parse_dict_recursive(subdict_str)
                pos_tmp = pos + len(subdict_str)
            else:
                raise ValueError('error when parsing dict: unknown elem')

            key = key.strip('"')
            if len(key) > 0:
                dict_out[key] = value

            pos_last = dict_str.find(',', pos_tmp)
            if pos_last < 0:
                break
            pos_last += 1
            pos = dict_str.find(':', pos_last)

        return dict_out

    # Make sure searchindex uses UTF-8 encoding
    if hasattr(searchindex, 'decode'):
        searchindex = searchindex.decode('UTF-8')

    # parse objects
    query = 'objects:'
    pos = searchindex.find(query)
    if pos < 0:
        raise ValueError('"objects:" not found in search index')

    sel = _select_block(searchindex[pos:], '{', '}')
    objects = _parse_dict_recursive(sel)

    # parse filenames
    query = 'filenames:'
    pos = searchindex.find(query)
    if pos < 0:
        raise ValueError('"filenames:" not found in search index')
    filenames = searchindex[pos + len(query) + 1:]
    filenames = filenames[:filenames.find(']')]
    filenames = [f.strip('"') for f in filenames.split(',')]

    return filenames, objects


class SphinxDocLinkResolver(object):
    """ Resolve documentation links using searchindex.js generated by Sphinx

    Parameters
    ----------
    doc_url : str
        The base URL of the project website.
    searchindex : str
        Filename of searchindex, relative to doc_url.
    extra_modules_test : list of str
        List of extra module names to test.
    relative : bool
        Return relative links (only useful for links to documentation of this
        package).
    """

    def __init__(self, doc_url, searchindex='searchindex.js',
                 extra_modules_test=None, relative=False):
        self.doc_url = doc_url
        self.relative = relative
        self._link_cache = {}

        self.extra_modules_test = extra_modules_test
        self._page_cache = {}
        if doc_url.startswith('http://'):
            if relative:
                raise ValueError('Relative links are only supported for local '
                                 'URLs (doc_url cannot start with "http://)"')
            searchindex_url = doc_url + '/' + searchindex
        else:
            searchindex_url = os.path.join(doc_url, searchindex)

        # detect if we are using relative links on a Windows system
        if os.name.lower() == 'nt' and not doc_url.startswith('http://'):
            if not relative:
                raise ValueError('You have to use relative=True for the local'
                                 ' package on a Windows system.')
            self._is_windows = True
        else:
            self._is_windows = False

        # download and initialize the search index
        sindex = get_data(searchindex_url)
        filenames, objects = parse_sphinx_searchindex(sindex)

        self._searchindex = dict(filenames=filenames, objects=objects)

    def _get_link(self, cobj):
        """Get a valid link, False if not found"""

        fname_idx = None
        full_name = cobj['module_short'] + '.' + cobj['name']
        if full_name in self._searchindex['objects']:
            value = self._searchindex['objects'][full_name]
            if isinstance(value, dict):
                value = value[next(iter(value.keys()))]
            fname_idx = value[0]
        elif cobj['module_short'] in self._searchindex['objects']:
            value = self._searchindex['objects'][cobj['module_short']]
            if cobj['name'] in value.keys():
                fname_idx = value[cobj['name']][0]

        if fname_idx is not None:
            fname = self._searchindex['filenames'][fname_idx] + '.html'

            if self._is_windows:
                fname = fname.replace('/', '\\')
                link = os.path.join(self.doc_url, fname)
            else:
                link = posixpath.join(self.doc_url, fname)

            if hasattr(link, 'decode'):
                link = link.decode('utf-8', 'replace')

            if link in self._page_cache:
                html = self._page_cache[link]
            else:
                html = get_data(link)
                self._page_cache[link] = html

            # test if cobj appears in page
            comb_names = [cobj['module_short'] + '.' + cobj['name']]
            if self.extra_modules_test is not None:
                for mod in self.extra_modules_test:
                    comb_names.append(mod + '.' + cobj['name'])
            url = False
            if hasattr(html, 'decode'):
                # Decode bytes under Python 3
                html = html.decode('utf-8', 'replace')

            for comb_name in comb_names:
                if hasattr(comb_name, 'decode'):
                    # Decode bytes under Python 3
                    comb_name = comb_name.decode('utf-8', 'replace')
                if comb_name in html:
                    url = link + u'#' + comb_name
            link = url
        else:
            link = False

        return link

    def resolve(self, cobj, this_url):
        """Resolve the link to the documentation, returns None if not found

        Parameters
        ----------
        cobj : dict
            Dict with information about the "code object" for which we are
            resolving a link.
            cobi['name'] : function or class name (str)
            cobj['module_short'] : shortened module name (str)
            cobj['module'] : module name (str)
        this_url: str
            URL of the current page. Needed to construct relative URLs
            (only used if relative=True in constructor).

        Returns
        -------
        link : str | None
            The link (URL) to the documentation.
        """
        full_name = cobj['module_short'] + '.' + cobj['name']
        link = self._link_cache.get(full_name, None)
        if link is None:
            # we don't have it cached
            link = self._get_link(cobj)
            # cache it for the future
            self._link_cache[full_name] = link

        if link is False or link is None:
            # failed to resolve
            return None

        if self.relative:
            link = os.path.relpath(link, start=this_url)
            if self._is_windows:
                # replace '\' with '/' so it on the web
                link = link.replace('\\', '/')

            # for some reason, the relative link goes one directory too high up
            link = link[3:]

        return link


###############################################################################
rst_template = """

.. _example_%(short_fname)s:

%(docstring)s

**Python source code:** :download:`%(fname)s <%(fname)s>`

.. literalinclude:: %(fname)s
    :lines: %(end_row)s-
    """

plot_rst_template = """

.. _example_%(short_fname)s:

%(docstring)s

%(image_list)s

%(stdout)s

**Python source code:** :download:`%(fname)s <%(fname)s>`

.. literalinclude:: %(fname)s
    :lines: %(end_row)s-

**Total running time of the example:** %(time_elapsed) .2f seconds
(%(time_m) .0f minutes %(time_s) .2f seconds)
    """

# The following strings are used when we have several pictures: we use
# an html div tag that our CSS uses to turn the lists into horizontal
# lists.
HLIST_HEADER = """
.. rst-class:: horizontal

"""

HLIST_IMAGE_TEMPLATE = """
    *

      .. image:: images/%s
            :scale: 47
"""

SINGLE_IMAGE = """
.. image:: images/%s
    :align: center
"""

# The following dictionary contains the information used to create the
# thumbnails for the front page of the scikit-learn home page.
# key: first image in set
# values: (number of plot in set, height of thumbnail)
carousel_thumbs = {'plot_classifier_comparison_001.png': (1, 600),
                   'plot_outlier_detection_001.png': (3, 372),
                   'plot_gpr_co2_001.png': (1, 350),
                   'plot_adaboost_twoclass_001.png': (1, 372),
                   'plot_compare_methods_001.png': (1, 349)}


def extract_docstring(filename, ignore_heading=False):
    """ Extract a module-level docstring, if any
    """
    if six.PY2:
        lines = open(filename).readlines()
    else:
        lines = open(filename, encoding='utf-8').readlines()
    start_row = 0
    if lines[0].startswith('#!'):
        lines.pop(0)
        start_row = 1
    docstring = ''
    first_par = ''
    line_iterator = iter(lines)
    tokens = tokenize.generate_tokens(lambda: next(line_iterator))
    for tok_type, tok_content, _, (erow, _), _ in tokens:
        tok_type = token.tok_name[tok_type]
        if tok_type in ('NEWLINE', 'COMMENT', 'NL', 'INDENT', 'DEDENT'):
            continue
        elif tok_type == 'STRING':
            docstring = eval(tok_content)
            # If the docstring is formatted with several paragraphs, extract
            # the first one:
            paragraphs = '\n'.join(
                line.rstrip() for line
                in docstring.split('\n')).split('\n\n')
            if paragraphs:
                if ignore_heading:
                    if len(paragraphs) > 1:
                        first_par = re.sub('\n', ' ', paragraphs[1])
                        first_par = ((first_par[:95] + '...')
                                     if len(first_par) > 95 else first_par)
                    else:
                        raise ValueError("Docstring not found by gallery.\n"
                                         "Please check the layout of your"
                                         " example file:\n {}\n and make sure"
                                         " it's correct".format(filename))
                else:
                    first_par = paragraphs[0]

        break
    return docstring, first_par, erow + 1 + start_row


def generate_example_rst(app):
    """ Generate the list of examples, as well as the contents of
        examples.
    """
    root_dir = os.path.join(app.builder.srcdir, 'auto_examples')
    example_dir = os.path.abspath(os.path.join(app.builder.srcdir, '..',
                                               'examples'))
    generated_dir = os.path.abspath(os.path.join(app.builder.srcdir,
                                                 'modules', 'generated'))

    try:
        plot_gallery = eval(app.builder.config.plot_gallery)
    except TypeError:
        plot_gallery = bool(app.builder.config.plot_gallery)
    if not os.path.exists(example_dir):
        os.makedirs(example_dir)
    if not os.path.exists(root_dir):
        os.makedirs(root_dir)
    if not os.path.exists(generated_dir):
        os.makedirs(generated_dir)

    # we create an index.rst with all examples
    fhindex = open(os.path.join(root_dir, 'index.rst'), 'w')
    # Note: The sidebar button has been removed from the examples page for now
    #      due to how it messes up the layout. Will be fixed at a later point
    fhindex.write("""\



.. raw:: html


    <style type="text/css">
    div#sidebarbutton {
        /* hide the sidebar collapser, while ensuring vertical arrangement */
        display: none;
    }
    </style>

.. _examples-index:

Examples
========

""")
    # Here we don't use an os.walk, but we recurse only twice: flat is
    # better than nested.
    seen_backrefs = set()
    generate_dir_rst('.', fhindex, example_dir, root_dir, plot_gallery, seen_backrefs)
    for directory in sorted(os.listdir(example_dir)):
        if os.path.isdir(os.path.join(example_dir, directory)):
            generate_dir_rst(directory, fhindex, example_dir, root_dir, plot_gallery, seen_backrefs)
    fhindex.flush()


def extract_line_count(filename, target_dir):
    # Extract the line count of a file
    example_file = os.path.join(target_dir, filename)
    if six.PY2:
        lines = open(example_file).readlines()
    else:
        lines = open(example_file, encoding='utf-8').readlines()
    start_row = 0
    if lines and lines[0].startswith('#!'):
        lines.pop(0)
        start_row = 1
    line_iterator = iter(lines)
    tokens = tokenize.generate_tokens(lambda: next(line_iterator))
    check_docstring = True
    erow_docstring = 0
    for tok_type, _, _, (erow, _), _ in tokens:
        tok_type = token.tok_name[tok_type]
        if tok_type in ('NEWLINE', 'COMMENT', 'NL', 'INDENT', 'DEDENT'):
            continue
        elif (tok_type == 'STRING') and check_docstring:
            erow_docstring = erow
            check_docstring = False
    return erow_docstring+1+start_row, erow+1+start_row


def line_count_sort(file_list, target_dir):
    # Sort the list of examples by line-count
    new_list = [x for x in file_list if x.endswith('.py')]
    unsorted = np.zeros(shape=(len(new_list), 2))
    unsorted = unsorted.astype(np.object)
    for count, exmpl in enumerate(new_list):
        docstr_lines, total_lines = extract_line_count(exmpl, target_dir)
        unsorted[count][1] = total_lines - docstr_lines
        unsorted[count][0] = exmpl
    index = np.lexsort((unsorted[:, 0].astype(np.str),
                        unsorted[:, 1].astype(np.float)))
    if not len(unsorted):
        return []
    return np.array(unsorted[index][:, 0]).tolist()


def _thumbnail_div(subdir, full_dir, fname, snippet, is_backref=False):
    """Generates RST to place a thumbnail in a gallery"""
    thumb = os.path.join(full_dir, 'images', 'thumb', fname[:-3] + '.png')
    link_name = os.path.join(full_dir, fname).replace(os.path.sep, '_')
    ref_name = os.path.join(subdir, fname).replace(os.path.sep, '_')
    if ref_name.startswith('._'):
        ref_name = ref_name[2:]
    out = []
    out.append("""

.. raw:: html

    <div class="thumbnailContainer" tooltip="{}">

""".format(snippet))

    out.append('.. only:: html\n\n')
    out.append('  .. figure:: %s\n' % thumb)
    if link_name.startswith('._'):
        link_name = link_name[2:]
    if full_dir != '.':
        out.append('    :target: ./%s/%s.html\n\n' % (full_dir, fname[:-3]))
    else:
        out.append('    :target: ./%s.html\n\n' % link_name[:-3])
    out.append("""    :ref:`example_%s`


.. raw:: html

    </div>

""" % (ref_name))
    if is_backref:
        out.append('.. only:: not html\n\n  * :ref:`example_%s`' % ref_name)
    return ''.join(out)


def generate_dir_rst(directory, fhindex, example_dir, root_dir, plot_gallery, seen_backrefs):
    """ Generate the rst file for an example directory.
    """
    if not directory == '.':
        target_dir = os.path.join(root_dir, directory)
        src_dir = os.path.join(example_dir, directory)
    else:
        target_dir = root_dir
        src_dir = example_dir
    if not os.path.exists(os.path.join(src_dir, 'README.txt')):
        raise ValueError('Example directory %s does not have a README.txt' %
                         src_dir)

    fhindex.write("""


%s


""" % open(os.path.join(src_dir, 'README.txt')).read())
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    sorted_listdir = line_count_sort(os.listdir(src_dir),
                                     src_dir)
    if not os.path.exists(os.path.join(directory, 'images', 'thumb')):
        os.makedirs(os.path.join(directory, 'images', 'thumb'))
    for fname in sorted_listdir:
        if fname.endswith('py'):
            backrefs = generate_file_rst(fname, target_dir, src_dir, root_dir, plot_gallery)
            new_fname = os.path.join(src_dir, fname)
            _, snippet, _ = extract_docstring(new_fname, True)
            fhindex.write(_thumbnail_div(directory, directory, fname, snippet))
            fhindex.write("""

.. toctree::
   :hidden:

   %s/%s

""" % (directory, fname[:-3]))
            for backref in backrefs:
                include_path = os.path.join(root_dir, '../modules/generated/%s.examples' % backref)
                seen = backref in seen_backrefs
                with open(include_path, 'a' if seen else 'w') as ex_file:
                    if not seen:
                        # heading
                        print(file=ex_file)
                        print('Examples using ``%s``' % backref, file=ex_file)
                        print('-----------------%s--' % ('-' * len(backref)),
                              file=ex_file)
                        print(file=ex_file)
                    rel_dir = os.path.join('../../auto_examples', directory)
                    ex_file.write(_thumbnail_div(directory, rel_dir, fname, snippet, is_backref=True))
                    seen_backrefs.add(backref)
    fhindex.write("""
.. raw:: html

    <div class="clearer"></div>
    """)  # clear at the end of the section

# modules for which we embed links into example code
DOCMODULES = ['sklearn', 'matplotlib', 'numpy', 'scipy']


def make_thumbnail(in_fname, out_fname, width, height):
    """Make a thumbnail with the same aspect ratio centered in an
       image with a given width and height
    """
    # local import to avoid testing dependency on PIL:
    try:
        from PIL import Image
    except ImportError:
        import Image
    img = Image.open(in_fname)
    width_in, height_in = img.size
    scale_w = width / float(width_in)
    scale_h = height / float(height_in)

    if height_in * scale_w <= height:
        scale = scale_w
    else:
        scale = scale_h

    width_sc = int(round(scale * width_in))
    height_sc = int(round(scale * height_in))

    # resize the image
    img.thumbnail((width_sc, height_sc), Image.ANTIALIAS)

    # insert centered
    thumb = Image.new('RGB', (width, height), (255, 255, 255))
    pos_insert = ((width - width_sc) // 2, (height - height_sc) // 2)
    thumb.paste(img, pos_insert)

    thumb.save(out_fname)
    # Use optipng to perform lossless compression on the resized image if
    # software is installed
    if os.environ.get('SKLEARN_DOC_OPTIPNG', False):
        try:
            subprocess.call(["optipng", "-quiet", "-o", "9", out_fname])
        except Exception:
            warnings.warn('Install optipng to reduce the size of the generated images')


def get_short_module_name(module_name, obj_name):
    """ Get the shortest possible module name """
    parts = module_name.split('.')
    short_name = module_name
    for i in range(len(parts) - 1, 0, -1):
        short_name = '.'.join(parts[:i])
        try:
            exec('from %s import %s' % (short_name, obj_name))
        except ImportError:
            # get the last working module name
            short_name = '.'.join(parts[:(i + 1)])
            break
    return short_name


class NameFinder(ast.NodeVisitor):
    """Finds the longest form of variable names and their imports in code

    Only retains names from imported modules.
    """

    def __init__(self):
        super(NameFinder, self).__init__()
        self.imported_names = {}
        self.accessed_names = set()

    def visit_Import(self, node, prefix=''):
        for alias in node.names:
            local_name = alias.asname or alias.name
            self.imported_names[local_name] = prefix + alias.name

    def visit_ImportFrom(self, node):
        self.visit_Import(node, node.module + '.')

    def visit_Name(self, node):
        self.accessed_names.add(node.id)

    def visit_Attribute(self, node):
        attrs = []
        while isinstance(node, ast.Attribute):
            attrs.append(node.attr)
            node = node.value

        if isinstance(node, ast.Name):
            # This is a.b, not e.g. a().b
            attrs.append(node.id)
            self.accessed_names.add('.'.join(reversed(attrs)))
        else:
            # need to get a in a().b
            self.visit(node)

    def get_mapping(self):
        for name in self.accessed_names:
            local_name = name.split('.', 1)[0]
            remainder = name[len(local_name):]
            if local_name in self.imported_names:
                # Join import path to relative path
                full_name = self.imported_names[local_name] + remainder
                yield name, full_name


def identify_names(code):
    """Builds a codeobj summary by identifying and resovles used names

    >>> code = '''
    ... from a.b import c
    ... import d as e
    ... print(c)
    ... e.HelloWorld().f.g
    ... '''
    >>> for name, o in sorted(identify_names(code).items()):
    ...     print(name, o['name'], o['module'], o['module_short'])
    c c a.b a.b
    e.HelloWorld HelloWorld d d
    """
    finder = NameFinder()
    finder.visit(ast.parse(code))

    example_code_obj = {}
    for name, full_name in finder.get_mapping():
        # name is as written in file (e.g. np.asarray)
        # full_name includes resolved import path (e.g. numpy.asarray)
        module, attribute = full_name.rsplit('.', 1)
        # get shortened module name
        module_short = get_short_module_name(module, attribute)
        cobj = {'name': attribute, 'module': module,
                'module_short': module_short}
        example_code_obj[name] = cobj
    return example_code_obj


def generate_file_rst(fname, target_dir, src_dir, root_dir, plot_gallery):
    """ Generate the rst file for a given example.

    Returns the set of sklearn functions/classes imported in the example.
    """
    base_image_name = os.path.splitext(fname)[0]
    image_fname = '%s_%%03d.png' % base_image_name

    this_template = rst_template
    last_dir = os.path.split(src_dir)[-1]
    # to avoid leading . in file names, and wrong names in links
    if last_dir == '.' or last_dir == 'examples':
        last_dir = ''
    else:
        last_dir += '_'
    short_fname = last_dir + fname
    src_file = os.path.join(src_dir, fname)
    example_file = os.path.join(target_dir, fname)
    shutil.copyfile(src_file, example_file)

    # The following is a list containing all the figure names
    figure_list = []

    image_dir = os.path.join(target_dir, 'images')
    thumb_dir = os.path.join(image_dir, 'thumb')
    if not os.path.exists(image_dir):
        os.makedirs(image_dir)
    if not os.path.exists(thumb_dir):
        os.makedirs(thumb_dir)
    image_path = os.path.join(image_dir, image_fname)
    stdout_path = os.path.join(image_dir,
                               'stdout_%s.txt' % base_image_name)
    time_path = os.path.join(image_dir,
                             'time_%s.txt' % base_image_name)
    thumb_file = os.path.join(thumb_dir, base_image_name + '.png')
    time_elapsed = 0
    if plot_gallery and fname.startswith('plot'):
        # generate the plot as png image if file name
        # starts with plot and if it is more recent than an
        # existing image.
        first_image_file = image_path % 1
        if os.path.exists(stdout_path):
            stdout = open(stdout_path).read()
        else:
            stdout = ''
        if os.path.exists(time_path):
            time_elapsed = float(open(time_path).read())

        if not os.path.exists(first_image_file) or \
           os.stat(first_image_file).st_mtime <= os.stat(src_file).st_mtime:
            # We need to execute the code
            print('plotting %s' % fname)
            t0 = time()
            import matplotlib.pyplot as plt
            plt.close('all')
            cwd = os.getcwd()
            try:
                # First CD in the original example dir, so that any file
                # created by the example get created in this directory
                orig_stdout = sys.stdout
                os.chdir(os.path.dirname(src_file))
                my_buffer = StringIO()
                my_stdout = Tee(sys.stdout, my_buffer)
                sys.stdout = my_stdout
                my_globals = {'pl': plt}
                execfile(os.path.basename(src_file), my_globals)
                time_elapsed = time() - t0
                sys.stdout = orig_stdout
                my_stdout = my_buffer.getvalue()

                if '__doc__' in my_globals:
                    # The __doc__ is often printed in the example, we
                    # don't with to echo it
                    my_stdout = my_stdout.replace(
                        my_globals['__doc__'],
                        '')
                my_stdout = my_stdout.strip().expandtabs()
                if my_stdout:
                    stdout = '**Script output**::\n\n  %s\n\n' % (
                        '\n  '.join(my_stdout.split('\n')))
                open(stdout_path, 'w').write(stdout)
                open(time_path, 'w').write('%f' % time_elapsed)
                os.chdir(cwd)

                # In order to save every figure we have two solutions :
                # * iterate from 1 to infinity and call plt.fignum_exists(n)
                #   (this requires the figures to be numbered
                #    incrementally: 1, 2, 3 and not 1, 2, 5)
                # * iterate over [fig_mngr.num for fig_mngr in
                #   matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
                fig_managers = matplotlib._pylab_helpers.Gcf.get_all_fig_managers()
                for fig_mngr in fig_managers:
                    # Set the fig_num figure as the current figure as we can't
                    # save a figure that's not the current figure.
                    fig = plt.figure(fig_mngr.num)
                    kwargs = {}
                    to_rgba = matplotlib.colors.colorConverter.to_rgba
                    for attr in ['facecolor', 'edgecolor']:
                        fig_attr = getattr(fig, 'get_' + attr)()
                        default_attr = matplotlib.rcParams['figure.' + attr]
                        if to_rgba(fig_attr) != to_rgba(default_attr):
                            kwargs[attr] = fig_attr

                    fig.savefig(image_path % fig_mngr.num, **kwargs)
                    figure_list.append(image_fname % fig_mngr.num)
            except:
                print(80 * '_')
                print('%s is not compiling:' % fname)
                traceback.print_exc()
                print(80 * '_')
            finally:
                os.chdir(cwd)
                sys.stdout = orig_stdout

            print(" - time elapsed : %.2g sec" % time_elapsed)
        else:
            figure_list = [f[len(image_dir):]
                           for f in glob.glob(image_path.replace("%03d",
                                                '[0-9][0-9][0-9]'))]
        figure_list.sort()

        # generate thumb file
        this_template = plot_rst_template
        car_thumb_path = os.path.join(os.path.split(root_dir)[0], '_build/html/stable/_images/')
        # Note: normaly, make_thumbnail is used to write to the path contained in `thumb_file`
        # which is within `auto_examples/../images/thumbs` depending on the example.
        # Because the carousel has different dimensions than those of the examples gallery,
        # I did not simply reuse them all as some contained whitespace due to their default gallery
        # thumbnail size. Below, for a few cases, seperate thumbnails are created (the originals can't
        # just be overwritten with the carousel dimensions as it messes up the examples gallery layout).
        # The special carousel thumbnails are written directly to _build/html/stable/_images/,
        # as for some reason unknown to me, Sphinx refuses to copy my 'extra' thumbnails from the
        # auto examples gallery to the _build folder. This works fine as is, but it would be cleaner to
        # have it happen with the rest. Ideally the should be written to 'thumb_file' as well, and then
        # copied to the _images folder during the `Copying Downloadable Files` step like the rest.
        if not os.path.exists(car_thumb_path):
            os.makedirs(car_thumb_path)
        if os.path.exists(first_image_file):
            # We generate extra special thumbnails for the carousel
            carousel_tfile = os.path.join(car_thumb_path, base_image_name + '_carousel.png')
            first_img = image_fname % 1
            if first_img in carousel_thumbs:
                make_thumbnail((image_path % carousel_thumbs[first_img][0]),
                               carousel_tfile, carousel_thumbs[first_img][1], 190)
            make_thumbnail(first_image_file, thumb_file, 400, 280)

    if not os.path.exists(thumb_file):
        # create something to replace the thumbnail
        make_thumbnail('images/no_image.png', thumb_file, 200, 140)

    docstring, short_desc, end_row = extract_docstring(example_file)

    # Depending on whether we have one or more figures, we're using a
    # horizontal list or a single rst call to 'image'.
    if len(figure_list) == 1:
        figure_name = figure_list[0]
        image_list = SINGLE_IMAGE % figure_name.lstrip('/')
    else:
        image_list = HLIST_HEADER
        for figure_name in figure_list:
            image_list += HLIST_IMAGE_TEMPLATE % figure_name.lstrip('/')

    time_m, time_s = divmod(time_elapsed, 60)
    f = open(os.path.join(target_dir, base_image_name + '.rst'), 'w')
    f.write(this_template % locals())
    f.flush()

    # save variables so we can later add links to the documentation
    if six.PY2:
        example_code_obj = identify_names(open(example_file).read())
    else:
        example_code_obj = \
            identify_names(open(example_file, encoding='utf-8').read())
    if example_code_obj:
        codeobj_fname = example_file[:-3] + '_codeobj.pickle'
        with open(codeobj_fname, 'wb') as fid:
            pickle.dump(example_code_obj, fid, pickle.HIGHEST_PROTOCOL)

    backrefs = set('{module_short}.{name}'.format(**entry)
                   for entry in example_code_obj.values()
                   if entry['module'].startswith('sklearn'))
    return backrefs


def embed_code_links(app, exception):
    """Embed hyperlinks to documentation into example code"""
    if exception is not None:
        return
    print('Embedding documentation hyperlinks in examples..')

    if app.builder.name == 'latex':
        # Don't embed hyperlinks when a latex builder is used.
        return

    # Add resolvers for the packages for which we want to show links
    doc_resolvers = {}
    doc_resolvers['sklearn'] = SphinxDocLinkResolver(app.builder.outdir,
                                                     relative=True)

    resolver_urls = {
        'matplotlib': 'http://matplotlib.org',
        'numpy': 'http://docs.scipy.org/doc/numpy-1.6.0',
        'scipy': 'http://docs.scipy.org/doc/scipy-0.11.0/reference',
    }
    for this_module, url in resolver_urls.items():
        try:
            doc_resolvers[this_module] = SphinxDocLinkResolver(url)
        except HTTPError as e:
            print("The following HTTP Error has occurred:\n")
            print(e.code)
        except URLError as e:
            print("\n...\n"
                  "Warning: Embedding the documentation hyperlinks requires "
                  "internet access.\nPlease check your network connection.\n"
                  "Unable to continue embedding `{0}` links due to a URL "
                  "Error:\n".format(this_module))
            print(e.args)

    example_dir = os.path.join(app.builder.srcdir, 'auto_examples')
    html_example_dir = os.path.abspath(os.path.join(app.builder.outdir,
                                                    'auto_examples'))

    # patterns for replacement
    link_pattern = '<a href="%s">%s</a>'
    orig_pattern = '<span class="n">%s</span>'
    period = '<span class="o">.</span>'

    for dirpath, _, filenames in os.walk(html_example_dir):
        for fname in filenames:
            print('\tprocessing: %s' % fname)
            full_fname = os.path.join(html_example_dir, dirpath, fname)
            subpath = dirpath[len(html_example_dir) + 1:]
            pickle_fname = os.path.join(example_dir, subpath,
                                        fname[:-5] + '_codeobj.pickle')

            if os.path.exists(pickle_fname):
                # we have a pickle file with the objects to embed links for
                with open(pickle_fname, 'rb') as fid:
                    example_code_obj = pickle.load(fid)
                fid.close()
                str_repl = {}
                # generate replacement strings with the links
                for name, cobj in example_code_obj.items():
                    this_module = cobj['module'].split('.')[0]

                    if this_module not in doc_resolvers:
                        continue

                    try:
                        link = doc_resolvers[this_module].resolve(cobj,
                                                                  full_fname)
                    except (HTTPError, URLError) as e:
                        print("The following error has occurred:\n")
                        print(repr(e))
                        continue

                    if link is not None:
                        parts = name.split('.')
                        name_html = period.join(orig_pattern % part
                                                for part in parts)
                        str_repl[name_html] = link_pattern % (link, name_html)
                # do the replacement in the html file

                # ensure greediness
                names = sorted(str_repl, key=len, reverse=True)
                expr = re.compile(r'(?<!\.)\b' +  # don't follow . or word
                                  '|'.join(re.escape(name)
                                           for name in names))

                def substitute_link(match):
                    return str_repl[match.group()]

                if len(str_repl) > 0:
                    with open(full_fname, 'rb') as fid:
                        lines_in = fid.readlines()
                    with open(full_fname, 'wb') as fid:
                        for line in lines_in:
                            line = line.decode('utf-8')
                            line = expr.sub(substitute_link, line)
                            fid.write(line.encode('utf-8'))
    print('[done]')


def setup(app):
    app.connect('builder-inited', generate_example_rst)
    app.add_config_value('plot_gallery', True, 'html')

    # embed links after build is finished
    app.connect('build-finished', embed_code_links)

    # Sphinx hack: sphinx copies generated images to the build directory
    #  each time the docs are made.  If the desired image name already
    #  exists, it appends a digit to prevent overwrites.  The problem is,
    #  the directory is never cleared.  This means that each time you build
    #  the docs, the number of images in the directory grows.
    #
    # This question has been asked on the sphinx development list, but there
    #  was no response: http://osdir.com/ml/sphinx-dev/2011-02/msg00123.html
    #
    # The following is a hack that prevents this behavior by clearing the
    #  image build directory each time the docs are built.  If sphinx
    #  changes their layout between versions, this will not work (though
    #  it should probably not cause a crash).  Tested successfully
    #  on Sphinx 1.0.7
    build_image_dir = '_build/html/_images'
    if os.path.exists(build_image_dir):
        filelist = os.listdir(build_image_dir)
        for filename in filelist:
            if filename.endswith('png'):
                os.remove(os.path.join(build_image_dir, filename))


def setup_module():
    # HACK: Stop nosetests running setup() above
    pass
