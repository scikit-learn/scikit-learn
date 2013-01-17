"""
Example generation for the scikit learn

Generate the rst files for the examples by iterating over the python
example files.

Files that generate images should start with 'plot'

"""
from time import time
import os
import shutil
import traceback
import glob
import sys
from StringIO import StringIO
import cPickle
import re
import urllib2
import gzip

try:
    from PIL import Image
except:
    import Image

import matplotlib
matplotlib.use('Agg')

import token
import tokenize

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


class DocLinkResolver(object):
    """ Resolver for packages that document each object/function on an
        individual page """
    def __init__(self, doc_base_url, relative=False):
        if relative and doc_base_url.startswith('http'):
            raise ValueError('relative and http url not supported')
        self.doc_base_url = doc_base_url
        self.relative = relative
        self._link_cache = {}

    def get_link(self, cobj):
        """Get a valid link, False if not found"""
        link = (self.doc_base_url + '/' + cobj['module_short'] + '.'
                + cobj['name'] + '.html')
        if link.startswith('http://'):
            try:
                resp = urllib2.urlopen(link)
                if resp.code != 200:
                    link = False
            except urllib2.HTTPError:
                link = False
        else:
            # assume it is a local file
            if not os.path.exists(link):
                link = False

        if link is not False:
            # make a link that highlights the occurence
            link = link + '#' + cobj['module_short'] + '.' + cobj['name']

        return link

    def resolve(self, cobj, this_url):
        """Resolve the link to the documentation, returns None if not found"""
        link = self._link_cache.get(cobj['name'], None)
        if link is None:
            # we don't have it cached
            link = self.get_link(cobj)
            # cache it for the future
            self._link_cache[cobj['name']] = link

        if link is False or link is None:
            # failed to resolve
            return None

        if self.relative:
            link = os.path.relpath(link, start=this_url)
            # for some reason, the relative link goes one directory too high up
            link = link[3:]

        return link


class SinglePageDocLinkResolver(DocLinkResolver):
    """ Resolver for packages that use pages that document several items """

    def __init__(self,  doc_pages_url, extra_modules_test=None):
        super(SinglePageDocLinkResolver, self).__init__(self, None)
        # get all the pages
        doc_pages_html = []
        for url in doc_pages_url:
            try:
                resp = urllib2.urlopen(url)
                encoding = resp.headers.dict.get('content-encoding', 'plain')
                data = resp.read()
                if encoding == 'plain':
                    html = resp.read()
                elif encoding == 'gzip':
                    data = StringIO(data)
                    html = gzip.GzipFile(fileobj=data).read()
                else:
                    print '%s has unknown encoding %s' % (url, encoding)
                    html = ''

                doc_pages_html.append(html)
            except urllib2.HTTPError:
                print 'error when retrieving %s' % url
                doc_pages_html.append('')
        self.doc_pages_html = doc_pages_html
        self.doc_pages_url = doc_pages_url
        self.extra_modules_test = extra_modules_test

    def get_link(self, cobj):
        comb_names = [cobj['module_short'] + '.' + cobj['name']]
        if self.extra_modules_test is not None:
            for mod in self.extra_modules_test:
                comb_names.append(mod + '.' + cobj['name'])
        for comb_name in comb_names:
            # find if the item we search for appears in any of the pages
            for html, url in zip(self.doc_pages_html, self.doc_pages_url):
                if html.find(comb_name) >= 0:
                    return url + '#' + comb_name
        return False

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


def extract_docstring(filename):
    """ Extract a module-level docstring, if any
    """
    lines = file(filename).readlines()
    start_row = 0
    if lines[0].startswith('#!'):
        lines.pop(0)
        start_row = 1

    docstring = ''
    first_par = ''
    tokens = tokenize.generate_tokens(iter(lines).next)
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
    return docstring, first_par, erow + 1 + start_row


def generate_example_rst(app):
    """ Generate the list of examples, as well as the contents of
        examples.
    """
    root_dir = os.path.join(app.builder.srcdir, 'auto_examples')
    example_dir = os.path.abspath(app.builder.srcdir + '/../' + 'examples')
    try:
        plot_gallery = eval(app.builder.config.plot_gallery)
    except TypeError:
        plot_gallery = bool(app.builder.config.plot_gallery)
    if not os.path.exists(example_dir):
        os.makedirs(example_dir)
    if not os.path.exists(root_dir):
        os.makedirs(root_dir)

    # we create an index.rst with all examples
    fhindex = file(os.path.join(root_dir, 'index.rst'), 'w')
    #Note: The sidebar button has been removed from the examples page for now
    #      due to how it messes up the layout. Will be fixed at a later point
    fhindex.write("""\

.. raw:: html


    <style type="text/css">

    div#sidebarbutton {
        display: none;
    }

    .figure {
        float: left;
        margin: 10px;
        width: auto;
        height: 200px;
        width: 180px;
    }

    .figure img {
        display: inline;
        }

    .figure .caption {
        width: 170px;
        text-align: center !important;
    }
    </style>

Examples
========

.. _examples-index:
""")
    # Here we don't use an os.walk, but we recurse only twice: flat is
    # better than nested.
    generate_dir_rst('.', fhindex, example_dir, root_dir, plot_gallery)
    for dir in sorted(os.listdir(example_dir)):
        if os.path.isdir(os.path.join(example_dir, dir)):
            generate_dir_rst(dir, fhindex, example_dir, root_dir, plot_gallery)
    fhindex.flush()


def generate_dir_rst(dir, fhindex, example_dir, root_dir, plot_gallery):
    """ Generate the rst file for an example directory.
    """
    if not dir == '.':
        target_dir = os.path.join(root_dir, dir)
        src_dir = os.path.join(example_dir, dir)
    else:
        target_dir = root_dir
        src_dir = example_dir
    if not os.path.exists(os.path.join(src_dir, 'README.txt')):
        print 80 * '_'
        print ('Example directory %s does not have a README.txt file'
                        % src_dir)
        print 'Skipping this directory'
        print 80 * '_'
        return
    fhindex.write("""


%s


""" % file(os.path.join(src_dir, 'README.txt')).read())
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    def sort_key(a):
        # put last elements without a plot
        if not a.startswith('plot') and a.endswith('.py'):
            return 'zz' + a
        return a
    for fname in sorted(os.listdir(src_dir), key=sort_key):
        if fname.endswith('py'):
            generate_file_rst(fname, target_dir, src_dir, plot_gallery)
            thumb = os.path.join(dir, 'images', 'thumb', fname[:-3] + '.png')
            link_name = os.path.join(dir, fname).replace(os.path.sep, '_')
            fhindex.write('.. figure:: %s\n' % thumb)
            if link_name.startswith('._'):
                link_name = link_name[2:]
            if dir != '.':
                fhindex.write('   :target: ./%s/%s.html\n\n' % (dir,
                                                               fname[:-3]))
            else:
                fhindex.write('   :target: ./%s.html\n\n' % link_name[:-3])
            fhindex.write("""   :ref:`example_%s`

.. toctree::
   :hidden:

   %s/%s

""" % (link_name, dir, fname[:-3]))
    fhindex.write("""
.. raw:: html

    <div style="clear: both"></div>
    """)  # clear at the end of the section

# modules for which we embed links into example code
DOCMODULES = ['sklearn', 'matplotlib', 'numpy', 'mayavi']


def make_thumbnail(in_fname, out_fname, width, height):
    """Make a thumbnail with the same aspect ratio centered in an
       image with a given width and height
    """
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
    pos_insert = ((width - width_sc) / 2, (height - height_sc) / 2)
    thumb.paste(img, pos_insert)

    thumb.save(out_fname)


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


def generate_file_rst(fname, target_dir, src_dir, plot_gallery):
    """ Generate the rst file for a given example.
    """
    base_image_name = os.path.splitext(fname)[0]
    image_fname = '%s_%%s.png' % base_image_name

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
    thumb_file = os.path.join(thumb_dir, fname[:-3] + '.png')
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

        if (not os.path.exists(first_image_file) or
                os.stat(first_image_file).st_mtime <=
                                    os.stat(src_file).st_mtime):
            # We need to execute the code
            print 'plotting %s' % fname
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

                # get variables so we can later add links to the documentation
                example_code_obj = {}
                for var_name, var in my_globals.iteritems():
                    if not hasattr(var, '__module__'):
                        continue
                    if not isinstance(var.__module__, basestring):
                        continue
                    if var.__module__.split('.')[0] not in DOCMODULES:
                        continue

                    # get the type as a string with other things stripped
                    tstr = str(type(var))
                    tstr = (tstr[tstr.find('\'')
                            + 1:tstr.rfind('\'')].split('.')[-1])
                    # get shortened module name
                    module_short = get_short_module_name(var.__module__,
                                                         tstr)
                    cobj = {'name': tstr, 'module': var.__module__,
                            'module_short': module_short,
                            'obj_type': 'object'}
                    example_code_obj[var_name] = cobj

                # find functions so we can later add links to the documentation
                funregex = re.compile('[\w.]+\(')
                fid = open(src_file, 'rt')
                for line in fid.readlines():
                    if line.startswith('#'):
                        continue
                    for match in funregex.findall(line):
                        fun_name = match[:-1]
                        try:
                            exec('this_fun = %s' % fun_name, my_globals)
                        except Exception as err:
                            print 'extracting function failed'
                            print err
                            continue
                        this_fun = my_globals['this_fun']
                        if not callable(this_fun):
                            continue
                        if not hasattr(this_fun, '__module__'):
                            continue
                        if not isinstance(this_fun.__module__, basestring):
                            continue
                        if this_fun.__module__.split('.')[0] not in DOCMODULES:
                            continue

                        # get shortened module name
                        fun_name_short = fun_name.split('.')[-1]
                        module_short = get_short_module_name(
                            this_fun.__module__, fun_name_short)
                        cobj = {'name': fun_name_short,
                                'module': this_fun.__module__,
                                'module_short': module_short,
                                'obj_type': 'function'}
                        example_code_obj[fun_name] = cobj

                fid.close()
                if len(example_code_obj) > 0:
                    # save the dictionary, so we can later add hyperlinks
                    codeobj_fname = example_file[:-3] + '_codeobj.pickle'
                    fid = open(codeobj_fname, 'wb')
                    cPickle.dump(example_code_obj, fid,
                                 cPickle.HIGHEST_PROTOCOL)

                if '__doc__' in my_globals:
                    # The __doc__ is often printed in the example, we
                    # don't with to echo it
                    my_stdout = my_stdout.replace(
                                            my_globals['__doc__'],
                                            '')
                my_stdout = my_stdout.strip()
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
                for fig_num in (fig_mngr.num for fig_mngr in
                        matplotlib._pylab_helpers.Gcf.get_all_fig_managers()):
                    # Set the fig_num figure as the current figure as we can't
                    # save a figure that's not the current figure.
                    plt.figure(fig_num)
                    plt.savefig(image_path % fig_num)
                    figure_list.append(image_fname % fig_num)
            except:
                print 80 * '_'
                print '%s is not compiling:' % fname
                traceback.print_exc()
                print 80 * '_'
            finally:
                os.chdir(cwd)
                sys.stdout = orig_stdout

            print " - time elapsed : %.2g sec" % time_elapsed
        else:
            figure_list = [f[len(image_dir):]
                            for f in glob.glob(image_path % '[1-9]')]
                            #for f in glob.glob(image_path % '*')]

        # generate thumb file
        this_template = plot_rst_template
        from matplotlib import image
        if os.path.exists(first_image_file):
            make_thumbnail(first_image_file, thumb_file, 180, 120)

    if not os.path.exists(thumb_file):
        # create something not to replace the thumbnail
        shutil.copy('images/blank_image.png', thumb_file)

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

    f = open(os.path.join(target_dir, fname[:-2] + 'rst'), 'w')
    f.write(this_template % locals())
    f.flush()


def embed_code_links(app, exception):
    """Embed hyperlinks to documentation into example code"""
    if exception is not None:
        return
    print 'Embedding documentation hyperlinks in examples..'

    # Add resolvers for the packages for which we want to show links
    doc_resolvers = {}
    doc_resolvers['sklearn'] = DocLinkResolver(
        os.path.abspath(os.path.join(app.builder.outdir, 'modules',
                                     'generated')),
        relative=True)

    doc_resolvers['numpy'] = DocLinkResolver(
        'http://docs.scipy.org/doc/numpy-1.6.0/reference/generated')

    # matplotlib and mayavi document several items on the same page
    doc_resolvers['matplotlib'] = SinglePageDocLinkResolver(
        ['http://matplotlib.org/api/animation_api.html',
         'http://matplotlib.org/api/artist_api.html',
         'http://matplotlib.org/api/axes_api.html',
         'http://matplotlib.org/api/axis_api.html',
         'http://matplotlib.org/api/index_backend_api.html',
         'http://matplotlib.org/api/cbook_api.html',
         'http://matplotlib.org/api/cm_api.html',
         'http://matplotlib.org/api/collections_api.html',
         'http://matplotlib.org/api/colorbar_api.html',
         'http://matplotlib.org/api/colors_api.html',
         'http://matplotlib.org/api/dates_api.html',
         'http://matplotlib.org/api/figure_api.html',
         'http://matplotlib.org/api/font_manager_api.html',
         'http://matplotlib.org/api/gridspec_api.html',
         'http://matplotlib.org/api/legend_api.html',
         'http://matplotlib.org/api/mathtext_api.html',
         'http://matplotlib.org/api/mlab_api.html',
         'http://matplotlib.org/api/nxutils_api.html',
         'http://matplotlib.org/api/path_api.html',
         'http://matplotlib.org/api/pyplot_api.html',
         'http://matplotlib.org/api/sankey_api.html',
         'http://matplotlib.org/api/spines_api.html',
         'http://matplotlib.org/api/ticker_api.html',
         'http://matplotlib.org/api/tight_layout_api.html',
         'http://matplotlib.org/api/units_api.html',
         'http://matplotlib.org/api/widgets_api.html',
          'http://matplotlib.org/api/pyplot_api.html'])

    mayavi_base = 'http://docs.enthought.com/mayavi/mayavi/auto/'
    doc_resolvers['mayavi'] = SinglePageDocLinkResolver(
        [mayavi_base + 'mlab_helper_functions.html',
         mayavi_base + 'mlab_figure.html',
         mayavi_base + 'mlab_decorations.html',
         mayavi_base + 'mlab_camera.html',
         mayavi_base + 'mlab_other_functions.html'],
        extra_modules_test=['mayavi.mlab'])

    example_dir = os.path.join(app.builder.srcdir, 'auto_examples')
    html_example_dir = os.path.abspath(app.builder.outdir + '/auto_examples')

    # patterns for replacement
    link_pattern = '<a href="%s">%s</a>'
    orig_pattern = '<span class="n">%s</span>'
    period = '<span class="o">.</span>'

    for dirpath, _, filenames in os.walk(html_example_dir):
        for fname in filenames:
            print '\tprocessing: %s' % fname
            full_fname = os.path.join(html_example_dir, dirpath, fname)
            subpath = dirpath[len(html_example_dir):]
            pickle_fname = (example_dir + '/' + subpath + '/' + fname[:-5]
                            + '_codeobj.pickle')
            if os.path.exists(pickle_fname):
                # we have a pickle file with the objects to embed links for
                fid = open(pickle_fname, 'rb')
                example_code_obj = cPickle.load(fid)
                fid.close()
                str_repl = {}
                # generate replacement strings with the links
                for name, cobj in example_code_obj.iteritems():
                    this_module = cobj['module'].split('.')[0]
                    if this_module not in doc_resolvers:
                        continue
                    link = doc_resolvers[this_module].resolve(cobj,
                                                              full_fname)
                    if link is not None:
                        parts = name.split('.')
                        name_html = orig_pattern % parts[0]
                        for part in parts[1:]:
                            name_html += period + orig_pattern % part
                        str_repl[name_html] = link_pattern % (link, name_html)
                # do the replacement in the html file
                if len(str_repl) > 0:
                    fid = open(full_fname, 'rt')
                    lines_in = fid.readlines()
                    fid.close()
                    fid = open(full_fname, 'wt')
                    for line in lines_in:
                        for name, link in str_repl.iteritems():
                            line = line.replace(name, link)
                        fid.write(line)
                    fid.close()
    print '[done]'


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
