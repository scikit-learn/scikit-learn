"""
A directive for including a matplotlib plot in a Sphinx document.

By default, in HTML output, `plot` will include a .png file with a
link to a high-res .png and .pdf.  In LaTeX output, it will include a
.pdf.

The source code for the plot may be included in one of three ways:

  1. **A path to a source file** as the argument to the directive::

       .. plot:: path/to/plot.py

     When a path to a source file is given, the content of the
     directive may optionally contain a caption for the plot::

       .. plot:: path/to/plot.py

          This is the caption for the plot

     Additionally, one may specify the name of a function to call (with
     no arguments) immediately after importing the module::

       .. plot:: path/to/plot.py plot_function1

  2. Included as **inline content** to the directive::

       .. plot::

          import matplotlib.pyplot as plt
          import matplotlib.image as mpimg
          import numpy as np
          img = mpimg.imread('_static/stinkbug.png')
          imgplot = plt.imshow(img)

  3. Using **doctest** syntax::

       .. plot::
          A plotting example:
          >>> import matplotlib.pyplot as plt
          >>> plt.plot([1,2,3], [4,5,6])

Options
-------

The ``plot`` directive supports the following options:

    format : {'python', 'doctest'}
        Specify the format of the input

    include-source : bool
        Whether to display the source code. The default can be changed
        using the `plot_include_source` variable in conf.py

    encoding : str
        If this source file is in a non-UTF8 or non-ASCII encoding,
        the encoding must be specified using the `:encoding:` option.
        The encoding will not be inferred using the ``-*- coding -*-``
        metacomment.

    context : bool or str
        If provided, the code will be run in the context of all
        previous plot directives for which the `:context:` option was
        specified.  This only applies to inline code plot directives,
        not those run from files. If the ``:context: reset`` option is
        specified, the context is reset for this and future plots, and
        previous figures are closed prior to running the code.
        ``:context:close-figs`` keeps the context but closes previous figures
        before running the code.

    nofigs : bool
        If specified, the code block will be run, but no figures will
        be inserted.  This is usually useful with the ``:context:``
        option.

Additionally, this directive supports all of the options of the
`image` directive, except for `target` (since plot will add its own
target).  These include `alt`, `height`, `width`, `scale`, `align` and
`class`.

Configuration options
---------------------

The plot directive has the following configuration options:

    plot_include_source
        Default value for the include-source option

    plot_html_show_source_link
        Whether to show a link to the source in HTML.

    plot_pre_code
        Code that should be executed before each plot. If not specified or None
        it will default to a string containing::

            import numpy as np
            from matplotlib import pyplot as plt

    plot_basedir
        Base directory, to which ``plot::`` file names are relative
        to.  (If None or empty, file names are relative to the
        directory where the file containing the directive is.)

    plot_formats
        File formats to generate. List of tuples or strings::

            [(suffix, dpi), suffix, ...]

        that determine the file format and the DPI. For entries whose
        DPI was omitted, sensible defaults are chosen. When passing from
        the command line through sphinx_build the list should be passed as
        suffix:dpi,suffix:dpi, ....

    plot_html_show_formats
        Whether to show links to the files in HTML.

    plot_rcparams
        A dictionary containing any non-standard rcParams that should
        be applied before each plot.

    plot_apply_rcparams
        By default, rcParams are applied when `context` option is not used in
        a plot directive.  This configuration option overrides this behavior
        and applies rcParams before each plot.

    plot_working_directory
        By default, the working directory will be changed to the directory of
        the example, so the code can get at its data files, if any.  Also its
        path will be added to `sys.path` so it can import any helper modules
        sitting beside it.  This configuration option can be used to specify
        a central directory (also added to `sys.path`) where data files and
        helper modules for all code are located.

    plot_template
        Provide a customized template for preparing restructured text.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange

import sys, os, shutil, io, re, textwrap
from os.path import relpath
import traceback
import warnings

if not six.PY3:
    import cStringIO

from docutils.parsers.rst import directives
from docutils.parsers.rst.directives.images import Image
align = Image.align
import sphinx

sphinx_version = sphinx.__version__.split(".")
# The split is necessary for sphinx beta versions where the string is
# '6b1'
sphinx_version = tuple([int(re.split('[^0-9]', x)[0])
                        for x in sphinx_version[:2]])

import jinja2  # Sphinx dependency.

import matplotlib
import matplotlib.cbook as cbook
try:
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("error", UserWarning)
        matplotlib.use('Agg')
except UserWarning:
    import matplotlib.pyplot as plt
    plt.switch_backend("Agg")
else:
    import matplotlib.pyplot as plt
from matplotlib import _pylab_helpers

__version__ = 2

#------------------------------------------------------------------------------
# Registration hook
#------------------------------------------------------------------------------

def plot_directive(name, arguments, options, content, lineno,
                   content_offset, block_text, state, state_machine):
    """Implementation of the ``.. plot::`` directive.

    See the module docstring for details.
    """
    return run(arguments, content, options, state_machine, state, lineno)


def _option_boolean(arg):
    if not arg or not arg.strip():
        # no argument given, assume used as a flag
        return True
    elif arg.strip().lower() in ('no', '0', 'false'):
        return False
    elif arg.strip().lower() in ('yes', '1', 'true'):
        return True
    else:
        raise ValueError('"%s" unknown boolean' % arg)


def _option_context(arg):
    if arg in [None, 'reset', 'close-figs']:
        return arg
    raise ValueError("argument should be None or 'reset' or 'close-figs'")


def _option_format(arg):
    return directives.choice(arg, ('python', 'doctest'))


def _option_align(arg):
    return directives.choice(arg, ("top", "middle", "bottom", "left", "center",
                                   "right"))


def mark_plot_labels(app, document):
    """
    To make plots referenceable, we need to move the reference from
    the "htmlonly" (or "latexonly") node to the actual figure node
    itself.
    """
    for name, explicit in six.iteritems(document.nametypes):
        if not explicit:
            continue
        labelid = document.nameids[name]
        if labelid is None:
            continue
        node = document.ids[labelid]
        if node.tagname in ('html_only', 'latex_only'):
            for n in node:
                if n.tagname == 'figure':
                    sectname = name
                    for c in n:
                        if c.tagname == 'caption':
                            sectname = c.astext()
                            break

                    node['ids'].remove(labelid)
                    node['names'].remove(name)
                    n['ids'].append(labelid)
                    n['names'].append(name)
                    document.settings.env.labels[name] = \
                        document.settings.env.docname, labelid, sectname
                    break


def setup(app):
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    options = {'alt': directives.unchanged,
               'height': directives.length_or_unitless,
               'width': directives.length_or_percentage_or_unitless,
               'scale': directives.nonnegative_int,
               'align': _option_align,
               'class': directives.class_option,
               'include-source': _option_boolean,
               'format': _option_format,
               'context': _option_context,
               'nofigs': directives.flag,
               'encoding': directives.encoding
               }

    app.add_directive('plot', plot_directive, True, (0, 2, False), **options)
    app.add_config_value('plot_pre_code', None, True)
    app.add_config_value('plot_include_source', False, True)
    app.add_config_value('plot_html_show_source_link', True, True)
    app.add_config_value('plot_formats', ['png', 'hires.png', 'pdf'], True)
    app.add_config_value('plot_basedir', None, True)
    app.add_config_value('plot_html_show_formats', True, True)
    app.add_config_value('plot_rcparams', {}, True)
    app.add_config_value('plot_apply_rcparams', False, True)
    app.add_config_value('plot_working_directory', None, True)
    app.add_config_value('plot_template', None, True)

    app.connect(str('doctree-read'), mark_plot_labels)

    metadata = {'parallel_read_safe': True, 'parallel_write_safe': True}
    return metadata

#------------------------------------------------------------------------------
# Doctest handling
#------------------------------------------------------------------------------

def contains_doctest(text):
    try:
        # check if it's valid Python as-is
        compile(text, '<string>', 'exec')
        return False
    except SyntaxError:
        pass
    r = re.compile(r'^\s*>>>', re.M)
    m = r.search(text)
    return bool(m)


def unescape_doctest(text):
    """
    Extract code from a piece of text, which contains either Python code
    or doctests.

    """
    if not contains_doctest(text):
        return text

    code = ""
    for line in text.split("\n"):
        m = re.match(r'^\s*(>>>|\.\.\.) (.*)$', line)
        if m:
            code += m.group(2) + "\n"
        elif line.strip():
            code += "# " + line.strip() + "\n"
        else:
            code += "\n"
    return code


def split_code_at_show(text):
    """
    Split code at plt.show()

    """

    parts = []
    is_doctest = contains_doctest(text)

    part = []
    for line in text.split("\n"):
        if (not is_doctest and line.strip() == 'plt.show()') or \
               (is_doctest and line.strip() == '>>> plt.show()'):
            part.append(line)
            parts.append("\n".join(part))
            part = []
        else:
            part.append(line)
    if "\n".join(part).strip():
        parts.append("\n".join(part))
    return parts


def remove_coding(text):
    r"""
    Remove the coding comment, which six.exec\_ doesn't like.
    """
    sub_re = re.compile("^#\s*-\*-\s*coding:\s*.*-\*-$", flags=re.MULTILINE)
    return sub_re.sub("", text)

#------------------------------------------------------------------------------
# Template
#------------------------------------------------------------------------------


TEMPLATE = """
{{ source_code }}

{{ only_html }}

   {% if source_link or (html_show_formats and not multi_image) %}
   (
   {%- if source_link -%}
   `Source code <{{ source_link }}>`__
   {%- endif -%}
   {%- if html_show_formats and not multi_image -%}
     {%- for img in images -%}
       {%- for fmt in img.formats -%}
         {%- if source_link or not loop.first -%}, {% endif -%}
         `{{ fmt }} <{{ dest_dir }}/{{ img.basename }}.{{ fmt }}>`__
       {%- endfor -%}
     {%- endfor -%}
   {%- endif -%}
   )
   {% endif %}

   {% for img in images %}
   .. figure:: {{ build_dir }}/{{ img.basename }}.{{ default_fmt }}
      {% for option in options -%}
      {{ option }}
      {% endfor %}

      {% if html_show_formats and multi_image -%}
        (
        {%- for fmt in img.formats -%}
        {%- if not loop.first -%}, {% endif -%}
        `{{ fmt }} <{{ dest_dir }}/{{ img.basename }}.{{ fmt }}>`__
        {%- endfor -%}
        )
      {%- endif -%}

      {{ caption }}
   {% endfor %}

{{ only_latex }}

   {% for img in images %}
   {% if 'pdf' in img.formats -%}
   .. figure:: {{ build_dir }}/{{ img.basename }}.pdf
      {% for option in options -%}
      {{ option }}
      {% endfor %}

      {{ caption }}
   {% endif -%}
   {% endfor %}

{{ only_texinfo }}

   {% for img in images %}
   .. image:: {{ build_dir }}/{{ img.basename }}.png
      {% for option in options -%}
      {{ option }}
      {% endfor %}

   {% endfor %}

"""

exception_template = """
.. htmlonly::

   [`source code <%(linkdir)s/%(basename)s.py>`__]

Exception occurred rendering plot.

"""

# the context of the plot for all directives specified with the
# :context: option
plot_context = dict()

class ImageFile(object):
    def __init__(self, basename, dirname):
        self.basename = basename
        self.dirname = dirname
        self.formats = []

    def filename(self, format):
        return os.path.join(self.dirname, "%s.%s" % (self.basename, format))

    def filenames(self):
        return [self.filename(fmt) for fmt in self.formats]


def out_of_date(original, derived):
    """
    Returns True if derivative is out-of-date wrt original,
    both of which are full file paths.
    """
    return (not os.path.exists(derived) or
            (os.path.exists(original) and
             os.stat(derived).st_mtime < os.stat(original).st_mtime))


class PlotError(RuntimeError):
    pass


def run_code(code, code_path, ns=None, function_name=None):
    """
    Import a Python module from a path, and run the function given by
    name, if function_name is not None.
    """

    # Change the working directory to the directory of the example, so
    # it can get at its data files, if any.  Add its path to sys.path
    # so it can import any helper modules sitting beside it.
    if six.PY2:
        pwd = os.getcwdu()
    else:
        pwd = os.getcwd()
    old_sys_path = list(sys.path)
    if setup.config.plot_working_directory is not None:
        try:
            os.chdir(setup.config.plot_working_directory)
        except OSError as err:
            raise OSError(str(err) + '\n`plot_working_directory` option in'
                          'Sphinx configuration file must be a valid '
                          'directory path')
        except TypeError as err:
            raise TypeError(str(err) + '\n`plot_working_directory` option in '
                            'Sphinx configuration file must be a string or '
                            'None')
        sys.path.insert(0, setup.config.plot_working_directory)
    elif code_path is not None:
        dirname = os.path.abspath(os.path.dirname(code_path))
        os.chdir(dirname)
        sys.path.insert(0, dirname)

    # Reset sys.argv
    old_sys_argv = sys.argv
    sys.argv = [code_path]

    # Redirect stdout
    stdout = sys.stdout
    if six.PY3:
        sys.stdout = io.StringIO()
    else:
        sys.stdout = cStringIO.StringIO()

    # Assign a do-nothing print function to the namespace.  There
    # doesn't seem to be any other way to provide a way to (not) print
    # that works correctly across Python 2 and 3.
    def _dummy_print(*arg, **kwarg):
        pass

    try:
        try:
            code = unescape_doctest(code)
            if ns is None:
                ns = {}
            if not ns:
                if setup.config.plot_pre_code is None:
                    six.exec_(six.text_type("import numpy as np\n" +
                    "from matplotlib import pyplot as plt\n"), ns)
                else:
                    six.exec_(six.text_type(setup.config.plot_pre_code), ns)
            ns['print'] = _dummy_print
            if "__main__" in code:
                six.exec_("__name__ = '__main__'", ns)
            code = remove_coding(code)
            six.exec_(code, ns)
            if function_name is not None:
                six.exec_(function_name + "()", ns)
        except (Exception, SystemExit) as err:
            raise PlotError(traceback.format_exc())
    finally:
        os.chdir(pwd)
        sys.argv = old_sys_argv
        sys.path[:] = old_sys_path
        sys.stdout = stdout
    return ns


def clear_state(plot_rcparams, close=True):
    if close:
        plt.close('all')
    matplotlib.rc_file_defaults()
    matplotlib.rcParams.update(plot_rcparams)


def get_plot_formats(config):
    default_dpi = {'png': 80, 'hires.png': 200, 'pdf': 200}
    formats = []
    plot_formats = config.plot_formats
    if isinstance(plot_formats, six.string_types):
        # String Sphinx < 1.3, Split on , to mimic
        # Sphinx 1.3 and later. Sphinx 1.3 always
        # returns a list.
        plot_formats = plot_formats.split(',')
    for fmt in plot_formats:
        if isinstance(fmt, six.string_types):
            if ':' in fmt:
                suffix, dpi = fmt.split(':')
                formats.append((str(suffix), int(dpi)))
            else:
                formats.append((fmt, default_dpi.get(fmt, 80)))
        elif type(fmt) in (tuple, list) and len(fmt) == 2:
            formats.append((str(fmt[0]), int(fmt[1])))
        else:
            raise PlotError('invalid image format "%r" in plot_formats' % fmt)
    return formats


def render_figures(code, code_path, output_dir, output_base, context,
                   function_name, config, context_reset=False,
                   close_figs=False):
    """
    Run a pyplot script and save the images in *output_dir*.

    Save the images under *output_dir* with file names derived from
    *output_base*
    """
    formats = get_plot_formats(config)

    # -- Try to determine if all images already exist

    code_pieces = split_code_at_show(code)

    # Look for single-figure output files first
    all_exists = True
    img = ImageFile(output_base, output_dir)
    for format, dpi in formats:
        if out_of_date(code_path, img.filename(format)):
            all_exists = False
            break
        img.formats.append(format)

    if all_exists:
        return [(code, [img])]

    # Then look for multi-figure output files
    results = []
    all_exists = True
    for i, code_piece in enumerate(code_pieces):
        images = []
        for j in xrange(1000):
            if len(code_pieces) > 1:
                img = ImageFile('%s_%02d_%02d' % (output_base, i, j),
                                output_dir)
            else:
                img = ImageFile('%s_%02d' % (output_base, j), output_dir)
            for format, dpi in formats:
                if out_of_date(code_path, img.filename(format)):
                    all_exists = False
                    break
                img.formats.append(format)

            # assume that if we have one, we have them all
            if not all_exists:
                all_exists = (j > 0)
                break
            images.append(img)
        if not all_exists:
            break
        results.append((code_piece, images))

    if all_exists:
        return results

    # We didn't find the files, so build them

    results = []
    if context:
        ns = plot_context
    else:
        ns = {}

    if context_reset:
        clear_state(config.plot_rcparams)
        plot_context.clear()

    close_figs = not context or close_figs

    for i, code_piece in enumerate(code_pieces):

        if not context or config.plot_apply_rcparams:
            clear_state(config.plot_rcparams, close_figs)
        elif close_figs:
            plt.close('all')

        run_code(code_piece, code_path, ns, function_name)

        images = []
        fig_managers = _pylab_helpers.Gcf.get_all_fig_managers()
        for j, figman in enumerate(fig_managers):
            if len(fig_managers) == 1 and len(code_pieces) == 1:
                img = ImageFile(output_base, output_dir)
            elif len(code_pieces) == 1:
                img = ImageFile("%s_%02d" % (output_base, j), output_dir)
            else:
                img = ImageFile("%s_%02d_%02d" % (output_base, i, j),
                                output_dir)
            images.append(img)
            for format, dpi in formats:
                try:
                    figman.canvas.figure.savefig(img.filename(format), dpi=dpi)
                except Exception as err:
                    raise PlotError(traceback.format_exc())
                img.formats.append(format)

        results.append((code_piece, images))

    if not context or config.plot_apply_rcparams:
        clear_state(config.plot_rcparams, close=not context)

    return results


def run(arguments, content, options, state_machine, state, lineno):
    document = state_machine.document
    config = document.settings.env.config
    nofigs = 'nofigs' in options

    formats = get_plot_formats(config)
    default_fmt = formats[0][0]

    options.setdefault('include-source', config.plot_include_source)
    keep_context = 'context' in options
    context_opt = None if not keep_context else options['context']

    rst_file = document.attributes['source']
    rst_dir = os.path.dirname(rst_file)

    if len(arguments):
        if not config.plot_basedir:
            source_file_name = os.path.join(setup.app.builder.srcdir,
                                            directives.uri(arguments[0]))
        else:
            source_file_name = os.path.join(setup.confdir, config.plot_basedir,
                                            directives.uri(arguments[0]))

        # If there is content, it will be passed as a caption.
        caption = '\n'.join(content)

        # If the optional function name is provided, use it
        if len(arguments) == 2:
            function_name = arguments[1]
        else:
            function_name = None

        with io.open(source_file_name, 'r', encoding='utf-8') as fd:
            code = fd.read()
        output_base = os.path.basename(source_file_name)
    else:
        source_file_name = rst_file
        code = textwrap.dedent("\n".join(map(six.text_type, content)))
        counter = document.attributes.get('_plot_counter', 0) + 1
        document.attributes['_plot_counter'] = counter
        base, ext = os.path.splitext(os.path.basename(source_file_name))
        output_base = '%s-%d.py' % (base, counter)
        function_name = None
        caption = ''

    base, source_ext = os.path.splitext(output_base)
    if source_ext in ('.py', '.rst', '.txt'):
        output_base = base
    else:
        source_ext = ''

    # ensure that LaTeX includegraphics doesn't choke in foo.bar.pdf filenames
    output_base = output_base.replace('.', '-')

    # is it in doctest format?
    is_doctest = contains_doctest(code)
    if 'format' in options:
        if options['format'] == 'python':
            is_doctest = False
        else:
            is_doctest = True

    # determine output directory name fragment
    source_rel_name = relpath(source_file_name, setup.confdir)
    source_rel_dir = os.path.dirname(source_rel_name)
    while source_rel_dir.startswith(os.path.sep):
        source_rel_dir = source_rel_dir[1:]

    # build_dir: where to place output files (temporarily)
    build_dir = os.path.join(os.path.dirname(setup.app.doctreedir),
                             'plot_directive',
                             source_rel_dir)
    # get rid of .. in paths, also changes pathsep
    # see note in Python docs for warning about symbolic links on Windows.
    # need to compare source and dest paths at end
    build_dir = os.path.normpath(build_dir)

    if not os.path.exists(build_dir):
        os.makedirs(build_dir)

    # output_dir: final location in the builder's directory
    dest_dir = os.path.abspath(os.path.join(setup.app.builder.outdir,
                                            source_rel_dir))
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir) # no problem here for me, but just use built-ins

    # how to link to files from the RST file
    dest_dir_link = os.path.join(relpath(setup.confdir, rst_dir),
                                 source_rel_dir).replace(os.path.sep, '/')
    try:
        build_dir_link = relpath(build_dir, rst_dir).replace(os.path.sep, '/')
    except ValueError:
        # on Windows, relpath raises ValueError when path and start are on
        # different mounts/drives
        build_dir_link = build_dir
    source_link = dest_dir_link + '/' + output_base + source_ext

    # make figures
    try:
        results = render_figures(code,
                                 source_file_name,
                                 build_dir,
                                 output_base,
                                 keep_context,
                                 function_name,
                                 config,
                                 context_reset=context_opt == 'reset',
                                 close_figs=context_opt == 'close-figs')
        errors = []
    except PlotError as err:
        reporter = state.memo.reporter
        sm = reporter.system_message(
            2, "Exception occurred in plotting {}\n from {}:\n{}".format(
                output_base, source_file_name, err),
            line=lineno)
        results = [(code, [])]
        errors = [sm]

    # Properly indent the caption
    caption = '\n'.join('      ' + line.strip()
                        for line in caption.split('\n'))

    # generate output restructuredtext
    total_lines = []
    for j, (code_piece, images) in enumerate(results):
        if options['include-source']:
            if is_doctest:
                lines = ['']
                lines += [row.rstrip() for row in code_piece.split('\n')]
            else:
                lines = ['.. code-block:: python', '']
                lines += ['    %s' % row.rstrip()
                          for row in code_piece.split('\n')]
            source_code = "\n".join(lines)
        else:
            source_code = ""

        if nofigs:
            images = []

        opts = [
            ':%s: %s' % (key, val) for key, val in six.iteritems(options)
            if key in ('alt', 'height', 'width', 'scale', 'align', 'class')]

        only_html = ".. only:: html"
        only_latex = ".. only:: latex"
        only_texinfo = ".. only:: texinfo"

        # Not-None src_link signals the need for a source link in the generated
        # html
        if j == 0 and config.plot_html_show_source_link:
            src_link = source_link
        else:
            src_link = None

        result = jinja2.Template(config.plot_template or TEMPLATE).render(
            default_fmt=default_fmt,
            dest_dir=dest_dir_link,
            build_dir=build_dir_link,
            source_link=src_link,
            multi_image=len(images) > 1,
            only_html=only_html,
            only_latex=only_latex,
            only_texinfo=only_texinfo,
            options=opts,
            images=images,
            source_code=source_code,
            html_show_formats=config.plot_html_show_formats and len(images),
            caption=caption)

        total_lines.extend(result.split("\n"))
        total_lines.extend("\n")

    if total_lines:
        state_machine.insert_input(total_lines, source=source_file_name)

    # copy image files to builder's output directory, if necessary
    if not os.path.exists(dest_dir):
        cbook.mkdirs(dest_dir)

    for code_piece, images in results:
        for img in images:
            for fn in img.filenames():
                destimg = os.path.join(dest_dir, os.path.basename(fn))
                if fn != destimg:
                    shutil.copyfile(fn, destimg)

    # copy script (if necessary)
    target_name = os.path.join(dest_dir, output_base + source_ext)
    with io.open(target_name, 'w', encoding="utf-8") as f:
        if source_file_name == rst_file:
            code_escaped = unescape_doctest(code)
        else:
            code_escaped = code
        f.write(code_escaped)

    return errors
