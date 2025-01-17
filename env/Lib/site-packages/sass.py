""":mod:`sass` --- Binding of ``libsass``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This simple C extension module provides a very simple binding of ``libsass``,
which is written in C/C++.  It contains only one function and one exception
type.

>>> import sass
>>> sass.compile(string='a { b { color: blue; } }')
'a b {\n  color: blue; }\n'

"""
import collections.abc
import inspect
import os.path
import re
import sys
import warnings

import _sass

__all__ = (
    'MODES', 'OUTPUT_STYLES', 'SOURCE_COMMENTS', 'CompileError', 'SassColor',
    'SassError', 'SassFunction', 'SassList', 'SassMap', 'SassNumber',
    'SassWarning', 'and_join', 'compile', 'libsass_version',
)
__version__ = '0.23.0'
libsass_version = _sass.libsass_version


#: (:class:`collections.abc.Mapping`) The dictionary of output styles.
#: Keys are output name strings, and values are flag integers.
OUTPUT_STYLES = _sass.OUTPUT_STYLES

#: (:class:`collections.abc.Mapping`) The dictionary of source comments styles.
#: Keys are mode names, and values are corresponding flag integers.
#:
#: .. versionadded:: 0.4.0
#:
#: .. deprecated:: 0.6.0
SOURCE_COMMENTS = {'none': 0, 'line_numbers': 1, 'default': 1, 'map': 2}

#: (:class:`frozenset`) The set of keywords :func:`compile()` can take.
MODES = frozenset(('string', 'filename', 'dirname'))


def to_native_s(s):
    if isinstance(s, bytes):
        return s.decode('UTF-8')
    else:
        return s


class CompileError(ValueError):
    """The exception type that is raised by :func:`compile()`.
    It is a subtype of :exc:`exceptions.ValueError`.
    """

    def __init__(self, msg):
        super().__init__(to_native_s(msg))


def mkdirp(path):
    try:
        os.makedirs(path)
    except OSError:
        if os.path.isdir(path):
            return
        raise


class SassFunction:
    """Custom function for Sass.  It can be instantiated using
    :meth:`from_lambda()` and :meth:`from_named_function()` as well.

    :param name: the function name
    :type name: :class:`str`
    :param arguments: the argument names
    :type arguments: :class:`collections.abc.Sequence`
    :param callable_: the actual function to be called
    :type callable_: :class:`collections.abc.Callable`

    .. versionadded:: 0.7.0

    """

    __slots__ = 'name', 'arguments', 'callable_'

    @classmethod
    def from_lambda(cls, name, lambda_):
        """Make a :class:`SassFunction` object from the given ``lambda_``
        function.  Since lambda functions don't have their name, it need
        its ``name`` as well.  Arguments are automatically inspected.

        :param name: the function name
        :type name: :class:`str`
        :param lambda_: the actual lambda function to be called
        :type lambda_: :class:`types.LambdaType`
        :returns: a custom function wrapper of the ``lambda_`` function
        :rtype: :class:`SassFunction`

        """
        a = inspect.getfullargspec(lambda_)
        varargs, varkw, defaults, kwonlyargs = (
            a.varargs, a.varkw, a.defaults, a.kwonlyargs,
        )

        if varargs or varkw or defaults or kwonlyargs:
            raise TypeError(
                'functions cannot have starargs or defaults: {} {}'.format(
                    name, lambda_,
                ),
            )
        return cls(name, a.args, lambda_)

    @classmethod
    def from_named_function(cls, function):
        """Make a :class:`SassFunction` object from the named ``function``.
        Function name and arguments are automatically inspected.

        :param function: the named function to be called
        :type function: :class:`types.FunctionType`
        :returns: a custom function wrapper of the ``function``
        :rtype: :class:`SassFunction`

        """
        if not getattr(function, '__name__', ''):
            raise TypeError('function must be named')
        return cls.from_lambda(function.__name__, function)

    def __init__(self, name, arguments, callable_):
        if not isinstance(name, str):
            raise TypeError('name must be a string, not ' + repr(name))
        elif not isinstance(arguments, collections.abc.Sequence):
            raise TypeError(
                'arguments must be a sequence, not ' +
                repr(arguments),
            )
        elif not callable(callable_):
            raise TypeError(repr(callable_) + ' is not callable')
        self.name = name
        self.arguments = tuple(
            arg if arg.startswith('$') else '$' + arg
            for arg in arguments
        )
        self.callable_ = callable_

    @property
    def signature(self):
        """Signature string of the function."""
        return '{}({})'.format(self.name, ', '.join(self.arguments))

    def __call__(self, *args, **kwargs):
        return self.callable_(*args, **kwargs)

    def __str__(self):
        return self.signature


def _normalize_importer_return_value(result):
    # An importer must return an iterable of iterables of 1-3 stringlike
    # objects
    if result is None:
        return result

    def _to_importer_result(single_result):
        single_result = tuple(single_result)
        if len(single_result) not in (1, 2, 3):
            raise ValueError(
                'Expected importer result to be a tuple of length (1, 2, 3) '
                'but got {}: {!r}'.format(len(single_result), single_result),
            )

        def _to_bytes(obj):
            if not isinstance(obj, bytes):
                return obj.encode('UTF-8')
            else:
                return obj

        return tuple(_to_bytes(s) for s in single_result)

    return tuple(_to_importer_result(x) for x in result)


def _importer_callback_wrapper(func):
    def inner(path, prev):
        path, prev = path.decode('UTF-8'), prev.decode('UTF-8')
        num_args = getattr(inner, '_num_args', None)
        if num_args is None:
            try:
                ret = func(path, prev)
            except TypeError:
                inner._num_args = 1
                ret = func(path)
            else:
                inner._num_args = 2
        elif num_args == 2:
            ret = func(path, prev)
        else:
            ret = func(path)
        return _normalize_importer_return_value(ret)
    return inner


def _validate_importers(importers):
    """Validates the importers and decorates the callables with our output
    formatter.
    """
    # They could have no importers, that's chill
    if importers is None:
        return None

    def _to_importer(priority, func):
        assert isinstance(priority, int), priority
        assert callable(func), func
        return (priority, _importer_callback_wrapper(func))

    # Our code assumes tuple of tuples
    return tuple(_to_importer(priority, func) for priority, func in importers)


def _raise(e):
    raise e


def compile_dirname(
    search_path, output_path, output_style, source_comments, include_paths,
    precision, custom_functions, importers, source_map_contents,
    source_map_embed, omit_source_map_url, source_map_root,
):
    fs_encoding = sys.getfilesystemencoding() or sys.getdefaultencoding()
    for dirpath, _, filenames in os.walk(search_path, onerror=_raise):
        filenames = [
            filename for filename in filenames
            if filename.endswith(('.scss', '.sass')) and
            not filename.startswith('_')
        ]
        for filename in filenames:
            input_filename = os.path.join(dirpath, filename)
            relpath_to_file = os.path.relpath(input_filename, search_path)
            output_filename = os.path.join(output_path, relpath_to_file)
            output_filename = re.sub('.s[ac]ss$', '.css', output_filename)
            input_filename = input_filename.encode(fs_encoding)
            s, v, _ = _sass.compile_filename(
                input_filename, output_style, source_comments, include_paths,
                precision, None, custom_functions, importers, None,
                source_map_contents, source_map_embed, omit_source_map_url,
                source_map_root,
            )
            if s:
                v = v.decode('UTF-8')
                mkdirp(os.path.dirname(output_filename))
                with open(
                    output_filename, 'w', encoding='UTF-8', newline='',
                ) as output_file:
                    output_file.write(v)
            else:
                return False, v
    return True, None


def _check_no_remaining_kwargs(func, kwargs):
    if kwargs:
        raise TypeError(
            '{}() got unexpected keyword argument(s) {}'.format(
                func.__name__,
                ', '.join(f"'{arg}'" for arg in sorted(kwargs)),
            ),
        )


def compile(**kwargs):
    r"""There are three modes of parameters :func:`compile()` can take:
    ``string``, ``filename``, and ``dirname``.

    The ``string`` parameter is the most basic way to compile Sass.
    It simply takes a string of Sass code, and then returns a compiled
    CSS string.

    :param string: Sass source code to compile.  it's exclusive to
                   ``filename`` and ``dirname`` parameters
    :type string: :class:`str`
    :param output_style: an optional coding style of the compiled result.
                         choose one of: ``'nested'`` (default), ``'expanded'``,
                         ``'compact'``, ``'compressed'``
    :type output_style: :class:`str`
    :param source_comments: whether to add comments about source lines.
                            :const:`False` by default
    :type source_comments: :class:`bool`
    :param source_map_contents: embed include contents in map
    :type source_map_contents: :class:`bool`
    :param source_map_embed: embed sourceMappingUrl as data URI
    :type source_map_embed: :class:`bool`
    :param omit_source_map_url: omit source map URL comment from output
    :type omit_source_map_url: :class:`bool`
    :param source_map_root: base path, will be emitted in source map as is
    :type source_map_root: :class:`str`
    :param include_paths: an optional list of paths to find ``@import``\ ed
                          Sass/CSS source files
    :type include_paths: :class:`collections.abc.Sequence`
    :param precision: optional precision for numbers. :const:`5` by default.
    :type precision: :class:`int`
    :param custom_functions: optional mapping of custom functions.
                             see also below `custom functions
                             <custom-functions_>`_ description
    :type custom_functions: :class:`set`,
                            :class:`collections.abc.Sequence`,
                            :class:`collections.abc.Mapping`
    :param custom_import_extensions: (ignored, for backward compatibility)
    :param indented: optional declaration that the string is Sass, not SCSS
                     formatted. :const:`False` by default
    :type indented: :class:`bool`
    :returns: the compiled CSS string
    :param importers: optional callback functions.
                     see also below `importer callbacks
                     <importer-callbacks_>`_ description
    :type importers: :class:`collections.abc.Callable`
    :rtype: :class:`str`
    :raises sass.CompileError: when it fails for any reason
                               (for example the given Sass has broken syntax)

    The ``filename`` is the most commonly used way.  It takes a string of
    Sass filename, and then returns a compiled CSS string.

    :param filename: the filename of Sass source code to compile.
                     it's exclusive to ``string`` and ``dirname`` parameters
    :type filename: :class:`str`
    :param output_style: an optional coding style of the compiled result.
                         choose one of: ``'nested'`` (default), ``'expanded'``,
                         ``'compact'``, ``'compressed'``
    :type output_style: :class:`str`
    :param source_comments: whether to add comments about source lines.
                            :const:`False` by default
    :type source_comments: :class:`bool`
    :param source_map_filename: use source maps and indicate the source map
                                output filename.  :const:`None` means not
                                using source maps.  :const:`None` by default.
    :type source_map_filename: :class:`str`
    :param source_map_contents: embed include contents in map
    :type source_map_contents: :class:`bool`
    :param source_map_embed: embed sourceMappingUrl as data URI
    :type source_map_embed: :class:`bool`
    :param omit_source_map_url: omit source map URL comment from output
    :type omit_source_map_url: :class:`bool`
    :param source_map_root: base path, will be emitted in source map as is
    :type source_map_root: :class:`str`
    :param include_paths: an optional list of paths to find ``@import``\ ed
                          Sass/CSS source files
    :type include_paths: :class:`collections.abc.Sequence`
    :param precision: optional precision for numbers. :const:`5` by default.
    :type precision: :class:`int`
    :param custom_functions: optional mapping of custom functions.
                             see also below `custom functions
                             <custom-functions_>`_ description
    :type custom_functions: :class:`set`,
                            :class:`collections.abc.Sequence`,
                            :class:`collections.abc.Mapping`
    :param custom_import_extensions: (ignored, for backward compatibility)
    :param importers: optional callback functions.
                     see also below `importer callbacks
                     <importer-callbacks_>`_ description
    :type importers: :class:`collections.abc.Callable`
    :returns: the compiled CSS string, or a pair of the compiled CSS string
              and the source map string if ``source_map_filename`` is set
    :rtype: :class:`str`, :class:`tuple`
    :raises sass.CompileError: when it fails for any reason
                               (for example the given Sass has broken syntax)
    :raises exceptions.IOError: when the ``filename`` doesn't exist or
                                cannot be read

    The ``dirname`` is useful for automation.  It takes a pair of paths.
    The first of the ``dirname`` pair refers the source directory, contains
    several Sass source files to compiled.  Sass source files can be nested
    in directories.  The second of the pair refers the output directory
    that compiled CSS files would be saved.  Directory tree structure of
    the source directory will be maintained in the output directory as well.
    If ``dirname`` parameter is used the function returns :const:`None`.

    :param dirname: a pair of ``(source_dir, output_dir)``.
                    it's exclusive to ``string`` and ``filename``
                    parameters
    :type dirname: :class:`tuple`
    :param output_style: an optional coding style of the compiled result.
                         choose one of: ``'nested'`` (default), ``'expanded'``,
                         ``'compact'``, ``'compressed'``
    :type output_style: :class:`str`
    :param source_comments: whether to add comments about source lines.
                            :const:`False` by default
    :type source_comments: :class:`bool`
    :param source_map_contents: embed include contents in map
    :type source_map_contents: :class:`bool`
    :param source_map_embed: embed sourceMappingUrl as data URI
    :type source_map_embed: :class:`bool`
    :param omit_source_map_url: omit source map URL comment from output
    :type omit_source_map_url: :class:`bool`
    :param source_map_root: base path, will be emitted in source map as is
    :type source_map_root: :class:`str`
    :param include_paths: an optional list of paths to find ``@import``\ ed
                          Sass/CSS source files
    :type include_paths: :class:`collections.abc.Sequence`
    :param precision: optional precision for numbers. :const:`5` by default.
    :type precision: :class:`int`
    :param custom_functions: optional mapping of custom functions.
                             see also below `custom functions
                             <custom-functions_>`_ description
    :type custom_functions: :class:`set`,
                            :class:`collections.abc.Sequence`,
                            :class:`collections.abc.Mapping`
    :param custom_import_extensions: (ignored, for backward compatibility)
    :raises sass.CompileError: when it fails for any reason
                               (for example the given Sass has broken syntax)

    .. _custom-functions:

    The ``custom_functions`` parameter can take three types of forms:

    :class:`~set`/:class:`~collections.abc.Sequence` of \
    :class:`SassFunction`\ s
       It is the most general form.  Although pretty verbose, it can take
       any kind of callables like type objects, unnamed functions,
       and user-defined callables.

       .. code-block:: python

          sass.compile(
              ...,
              custom_functions={
                  sass.SassFunction('func-name', ('$a', '$b'), some_callable),
                  ...
              }
          )

    :class:`~collections.abc.Mapping` of names to functions
       Less general, but easier-to-use form.  Although it's not it can take
       any kind of callables, it can take any kind of *functions* defined
       using :keyword:`def`/:keyword:`lambda` syntax.
       It cannot take callables other than them since inspecting arguments
       is not always available for every kind of callables.

       .. code-block:: python

          sass.compile(
              ...,
              custom_functions={
                  'func-name': lambda a, b: ...,
                  ...
              }
          )

    :class:`~set`/:class:`~collections.abc.Sequence` of \
    named functions
       Not general, but the easiest-to-use form for *named* functions.
       It can take only named functions, defined using :keyword:`def`.
       It cannot take lambdas sinc names are unavailable for them.

       .. code-block:: python

          def func_name(a, b):
              return ...

          sass.compile(
              ...,
              custom_functions={func_name}
          )

    .. _importer-callbacks:

    Newer versions of ``libsass`` allow developers to define callbacks to be
    called and given a chance to process ``@import`` directives. You can
    define yours by passing in a list of callables via the ``importers``
    parameter. The callables must be passed as 2-tuples in the form:

    .. code-block:: python

        (priority_int, callback_fn)

    A priority of zero is acceptable; priority determines the order callbacks
    are attempted.

    These callbacks can accept one or two string arguments. The first argument
    is the path that was passed to the ``@import`` directive; the second
    (optional) argument is the previous resolved path, where the ``@import``
    directive was found. The callbacks must either return ``None`` to
    indicate the path wasn't handled by that callback (to continue with others
    or fall back on internal ``libsass`` filesystem behaviour) or a list of
    one or more tuples, each in one of three forms:

    * A 1-tuple representing an alternate path to handle internally; or,
    * A 2-tuple representing an alternate path and the content that path
      represents; or,
    * A 3-tuple representing the same as the 2-tuple with the addition of a
      "sourcemap".

    All tuple return values must be strings. As a not overly realistic
    example:

    .. code-block:: python

        def my_importer(path, prev):
            return [(path, '#' + path + ' { color: red; }')]

        sass.compile(
                ...,
                importers=[(0, my_importer)]
            )

    Now, within the style source, attempting to ``@import 'button';`` will
    instead attach ``color: red`` as a property of an element with the
    imported name.

    .. versionadded:: 0.4.0
       Added ``source_comments`` and ``source_map_filename`` parameters.

    .. versionchanged:: 0.6.0
       The ``source_comments`` parameter becomes to take only :class:`bool`
       instead of :class:`str`.

    .. deprecated:: 0.6.0
       Values like ``'none'``, ``'line_numbers'``, and ``'map'`` for
       the ``source_comments`` parameter are deprecated.

    .. versionadded:: 0.7.0
       Added ``precision`` parameter.

    .. versionadded:: 0.7.0
       Added ``custom_functions`` parameter.

    .. versionadded:: 0.11.0
       ``source_map_filename`` no longer implies ``source_comments``.

    .. versionadded:: 0.17.0
       Added ``source_map_contents``, ``source_map_embed``,
       ``omit_source_map_url``, and ``source_map_root`` parameters.

    .. versionadded:: 0.18.0
        The importer callbacks can now take a second argument, the previously-
        resolved path, so that importers can do relative path resolution.

    """
    modes = set()
    for mode_name in MODES:
        if mode_name in kwargs:
            modes.add(mode_name)
    if not modes:
        raise TypeError('choose one at least in ' + and_join(MODES))
    elif len(modes) > 1:
        raise TypeError(
            and_join(modes) + ' are exclusive each other; '
            'cannot be used at a time',
        )
    precision = kwargs.pop('precision', 5)
    output_style = kwargs.pop('output_style', 'nested')
    if not isinstance(output_style, str):
        raise TypeError(
            'output_style must be a string, not ' +
            repr(output_style),
        )
    try:
        output_style = OUTPUT_STYLES[output_style]
    except KeyError:
        raise CompileError(
            '{} is unsupported output_style; choose one of {}'
            ''.format(output_style, and_join(OUTPUT_STYLES)),
        )
    source_comments = kwargs.pop('source_comments', False)
    if source_comments in SOURCE_COMMENTS:
        if source_comments == 'none':
            deprecation_message = (
                'you can simply pass False to '
                "source_comments instead of 'none'"
            )
            source_comments = False
        elif source_comments in ('line_numbers', 'default'):
            deprecation_message = (
                'you can simply pass True to '
                'source_comments instead of ' +
                repr(source_comments)
            )
            source_comments = True
        else:
            deprecation_message = (
                "you don't have to pass 'map' to "
                'source_comments but just need to '
                'specify source_map_filename'
            )
            source_comments = False
        warnings.warn(
            "values like 'none', 'line_numbers', and 'map' for "
            'the source_comments parameter are deprecated; ' +
            deprecation_message,
            FutureWarning,
        )
    if not isinstance(source_comments, bool):
        raise TypeError(
            'source_comments must be bool, not ' +
            repr(source_comments),
        )
    fs_encoding = sys.getfilesystemencoding() or sys.getdefaultencoding()

    def _get_file_arg(key):
        ret = kwargs.pop(key, None)
        if ret is not None and not isinstance(ret, str):
            raise TypeError(f'{key} must be a string, not {ret!r}')
        elif isinstance(ret, str):
            ret = ret.encode(fs_encoding)
        if ret and 'filename' not in modes:
            raise CompileError(
                '{} is only available with filename= keyword argument since '
                'has to be aware of it'.format(key),
            )
        return ret

    source_map_filename = _get_file_arg('source_map_filename')
    output_filename_hint = _get_file_arg('output_filename_hint')

    source_map_contents = kwargs.pop('source_map_contents', False)
    source_map_embed = kwargs.pop('source_map_embed', False)
    omit_source_map_url = kwargs.pop('omit_source_map_url', False)
    source_map_root = kwargs.pop('source_map_root', None)

    if isinstance(source_map_root, str):
        source_map_root = source_map_root.encode('utf-8')

    # #208: cwd is always included in include paths
    include_paths = (os.getcwd(),)
    include_paths += tuple(kwargs.pop('include_paths', ()) or ())
    include_paths = os.pathsep.join(include_paths)
    if isinstance(include_paths, str):
        include_paths = include_paths.encode(fs_encoding)

    custom_functions = kwargs.pop('custom_functions', ())
    if isinstance(custom_functions, collections.abc.Mapping):
        custom_functions = [
            SassFunction.from_lambda(name, lambda_)
            for name, lambda_ in custom_functions.items()
        ]
    elif isinstance(
            custom_functions,
            (collections.abc.Set, collections.abc.Sequence),
    ):
        custom_functions = [
            func if isinstance(func, SassFunction)
            else SassFunction.from_named_function(func)
            for func in custom_functions
        ]
    else:
        raise TypeError(
            'custom_functions must be one of:\n'
            '- a set/sequence of {0.__module__}.{0.__name__} objects,\n'
            '- a mapping of function name strings to lambda functions,\n'
            '- a set/sequence of named functions,\n'
            'not {1!r}'.format(SassFunction, custom_functions),
        )

    if kwargs.pop('custom_import_extensions', None) is not None:
        warnings.warn(
            '`custom_import_extensions` has no effect and will be removed in '
            'a future version.',
            FutureWarning,
        )

    importers = _validate_importers(kwargs.pop('importers', None))

    if 'string' in modes:
        string = kwargs.pop('string')
        if isinstance(string, str):
            string = string.encode('utf-8')
        indented = kwargs.pop('indented', False)
        if not isinstance(indented, bool):
            raise TypeError(
                'indented must be bool, not ' +
                repr(source_comments),
            )
        _check_no_remaining_kwargs(compile, kwargs)
        s, v = _sass.compile_string(
            string, output_style, source_comments, include_paths, precision,
            custom_functions, indented, importers,
            source_map_contents, source_map_embed, omit_source_map_url,
            source_map_root,
        )
        if s:
            return v.decode('utf-8')
    elif 'filename' in modes:
        filename = kwargs.pop('filename')
        if not isinstance(filename, str):
            raise TypeError('filename must be a string, not ' + repr(filename))
        elif not os.path.isfile(filename):
            raise OSError(f'{filename!r} seems not a file')
        elif isinstance(filename, str):
            filename = filename.encode(fs_encoding)
        _check_no_remaining_kwargs(compile, kwargs)
        s, v, source_map = _sass.compile_filename(
            filename, output_style, source_comments, include_paths, precision,
            source_map_filename, custom_functions, importers,
            output_filename_hint,
            source_map_contents, source_map_embed, omit_source_map_url,
            source_map_root,
        )
        if s:
            v = v.decode('utf-8')
            if source_map_filename:
                source_map = source_map.decode('utf-8')
                v = v, source_map
            return v
    elif 'dirname' in modes:
        try:
            search_path, output_path = kwargs.pop('dirname')
        except ValueError:
            raise ValueError(
                'dirname must be a pair of (source_dir, '
                'output_dir)',
            )
        _check_no_remaining_kwargs(compile, kwargs)
        s, v = compile_dirname(
            search_path, output_path, output_style, source_comments,
            include_paths, precision, custom_functions, importers,
            source_map_contents, source_map_embed, omit_source_map_url,
            source_map_root,
        )
        if s:
            return
    else:
        raise TypeError('something went wrong')
    assert not s
    raise CompileError(v)


def and_join(strings):
    """Join the given ``strings`` by commas with last `' and '` conjunction.

    >>> and_join(['Korea', 'Japan', 'China', 'Taiwan'])
    'Korea, Japan, China, and Taiwan'

    :param strings: a list of words to join
    :type string: :class:`collections.abc.Sequence`
    :returns: a joined string
    :rtype: :class:`str`, :class:`basestring`

    """
    last = len(strings) - 1
    if last == 0:
        return strings[0]
    elif last < 0:
        return ''
    iterator = enumerate(strings)
    return ', '.join('and ' + s if i == last else s for i, s in iterator)


"""
This module provides datatypes to be used in custom sass functions.

The following mappings from sass types to python types are used:

SASS_NULL: ``None``
SASS_BOOLEAN: ``True`` or ``False``
SASS_STRING: class:`str`
SASS_NUMBER: class:`SassNumber`
SASS_COLOR: class:`SassColor`
SASS_LIST: class:`SassList`
SASS_MAP: class:`dict` or class:`SassMap`
SASS_ERROR: class:`SassError`
SASS_WARNING: class:`SassWarning`
"""


class SassNumber(collections.namedtuple('SassNumber', ('value', 'unit'))):

    def __new__(cls, value, unit):
        value = float(value)
        if not isinstance(unit, str):
            unit = unit.decode('UTF-8')
        return super().__new__(cls, value, unit)


class SassColor(collections.namedtuple('SassColor', ('r', 'g', 'b', 'a'))):

    def __new__(cls, r, g, b, a):
        r = float(r)
        g = float(g)
        b = float(b)
        a = float(a)
        return super().__new__(cls, r, g, b, a)


SASS_SEPARATOR_COMMA = collections.namedtuple('SASS_SEPARATOR_COMMA', ())()
SASS_SEPARATOR_SPACE = collections.namedtuple('SASS_SEPARATOR_SPACE', ())()
SEPARATORS = frozenset((SASS_SEPARATOR_COMMA, SASS_SEPARATOR_SPACE))


class SassList(
    collections.namedtuple(
        'SassList', ('items', 'separator', 'bracketed'),
    ),
):

    def __new__(cls, items, separator, bracketed=False):
        items = tuple(items)
        assert separator in SEPARATORS, separator
        assert isinstance(bracketed, bool), bracketed
        return super().__new__(cls, items, separator, bracketed)


class SassError(collections.namedtuple('SassError', ('msg',))):

    def __new__(cls, msg):
        if not isinstance(msg, str):
            msg = msg.decode('UTF-8')
        return super().__new__(cls, msg)


class SassWarning(collections.namedtuple('SassWarning', ('msg',))):

    def __new__(cls, msg):
        if not isinstance(msg, str):
            msg = msg.decode('UTF-8')
        return super().__new__(cls, msg)


class SassMap(collections.abc.Mapping):
    """Because sass maps can have mapping types as keys, we need an immutable
    hashable mapping type.

    .. versionadded:: 0.7.0

    """

    __slots__ = '_dict', '_hash'

    def __init__(self, *args, **kwargs):
        self._dict = dict(*args, **kwargs)
        # An assertion that all things are hashable
        self._hash = hash(frozenset(self._dict.items()))

    # Mapping interface

    def __getitem__(self, key):
        return self._dict[key]

    def __iter__(self):
        return iter(self._dict)

    def __len__(self):
        return len(self._dict)

    # Our interface

    def __repr__(self):
        return f'{type(self).__name__}({frozenset(self.items())})'

    def __hash__(self):
        return self._hash

    def _immutable(self, *_):
        raise TypeError('SassMaps are immutable.')

    __setitem__ = __delitem__ = _immutable
