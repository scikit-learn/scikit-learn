""":mod:`sassutils.builder` --- Build the whole directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
import collections.abc
import os.path
import re
import warnings

from sass import compile

__all__ = 'SUFFIXES', 'SUFFIX_PATTERN', 'Manifest', 'build_directory'


#: (:class:`frozenset`) The set of supported filename suffixes.
SUFFIXES = frozenset(('sass', 'scss'))

#: (:class:`re.RegexObject`) The regular expression pattern which matches to
#: filenames of supported :const:`SUFFIXES`.
SUFFIX_PATTERN = re.compile(
    '[.](' + '|'.join(map(re.escape, sorted(SUFFIXES))) + ')$',
)


def build_directory(
    sass_path, css_path, output_style='nested',
    _root_sass=None, _root_css=None, strip_extension=False,
):
    """Compiles all Sass/SCSS files in ``path`` to CSS.

    :param sass_path: the path of the directory which contains source files
                      to compile
    :type sass_path: :class:`str`, :class:`basestring`
    :param css_path: the path of the directory compiled CSS files will go
    :type css_path: :class:`str`, :class:`basestring`
    :param output_style: an optional coding style of the compiled result.
                         choose one of: ``'nested'`` (default), ``'expanded'``,
                         ``'compact'``, ``'compressed'``
    :type output_style: :class:`str`
    :returns: a dictionary of source filenames to compiled CSS filenames
    :rtype: :class:`collections.abc.Mapping`

    .. versionadded:: 0.6.0
       The ``output_style`` parameter.

    """
    if _root_sass is None or _root_css is None:
        _root_sass = sass_path
        _root_css = css_path
    result = {}
    if not os.path.isdir(css_path):
        os.mkdir(css_path)
    for name in os.listdir(sass_path):
        sass_fullname = os.path.join(sass_path, name)
        if SUFFIX_PATTERN.search(name) and os.path.isfile(sass_fullname):
            if name[0] == '_':
                # Do not compile if it's partial
                continue
            if strip_extension:
                name, _ = os.path.splitext(name)
            css_fullname = os.path.join(css_path, name) + '.css'
            css = compile(
                filename=sass_fullname,
                output_style=output_style,
                include_paths=[_root_sass],
            )
            with open(
                css_fullname, 'w', encoding='utf-8', newline='',
            ) as css_file:
                css_file.write(css)
            result[os.path.relpath(sass_fullname, _root_sass)] = \
                os.path.relpath(css_fullname, _root_css)
        elif os.path.isdir(sass_fullname):
            css_fullname = os.path.join(css_path, name)
            subresult = build_directory(
                sass_fullname, css_fullname,
                output_style=output_style,
                _root_sass=_root_sass,
                _root_css=_root_css,
                strip_extension=strip_extension,
            )
            result.update(subresult)
    return result


class Manifest:
    """Building manifest of Sass/SCSS.

    :param sass_path: the path of the directory that contains Sass/SCSS
                      source files
    :type sass_path: :class:`str`, :class:`basestring`
    :param css_path: the path of the directory to store compiled CSS
                     files
    :type css_path: :class:`str`, :class:`basestring`
    :param strip_extension: whether to remove the original file extension
    :type strip_extension: :class:`bool`
    """

    @classmethod
    def normalize_manifests(cls, manifests):
        if manifests is None:
            manifests = {}
        elif isinstance(manifests, collections.abc.Mapping):
            manifests = dict(manifests)
        else:
            raise TypeError(
                'manifests must be a mapping object, not ' +
                repr(manifests),
            )
        for package_name, manifest in manifests.items():
            if not isinstance(package_name, str):
                raise TypeError(
                    'manifest keys must be a string of package '
                    'name, not ' + repr(package_name),
                )
            if isinstance(manifest, Manifest):
                continue
            elif isinstance(manifest, tuple):
                manifest = Manifest(*manifest)
            elif isinstance(manifest, collections.abc.Mapping):
                manifest = Manifest(**manifest)
            elif isinstance(manifest, str):
                manifest = Manifest(manifest)
            else:
                raise TypeError(
                    'manifest values must be a sassutils.builder.Manifest, '
                    'a pair of (sass_path, css_path), or a string of '
                    'sass_path, not ' + repr(manifest),
                )
            manifests[package_name] = manifest
        return manifests

    def __init__(
            self,
            sass_path,
            css_path=None,
            wsgi_path=None,
            strip_extension=None,
    ):
        if not isinstance(sass_path, str):
            raise TypeError(
                'sass_path must be a string, not ' +
                repr(sass_path),
            )
        if css_path is None:
            css_path = sass_path
        elif not isinstance(css_path, str):
            raise TypeError(
                'css_path must be a string, not ' +
                repr(css_path),
            )
        if wsgi_path is None:
            wsgi_path = css_path
        elif not isinstance(wsgi_path, str):
            raise TypeError(
                'wsgi_path must be a string, not ' +
                repr(wsgi_path),
            )
        if strip_extension is None:
            warnings.warn(
                '`strip_extension` was not specified, defaulting to `False`.\n'
                'In the future, `strip_extension` will default to `True`.',
                FutureWarning,
            )
            strip_extension = False
        elif not isinstance(strip_extension, bool):
            raise TypeError(
                'strip_extension must be bool not {!r}'.format(
                    strip_extension,
                ),
            )
        self.sass_path = sass_path
        self.css_path = css_path
        self.wsgi_path = wsgi_path
        self.strip_extension = strip_extension

    def resolve_filename(self, package_dir, filename):
        """Gets a proper full relative path of Sass source and
        CSS source that will be generated, according to ``package_dir``
        and ``filename``.

        :param package_dir: the path of package directory
        :type package_dir: :class:`str`, :class:`basestring`
        :param filename: the filename of Sass/SCSS source to compile
        :type filename: :class:`str`, :class:`basestring`
        :returns: a pair of (sass, css) path
        :rtype: :class:`tuple`

        """
        sass_path = os.path.join(package_dir, self.sass_path, filename)
        if self.strip_extension:
            filename, _ = os.path.splitext(filename)
        css_filename = filename + '.css'
        css_path = os.path.join(package_dir, self.css_path, css_filename)
        return sass_path, css_path

    def unresolve_filename(self, package_dir, filename):
        """Retrieves the probable source path from the output filename.  Pass
        in a .css path to get out a .scss path.

        :param package_dir: the path of the package directory
        :type package_dir: :class:`str`
        :param filename: the css filename
        :type filename: :class:`str`
        :returns: the scss filename
        :rtype: :class:`str`
        """
        filename, _ = os.path.splitext(filename)
        if self.strip_extension:
            for ext in ('.scss', '.sass'):
                test_path = os.path.join(
                    package_dir, self.sass_path, filename + ext,
                )
                if os.path.exists(test_path):
                    return filename + ext
            else:  # file not found, let it error with `.scss` extension
                return filename + '.scss'
        else:
            return filename

    def build(self, package_dir, output_style='nested'):
        """Builds the Sass/SCSS files in the specified :attr:`sass_path`.
        It finds :attr:`sass_path` and locates :attr:`css_path`
        as relative to the given ``package_dir``.

        :param package_dir: the path of package directory
        :type package_dir: :class:`str`, :class:`basestring`
        :param output_style: an optional coding style of the compiled result.
                             choose one of: ``'nested'`` (default),
                             ``'expanded'``, ``'compact'``, ``'compressed'``
        :type output_style: :class:`str`
        :returns: the set of compiled CSS filenames
        :rtype: :class:`frozenset`

        .. versionadded:: 0.6.0
           The ``output_style`` parameter.

        """
        sass_path = os.path.join(package_dir, self.sass_path)
        css_path = os.path.join(package_dir, self.css_path)
        css_files = build_directory(
            sass_path, css_path,
            output_style=output_style,
            strip_extension=self.strip_extension,
        ).values()
        return frozenset(
            os.path.join(self.css_path, filename)
            for filename in css_files
        )

    def build_one(self, package_dir, filename, source_map=False):
        """Builds one Sass/SCSS file.

        :param package_dir: the path of package directory
        :type package_dir: :class:`str`, :class:`basestring`
        :param filename: the filename of Sass/SCSS source to compile
        :type filename: :class:`str`, :class:`basestring`
        :param source_map: whether to use source maps.  if :const:`True`
                           it also write a source map to a ``filename``
                           followed by :file:`.map` suffix.
                           default is :const:`False`
        :type source_map: :class:`bool`
        :returns: the filename of compiled CSS
        :rtype: :class:`str`, :class:`basestring`

        .. versionadded:: 0.4.0
           Added optional ``source_map`` parameter.

        """
        sass_filename, css_filename = self.resolve_filename(
            package_dir, filename,
        )
        root_path = os.path.join(package_dir, self.sass_path)
        css_path = os.path.join(package_dir, self.css_path, css_filename)
        if source_map:
            source_map_path = css_filename + '.map'
            css, source_map = compile(
                filename=sass_filename,
                include_paths=[root_path],
                source_map_filename=source_map_path,  # FIXME
                output_filename_hint=css_path,
            )
        else:
            css = compile(filename=sass_filename, include_paths=[root_path])
            source_map_path = None
            source_map = None
        css_folder = os.path.dirname(css_path)
        if not os.path.exists(css_folder):
            os.makedirs(css_folder)
        with open(css_path, 'w', encoding='utf-8', newline='') as f:
            f.write(css)
        if source_map:
            # Source maps are JSON, and JSON has to be UTF-8 encoded
            with open(
                source_map_path, 'w', encoding='utf-8', newline='',
            ) as f:
                f.write(source_map)
        return css_filename
