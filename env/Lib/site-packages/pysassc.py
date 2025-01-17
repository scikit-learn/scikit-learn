#!/usr/bin/env python
r""":mod:`pysassc` --- SassC compliant command line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This provides SassC_ compliant CLI executable named :program:`pysassc`:

.. sourcecode:: console

   $ pysassc
   Usage: pysassc [options] SCSS_FILE [CSS_FILE]

There are options as well:

.. option:: -t <style>, --style <style>

   Coding style of the compiled result.  The same as :func:`sass.compile()`
   function's ``output_style`` keyword argument.  Default is ``nested``.

.. option:: -s <style>, --output-style <style>

    Alias for -t / --style.

    .. deprecated:: 0.11.0

.. option:: -I <dir>, --include-path <dir>

   Optional directory path to find ``@import``\ ed (S)CSS files.
   Can be multiply used.

.. option:: -m, -g, --sourcemap

   Emit source map.  Requires the second argument (output CSS filename).
   The filename of source map will be the output CSS filename followed by
   :file:`.map`.

   .. versionadded:: 0.4.0

.. option:: -p, --precision

   Set the precision for numbers. Default is 5.

   .. versionadded:: 0.7.0

.. option:: --source-comments

   Include debug info in output.

   .. versionadded:: 0.11.0

.. option:: --sourcemap-file

   Output file for source map

   .. versionadded:: 0.17.0

.. option:: --sourcemap-contents

   Embed sourcesContent in source map.

   .. versionadded:: 0.17.0

.. option:: --sourcemap-embed

   Embed sourceMappingUrl as data URI

   .. versionadded:: 0.17.0

.. option:: --omit-sourcemap-url

   Omit source map URL comment from output

   .. versionadded:: 0.17.0

.. option:: --sourcemap-root

   Base path, will be emitted to sourceRoot in source-map as is

   .. versionadded:: 0.17.0

.. option:: -v, --version

   Prints the program version.

.. option:: -h, --help

   Prints the help message.

.. _SassC: https://github.com/sass/sassc

"""
import functools
import optparse
import sys
import warnings

import sass


def main(argv=sys.argv, stdout=sys.stdout, stderr=sys.stderr):
    parser = optparse.OptionParser(
        usage='%prog [options] SCSS_FILE [OUT_CSS_FILE]',
        version='%prog {} (sass/libsass {})'.format(
            sass.__version__, sass.libsass_version,
        ),
    )
    output_styles = list(sass.OUTPUT_STYLES)
    output_styles = ', '.join(output_styles[:-1]) + ', or ' + output_styles[-1]
    parser.add_option(
        '-t', '--style', '-s', '--output-style', metavar='STYLE',
        type='choice', choices=list(sass.OUTPUT_STYLES), default='nested',
        help=(
            'Coding style of the compiled result.  Choose one of ' +
            output_styles + '. [default: %default]'
        ),
    )
    parser.add_option(
        '-m', '-g', '--sourcemap', dest='source_map',
        action='store_true', default=False,
        help='Emit source map.  Requires the second argument '
             '(output css filename).',
    )
    parser.add_option(
        '--sourcemap-file', dest='source_map_file', metavar='FILE',
        action='store',
        help='Output file for source map. If omitted, source map is based on '
             'the output css filename',
    )
    parser.add_option(
        '--sourcemap-contents', dest='source_map_contents',
        action='store_true', default=False,
        help='Embed sourcesContent in source map',
    )
    parser.add_option(
        '--sourcemap-embed', dest='source_map_embed',
        action='store_true', default=False,
        help='Embed sourceMappingUrl as data URI',
    )
    parser.add_option(
        '--omit-sourcemap-url', dest='omit_source_map_url',
        action='store_true', default=False,
        help='Omit source map URL comment from output',
    )
    parser.add_option(
        '--sourcemap-root', metavar='DIR',
        dest='source_map_root', action='store',
        help='Base path, will be emitted to sourceRoot in source-map as is',
    )
    parser.add_option(
        '-I', '--include-path', metavar='DIR',
        dest='include_paths', action='append',
        help='Path to find "@import"ed (S)CSS source files. '
             'Can be multiply used.',
    )
    parser.add_option(
        '-p', '--precision', action='store', type='int', default=5,
        help='Set the precision for numbers. [default: %default]',
    )
    parser.add_option(
        '--source-comments', action='store_true', default=False,
        help='Include debug info in output',
    )
    parser.add_option('--import-extensions', help=optparse.SUPPRESS_HELP)
    options, args = parser.parse_args(argv[1:])
    error = functools.partial(
        print,
        parser.get_prog_name() + ': error:',
        file=stderr,
    )
    if not args:
        parser.print_usage(stderr)
        error('too few arguments')
        return 2
    elif len(args) > 2:
        parser.print_usage(stderr)
        error('too many arguments')
        return 2
    filename = args[0]
    if options.source_map and len(args) < 2:
        parser.print_usage(stderr)
        error(
            '-m/-g/--sourcemap requires the second argument, the output '
            'css filename.',
        )
        return 2

    if options.import_extensions:
        warnings.warn(
            '`--import-extensions` has no effect and will be removed in '
            'a future version.',
            FutureWarning,
        )

    try:
        if options.source_map:
            source_map_filename = options.source_map_file or args[1] + '.map'
            css, source_map = sass.compile(
                filename=filename,
                output_style=options.style,
                source_comments=options.source_comments,
                source_map_filename=source_map_filename,
                source_map_contents=options.source_map_contents,
                source_map_embed=options.source_map_embed,
                omit_source_map_url=options.omit_source_map_url,
                source_map_root=options.source_map_root,
                output_filename_hint=args[1],
                include_paths=options.include_paths,
                precision=options.precision,
            )
        else:
            source_map_filename = None
            source_map = None
            css = sass.compile(
                filename=filename,
                output_style=options.style,
                source_comments=options.source_comments,
                include_paths=options.include_paths,
                precision=options.precision,
            )
    except OSError as e:
        error(e)
        return 3
    except sass.CompileError as e:
        error(e)
        return 1
    else:
        if len(args) < 2:
            print(css, file=stdout)
        else:
            with open(args[1], 'w', encoding='utf-8', newline='') as f:
                f.write(css)
        if source_map_filename:
            with open(
                source_map_filename, 'w', encoding='utf-8', newline='',
            ) as f:
                f.write(source_map)
    return 0


if __name__ == '__main__':
    sys.exit(main())
