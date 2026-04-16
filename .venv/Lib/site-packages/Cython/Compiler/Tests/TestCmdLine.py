import os
import sys
import re
from io import StringIO
from unittest import TestCase
from unittest.mock import patch, Mock

from .. import Options
from ..CmdLine import parse_command_line

from .Utils import backup_Options, restore_Options, check_global_options

unpatched_exists = os.path.exists

def patched_exists(path):
    # avoid the Cython command raising a file not found error
    if path in (
        'source.pyx',
        os.path.join('/work/dir', 'source.pyx'),
        os.path.join('my_working_path', 'source.pyx'),
        'file.pyx',
        'file1.pyx',
        'file2.pyx',
        'file3.pyx',
        'foo.pyx',
        'bar.pyx',
    ):
        return True
    return unpatched_exists(path)

@patch('os.path.exists', new=Mock(side_effect=patched_exists))
class CmdLineParserTest(TestCase):
    def setUp(self):
        self._options_backup = backup_Options()

    def tearDown(self):
        restore_Options(self._options_backup)

    def check_default_global_options(self, white_list=[]):
        self.assertEqual(check_global_options(self._options_backup, white_list), "")

    def check_default_options(self, options, white_list=[]):
        default_options = Options.CompilationOptions(Options.default_options)
        no_value = object()
        for name in default_options.__dict__.keys():
            if name not in white_list:
                self.assertEqual(getattr(options, name, no_value), getattr(default_options, name), msg="error in option " + name)

    def test_short_options(self):
        options, sources = parse_command_line([
            '-V', '-l', '-+', '-t', '-v', '-v', '-v', '-p', '-D', '-a', '-3',
        ])
        self.assertFalse(sources)
        self.assertTrue(options.show_version)
        self.assertTrue(options.use_listing_file)
        self.assertTrue(options.cplus)
        self.assertTrue(options.timestamps)
        self.assertTrue(options.verbose >= 3)
        self.assertTrue(Options.embed_pos_in_docstring)
        self.assertFalse(Options.docstrings)
        self.assertTrue(Options.annotate)
        self.assertEqual(options.language_level, 3)

        options, sources = parse_command_line([
            '-f', '-2', 'source.pyx',
        ])
        self.assertTrue(sources)
        self.assertTrue(len(sources) == 1)
        self.assertFalse(options.timestamps)
        self.assertEqual(options.language_level, 2)

    def test_long_options(self):
        options, sources = parse_command_line([
            '--version', '--create-listing', '--cplus', '--embed', '--timestamps',
            '--verbose', '--verbose', '--verbose',
            '--embed-positions', '--no-docstrings', '--annotate', '--lenient',
        ])
        self.assertFalse(sources)
        self.assertTrue(options.show_version)
        self.assertTrue(options.use_listing_file)
        self.assertTrue(options.cplus)
        self.assertEqual(Options.embed, 'main')
        self.assertTrue(options.timestamps)
        self.assertTrue(options.verbose >= 3)
        self.assertTrue(Options.embed_pos_in_docstring)
        self.assertFalse(Options.docstrings)
        self.assertTrue(Options.annotate)
        self.assertFalse(Options.error_on_unknown_names)
        self.assertFalse(Options.error_on_uninitialized)

        options, sources = parse_command_line([
            '--force', 'source.pyx',
        ])
        self.assertTrue(sources)
        self.assertTrue(len(sources) == 1)
        self.assertFalse(options.timestamps)

    def test_options_with_values(self):
        options, sources = parse_command_line([
            '--embed=huhu',
            '-I/test/include/dir1', '--include-dir=/test/include/dir2',
            '--include-dir', '/test/include/dir3',
            '--working=/work/dir',
            'source.pyx',
            '--output-file=/output/dir',
            '--pre-import=/pre/import',
            '--cleanup=3',
            '--annotate-coverage=cov.xml',
            '--gdb-outdir=/gdb/outdir',
            '--directive=wraparound=false',
            '--shared=foo.shared',
        ])
        self.assertEqual(sources, ['source.pyx'])
        self.assertEqual(Options.embed, 'huhu')
        self.assertEqual(options.include_path, ['/test/include/dir1', '/test/include/dir2', '/test/include/dir3'])
        self.assertEqual(options.working_path, '/work/dir')
        self.assertEqual(options.output_file, '/output/dir')
        self.assertEqual(Options.pre_import, '/pre/import')
        self.assertEqual(Options.generate_cleanup_code, 3)
        self.assertTrue(Options.annotate)
        self.assertEqual(Options.annotate_coverage_xml, 'cov.xml')
        self.assertTrue(options.gdb_debug)
        self.assertEqual(options.output_dir, '/gdb/outdir')
        self.assertEqual(options.compiler_directives['wraparound'], False)
        self.assertEqual(options.shared_utility_qualified_name, 'foo.shared')

    def test_embed_before_positional(self):
        options, sources = parse_command_line([
            '--embed',
            'source.pyx',
        ])
        self.assertEqual(sources, ['source.pyx'])
        self.assertEqual(Options.embed, 'main')

    def test_two_embeds(self):
        options, sources = parse_command_line([
            '--embed', '--embed=huhu',
            'source.pyx',
        ])
        self.assertEqual(sources, ['source.pyx'])
        self.assertEqual(Options.embed, 'huhu')

    def test_two_embeds2(self):
        options, sources = parse_command_line([
            '--embed=huhu', '--embed',
            'source.pyx',
        ])
        self.assertEqual(sources, ['source.pyx'])
        self.assertEqual(Options.embed, 'main')

    def test_no_annotate(self):
        options, sources = parse_command_line([
            '--embed=huhu', 'source.pyx'
        ])
        self.assertFalse(Options.annotate)

    def test_annotate_short(self):
        options, sources = parse_command_line([
            '-a',
            'source.pyx',
        ])
        self.assertEqual(Options.annotate, 'default')

    def test_annotate_long(self):
        options, sources = parse_command_line([
            '--annotate',
            'source.pyx',
        ])
        self.assertEqual(Options.annotate, 'default')

    def test_annotate_fullc(self):
        options, sources = parse_command_line([
            '--annotate-fullc',
            'source.pyx',
        ])
        self.assertEqual(Options.annotate, 'fullc')

    def test_short_w(self):
        options, sources = parse_command_line([
            '-w', 'my_working_path',
            'source.pyx'
        ])
        self.assertEqual(options.working_path, 'my_working_path')
        self.check_default_global_options()
        self.check_default_options(options, ['working_path'])

    def test_short_o(self):
        options, sources = parse_command_line([
            '-o', 'my_output',
            'source.pyx'
        ])
        self.assertEqual(options.output_file, 'my_output')
        self.check_default_global_options()
        self.check_default_options(options, ['output_file'])

    def test_short_z(self):
        options, sources = parse_command_line([
            '-z', 'my_preimport',
            'source.pyx'
        ])
        self.assertEqual(Options.pre_import, 'my_preimport')
        self.check_default_global_options(['pre_import'])
        self.check_default_options(options)

    def test_convert_range(self):
        options, sources = parse_command_line([
            '--convert-range',
            'source.pyx'
        ])
        self.assertEqual(Options.convert_range, True)
        self.check_default_global_options(['convert_range'])
        self.check_default_options(options)

    def test_line_directives(self):
        options, sources = parse_command_line([
            '--line-directives',
            'source.pyx'
        ])
        self.assertEqual(options.emit_linenums, True)
        self.check_default_global_options()
        self.check_default_options(options, ['emit_linenums'])

    def test_no_c_in_traceback(self):
        options, sources = parse_command_line([
            '--no-c-in-traceback',
            'source.pyx'
        ])
        self.assertEqual(options.c_line_in_traceback, False)
        self.check_default_global_options()
        self.check_default_options(options, ['c_line_in_traceback'])

    def test_gdb(self):
        options, sources = parse_command_line([
            '--gdb',
            'source.pyx'
        ])
        self.assertEqual(options.gdb_debug, True)
        self.assertEqual(options.output_dir, os.curdir)
        self.check_default_global_options()
        self.check_default_options(options, ['gdb_debug', 'output_dir'])

    def test_3str(self):
        options, sources = parse_command_line([
            '--3str',
            'source.pyx'
        ])
        self.assertEqual(options.language_level, '3')
        self.check_default_global_options()
        self.check_default_options(options, ['language_level'])

    def test_capi_reexport_cincludes(self):
        options, sources = parse_command_line([
            '--capi-reexport-cincludes',
            'source.pyx'
        ])
        self.assertEqual(options.capi_reexport_cincludes, True)
        self.check_default_global_options()
        self.check_default_options(options, ['capi_reexport_cincludes'])

    def test_fast_fail(self):
        options, sources = parse_command_line([
            '--fast-fail',
            'source.pyx'
        ])
        self.assertEqual(Options.fast_fail, True)
        self.check_default_global_options(['fast_fail'])
        self.check_default_options(options)

    def test_cimport_from_pyx(self):
        options, sources = parse_command_line([
            '--cimport-from-pyx',
            'source.pyx'
        ])
        self.assertEqual(Options.cimport_from_pyx, True)
        self.check_default_global_options(['cimport_from_pyx'])
        self.check_default_options(options)

    def test_Werror(self):
        options, sources = parse_command_line([
            '-Werror',
            'source.pyx'
        ])
        self.assertEqual(Options.warning_errors, True)
        self.check_default_global_options(['warning_errors'])
        self.check_default_options(options)

    def test_warning_errors(self):
        options, sources = parse_command_line([
            '--warning-errors',
            'source.pyx'
        ])
        self.assertEqual(Options.warning_errors, True)
        self.check_default_global_options(['warning_errors'])
        self.check_default_options(options)

    def test_Wextra(self):
        options, sources = parse_command_line([
            '-Wextra',
            'source.pyx'
        ])
        self.assertEqual(options.compiler_directives, Options.extra_warnings)
        self.check_default_global_options()
        self.check_default_options(options, ['compiler_directives'])

    def test_warning_extra(self):
        options, sources = parse_command_line([
            '--warning-extra',
            'source.pyx'
        ])
        self.assertEqual(options.compiler_directives, Options.extra_warnings)
        self.check_default_global_options()
        self.check_default_options(options, ['compiler_directives'])

    def test_old_style_globals(self):
        options, sources = parse_command_line([
            '--old-style-globals',
            'source.pyx'
        ])
        self.assertEqual(Options.old_style_globals, True)
        self.check_default_global_options(['old_style_globals'])
        self.check_default_options(options)

    def test_directive_multiple(self):
        options, source = parse_command_line([
               '-X', 'cdivision=True',
               '-X', 'c_string_type=bytes',
               'source.pyx'
        ])
        self.assertEqual(options.compiler_directives['cdivision'], True)
        self.assertEqual(options.compiler_directives['c_string_type'], 'bytes')
        self.check_default_global_options()
        self.check_default_options(options, ['compiler_directives'])

    def test_directive_multiple_v2(self):
        options, source = parse_command_line([
               '-X', 'cdivision=True,c_string_type=bytes',
               'source.pyx'
        ])
        self.assertEqual(options.compiler_directives['cdivision'], True)
        self.assertEqual(options.compiler_directives['c_string_type'], 'bytes')
        self.check_default_global_options()
        self.check_default_options(options, ['compiler_directives'])

    def test_directive_value_yes(self):
        options, source = parse_command_line([
               '-X', 'cdivision=YeS',
               'source.pyx'
        ])
        self.assertEqual(options.compiler_directives['cdivision'], True)
        self.check_default_global_options()
        self.check_default_options(options, ['compiler_directives'])

    def test_directive_value_no(self):
        options, source = parse_command_line([
               '-X', 'cdivision=no',
               'source.pyx'
        ])
        self.assertEqual(options.compiler_directives['cdivision'], False)
        self.check_default_global_options()
        self.check_default_options(options, ['compiler_directives'])

    def test_directive_value_invalid(self):
        self.assertRaises(ValueError, parse_command_line, [
               '-X', 'cdivision=sadfasd',
               'source.pyx'
        ])

    def test_directive_key_invalid(self):
        self.assertRaises(ValueError, parse_command_line, [
               '-X', 'abracadabra',
               'source.pyx'
        ])

    def test_directive_no_value(self):
        self.assertRaises(ValueError, parse_command_line, [
               '-X', 'cdivision',
               'source.pyx'
        ])

    def test_compile_time_env_short(self):
        options, source = parse_command_line([
               '-E', 'MYSIZE=10',
               'source.pyx'
        ])
        self.assertEqual(options.compile_time_env['MYSIZE'], 10)
        self.check_default_global_options()
        self.check_default_options(options, ['compile_time_env'])

    def test_compile_time_env_long(self):
        options, source = parse_command_line([
               '--compile-time-env', 'MYSIZE=10',
               'source.pyx'
        ])
        self.assertEqual(options.compile_time_env['MYSIZE'], 10)
        self.check_default_global_options()
        self.check_default_options(options, ['compile_time_env'])

    def test_compile_time_env_multiple(self):
        options, source = parse_command_line([
               '-E', 'MYSIZE=10', '-E', 'ARRSIZE=11',
               'source.pyx'
        ])
        self.assertEqual(options.compile_time_env['MYSIZE'], 10)
        self.assertEqual(options.compile_time_env['ARRSIZE'], 11)
        self.check_default_global_options()
        self.check_default_options(options, ['compile_time_env'])

    def test_compile_time_env_multiple_v2(self):
        options, source = parse_command_line([
               '-E', 'MYSIZE=10,ARRSIZE=11',
               'source.pyx'
        ])
        self.assertEqual(options.compile_time_env['MYSIZE'], 10)
        self.assertEqual(options.compile_time_env['ARRSIZE'], 11)
        self.check_default_global_options()
        self.check_default_options(options, ['compile_time_env'])

    def test_option_first(self):
        options, sources = parse_command_line(['-V', 'file.pyx'])
        self.assertEqual(sources, ['file.pyx'])

    def test_file_inbetween(self):
        options, sources = parse_command_line(['-V', 'file.pyx', '-a'])
        self.assertEqual(sources, ['file.pyx'])

    def test_option_trailing(self):
        options, sources = parse_command_line(['file.pyx', '-V'])
        self.assertEqual(sources, ['file.pyx'])

    def test_multiple_files(self):
        options, sources = parse_command_line([
             'file1.pyx', '-V',
             'file2.pyx', '-a',
             'file3.pyx'
        ])
        self.assertEqual(sources, ['file1.pyx', 'file2.pyx', 'file3.pyx'])

    def test_debug_flags(self):
        options, sources = parse_command_line([
             '--debug-disposal-code', '--debug-coercion',
             'file3.pyx'
        ])
        from Cython.Compiler import DebugFlags
        for name in ['debug_disposal_code', 'debug_temp_alloc', 'debug_coercion']:
            self.assertEqual(getattr(DebugFlags, name), name in ['debug_disposal_code', 'debug_coercion'])
            setattr(DebugFlags, name, 0)  # restore original value

    def test_gdb_overwrites_gdb_outdir(self):
        options, sources = parse_command_line([
            '--gdb-outdir=my_dir', '--gdb',
            'file3.pyx'
        ])
        self.assertEqual(options.gdb_debug, True)
        self.assertEqual(options.output_dir, os.curdir)
        self.check_default_global_options()
        self.check_default_options(options, ['gdb_debug', 'output_dir'])

    def test_gdb_first(self):
        options, sources = parse_command_line([
            '--gdb', '--gdb-outdir=my_dir',
            'file3.pyx'
        ])
        self.assertEqual(options.gdb_debug, True)
        self.assertEqual(options.output_dir, 'my_dir')
        self.check_default_global_options()
        self.check_default_options(options, ['gdb_debug', 'output_dir'])

    def test_coverage_overwrites_annotation(self):
        options, sources = parse_command_line([
            '--annotate-fullc', '--annotate-coverage=my.xml',
            'file3.pyx'
        ])
        self.assertEqual(Options.annotate, True)
        self.assertEqual(Options.annotate_coverage_xml, 'my.xml')
        self.check_default_global_options(['annotate', 'annotate_coverage_xml'])
        self.check_default_options(options)

    def test_coverage_first(self):
        options, sources = parse_command_line([
            '--annotate-coverage=my.xml', '--annotate-fullc',
            'file3.pyx'
        ])
        self.assertEqual(Options.annotate, 'fullc')
        self.assertEqual(Options.annotate_coverage_xml, 'my.xml')
        self.check_default_global_options(['annotate', 'annotate_coverage_xml'])
        self.check_default_options(options)

    def test_annotate_first_fullc_second(self):
        options, sources = parse_command_line([
            '--annotate', '--annotate-fullc',
            'file3.pyx'
        ])
        self.assertEqual(Options.annotate, 'fullc')
        self.check_default_global_options(['annotate'])
        self.check_default_options(options)

    def test_annotate_fullc_first(self):
        options, sources = parse_command_line([
            '--annotate-fullc', '--annotate',
            'file3.pyx'
        ])
        self.assertEqual(Options.annotate, 'default')
        self.check_default_global_options(['annotate'])
        self.check_default_options(options)

    def test_warning_extra_dont_overwrite(self):
        options, sources = parse_command_line([
            '-X', 'cdivision=True',
            '--warning-extra',
            '-X', 'c_string_type=bytes',
            'source.pyx'
        ])
        self.assertTrue(len(options.compiler_directives), len(Options.extra_warnings) + 1)
        self.check_default_global_options()
        self.check_default_options(options, ['compiler_directives'])

    def test_module_name(self):
        options, sources = parse_command_line([
            'source.pyx'
        ])
        self.assertEqual(options.module_name, None)
        self.check_default_global_options()
        self.check_default_options(options)
        options, sources = parse_command_line([
            '--module-name', 'foo.bar',
            'source.pyx'
        ])
        self.assertEqual(options.module_name, 'foo.bar')
        self.check_default_global_options()
        self.check_default_options(options, ['module_name'])

    def test_generate_shared(self):
        options, sources = parse_command_line([
            '--generate-shared=foo/shared.c',
        ])
        self.assertEqual(sources, [])
        self.assertEqual(options.shared_c_file_path, 'foo/shared.c')

    def test_errors(self):
        def error(args, regex=None):
            old_stderr = sys.stderr
            stderr = sys.stderr = StringIO()
            try:
                self.assertRaises(SystemExit, parse_command_line, list(args))
            finally:
                sys.stderr = old_stderr
            msg = stderr.getvalue()
            err_msg = 'Message "{}"'.format(msg.strip())
            self.assertTrue(msg.startswith('usage: '),
                            '%s does not start with "usage :"' % err_msg)
            self.assertTrue(': error: ' in msg,
                            '%s does not contain ": error :"' % err_msg)
            if regex:
                self.assertTrue(re.search(regex, msg),
                                '%s does not match search "%s"' %
                                (err_msg, regex))

        error(['-1'],
              'unknown option -1')
        error(['-I'],
              'argument -I/--include-dir: expected one argument')
        error(['--version=-a'],
              "argument -V/--version: ignored explicit argument '-a'")
        error(['--version=--annotate=true'],
              "argument -V/--version: ignored explicit argument "
              "'--annotate=true'")
        error(['--working'],
              "argument -w/--working: expected one argument")
        error(['--verbose=1'],
              "argument -v/--verbose: ignored explicit argument '1'")
        error(['--cleanup'],
              "argument --cleanup: expected one argument")
        error(['--debug-disposal-code-wrong-name', 'file3.pyx'],
              "unknown option --debug-disposal-code-wrong-name")
        error(['--module-name', 'foo.pyx'],
              "Need at least one source file")
        error(['--module-name', 'foo.bar'],
              "Need at least one source file")
        error(['--module-name', 'foo.bar', 'foo.pyx', 'bar.pyx'],
              "Only one source file allowed when using --module-name")
        error(['--module-name', 'foo.bar', '--timestamps', 'foo.pyx'],
              "Cannot use --module-name with --timestamps")
        error(['--generate-shared=shared.c', 'foo.pyx'],
              "Source file not allowed when using --generate-shared")
        error(['--generate-shared'],
              "argument --generate-shared: expected one argument")
