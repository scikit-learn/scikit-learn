import pathlib
import re
import unittest

from ...Utils import open_source_file
from ..Dependencies import strip_string_literals


class TestStripLiterals(unittest.TestCase):
    maxDiff = None

    @staticmethod
    def _rebuild_string(stripped, literals):
        def lookup(match):
            return literals[match.group()]

        return re.sub("__Pyx_L[0-9]+_", lookup, stripped)

    def test_strip_string_literals(self):
        def strip_equals(s, expected):
            stripped, literals = strip_string_literals(s)
            self.assertEqual(expected, stripped)

            recovered = self._rebuild_string(stripped, literals)
            self.assertEqual(s, recovered)

        unchanged = [
            "",
            """abc""",
            """123""",
            """func(123)""",
            """ '' """,
            """ '''''''''''' """,
            """ '''''''''''''' """,
        ]

        tests = [(code, code) for code in unchanged] + [
            # strings and quotes
            ('"x"',
             '"__Pyx_L1_"'),
            ("'x'",
             "'__Pyx_L1_'"),
            (""" '"' "'" """,
             """ '__Pyx_L1_' "__Pyx_L2_" """),
            (""" '''' ''' """,
             """ '''__Pyx_L1_''' """),
            (''' """" """ ''',
             ''' """__Pyx_L1_""" '''),
            (" '''a\n''' ",
             " '''__Pyx_L1_''' "),

            # escapes
            (r"'a\'b'",
             "'__Pyx_L1_'"),
            (r"'a\\'",
             "'__Pyx_L1_'"),
            (r"'a\\\'b'",
             "'__Pyx_L1_'"),

            # string prefixes
            ("u'abc'",
             "u'__Pyx_L1_'"),
            (r"r'abc\\'",
             "r'__Pyx_L1_'"),
            (r"ru'abc\\'",
             "ru'__Pyx_L1_'"),

            # comments
            ("abc # foo",
             "abc #__Pyx_L1_"),
            ("abc # 'x'",
             "abc #__Pyx_L1_"),
            ("'abc#'",
             "'__Pyx_L1_'"),

            # special commands
            ("include 'a.pxi' # something here",
             "include '__Pyx_L1_' #__Pyx_L2_"),
            ("cdef extern from 'a.h': # comment",
             "cdef extern from '__Pyx_L1_': #__Pyx_L2_"),

            # mixed strings
            (""" func('xyz') + " " + "" '' # '' | "" "123" 'xyz' "' """,
             """ func('__Pyx_L1_') + "__Pyx_L2_" + "" '' #__Pyx_L3_"""),

            (""" f'f' """,
             """ f'__Pyx_L1_' """),

            (""" f'a{123}b' """,
             """ f'__Pyx_L1_{123}__Pyx_L2_' """),

            (""" f'{1}{f'xyz'}' """,
             """ f'{1}{f'__Pyx_L1_'}' """),

            (""" f'{f'''xyz{f\"""abc\"""}'''}' """,
             """ f'{f'''__Pyx_L1_{f\"""__Pyx_L2_\"""}'''}' """),

            (""" f'{{{{{"abc"}}}}}{{}}{{' == '{{abc}}{}{' """,
             """ f'__Pyx_L1_{"__Pyx_L2_"}__Pyx_L3_' == '__Pyx_L4_' """),

            ("f'" + ('{x} ' * 250) + "{x:{width}} '",
             "f'" + ''.join([f'{{x}}__Pyx_L{n}_' for n in range(1, 251)]) + "{x:{width}}__Pyx_L251_'")
        ]

        for code, expected in tests:
            with self.subTest(code=code):
                strip_equals(code, expected)  # plain
            code = code.strip()
            expected = expected.strip()
            with self.subTest(code=code):
                strip_equals(code, expected)  # stripped
            code += "\n"
            expected += "\n"
            with self.subTest(code=code):
                strip_equals(code, expected)  # +EOL

        # GH-5977: unclosed string literal
        strip_equals(
            """ print("Say something: %s' % something) """,
            """ print("__Pyx_L1_"""
        )

    def _test_all_files(self, base_dir, file_paths):
        _find_leftover_string = re.compile(r"""[^_'"}](['"]+)[^_'"{]""").search
        for file_path in sorted(file_paths):
            with self.subTest(file=str(file_path.relative_to(base_dir))):
                with open_source_file(str(file_path)) as f:
                    code = f.read()
                stripped, literals = strip_string_literals(code)

                match = _find_leftover_string(stripped)
                if match and len(match.group(1)) != 2:
                    match_pos = match.start() + 1
                    self.fail(f"Leftover string found: {stripped[match_pos - 12 : match_pos + 12]!r}")

                recovered = self._rebuild_string(stripped, literals)
                self.assertEqual(code, recovered)

    def test_strip_string_literals_py_files(self):
        # process all .py files in the Cython package
        package_dir = pathlib.Path(__file__).absolute().parents[2]
        assert package_dir.name == 'Cython'
        base_dir = package_dir.parent
        self._test_all_files(base_dir, package_dir.rglob("*.py"))

    def test_strip_string_literals_test_files(self):
        # process all .py[x] files in the tests package
        base_dir = pathlib.Path(__file__).absolute().parents[3]
        tests_dir = base_dir / 'tests'
        test_files = []
        for test_subdir in tests_dir.iterdir():
            if test_subdir.is_dir() and test_subdir.name != 'errors':
                test_files.extend(test_subdir.rglob("*.py"))
                test_files.extend(test_subdir.rglob("*.pyx"))
        self._test_all_files(base_dir, test_files)
