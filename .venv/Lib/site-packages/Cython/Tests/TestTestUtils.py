import os.path
import unittest
import tempfile
import textwrap
import shutil

from ..TestUtils import write_file, write_newer_file, _parse_pattern


class TestTestUtils(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        if self.temp_dir and os.path.isdir(self.temp_dir):
            shutil.rmtree(self.temp_dir)
        super().tearDown()

    def _test_path(self, filename):
        return os.path.join(self.temp_dir, filename)

    def _test_write_file(self, content, expected, **kwargs):
        file_path = self._test_path("abcfile")
        write_file(file_path, content, **kwargs)
        assert os.path.isfile(file_path)

        with open(file_path, 'rb') as f:
            found = f.read()
        assert found == expected, (repr(expected), repr(found))

    def test_write_file_text(self):
        text = "abcüöä"
        self._test_write_file(text, text.encode('utf8'))

    def test_write_file_dedent(self):
        text = """
        A horse is a horse,
        of course, of course,
        And no one can talk to a horse
        of course
        """
        self._test_write_file(text, textwrap.dedent(text).encode('utf8'), dedent=True)

    def test_write_file_bytes(self):
        self._test_write_file(b"ab\0c", b"ab\0c")

    def test_write_newer_file(self):
        file_path_1 = self._test_path("abcfile1.txt")
        file_path_2 = self._test_path("abcfile2.txt")
        write_file(file_path_1, "abc")
        assert os.path.isfile(file_path_1)
        write_newer_file(file_path_2, file_path_1, "xyz")
        assert os.path.isfile(file_path_2)
        assert os.path.getmtime(file_path_2) > os.path.getmtime(file_path_1)

    def test_write_newer_file_same(self):
        file_path = self._test_path("abcfile.txt")
        write_file(file_path, "abc")
        mtime = os.path.getmtime(file_path)
        write_newer_file(file_path, file_path, "xyz")
        assert os.path.getmtime(file_path) > mtime

    def test_write_newer_file_fresh(self):
        file_path = self._test_path("abcfile.txt")
        assert not os.path.exists(file_path)
        write_newer_file(file_path, file_path, "xyz")
        assert os.path.isfile(file_path)

    def test_parse_pattern(self):
        self.assertEqual(
            _parse_pattern("pattern"),
            (None, None, 'pattern')
        )
        self.assertEqual(
            _parse_pattern("/start/:pattern"),
            ('start', None, 'pattern')
        )
        self.assertEqual(
            _parse_pattern(":/end/  pattern"),
            (None, 'end', 'pattern')
        )
        self.assertEqual(
            _parse_pattern("/start/:/end/  pattern"),
            ('start', 'end', 'pattern')
        )
        self.assertEqual(
            _parse_pattern("/start/:/end/pattern"),
            ('start', 'end', 'pattern')
        )
