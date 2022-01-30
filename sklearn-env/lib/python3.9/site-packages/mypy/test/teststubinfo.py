import unittest

from mypy.stubinfo import is_legacy_bundled_package


class TestStubInfo(unittest.TestCase):
    def test_is_legacy_bundled_packages(self) -> None:
        assert not is_legacy_bundled_package('foobar_asdf', 2)
        assert not is_legacy_bundled_package('foobar_asdf', 3)

        assert is_legacy_bundled_package('certifi', 2)
        assert is_legacy_bundled_package('certifi', 3)

        assert is_legacy_bundled_package('scribe', 2)
        assert not is_legacy_bundled_package('scribe', 3)

        assert not is_legacy_bundled_package('dataclasses', 2)
        assert is_legacy_bundled_package('dataclasses', 3)
