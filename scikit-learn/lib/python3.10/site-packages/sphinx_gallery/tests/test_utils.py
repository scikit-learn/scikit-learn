"""Test utility functions."""

from pathlib import Path

from sphinx_gallery.utils import _combine_backreferences, _read_json, _write_json


def test_combine_backreferences():
    """Check `_combine_backreferences` works as expected."""
    backrefs_a = {
        "a": [1, 2],
        "b": [11, 12],
    }
    assert _combine_backreferences({}, backrefs_a) == backrefs_a

    backrefs_b = {
        "a": [3, 4],
        "c": [21, 22],
    }
    assert _combine_backreferences(backrefs_a, backrefs_b) == {
        "a": [1, 2, 3, 4],
        "b": [11, 12],
        "c": [21, 22],
    }


def test_read_write_json(tmpdir):
    """Check `_read_json` and `_write_json` work as expected."""
    path = Path(tmpdir, "test")
    data = {
        "object1": ("path/file.py", "first intro", "first title"),
        "object2": ("path2/file2.py", "second intro", "second title"),
    }
    _write_json(path, data, "test_dict")
    # Writing converts tuples to lists
    assert _read_json(Path(path).with_name(path.stem + "test_dict.json")) == {
        key: list(value) for key, value in data.items()
    }
