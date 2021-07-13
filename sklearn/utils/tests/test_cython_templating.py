import pathlib
import pytest
import sklearn

TEMPITA_EXTENSION = "tp"


def test_files_generated_by_templates_are_git_ignored():
    gitignore_file = pathlib.Path(sklearn.__file__).parent.parent / ".gitignore"
    if not gitignore_file.exists():
        pytest.skip("Tests are not run from the source folder")
    base_dir = pathlib.Path(sklearn.__file__).parent
    for filename in base_dir.glob(f"**/*.{TEMPITA_EXTENSION}"):
        filename_wo_tempita_suffix = filename.with_suffix("")
        assert not filename_wo_tempita_suffix.exists()
