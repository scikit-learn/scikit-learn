import pathlib
import pytest
import sklearn

TEMPITA_EXTENSION = "tp"


def test_files_generated_by_templates_are_git_ignored():
    gitignore_file = pathlib.Path(sklearn.__file__).parent.parent / ".gitignore"
    if not gitignore_file.exists():
        pytest.skip("Tests are not run from the source folder")

    base_dir = pathlib.Path(sklearn.__file__).parent
    ignored_files = open(gitignore_file, "r").readlines()
    ignored_files = list(map(lambda line: line.strip("\n"), ignored_files))

    for filename in base_dir.glob(f"**/*.{TEMPITA_EXTENSION}"):
        filename = str(filename).split("scikit-learn/")[-1]
        filename_wo_tempita_suffix = filename.strip(f".{TEMPITA_EXTENSION}")
        assert filename_wo_tempita_suffix in ignored_files
