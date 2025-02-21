
"""
Module to expose more detailed version info for the installed `numpy`
"""
version = "2.2.2"
__version__ = version
full_version = version

git_revision = "fd8a68e4978defd094cffa23c71d7de7fb213e3f"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
