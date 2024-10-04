
"""
Module to expose more detailed version info for the installed `numpy`
"""
version = "2.1.1"
__version__ = version
full_version = version

git_revision = "48606ab22bfdb0e9d7ec4ed5eef984b873b7796d"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
