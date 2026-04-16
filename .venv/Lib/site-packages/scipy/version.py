
"""
Module to expose more detailed version info for the installed `scipy`
"""
version = "1.17.1"
full_version = version
short_version = version.split('.dev')[0]
git_revision = "527eb7fd7953a1de068f94bf8b322f249b9405ae"
release = 'dev' not in version and '+' not in version

if not release:
    version = full_version
