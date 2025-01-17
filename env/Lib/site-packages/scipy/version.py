
"""
Module to expose more detailed version info for the installed `scipy`
"""
version = "1.15.0"
full_version = version
short_version = version.split('.dev')[0]
git_revision = "6e246d0b54dd55dc69232a0caae6772228a7ac25"
release = 'dev' not in version and '+' not in version

if not release:
    version = full_version
