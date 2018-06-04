"""
Archive tools for wheel.
"""

import os
import os.path
import time
import zipfile
from distutils import log


def archive_wheelfile(base_name, base_dir):
    """Archive all files under `base_dir` in a whl file and name it like
    `base_name`.
    """
    olddir = os.path.abspath(os.curdir)
    base_name = os.path.abspath(base_name)
    try:
        os.chdir(base_dir)
        return make_wheelfile_inner(base_name)
    finally:
        os.chdir(olddir)


def make_wheelfile_inner(base_name, base_dir='.'):
    """Create a whl file from all the files under 'base_dir'.

    Places .dist-info at the end of the archive."""

    zip_filename = base_name + ".whl"

    log.info("creating '%s' and adding '%s' to it", zip_filename, base_dir)

    # Some applications need reproducible .whl files, but they can't do this
    # without forcing the timestamp of the individual ZipInfo objects.  See
    # issue #143.
    timestamp = os.environ.get('SOURCE_DATE_EPOCH')
    if timestamp is None:
        date_time = None
    else:
        date_time = time.gmtime(int(timestamp))[0:6]

    score = {'WHEEL': 1, 'METADATA': 2, 'RECORD': 3}

    def writefile(path, date_time):
        st = os.stat(path)
        if date_time is None:
            mtime = time.gmtime(st.st_mtime)
            date_time = mtime[0:6]
        zinfo = zipfile.ZipInfo(path, date_time)
        zinfo.external_attr = st.st_mode << 16
        zinfo.compress_type = zipfile.ZIP_DEFLATED
        with open(path, 'rb') as fp:
            zip.writestr(zinfo, fp.read())
        log.info("adding '%s'" % path)

    with zipfile.ZipFile(zip_filename, "w", compression=zipfile.ZIP_DEFLATED,
                         allowZip64=True) as zip:
        deferred = []
        for dirpath, dirnames, filenames in os.walk(base_dir):
            # Sort the directory names so that `os.walk` will walk them in a
            # defined order on the next iteration.
            dirnames.sort()
            for name in sorted(filenames):
                path = os.path.normpath(os.path.join(dirpath, name))

                if os.path.isfile(path):
                    if dirpath.endswith('.dist-info'):
                        deferred.append((score.get(name, 0), path))
                    else:
                        writefile(path, date_time)

        deferred.sort()
        for score, path in deferred:
            writefile(path, date_time)

    return zip_filename
