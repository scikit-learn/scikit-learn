import logging
import os
import tempfile
import shutil
import json
from subprocess import check_call
from tarfile import TarFile

from dateutil.zoneinfo import METADATA_FN, ZONEFILENAME


def rebuild(filename, tag=None, format="gz", zonegroups=[], metadata=None):
    """Rebuild the internal timezone info in dateutil/zoneinfo/zoneinfo*tar*

    filename is the timezone tarball from ``ftp.iana.org/tz``.

    """
    tmpdir = tempfile.mkdtemp()
    zonedir = os.path.join(tmpdir, "zoneinfo")
    moduledir = os.path.dirname(__file__)
    try:
        with TarFile.open(filename) as tf:
            for name in zonegroups:
                tf.extract(name, tmpdir)
            filepaths = [os.path.join(tmpdir, n) for n in zonegroups]
            try:
                check_call(["zic", "-d", zonedir] + filepaths)
            except OSError as e:
                _print_on_nosuchfile(e)
                raise
        # write metadata file
        with open(os.path.join(zonedir, METADATA_FN), 'w') as f:
            json.dump(metadata, f, indent=4, sort_keys=True)
        target = os.path.join(moduledir, ZONEFILENAME)
        with TarFile.open(target, "w:%s" % format) as tf:
            for entry in os.listdir(zonedir):
                entrypath = os.path.join(zonedir, entry)
                tf.add(entrypath, entry)
    finally:
        shutil.rmtree(tmpdir)


def _print_on_nosuchfile(e):
    """Print helpful troubleshooting message

    e is an exception raised by subprocess.check_call()

    """
    if e.errno == 2:
        logging.error(
            "Could not find zic. Perhaps you need to install "
            "libc-bin or some other package that provides it, "
            "or it's not in your PATH?")
