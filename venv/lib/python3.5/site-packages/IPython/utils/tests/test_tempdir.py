#-----------------------------------------------------------------------------
#  Copyright (C) 2012-  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

import os

from IPython.utils.tempdir import NamedFileInTemporaryDirectory
from IPython.utils.tempdir import TemporaryWorkingDirectory


def test_named_file_in_temporary_directory():
    with NamedFileInTemporaryDirectory('filename') as file:
        name = file.name
        assert not file.closed
        assert os.path.exists(name)
        file.write(b'test')
    assert file.closed
    assert not os.path.exists(name)

def test_temporary_working_directory():
    with TemporaryWorkingDirectory() as dir:
        assert os.path.exists(dir)
        assert os.path.realpath(os.curdir) == os.path.realpath(dir)
    assert not os.path.exists(dir)
    assert os.path.abspath(os.curdir) != dir
