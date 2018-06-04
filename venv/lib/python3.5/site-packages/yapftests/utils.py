# -*- coding: utf-8 -*-
# Copyright 2017 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Utilities for tests."""

import contextlib
import io
import os
import sys
import tempfile


@contextlib.contextmanager
def stdout_redirector(stream):  # pylint: disable=invalid-name
  old_stdout = sys.stdout
  sys.stdout = stream
  try:
    yield
  finally:
    sys.stdout = old_stdout


# NamedTemporaryFile is useless because on Windows the temporary file would be
# created with O_TEMPORARY, which would not allow the file to be opened a
# second time, even by the same process, unless the same flag is used.
# Thus we provide a simplified version ourselves.
#
# Note: returns a tuple of (io.file_obj, file_path), instead of a file_obj with
# a .name attribute
#
# Note: `buffering` is set to -1 despite documentation of NamedTemporaryFile
# says None. This is probably a problem with the python documentation.
@contextlib.contextmanager
def NamedTempFile(mode='w+b',
                  buffering=-1,
                  encoding=None,
                  errors=None,
                  newline=None,
                  suffix=None,
                  prefix=None,
                  dirname=None,
                  text=False):
  """Context manager creating a new temporary file in text mode."""
  if sys.version_info < (3, 5):  # covers also python 2
    if suffix is None:
      suffix = ''
    if prefix is None:
      prefix = 'tmp'
  (fd, fname) = tempfile.mkstemp(
      suffix=suffix, prefix=prefix, dir=dirname, text=text)
  f = io.open(
      fd,
      mode=mode,
      buffering=buffering,
      encoding=encoding,
      errors=errors,
      newline=newline)
  yield f, fname
  f.close()
  os.remove(fname)


@contextlib.contextmanager
def TempFileContents(dirname,
                     contents,
                     encoding='utf-8',
                     newline='',
                     suffix=None):
  # Note: NamedTempFile properly handles unicode encoding when using mode='w'
  with NamedTempFile(
      dirname=dirname,
      mode='w',
      encoding=encoding,
      newline=newline,
      suffix=suffix) as (f, fname):
    f.write(contents)
    f.flush()
    yield fname
