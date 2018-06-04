# Copyright 2015 Google Inc. All Rights Reserved.
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
"""Interface to file resources.

This module provides functions for interfacing with files: opening, writing, and
querying.
"""

import fnmatch
import os
import re

from lib2to3.pgen2 import tokenize

from yapf.yapflib import errors
from yapf.yapflib import py3compat
from yapf.yapflib import style

CR = '\r'
LF = '\n'
CRLF = '\r\n'


def GetDefaultStyleForDir(dirname):
  """Return default style name for a given directory.

  Looks for .style.yapf or setup.cfg in the parent directories.

  Arguments:
    dirname: (unicode) The name of the directory.

  Returns:
    The filename if found, otherwise return the global default (pep8).
  """
  dirname = os.path.abspath(dirname)
  while True:
    # See if we have a .style.yapf file.
    style_file = os.path.join(dirname, style.LOCAL_STYLE)
    if os.path.exists(style_file):
      return style_file

    # See if we have a setup.cfg file with a '[yapf]' section.
    config_file = os.path.join(dirname, style.SETUP_CONFIG)
    if os.path.exists(config_file):
      with open(config_file) as fd:
        config = py3compat.ConfigParser()
        config.read_file(fd)
        if config.has_section('yapf'):
          return config_file

    dirname = os.path.dirname(dirname)
    if (not dirname or not os.path.basename(dirname) or
        dirname == os.path.abspath(os.path.sep)):
      break

  global_file = os.path.expanduser(style.GLOBAL_STYLE)
  if os.path.exists(global_file):
    return global_file

  return style.DEFAULT_STYLE


def GetCommandLineFiles(command_line_file_list, recursive, exclude):
  """Return the list of files specified on the command line."""
  return _FindPythonFiles(command_line_file_list, recursive, exclude)


def WriteReformattedCode(filename,
                         reformatted_code,
                         encoding='',
                         in_place=False):
  """Emit the reformatted code.

  Write the reformatted code into the file, if in_place is True. Otherwise,
  write to stdout.

  Arguments:
    filename: (unicode) The name of the unformatted file.
    reformatted_code: (unicode) The reformatted code.
    encoding: (unicode) The encoding of the file.
    in_place: (bool) If True, then write the reformatted code to the file.
  """
  if in_place:
    with py3compat.open_with_encoding(
        filename, mode='w', encoding=encoding, newline='') as fd:
      fd.write(reformatted_code)
  else:
    py3compat.EncodeAndWriteToStdout(reformatted_code)


def LineEnding(lines):
  """Retrieve the line ending of the original source."""
  endings = {CRLF: 0, CR: 0, LF: 0}
  for line in lines:
    if line.endswith(CRLF):
      endings[CRLF] += 1
    elif line.endswith(CR):
      endings[CR] += 1
    elif line.endswith(LF):
      endings[LF] += 1
  return (sorted(endings, key=endings.get, reverse=True) or [LF])[0]


def _FindPythonFiles(filenames, recursive, exclude):
  """Find all Python files."""
  if exclude and any(e.startswith('./') for e in exclude):
    raise errors.YapfError("path in '--exclude' should not start with ./")

  python_files = []
  for filename in filenames:
    if filename != '.' and exclude and IsIgnored(filename, exclude):
      continue
    if os.path.isdir(filename):
      if recursive:
        # TODO(morbo): Look into a version of os.walk that can handle recursion.
        excluded_dirs = []
        for dirpath, _, filelist in os.walk(filename):
          if dirpath != '.' and exclude and IsIgnored(dirpath, exclude):
            excluded_dirs.append(dirpath)
            continue
          elif any(dirpath.startswith(e) for e in excluded_dirs):
            continue
          for f in filelist:
            filepath = os.path.join(dirpath, f)
            if exclude and IsIgnored(filepath, exclude):
              continue
            if IsPythonFile(filepath):
              python_files.append(filepath)
      else:
        raise errors.YapfError(
            "directory specified without '--recursive' flag: %s" % filename)
    elif os.path.isfile(filename):
      python_files.append(filename)

  return python_files


def IsIgnored(path, exclude):
  """Return True if filename matches any patterns in exclude."""
  path = path.lstrip("/")
  while path.startswith("./"):
    path = path[2:]
  return any(fnmatch.fnmatch(path, e.rstrip('/')) for e in exclude)


def IsPythonFile(filename):
  """Return True if filename is a Python file."""
  if os.path.splitext(filename)[1] == '.py':
    return True

  try:
    with open(filename, 'rb') as fd:
      encoding = tokenize.detect_encoding(fd.readline)[0]

    # Check for correctness of encoding.
    with py3compat.open_with_encoding(
        filename, mode='r', encoding=encoding) as fd:
      fd.read()
  except UnicodeDecodeError:
    encoding = 'latin-1'
  except (IOError, SyntaxError):
    # If we fail to detect encoding (or the encoding cookie is incorrect - which
    # will make detect_encoding raise SyntaxError), assume it's not a Python
    # file.
    return False

  try:
    with py3compat.open_with_encoding(
        filename, mode='r', encoding=encoding) as fd:
      first_line = fd.readlines()[0]
  except (IOError, IndexError):
    return False

  return re.match(r'^#!.*\bpython[23]?\b', first_line)


def FileEncoding(filename):
  """Return the file's encoding."""
  with open(filename, 'rb') as fd:
    return tokenize.detect_encoding(fd.readline)[0]
