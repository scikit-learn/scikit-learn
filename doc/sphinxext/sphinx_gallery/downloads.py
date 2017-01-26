# -*- coding: utf-8 -*-
r"""
Utilities for downloadable items
================================

"""
# Author: Óscar Nájera
# License: 3-clause BSD

from __future__ import absolute_import, division, print_function

import os
import zipfile

CODE_DOWNLOAD = """
\n.. container:: sphx-glr-footer

\n  .. container:: sphx-glr-download

     :download:`Download Python source code: {0} <{0}>`\n

\n  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: {1} <{1}>`\n"""

CODE_ZIP_DOWNLOAD = """
\n.. container:: sphx-glr-footer

\n  .. container:: sphx-glr-download

    :download:`Download all examples in Python source code: {0} </{1}>`\n

\n  .. container:: sphx-glr-download

    :download:`Download all examples in Jupyter notebooks: {2} </{3}>`\n"""


def python_zip(file_list, gallery_path, extension='.py'):
    """Stores all files in file_list into an zip file

    Parameters
    ----------
    file_list : list of strings
        Holds all the file names to be included in zip file
    gallery_path : string
        path to where the zipfile is stored
    extension : str
        '.py' or '.ipynb' In order to deal with downloads of python
        sources and jupyter notebooks the file extension from files in
        file_list will be removed and replace with the value of this
        variable while generating the zip file
    Returns
    -------
    zipname : string
        zip file name, written as `target_dir_{python,jupyter}.zip`
        depending on the extension
    """
    zipname = gallery_path.replace(os.path.sep, '_')
    zipname += '_python' if extension == '.py' else '_jupyter'
    zipname = os.path.join(gallery_path, zipname + '.zip')

    zipf = zipfile.ZipFile(zipname, mode='w')
    for fname in file_list:
        file_src = os.path.splitext(fname)[0] + extension
        zipf.write(file_src)
    zipf.close()

    return zipname


def list_downloadable_sources(target_dir):
    """Returns a list of python source files is target_dir

    Parameters
    ----------
    target_dir : string
        path to the directory where python source file are
    Returns
    -------
    list
        list of paths to all Python source files in `target_dir`
    """
    return [os.path.join(target_dir, fname)
            for fname in os.listdir(target_dir)
            if fname.endswith('.py')]


def generate_zipfiles(gallery_dir):
    """
    Collects all Python source files and Jupyter notebooks in
    gallery_dir and makes zipfiles of them

    Parameters
    ----------
    gallery_dir : string
        path of the gallery to collect downloadable sources

    Return
    ------
    download_rst: string
        RestructuredText to include download buttons to the generated files
    """

    listdir = list_downloadable_sources(gallery_dir)
    for directory in sorted(os.listdir(gallery_dir)):
        if os.path.isdir(os.path.join(gallery_dir, directory)):
            target_dir = os.path.join(gallery_dir, directory)
            listdir.extend(list_downloadable_sources(target_dir))

    py_zipfile = python_zip(listdir, gallery_dir)
    jy_zipfile = python_zip(listdir, gallery_dir, ".ipynb")

    dw_rst = CODE_ZIP_DOWNLOAD.format(os.path.basename(py_zipfile),
                                      py_zipfile,
                                      os.path.basename(jy_zipfile),
                                      jy_zipfile)
    return dw_rst
