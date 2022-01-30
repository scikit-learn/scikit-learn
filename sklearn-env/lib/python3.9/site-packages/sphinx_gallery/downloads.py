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

from .utils import _replace_md5

CODE_DOWNLOAD = """
.. _sphx_glr_download_{3}:

\n.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example

{2}
\n  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: {0} <{0}>`\n

\n  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: {1} <{1}>`\n"""

CODE_ZIP_DOWNLOAD = """
\n.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-gallery

\n  .. container:: sphx-glr-download sphx-glr-download-python

    :download:`Download all examples in Python source code: {0} </{1}>`\n

\n  .. container:: sphx-glr-download sphx-glr-download-jupyter

    :download:`Download all examples in Jupyter notebooks: {2} </{3}>`\n"""


def python_zip(file_list, gallery_path, extension='.py'):
    """Stores all files in file_list into an zip file

    Parameters
    ----------
    file_list : list
        Holds all the file names to be included in zip file
    gallery_path : str
        path to where the zipfile is stored
    extension : str
        '.py' or '.ipynb' In order to deal with downloads of python
        sources and jupyter notebooks the file extension from files in
        file_list will be removed and replace with the value of this
        variable while generating the zip file
    Returns
    -------
    zipname : str
        zip file name, written as `target_dir_{python,jupyter}.zip`
        depending on the extension
    """
    zipname = os.path.basename(os.path.normpath(gallery_path))
    zipname += '_python' if extension == '.py' else '_jupyter'
    zipname = os.path.join(gallery_path, zipname + '.zip')
    zipname_new = zipname + '.new'
    with zipfile.ZipFile(zipname_new, mode='w') as zipf:
        for fname in file_list:
            file_src = os.path.splitext(fname)[0] + extension
            zipf.write(file_src, os.path.relpath(file_src, gallery_path))
    _replace_md5(zipname_new)
    return zipname


def list_downloadable_sources(target_dir):
    """Returns a list of python source files is target_dir

    Parameters
    ----------
    target_dir : str
        path to the directory where python source file are
    Returns
    -------
    list
        list of paths to all Python source files in `target_dir`
    """
    return [os.path.join(target_dir, fname)
            for fname in os.listdir(target_dir)
            if fname.endswith('.py')]


def generate_zipfiles(gallery_dir, src_dir):
    """
    Collects all Python source files and Jupyter notebooks in
    gallery_dir and makes zipfiles of them

    Parameters
    ----------
    gallery_dir : str
        path of the gallery to collect downloadable sources
    src_dir : str
        The build source directory. Needed to make the RST paths relative.

    Return
    ------
    download_rst: str
        RestructuredText to include download buttons to the generated files
    """

    listdir = list_downloadable_sources(gallery_dir)
    for directory in sorted(os.listdir(gallery_dir)):
        if os.path.isdir(os.path.join(gallery_dir, directory)):
            target_dir = os.path.join(gallery_dir, directory)
            listdir.extend(list_downloadable_sources(target_dir))

    py_zipfile = python_zip(listdir, gallery_dir)
    jy_zipfile = python_zip(listdir, gallery_dir, ".ipynb")

    def rst_path(filepath):
        filepath = os.path.relpath(filepath, os.path.normpath(src_dir))
        return filepath.replace(os.sep, '/')

    dw_rst = CODE_ZIP_DOWNLOAD.format(os.path.basename(py_zipfile),
                                      rst_path(py_zipfile),
                                      os.path.basename(jy_zipfile),
                                      rst_path(jy_zipfile))
    return dw_rst
