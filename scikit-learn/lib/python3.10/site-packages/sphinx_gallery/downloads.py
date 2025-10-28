r"""Utilities for downloadable items."""

# Author: Óscar Nájera
# License: 3-clause BSD

import os

from .utils import zip_files

CODE_ZIP_DOWNLOAD = """
.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-gallery

    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download all examples {0} source code: {1} </{2}>`
"""

NOTEBOOK_ZIP_DOWNLOAD = """
    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download all examples in Jupyter notebooks: {0} </{1}>`
"""


def python_zip(file_list, gallery_path, extension=".py"):
    """Store all files in file_list into an zip file.

    Parameters
    ----------
    file_list : list
        Holds all the file names to be included in zip file
    gallery_path : str
        path to where the zipfile is stored
    extension : str | None
        In order to deal with downloads of plain source files and jupyter notebooks, if
        this value is not None, the file extension from files in file_list will be
        removed and replace with the value of this variable while generating the zip
        file

    Returns
    -------
    zipname : str
        zip file name, written as `target_dir_python.zip`, `target_dir_jupyter.zip`,
        or `target_dir.zip` depending on the extension
    """
    zipname = os.path.basename(os.path.normpath(gallery_path))
    if extension == ".py":
        zipname += "_python"
    elif extension == ".ipynb":
        zipname += "_jupyter"
    zipname = os.path.join(gallery_path, zipname + ".zip")
    return zip_files(file_list, zipname, gallery_path, extension)


def list_downloadable_sources(target_dir, extensions=(".py",)):
    """Return a list of source files in target_dir.

    Parameters
    ----------
    target_dir : str
        path to the directory where source file are
    extensions : tuple[str]
        tuple of file extensions to include

    Returns
    -------
    list
        list of paths to all source files in `target_dir` ending with one of the
        specified extensions
    """
    return [
        os.path.join(target_dir, fname)
        for fname in os.listdir(target_dir)
        if fname.endswith(extensions)
    ]


def generate_zipfiles(gallery_dir, src_dir, gallery_conf):
    """Collects downloadable sources and makes zipfiles of them.

    Collects all source files and Jupyter notebooks in gallery_dir.

    Parameters
    ----------
    gallery_dir : str
        path of the gallery to collect downloadable sources
    src_dir : str
        The build source directory. Needed to make the reST paths relative.
    gallery_conf : dict[str, Any]
        Sphinx-Gallery configuration dictionary

    Return
    ------
    download_rst: str
        RestructuredText to include download buttons to the generated files
    """
    src_ext = tuple(gallery_conf["example_extensions"])
    notebook_ext = tuple(gallery_conf["notebook_extensions"])
    source_files = list_downloadable_sources(gallery_dir, src_ext)
    notebook_files = list_downloadable_sources(gallery_dir, notebook_ext)
    for directory in sorted(os.listdir(gallery_dir)):
        if os.path.isdir(os.path.join(gallery_dir, directory)):
            target_dir = os.path.join(gallery_dir, directory)
            source_files.extend(list_downloadable_sources(target_dir, src_ext))
            notebook_files.extend(list_downloadable_sources(target_dir, notebook_ext))

    def rst_path(filepath):
        filepath = os.path.relpath(filepath, os.path.normpath(src_dir))
        return filepath.replace(os.sep, "/")

    all_python = all(f.endswith(".py") for f in source_files)
    py_zipfile = python_zip(source_files, gallery_dir, ".py" if all_python else None)
    dw_rst = CODE_ZIP_DOWNLOAD.format(
        "in Python" if all_python else "as",
        os.path.basename(py_zipfile),
        rst_path(py_zipfile),
    )

    if notebook_files:
        jy_zipfile = python_zip(notebook_files, gallery_dir, ".ipynb")
        dw_rst += NOTEBOOK_ZIP_DOWNLOAD.format(
            os.path.basename(jy_zipfile),
            rst_path(jy_zipfile),
        )

    return dw_rst
