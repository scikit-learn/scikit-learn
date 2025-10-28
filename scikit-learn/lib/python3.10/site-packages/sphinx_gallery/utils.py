"""Utilities.

Miscellaneous utilities.
"""

# Author: Eric Larson
# License: 3-clause BSD

import hashlib
import json
import os
import re
import subprocess
import zipfile
from collections import defaultdict
from functools import partial
from pathlib import Path
from shutil import copyfile, move

import sphinx.util
from sphinx.errors import ExtensionError

try:
    from sphinx.util.display import status_iterator  # noqa: F401
except Exception:  # Sphinx < 6
    from sphinx.util import status_iterator  # noqa: F401


logger = sphinx.util.logging.getLogger("sphinx-gallery")


# Text writing kwargs for builtins.open
_W_KW = dict(encoding="utf-8", newline="\n")


def _get_image():
    try:
        from PIL import Image
    except ImportError as exc:  # capture the error for the modern way
        try:
            import Image
        except ImportError:
            raise ExtensionError(
                "Could not import pillow, which is required "
                f"to rescale images (e.g., for thumbnails): {exc}"
            )
    return Image


def scale_image(in_fname, out_fname, max_width, max_height):
    """Scales image centered in image box using `max_width` and `max_height`.

    The same aspect ratio is retained. If `in_fname` == `out_fname` the image can only
    be scaled down.
    """
    # local import to avoid testing dependency on PIL:
    Image = _get_image()
    img = Image.open(in_fname)
    # XXX someday we should just try img.thumbnail((max_width, max_height)) ...
    width_in, height_in = img.size
    scale_w = max_width / float(width_in)
    scale_h = max_height / float(height_in)

    if height_in * scale_w <= max_height:
        scale = scale_w
    else:
        scale = scale_h

    if scale >= 1.0 and in_fname == out_fname:
        return

    width_sc = int(round(scale * width_in))
    height_sc = int(round(scale * height_in))

    # resize the image using resize; if using .thumbnail and the image is
    # already smaller than max_width, max_height, then this won't scale up
    # at all (maybe could be an option someday...)
    try:  # Pillow 9+
        bicubic = Image.Resampling.BICUBIC
    except Exception:
        bicubic = Image.BICUBIC
    img = img.resize((width_sc, height_sc), bicubic)
    # img.thumbnail((width_sc, height_sc), Image.BICUBIC)
    # width_sc, height_sc = img.size  # necessary if using thumbnail

    # insert centered
    thumb = Image.new("RGBA", (max_width, max_height), (255, 255, 255, 0))
    pos_insert = ((max_width - width_sc) // 2, (max_height - height_sc) // 2)
    thumb.paste(img, pos_insert)

    try:
        thumb.save(out_fname)
    except OSError:
        # try again, without the alpha channel (e.g., for JPEG)
        thumb.convert("RGB").save(out_fname)


def optipng(fname, args=()):
    """Optimize a PNG in place.

    Parameters
    ----------
    fname : str
        The filename. If it ends with '.png', ``optipng -o7 fname`` will
        be run. If it fails because the ``optipng`` executable is not found
        or optipng fails, the function returns.
    args : tuple
        Extra command-line arguments, such as ``['-o7']``.
    """
    fname = str(fname)
    if fname.endswith(".png"):
        # -o7 because this is what CPython used
        # https://github.com/python/cpython/pull/8032
        try:
            subprocess.check_call(
                ["optipng"] + list(args) + [fname],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except (subprocess.CalledProcessError, OSError):  # FileNotFoundError
            pass


def _has_optipng():
    try:
        subprocess.check_call(
            ["optipng", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except OSError:  # FileNotFoundError
        return False
    else:
        return True


def get_md5sum(src_file, mode="b"):
    """Returns md5sum of file.

    Parameters
    ----------
    src_file : str
        Filename to get md5sum for.
    mode : 't' or 'b'
        File mode to open file with. When in text mode, universal line endings
        are used to ensure consistency in hashes between platforms.
    """
    if mode == "t":
        kwargs = {"errors": "surrogateescape", "encoding": "utf-8"}
    else:
        kwargs = {}
    # Universal newline mode is intentional here
    with open(src_file, f"r{mode}", **kwargs) as src_data:
        src_content = src_data.read()
        if mode == "t":
            src_content = src_content.encode(**kwargs)
        return hashlib.md5(src_content).hexdigest()


def _replace_md5(fname_new, fname_old=None, *, method="move", mode="b", check="md5"):
    fname_new = str(fname_new)  # convert possible Path
    assert method in ("move", "copy")
    if fname_old is None:
        assert fname_new.endswith(".new")
        fname_old = os.path.splitext(fname_new)[0]
    replace = True
    if os.path.isfile(fname_old):
        if check == "md5":  # default
            func = partial(get_md5sum, mode=mode)
        else:
            assert check == "json"

            def func(x):
                return json.loads(Path(x).read_text("utf-8"))

        try:
            equiv = func(fname_old) == func(fname_new)
        except Exception:  # e.g., old JSON file is a problem
            equiv = False
        if equiv:
            replace = False
            if method == "move":
                os.remove(fname_new)
        else:
            logger.debug(f"Replacing stale {fname_old} with {fname_new}")
    if replace:
        if method == "move":
            move(fname_new, fname_old)
        else:
            copyfile(fname_new, fname_old)
    assert os.path.isfile(fname_old)


def check_duplicate_filenames(files):
    """Check for duplicate filenames across gallery directories."""
    # Check whether we'll have duplicates
    used_names = set()
    dup_names = list()

    for this_file in files:
        this_fname = os.path.basename(this_file)
        if this_fname in used_names:
            dup_names.append(this_file)
        else:
            used_names.add(this_fname)

    if len(dup_names) > 0:
        logger.warning(
            "Duplicate example file name(s) found. Having duplicate file "
            "names will break some links. "
            "List of files: %s",
            sorted(dup_names),
        )


def check_spaces_in_filenames(files):
    """Check for spaces in filenames across example directories."""
    regex = re.compile(r"[\s]")
    files_with_space = list(filter(regex.search, files))
    if files_with_space:
        logger.warning(
            "Example file name(s) with spaces found. Having spaces in "
            "file names will break some links. "
            "List of files: %s",
            sorted(files_with_space),
        )


def _collect_gallery_files(examples_dirs, gallery_conf, check_filenames=False):
    """Collect files with `example_extensions`, accounting for `ignore_pattern`.

    If `check_filenames` we check one level of sub-folders as well as root
    `example_dirs` for gallery example files. We then check for duplicate and
    spaces in full file paths.
    """
    exts = gallery_conf["example_extensions"]
    max_depth = 1 if check_filenames else 0
    files = []
    for example_dir in examples_dirs:
        example_depth = os.path.abspath(example_dir).count(os.sep)
        for root, _, filenames in os.walk(example_dir):
            root = os.path.normpath(root)
            if (root.count(os.sep) - example_depth) > max_depth:
                break
            for filename in filenames:
                if (s := Path(filename).suffix) and s in exts:
                    if re.search(gallery_conf["ignore_pattern"], filename) is None:
                        file = filename
                        if check_filenames:
                            file = os.path.join(root, filename)
                        files.append(file)
    if check_filenames:
        check_duplicate_filenames(files)
        check_spaces_in_filenames(files)
    return files


def zip_files(file_list, zipname, relative_to, extension=None):
    """
    Creates a zip file with the given files.

    A zip file named `zipname` will be created containing the files listed in
    `file_list`. The zip file contents will be stored with their paths stripped to be
    relative to `relative_to`.
    """
    zipname_new = str(zipname) + ".new"
    with zipfile.ZipFile(zipname_new, mode="w") as zipf:
        for fname in file_list:
            if extension is not None:
                fname = os.path.splitext(fname)[0] + extension
            zipf.write(fname, os.path.relpath(fname, relative_to))
    _replace_md5(zipname_new)
    return zipname


def _has_pypandoc():
    """Check if pypandoc package available."""
    try:
        import pypandoc  # noqa

        # Import error raised only when function called
        version = pypandoc.get_pandoc_version()
    except (ImportError, OSError):
        return None, None
    else:
        return True, version


def _has_graphviz():
    try:
        import graphviz  # noqa F401
    except ImportError as exc:
        logger.info(
            "`graphviz` required for graphical visualization "
            f"but could not be imported, got: {exc}"
        )
        return False
    return True


def _escape_ansi(s):
    """Remove ANSI terminal formatting characters from a string."""
    return re.sub(r"(?:\x1B[@-_]|[\x80-\x9F])[0-?]*[ -/]*[@-~]", "", s)


def _format_toctree(items, includehidden=False):
    """Format a toc tree."""
    st = """
.. toctree::
   :hidden:"""
    if includehidden:
        st += """
   :includehidden:
"""
    st += """

   {}\n""".format("\n   ".join(items))

    st += "\n"

    return st


_CUSTOM_EXAMPLE_ORDER = [
    "plot_1.py",
    "plot_3.py",
    "plot_2.py",
    "plot_5.py",
    "plot_6.py",
    "plot_4.py",
    "plot_8.py",
    "plot_7.py",
    "plot_9.py",
]


def _custom_example_sorter(filename):
    """Importable custom sorter func, used in our test suite."""
    return _CUSTOM_EXAMPLE_ORDER.index(filename)


def _custom_subsection_sorter(foldername):
    """Importable custom sorter func for subsection folders, used in our test suite."""
    return foldername[::-1]


def custom_minigallery_sort_order_sorter(file):
    """Importable custom sorter for minigallery_sort_order, used in our test suite."""
    ORDER = [
        "plot_3.py",
        "plot_2.py",
        "plot_1.py",
    ]
    return ORDER.index(Path(file).name)


# Should be matched with `_read_json`
def _write_json(target_file, to_save, name=""):
    """Write dictionary to JSON file."""
    codeobj_fname = Path(target_file).with_name(target_file.stem + f"{name}.json.new")
    with open(codeobj_fname, "w", **_W_KW) as fid:
        json.dump(
            to_save,
            fid,
            sort_keys=True,
            ensure_ascii=False,
            indent=1,
            check_circular=False,
        )
    _replace_md5(codeobj_fname, check="json")


def _read_json(json_fname):
    """Read JSON dictionary from file."""
    with open(json_fname, "r", encoding="utf-8") as fid:
        json_dict = json.load(fid)
    return json_dict


def _combine_backreferences(dict_a, dict_b):
    """Combine backreferences dictionaries, joining lists when keys are the same."""
    # `dict_b` is None when `backreferences_dir` config not set
    if isinstance(dict_b, dict):
        for key, value in dict_b.items():
            dict_a.setdefault(key, []).extend(value)
    return dict_a
