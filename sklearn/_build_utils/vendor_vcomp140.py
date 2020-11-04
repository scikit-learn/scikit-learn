"""Embed vcomp140.dll after generating the scikit-learn Windows wheel."""


import os
import os.path as op
import zipfile
import textwrap
from glob import glob


GLOB_PATTERN = "*.whl"
DISTRIBUTOR_INIT = op.join("sklearn", "_distributor_init.py")
VCOMP140_SRC_PATH = "C:\\Windows\System32\\vcomp140.dll"  # noqa


def _delete_file(zipf, filename):
    """Delete a file from the ZIP file."""
    zipf_copy = zipfile.ZipFile(zipf.filename + ".copy", "a")

    for item in zipf.infolist():
        # The buffer hold the read file
        buffer = zipf.read(item.filename)

        if item.filename != filename:
            zipf_copy.writestr(item, buffer)

    return zipf_copy


def make_distributor_init(zipf, dll_filename):
    """Create a _distributor_init.py file for the vcomp140.dll.

    This file is imported first when importing the
    sklearn package so as to pre-load the vendored
    vcomp140.dll.
    """
    zipf.writestr(DISTRIBUTOR_INIT, textwrap.dedent("""
        '''Helper to preload vcomp140.dll to prevent "not found" errors.

        Once the vcomp140.dll is preloaded, the namespace is made
        available to any subsequent vcomp140.dll. This is created
        as part of the scripts that build the wheel.
        '''

        import os
        import os.path as op
        from ctypes import WinDLL


        if os.name == "nt":
            # Load the vcomp140.dll in sklearn/.libs by convention
            dll_path = op.join(op.dirname(__file__), ".libs", "{0}")
            WinDLL(op.abspath(dll_path))
        """.format(dll_filename)))


def embed_vcomp140(wheel_dirname):
    """Embed the vcomp140.dll in the wheel."""
    if not op.exists(VCOMP140_SRC_PATH):
        raise ValueError(f"Could not find {VCOMP140_SRC_PATH}.")

    if not op.isdir(wheel_dirname):
        raise RuntimeError(f"Could not find {wheel_dirname} folder.")

    wheel_file_glob_pattern = op.join(wheel_dirname, GLOB_PATTERN)
    wheel_files = glob(wheel_file_glob_pattern)

    if len(wheel_file_glob_pattern) == 0:
        raise RuntimeError(f"Could not find wheels matching "
                           f"{wheel_file_glob_pattern}.")

    if len(wheel_files) > 1:
        raise RuntimeError(f"Found too many wheels: "
                           f"{', '.join(wheel_files)}.")

    wheel_file = wheel_files[0]
    dll_filename = op.basename(VCOMP140_SRC_PATH)
    target_folder = op.join("sklearn", ".libs", "vcomp140.dll")

    # Remove the _distributor_init file to avoid duplicates
    with zipfile.ZipFile(wheel_file, "a") as zipf:
        zipf_copy = _delete_file(zipf, DISTRIBUTOR_INIT)

    # Remove temporary artifacts
    os.remove(zipf.filename)

    # Rename to the original filename
    os.rename(zipf_copy.filename, zipf.filename)

    with zipfile.ZipFile(wheel_file, "a") as zipf:
        print(f"Copying {VCOMP140_SRC_PATH} to {target_folder}.")
        zipf.write(VCOMP140_SRC_PATH, target_folder)

        # Generate the _distributor_init file in the source tree
        print("Generating the '_distributor_init.py' file.")
        make_distributor_init(zipf, dll_filename)
