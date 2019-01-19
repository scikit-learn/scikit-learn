"""Patches process_generate_options in sphinx.ext.autosummary to add a
``-lowercase`` to function names that have the same case-insenstive name as a
class name.

Assumptions:
1. Class names always begins with a capitalized letter.
2. If a function name matches a class name (case-insensitive), then it is the
only function that matches.
"""
import os
import inspect
from contextlib import contextmanager

from sphinx.ext.autosummary import get_rst_suffix
from sphinx.ext.autosummary import process_generate_options


@contextmanager
def patch_os_path_join(generated_dirname, orig_suffix, new_suffix):
    orig_os_path_join = os.path.join
    all_class_names = set()

    def custom_os_path_join(a, *p):
        path = orig_os_path_join(a, *p)
        if len(p) != 1:
            return path

        dirname, name_with_suffix = os.path.split(path)
        if not (dirname.endswith(generated_dirname) and
                name_with_suffix.startswith("sklearn") and
                name_with_suffix.endswith(orig_suffix)):
            return path

        name_without_suffix = name_with_suffix[:-len(orig_suffix)]
        name_parts = name_without_suffix.split(".")
        if not name_parts:
            return path

        py_name = name_parts[-1]
        if not py_name:
            return path

        name_without_suffix_lower = name_without_suffix.lower()
        # First letter is capitalized, thus a class
        if py_name[0] == py_name[0].upper():
            all_class_names.add(name_without_suffix_lower)
            return path

        if name_without_suffix_lower not in all_class_names:
            return path

        return orig_os_path_join(dirname, name_without_suffix + new_suffix)

    os.path.join = custom_os_path_join
    yield
    os.path.join = orig_os_path_join


def process_generate_options_custom_files(app):
    org_suffix = get_rst_suffix(app)
    new_suffix = "-lowercase" + org_suffix

    with patch_os_path_join("generated", org_suffix, new_suffix):
        process_generate_options(app)


def setup(app):
    app.setup_extension('sphinx.ext.autosummary')

    # Override process_generate_options added by sphinx.ext.autosummary
    builder_inited_listeners = app.events.listeners["builder-inited"]

    for listener_id, obj in builder_inited_listeners.items():
        if (inspect.isfunction(obj)
                and obj.__name__ == "process_generate_options"):
            builder_inited_listeners[listener_id] = \
                process_generate_options_custom_files
            break
