"""Patches process_generate_options in sphinx.ext.autosummary to add a hash
to filenames that have the same case-insensitive name as another generated
file.
"""
import os
import inspect
from contextlib import contextmanager
from collections import defaultdict
from hashlib import sha1

from sphinx.ext.autosummary import get_rst_suffix
from sphinx.ext.autosummary import process_generate_options


@contextmanager
def patch_os_path_join(generated_dirname, suffix):
    orig_os_path_join = os.path.join
    all_names = defaultdict(set)

    def custom_os_path_join(a, *p):
        path = orig_os_path_join(a, *p)
        if len(p) != 1:
            return path

        dirname, name_with_suffix = os.path.split(path)
        if not dirname.endswith(generated_dirname):
            return path

        name_with_suffix_lower = name_with_suffix.lower()
        if name_with_suffix_lower not in all_names:
            all_names[name_with_suffix_lower].add(name_with_suffix)
            return path

        if name_with_suffix in all_names[name_with_suffix_lower]:
            return path

        all_names[name_with_suffix_lower].add(name_with_suffix)

        # Include hash with name
        name_wo_suffix = name_with_suffix.rstrip(suffix)
        new_name_w_suffix = "{}{}{}".format(
            name_wo_suffix,
            sha1(name_wo_suffix.encode()).hexdigest(),
            suffix
        )
        return orig_os_path_join(dirname, new_name_w_suffix)

    os.path.join = custom_os_path_join
    yield
    os.path.join = orig_os_path_join


def process_generate_options_custom_files(app):
    """Patches os.path.join to add a hash to files that have the same
    case-insensitive name as another generated file.
    """
    with patch_os_path_join("generated", get_rst_suffix(app)):
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
