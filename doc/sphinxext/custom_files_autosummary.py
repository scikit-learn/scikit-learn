"""
Patches os.path.join when process_generate_options to map filenames.

The custom_autosummary_file_map sphinx configuration option is used to map
filenames to custom filenames:

custom_autosummary_include_hash_in_suffix = {
    "sklearn.cluster.dbscan", "sklearn.cluster.optics",
    "sklearn.covariance.oas"
}
"""
import os
import inspect
from contextlib import contextmanager

from sphinx.ext.autosummary import get_rst_suffix
from sphinx.ext.autosummary import process_generate_options


@contextmanager
def patch_os_path_join(hashed_set_with_suffix, suffix):
    orig_os_path_join = os.path.join

    def custom_os_path_join(a, *p):
        path = orig_os_path_join(a, *p)
        if len(p) != 1:
            return path
        name_with_suffix = p[0]

        if name_with_suffix not in hashed_set_with_suffix:
            return path

        name_wo_suffix = name_with_suffix.rstrip(suffix)
        return path.replace(
            name_with_suffix,
            name_wo_suffix + str(hash(name_wo_suffix)) + suffix
        )

    os.path.join = custom_os_path_join
    yield
    os.path.join = orig_os_path_join


def process_generate_options_custom_files(app):
    """Patches os.path.join to replace filenames with ones that are in
    custom_autosummary_file_map
    """

    suffix = get_rst_suffix(app)

    hashed_set_with_suffix = {
        name + suffix for name in
        app.config.custom_autosummary_include_hash_in_suffix
    }

    with patch_os_path_join(hashed_set_with_suffix, suffix):
        process_generate_options(app)


def setup(app):
    app.setup_extension('sphinx.ext.autosummary')
    app.add_config_value(
        "custom_autosummary_include_hash_in_suffix", set(), None)

    # Override process_generate_options added by sphinx.ext.autosummary
    builder_inited_listeners = app.events.listeners["builder-inited"]

    for listener_id, obj in builder_inited_listeners.items():
        if (inspect.isfunction(obj)
                and obj.__name__ == "process_generate_options"):
            builder_inited_listeners[listener_id] = \
                process_generate_options_custom_files
            break
