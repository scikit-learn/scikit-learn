"""
Patches os.path.join when process_generate_options to map filenames.

The custom_autosummary_file_map sphinx configuration option is used to map
filenames to custom filenames:

custom_autosummary_file_map = {
    "sklearn.cluster.dbscan": "sklearn.cluster.dbscan_lowercase",
    "sklearn.cluster.optics": "sklearn.cluster.optics_lowercase"
}
"""
import os
import inspect
from contextlib import suppress

from sphinx.ext.autosummary import get_rst_suffix
from sphinx.ext.autosummary import process_generate_options


def process_generate_options_custom_files(app):
    """Patches os.path.join to replace filenames with ones that are in
    custom_autosummary_file_map
    """

    orig_os_path_join = os.path.join
    suffix = get_rst_suffix(app)

    custom_autosummary_file_map_suffix = {
        k + suffix: v + suffix
        for k, v in app.config.custom_autosummary_file_map.items()
    }

    def custom_os_path_join(a, *p):
        if len(p) != 1:
            return orig_os_path_join(a, *p)
        name_with_suffix = p[0]
        with suppress(KeyError):
            name_with_suffix = custom_autosummary_file_map_suffix[
                name_with_suffix]
        return orig_os_path_join(a, name_with_suffix)

    os.path.join = custom_os_path_join
    process_generate_options(app)
    os.path.join = orig_os_path_join


def setup(app):
    app.setup_extension('sphinx.ext.autosummary')
    app.add_config_value("custom_autosummary_file_map", {}, None)

    # Override process_generate_options added by sphinx.ext.autosummary
    builder_inited_listeners = app.events.listeners["builder-inited"]

    for listener_id, obj in builder_inited_listeners.items():
        if (inspect.isfunction(obj)
                and obj.__name__ == "process_generate_options"):
            builder_inited_listeners[listener_id] = \
                process_generate_options_custom_files
            break
