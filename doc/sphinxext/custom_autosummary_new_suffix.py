"""Patches process_generate_options in sphinx.ext.autosummary to add change
the suffix of items in configuration option
custom_autosummary_names_with_new_suffix to
custom_autosummary_new_suffix. For example:

```python
custom_autosummary_names_with_new_suffix = {
    "sklearn.cluster.dbscan",
    "sklearn.cluster.optics",
    "sklearn.covariance.oas",
    "sklearn.decomposition.fastica"
}
custom_autosummary_new_suffix = "-lowercase.rst"
custom_autosummary_generated_dirname = os.path.join('modules', 'generated')

This extension monkeypatches `os.path.join` used by
https://github.com/sphinx-doc/sphinx/blob/7ffd6ccee8b0c6316159c4295e2f44f8c57b90d6/sphinx/ext/autosummary/generate.py#L168
to return a filename with a new suffix.
```
"""
import os
import inspect
from contextlib import contextmanager

from sphinx.ext.autosummary import get_rst_suffix
from sphinx.ext.autosummary import process_generate_options


@contextmanager
def patch_os_path_join(generated_dirname, filename_map):
    orig_os_path_join = os.path.join

    def custom_os_path_join(a, *p):
        path = orig_os_path_join(a, *p)
        if len(p) != 1:
            return path

        dirname, name_with_suffix = os.path.split(path)
        if not (dirname.endswith(generated_dirname) and
                name_with_suffix in filename_map):
            return path

        return orig_os_path_join(dirname, filename_map[name_with_suffix])

    os.path.join = custom_os_path_join
    yield
    os.path.join = orig_os_path_join


def process_generate_options_custom_files(app):

    orig_suffix = get_rst_suffix(app)
    new_suffix = app.config.custom_autosummary_new_suffix
    generated_dirname = app.config.custom_autosummary_generated_dirname
    filename_map = {
        name + orig_suffix: name + new_suffix
        for name in app.config.custom_autosummary_names_with_new_suffix
    }

    with patch_os_path_join(generated_dirname, filename_map):
        process_generate_options(app)


def setup(app):
    app.setup_extension('sphinx.ext.autosummary')
    app.add_config_value(
        'custom_autosummary_names_with_new_suffix', set(), None)
    app.add_config_value(
        'custom_autosummary_new_suffix', '-lowercase.rst', None)
    app.add_config_value(
        'custom_autosummary_generated_dirname', '', None)

    # Override process_generate_options added by sphinx.ext.autosummary
    builder_inited_listeners = app.events.listeners["builder-inited"]

    for listener_id, obj in builder_inited_listeners.items():
        if (inspect.isfunction(obj)
                and obj.__name__ == "process_generate_options"):
            builder_inited_listeners[listener_id] = \
                process_generate_options_custom_files
            break
