# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
import pprint
from collections import UserDict

from sklearn._config import get_config


class ParamsDict(UserDict):
    """Dictionary-like class to visualize parameters in HTML.

    It builds an HTML structure to be used with Jupyter notebooks or similar
    environments. It allows storing metadata to track non-default parameters.
    """

    def __init__(self, params=None, non_default=tuple()):
        super().__init__(params or {})
        self.non_default = non_default

    @property
    def _repr_html_(self):
        # Taken from sklearn.base.BaseEstimator
        """HTML representation of estimator.

        This is redundant with the logic of `_repr_mimebundle_`. The latter
        should be favored in the long term, `_repr_html_` is only
        implemented for consumers who do not interpret `_repr_mimbundle_`.
        """
        if get_config()["display"] != "diagram":
            raise AttributeError(
                "_repr_html_ is only defined when the "
                "'display' configuration option is set to "
                "'diagram'"
            )
        return self._repr_html_inner

    def _repr_html_inner(self):
        return _html_template(self)

    def _repr_mimebundle_(self, **kwargs):
        """Mime bundle used by jupyter kernels to display estimator"""
        output = {"text/plain": repr(self)}
        if get_config()["display"] == "diagram":
            output["text/html"] = _html_template(self)


def _read_params(param, value, non_default_params, row_template):

    if value != "deprecated" and isinstance(value, str):
        cleaned_value = f'"{value}"'
    else:
        cleaned_value = html.escape(str(value))
    if len(cleaned_value) > 50:
        if param == "param_distributions":
            formatted_value = pprint.pformat(cleaned_value)
            cleaned_value = f"<pre>{formatted_value}</pre>"
        else:
            cleaned_value = "(...)"
    param_type = "user-set" if param in non_default_params else "default"

    return row_template.format(
        param_type=param_type, param_name=param, param_value=cleaned_value
    )


def _html_template(params):

    HTML_TEMPLATE = """
        <div class="estimator-table">
            <details>
                <summary>Parameters</summary>
                <table>
                  <tbody>
                    {rows}
                  </tbody>
                </table>
            </details>
        </div>
    """
    ROW_TEMPLATE = """
        <tr class="{param_type}">
            <td><i class="fa-regular fa-copy"
                 onclick="copyToClipboard('{param_name}',
                          this.parentElement.nextElementSibling)"
            ></i></td>
            <td class="param">{param_name}&nbsp;</td>
            <td class="value">{param_value}</td>
        </tr>
    """

    rows = [
        _read_params(param, value, params.non_default, ROW_TEMPLATE)
        for param, value in params.items()
    ]

    return HTML_TEMPLATE.format(rows="\n".join(rows))
