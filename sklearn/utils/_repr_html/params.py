# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
import reprlib
from collections import UserDict

from sklearn.utils._repr_html.base import ReprHTMLMixin


def _read_params(name, value, non_default_params):
    """Categorizes parameters as 'default' or 'user-set' and formats their values.
    Escapes or truncates parameter values for display safety and readability.
    """
    r = reprlib.Repr()
    r.maxlist = 2  # Show only first 2 items of lists
    r.maxtuple = 1  # Show only first item of tuples
    r.maxstring = 50  # Limit string length
    cleaned_value = html.escape(r.repr(value))

    param_type = "user-set" if name in non_default_params else "default"

    return {"param_type": param_type, "param_name": name, "param_value": cleaned_value}


def _params_html_repr(params):
    """Generate HTML representation of estimator parameters.

    Creates an HTML table with parameter names and values, wrapped in a
    collapsible details element. Parameters are styled differently based
    on whether they are default or user-set values.
    """
    HTML_TEMPLATE = """
        <div class="estimator-table">
            <details>
                <summary>Parameters</summary>
                <table class="parameters-table">
                  <tbody>
                    {rows}
                  </tbody>
                </table>
            </details>
        </div>
    """
    ROW_TEMPLATE = """
        <tr class="{param_type}">
            <td><i class="copy-paste-icon"
                 onclick="copyToClipboard('{param_name}',
                          this.parentElement.nextElementSibling)"
            ></i></td>
            <td class="param">{param_name}&nbsp;</td>
            <td class="value">{param_value}</td>
        </tr>
    """

    rows = [
        ROW_TEMPLATE.format(**_read_params(name, value, params.non_default))
        for name, value in params.items()
    ]

    return HTML_TEMPLATE.format(rows="\n".join(rows))


class ParamsDict(ReprHTMLMixin, UserDict):
    """Dictionary-like class to store and provide an HTML representation.

    It builds an HTML structure to be used with Jupyter notebooks or similar
    environments. It allows storing metadata to track non-default parameters.

    Parameters
    ----------
    params : dict, default=None
        The original dictionary of parameters and their values.

    non_default : tuple
        The list of non-default parameters.
    """

    _html_repr = _params_html_repr

    def __init__(self, params=None, non_default=tuple()):
        super().__init__(params or {})
        self.non_default = non_default
