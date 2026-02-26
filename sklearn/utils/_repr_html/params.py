# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
import reprlib
from collections import UserDict

from sklearn.utils._repr_html.base import ReprHTMLMixin
from sklearn.utils._repr_html.common import (
    generate_link_to_param_doc,
    get_docstring,
)


def _read_params(name, value, non_default_params):
    """Categorizes parameters as 'default' or 'user-set' and formats their values.
    Escapes or truncates parameter values for display safety and readability.
    """
    name = html.escape(name)
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
    PARAMS_TABLE_TEMPLATE = """
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

    PARAM_ROW_TEMPLATE = """
        <tr class="{param_type}">
            <td><i class="copy-paste-icon"
                 onclick="copyToClipboard('{param_name}',
                          this.parentElement.nextElementSibling)"
            ></i></td>
            <td class="param">{param_display}</td>
            <td class="value">{param_value}</td>
        </tr>
    """

    PARAM_AVAILABLE_DOC_LINK_TEMPLATE = """
        <a class="param-doc-link"
            style="anchor-name: --doc-link-{param_name};"
            rel="noreferrer" target="_blank" href="{link}">
            {param_name}
            <span class="param-doc-description"
            style="position-anchor: --doc-link-{param_name};">
            {param_description}</span>
        </a>
    """

    rows = []
    for row in params:
        param = _read_params(row, params[row], params.non_default)
        link = generate_link_to_param_doc(params.estimator_class, row, params.doc_link)

        param_description = get_docstring(params.estimator_class, "Parameters", row)

        if params.doc_link and link and param_description:
            # Create clickable parameter name with documentation link
            param_display = PARAM_AVAILABLE_DOC_LINK_TEMPLATE.format(
                link=link,
                param_name=param["param_name"],
                param_description=param_description,
            )
        else:
            # Just show the parameter name without link
            param_display = param["param_name"]

        rows.append(PARAM_ROW_TEMPLATE.format(**param, param_display=param_display))

    return PARAMS_TABLE_TEMPLATE.format(rows="\n".join(rows))


class ParamsDict(ReprHTMLMixin, UserDict):
    """Dictionary-like class to store and provide an HTML representation.

    It builds an HTML structure to be used with Jupyter notebooks or similar
    environments. It allows storing metadata to track non-default parameters.

    Parameters
    ----------
    params : dict, default=None
        The original dictionary of parameters and their values.

    non_default : tuple, default=(,)
        The list of non-default parameters.

    estimator_class : type, default=None
        The class of the estimator. It allows to find the online documentation
        link for each parameter.

    doc_link : str, default=""
        The base URL to the online documentation for the estimator class.
        Used to generate parameter-specific documentation links in the HTML
        representation. If empty, documentation links will not be generated.
    """

    _html_repr = _params_html_repr

    def __init__(
        self, *, params=None, non_default=tuple(), estimator_class=None, doc_link=""
    ):
        super().__init__(params or {})
        self.non_default = non_default
        self.estimator_class = estimator_class
        self.doc_link = doc_link
