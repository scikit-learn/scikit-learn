# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
import re
import reprlib
from collections import UserDict
from urllib.parse import quote

from sklearn.utils._repr_html.base import ReprHTMLMixin

CLASS_DOC_URL_PREFIX = "https://scikit-learn.org/{doc_version}/modules/generated/"


def link_to_param_doc(estimator_class, param_name, doc_link):
    """URL to the relevant section of the docstring using a Text Fragment

    https://developer.mozilla.org/en-US/docs/Web/URI/Reference/Fragment/Text_fragments
    """
    docstring = estimator_class.__doc__
    m = re.search(f"{param_name} *: *(.+)", docstring)
    if m is None:
        # No match found in the docstring, return None to indicate that we
        # cannot link.
        return None

    # Extract the whole line of the type information, up to the line break as
    # disambiguation suffix to build the fragment
    param_type = m.group(1)

    text_fragment = f"{quote(param_name)}{quote(param_type)}"

    return f"{doc_link}#:~:text={text_fragment}"


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


def _doc_row(estimator_class, row, param_name, doc_link):
    """
    Generate an HTML table row containing a link to the online
    documentation for a specific parameter of an estimator.
    If the link cannot be generated, an empty string is returned.
    """

    link = link_to_param_doc(estimator_class, row, doc_link)

    if link:
        link_string = (
            f'rel="noreferrer" target="_blank" href='
            f"{link} "
            f'style="color: white; background: black;">?<span>Documentation'
            f" for `{param_name}`</span>"
        )
    else:
        link_string = (
            f'style="color: white; background: black;">?<span>Documentation'
            f" for `{param_name}` not found </span>"
        )

    return link_string


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
            <td class="doc_link"><a class="sk-estimator-doc-link"
                                 {doc_link}
                                </a>
            </td>
        </tr>
    """
    rows = []
    for row in params:
        param = _read_params(row, params[row], params.non_default)
        rows.append(
            ROW_TEMPLATE.format(
                **param,
                doc_link=_doc_row(
                    params.estimator_class, row, param["param_name"], params.doc_link
                ),
            )
        )

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

    estimator_class : type
        The class of the estimator. It allows to find the online documentation
        link for each paramter.
    """

    _html_repr = _params_html_repr

    def __init__(
        self, params=None, non_default=tuple(), estimator_class=None, doc_link=""
    ):
        super().__init__(params or {})
        self.non_default = non_default
        self.estimator_class = estimator_class
        self.doc_link = doc_link
