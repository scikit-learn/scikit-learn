# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
import inspect
import re
import reprlib
from collections import UserDict
from functools import lru_cache
from urllib.parse import quote

from sklearn.externals._numpydoc import docscrape
from sklearn.utils._repr_html.base import ReprHTMLMixin


def _generate_link_to_param_doc(estimator_class, param_name, doc_link):
    """URL to the relevant section of the docstring using a Text Fragment

    https://developer.mozilla.org/en-US/docs/Web/URI/Reference/Fragment/Text_fragments
    """
    docstring = estimator_class.__doc__

    m = re.search(f"{param_name} : (.+)\\n", docstring or "")

    if m is None:
        # No match found in the docstring, return None to indicate that we
        # cannot link.
        return None

    # Extract the whole line of the type information, up to the line break as
    # disambiguation suffix to build the fragment
    param_type = m.group(1)
    text_fragment = f"{quote(param_name)},-{quote(param_type)}"

    return f"{doc_link}#:~:text={text_fragment}"


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


@lru_cache
def _scrape_estimator_docstring(docstring):
    return docscrape.NumpyDocString(docstring)


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
            rel="noreferrer" target="_blank" href="{link}">
            {param_name}
            <span class="param-doc-description">{param_description}</span>
        </a>
    """
    estimator_class_docs = inspect.getdoc(params.estimator_class)
    if estimator_class_docs and (
        structured_docstring := _scrape_estimator_docstring(estimator_class_docs)
    ):
        param_map = {
            param_docstring.name: param_docstring
            for param_docstring in structured_docstring["Parameters"]
        }
    else:
        param_map = {}
    rows = []
    for row in params:
        param = _read_params(row, params[row], params.non_default)
        link = _generate_link_to_param_doc(params.estimator_class, row, params.doc_link)
        if param_numpydoc := param_map.get(row, None):
            param_description = (
                f"{param_numpydoc.name}: {param_numpydoc.type}<br><br>"
                f"{'<br>'.join(param_numpydoc.desc)}"
            )
        else:
            param_description = None

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
