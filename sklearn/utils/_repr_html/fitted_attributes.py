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


def _read_fitted_attr(value):
    r = reprlib.Repr()
    r.maxlist = 2
    r.maxdict = 1
    r.maxstring = 50

    if isinstance(value, float):
        value = round(value, 3)

    return html.escape(r.repr(value))


@lru_cache
def _scrape_estimator_docstring(docstring):
    return docscrape.NumpyDocString(docstring)


def _fitted_attr_html_repr(fitted_attributes):
    """Generate HTML representation of estimator fitted attributes.

    Creates an HTML table with fitted attribute names and values
    wrapped in a collapsible details element. When attributes are arrays,
    shape is shown.
    """

    FITTED_ATTR_TEMPLATE = """
        <div class="estimator-table">
            <details>
                <summary>Fitted attributes</summary>
                <table class="parameters-table">
                    <tbody>
                        <tr>
                        <th>Class</th>
                        <th>Type</th>
                        <th>Shape</th>
                        <th>dtype</th>
                        <th>size</th>
                        <th>value</th>
                        </tr>
                        {rows}
                    </tbody>
                </table>
            </details>
        </div>
    """
    FITTED_ATTR_ROW_TEMPLATE = """
       <tr class="default">
           <td class="param">{fitted_attr_display}</td>
           <td>{0}</td>
           <td>{1}</td>
           <td>{2}</td>
           <td>{3}</td>
           <td>{4}</td>
       </tr>
    """

    FITTED_ATTR_AVAILABLE_DOC_LINK_TEMPLATE = """
        <a class="param-doc-link"
            style="anchor-name: --doc-link-{fitted_attr_name};"
            rel="noreferrer" target="_blank" href="{link}">
            {fitted_attr_name}
            <span class="param-doc-description"
            style="position-anchor: --doc-link-{fitted_attr_name};">
            {fitted_attr_description}</span>
        </a>
    """
    estimator_class_docs = inspect.getdoc(fitted_attributes.estimator_class)
    if estimator_class_docs and (
        structured_docstring := _scrape_estimator_docstring(estimator_class_docs)
    ):
        fitted_attr_map = {
            fitted_attr_docstring.name: fitted_attr_docstring
            for fitted_attr_docstring in structured_docstring["Attributes"]
        }
    else:
        fitted_attr_map = {}

    rows = []
    # for fitted_attr_name, attr_info in fitted_attributes.items():
    for name, value in fitted_attributes.items():
        link = _generate_link_to_param_doc(
            fitted_attributes.estimator_class,
            name,
            fitted_attributes.doc_link,
        )

        if fitted_attr_numpydoc := fitted_attr_map.get(name, None):
            escaped_lines = (html.escape(line) for line in fitted_attr_numpydoc.desc)
            fitted_attr_description = (
                f"{html.escape(fitted_attr_numpydoc.name)}:"
                f"{html.escape(fitted_attr_numpydoc.type)}<br><br>"
                f"{'<br>'.join(escaped_lines)}"
            )
        else:
            fitted_attr_description = None

        if fitted_attributes.doc_link and link and fitted_attr_description:
            # Create clickable parameter name with documentation link
            fitted_attr_display = FITTED_ATTR_AVAILABLE_DOC_LINK_TEMPLATE.format(
                link=link,
                fitted_attr_name=html.escape(name),
                fitted_attr_description=fitted_attr_description,
            )
        else:
            # Just show the parameter name without link
            fitted_attr_display = html.escape(name)
        if len(value) == 2:
            html_row_values = (value[0], "", "", "", _read_fitted_attr(value[1]))
        else:
            html_row_values = (
                value[0],
                _read_fitted_attr(value[1]),
                value[2],
                value[3],
                "",
            )

        rows.append(
            FITTED_ATTR_ROW_TEMPLATE.format(
                *html_row_values,
                fitted_attr_display=fitted_attr_display,
            )
        )

    return FITTED_ATTR_TEMPLATE.format(rows="\n".join(rows))


class AttrsDict(ReprHTMLMixin, UserDict):
    """Dictionary-like class to store and provide an HTML representation.

    It builds an HTML structure to be used with Jupyter notebooks or similar
    environments.

    Parameters
    ----------
    fitted_attributes : dict, default=None
        Dictionary of fitted attributes and their values. When this is
        an array, it includes its size.
    """

    _html_repr = _fitted_attr_html_repr

    def __init__(self, *, fitted_attrs=None, estimator_class=None, doc_link=""):
        super().__init__(fitted_attrs or {})
        self.estimator_class = estimator_class
        self.doc_link = doc_link
