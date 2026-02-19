# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
import reprlib
from collections import UserDict

import numpy as np

from sklearn.utils._repr_html.base import ReprHTMLMixin
from sklearn.utils._repr_html.common import (
    generate_link_to_param_doc,
    get_docstring,
)


def _read_fitted_attr(value):
    r = reprlib.Repr()
    for attr in (
        "maxlist",
        "maxdict",
        "maxtuple",
        "maxset",
        "maxfrozenset",
        "maxdeque",
        "maxarray",
    ):
        setattr(r, attr, 4)
    r.maxstring = 9

    if isinstance(value, float):
        value = round(value, 3)

    if isinstance(value, np.ndarray):
        value = np.array2string(
            value, precision=2, separator=",", suppress_small=True, threshold=4
        )
        r.maxstring = 18

        return html.escape(value)

    return html.escape(r.repr(value))


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
                        <th>Name</th>
                        <th>Type</th>
                        <th>Shape</th>
                        <th>dtype</th>
                        <th>Value</th>
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

    rows = []
    # for fitted_attr_name, attr_info in fitted_attributes.items():
    for name, value in fitted_attributes.items():
        link = generate_link_to_param_doc(
            fitted_attributes.estimator_class,
            name,
            fitted_attributes.doc_link,
        )
        fitted_attr_description = get_docstring(
            fitted_attributes.estimator_class, "Attributes", name
        )

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
            html_row_values = (value[0], "", "", _read_fitted_attr(value[1]))
        else:
            html_row_values = (
                value[0],
                value[1],
                value[2],
                _read_fitted_attr(value[3]),
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
    fitted_attrs : dict, default=None
        Dictionary of fitted attributes and their values.

    estimator_class : type, default=None
        The class of the estimator. It allows to find the online documentation
        link for each parameter.

    doc_link : str, default=""
        The base URL to the online documentation for the estimator class.
        Used to generate parameter-specific documentation links in the HTML
        representation. If empty, documentation links will not be generated.
    """

    _html_repr = _fitted_attr_html_repr

    def __init__(self, *, fitted_attrs=None, estimator_class=None, doc_link=""):
        super().__init__(fitted_attrs or {})
        self.estimator_class = estimator_class
        self.doc_link = doc_link
