# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

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


def _read_params(value):
    r = reprlib.Repr()
    r.maxlist = 2
    r.maxdict = 1
    r.maxstring = 50
    cleaned_value = r.repr(value)

    return cleaned_value


@lru_cache
def _scrape_estimator_docstring(docstring):
    return docscrape.NumpyDocString(docstring)


def _fitted_attr_html_repr(fitted_attributes):
    """Generate HTML representation of estimator fitted attributes.

    Creates an HTML table with fitted attribute names and values
    wrapped in a collapsible details element. When attributes are arrays,
    shape is shown.
    """

    HTML_TEMPLATE = """
       <div class="estimator-table">
           <details>
               <summary>Fitted attributes</summary>
               <table class="body-table">
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
    ROW_TEMPLATE = """
       <tr class="default">
           <td>{name}&nbsp;</td>
           <td>{type}</td>
           <td>{shape}</td>
           <td>{dtype}</td>
           <td>{attr_size}</td>
           <td>{attr_value}</td>
       </tr>
    """

    FITTED_ATTR_AVAILABLE_DOC_LINK_TEMPLATE = """
        <a class="param-doc-link"
            rel="noreferrer" target="_blank" href="{link}">
            {param_name}
            <span class="param-doc-description">{param_description}</span>
        </a>
    """
    estimator_class_docs = inspect.getdoc(fitted_attributes.estimator_class)

    rows = []
    for fitted_attr_name, attr_info in fitted_attributes.items():
        link = _generate_link_to_param_doc(
            fitted_attributes.estimator_class,
            fitted_attr_name,
            fitted_attributes.doc_link,
        )
        formated_attr_value = _read_params(attr_info[1])

        if len(attr_info) == 2:
            rows.append(
                ROW_TEMPLATE.format(
                    name=fitted_attr_name,
                    type=attr_info[0],
                    shape="",
                    dtype="",
                    attr_size="",
                    attr_value=formated_attr_value,
                )
            )
        else:
            rows.append(
                ROW_TEMPLATE.format(
                    name=fitted_attr_name,
                    type=attr_info[0],
                    shape=attr_info[1],
                    dtype=attr_info[2],
                    attr_size=attr_info[3],
                    attr_value="",
                )
            )

    return HTML_TEMPLATE.format(rows="\n".join(rows))


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
