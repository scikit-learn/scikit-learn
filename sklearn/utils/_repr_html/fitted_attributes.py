# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import re
from collections import UserDict
from urllib.parse import quote

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
                   {rows}
                 </tbody>
               </table>
           </details>
       </div>
    """
    ROW_TEMPLATE = """
       <tr class="default">
           <td>{name}&nbsp;</td>
           <td>{value}</td>
       </tr>
    """

    rows = [
        ROW_TEMPLATE.format(name=name, value=value)
        for name, value in fitted_attributes.items()
    ]
    # for row in fitted_attributes:
    # link = _generate_link_to_param_doc(row.estimator_class, row, row.doc_link)

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
