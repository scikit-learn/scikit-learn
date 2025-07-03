# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.utils._repr_html.base import ReprHTMLMixin


def _methods_html_repr(methods_dict):
    """Generate HTML representation of estimator methods.
    Creates an HTML table with methods names, signature
    and link to its documentation. It is wrapped in a collapsible
    details element.
    """

    HTML_TEMPLATE = """
       <div class="estimator-table">
           <details>
               <summary>Methods</summary>
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
           <td>{signature}</td>

       </tr>
    """
    # Add this <td>{doclink}</td>
    rows = [
        ROW_TEMPLATE.format(name=name, signature=signature)  # Fix me: doclink=doclink)
        for name, signature in methods_dict.items()
    ]

    return HTML_TEMPLATE.format(rows="\n".join(rows))


class MethodsDict(ReprHTMLMixin, dict):
    _html_repr = _methods_html_repr
