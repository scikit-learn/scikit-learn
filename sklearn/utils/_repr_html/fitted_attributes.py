# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.utils._repr_html.base import ReprHTMLMixin


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

    return HTML_TEMPLATE.format(rows="\n".join(rows))


class AttrsDict(ReprHTMLMixin, dict):
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
