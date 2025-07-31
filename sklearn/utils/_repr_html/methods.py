# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import reprlib
from collections import UserDict
from urllib.parse import quote

from sklearn.utils._repr_html.base import ReprHTMLMixin

CLASS_DOC_URL_PREFIX = "https://scikit-learn.org/{doc_version}/modules/generated/"


def get_module_doc_link(estimator_class, method_name):
    """URL to the relevant section of the docstring using a Text Fragment
    https://developer.mozilla.org/en-US/docs/Web/URI/Reference/Fragment/Text_fragments
    """
    import sklearn

    if hasattr(estimator_class, "__module__"):
        module_name = estimator_class.__module__
    else:
        module_name = None

    if module_name is None or not module_name.startswith("sklearn."):
        # Not a scikit-learn estimator. Do not link to the scikit-learn
        # documentation.
        return None

    if ".dev" in sklearn.__version__:
        doc_version = "dev"
    else:
        doc_version = ".".join(sklearn.__version__.split(".")[:2])
    class_doc_base_url = CLASS_DOC_URL_PREFIX.format(doc_version=doc_version)

    # Strip private submodule component if any:
    if "._" in module_name:
        module_name = module_name.split("._")[0]

    class_name = estimator_class.__qualname__

    doc = getattr(estimator_class, method_name).__doc__
    if doc:
        docstring = doc[:50]
    else:
        docstring = ""

    base_url = f"{class_doc_base_url}{quote(module_name)}.{quote(class_name)}.html"
    method_fragment = f"{quote(module_name)}.{class_name}.{method_name}"

    return f"{base_url}#{method_fragment}", docstring


def _doc_row(estimator_class, method_name):
    """
    Generate an HTML table row containing a link to the online
    documentation for a specific parameter of an estimator.
    If the link cannot be generated, an empty string is returned.
    """

    link, docstring = get_module_doc_link(estimator_class, method_name)

    if link:
        link_string = (
            f'rel="noreferrer" target="_blank" href='
            f"{link} "
            f'style="color: white; background: black;">?'
            # style adds bottom space to hover box
            f'<span style="bottom: 2rem; background-color:white;">{docstring}... '
            f"<br> click to see online documentation</span>"
        )
    else:
        link_string = (
            f'style="color: white; background: black;">?<span>Online documentation'
            f"for `{method_name}`not found </span>"
        )

    return link_string


def _methods_html_repr(methods):
    """Generate HTML representation of estimator methods.
    Creates an HTML table with methods names, signature
    and link to its documentation. It is wrapped in a collapsible
    details element.
    """

    HTML_TEMPLATE = """
       <div class="method-table">
           <details>
               <summary>Methods</summary>
               <table class="methods-table">
                 <tbody>
                   {rows}
                 </tbody>
               </table>
           </details>
       </div>
    """
    ROW_TEMPLATE = """
       <tr class="methods">
           <td class="methods-name">{name}&nbsp; </td>
           <td class="methods-signature">{signature}</td>
           <td><a class="sk-estimator-doc-link
                                 {doc_link}
                                </a>
            </td>
       </tr>
    """
    r = reprlib.Repr()
    r.maxlist = 2  # Show only first 2 items of lists
    r.maxtuple = 1  # Show only first item of tuples
    r.maxstring = 50  # Limit string length
    rows = [
        ROW_TEMPLATE.format(
            name=name,
            signature=r.repr(signature).replace("'", ""),
            doc_link=_doc_row(methods.estimator, name),
        )
        for name, signature in methods.items()
    ]

    return HTML_TEMPLATE.format(rows="\n".join(rows))


class MethodsDict(ReprHTMLMixin, UserDict):
    """
    A dictionary-like object for storing methods of an estimator,
    with an HTML representation.

    This class extends `UserDict` to store methods and their signatures,
    and provides an HTML representation for visualizing the methods in a
    collapsible table format.

    Parameters
    ----------
    methods : dict
        A dictionary where keys are method names and values are their
        signatures.
    estimator_class : class
        The class of the estimator whose methods are being represented.
    """

    _html_repr = _methods_html_repr

    def __init__(self, methods, estimator_class):
        super().__init__(methods)
        self.estimator = estimator_class
