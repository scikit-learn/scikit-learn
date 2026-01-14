# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
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

    # doc = inspect.getdoc(getattr(estimator_class, method_name))
    doc = getattr(estimator_class, method_name).__doc__
    if doc:
        docstring = doc[:50]
    else:
        docstring = ""

    base_url = f"{class_doc_base_url}{quote(module_name)}.{quote(class_name)}.html"
    method_fragment = f"{quote(module_name)}.{class_name}.{method_name}"

    return f"{base_url}#{method_fragment}", docstring


def _read_params(method, signature):
    method = html.escape(method)
    r = reprlib.Repr()
    r.maxlist = 2
    r.maxtuple = 1  # Show only first item of tuples
    r.maxstring = 20  # Limit string length
    cleaned_signature = html.escape(r.repr(signature).replace("'", ""))

    return method, cleaned_signature


def _methods_html_repr(methods):
    """Build HTML representation of estimator methods.
    Creates an HTML table with methods names, signature
    and link to its documentation. It is wrapped in a collapsible
    details element.
    """

    HTML_TEMPLATE = """
       <div class="estimator-table">
           <details>
               <summary>Methods</summary>
               <table class="parameters-table">
                 <tbody>
                   {rows}
                 </tbody>
               </table>
           </details>
       </div>
    """
    ROW_TEMPLATE = """
       <tr class="methods">
           <td class="param">{method_display}</td>
           <td class="method-signature">{signature}</td>
       </tr>
    """
    METHOD_AVAILABLE_DOC_LINK_TEMPLATE = """
        <a class="param-doc-link"
            style="anchor-name: --doc-link-{name};"
            rel="noreferrer" target="_blank" href="{link}">
            {name}
            <span class="param-doc-description"
            style="position-anchor: --doc-link-{name};">
            {method_display}</span>
        </a>
    """

    rows = []
    for row in methods:
        link, docstring = get_module_doc_link(methods.estimator_class, row)
        name, signature = _read_params(row, methods[row])

        if link and docstring:
            method_display = METHOD_AVAILABLE_DOC_LINK_TEMPLATE.format(
                link=link,
                name=name,
                method_display=docstring,
            )
        else:
            method_display = row

        rows.append(
            ROW_TEMPLATE.format(signature=signature, method_display=method_display)
        )
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
    doc_link : str, default=""
        The base URL to the online documentation for the estimator class.
        Used to generate parameter-specific documentation links in the HTML
        representation. If empty, documentation links will not be generated.
    """

    _html_repr = _methods_html_repr

    def __init__(self, *, methods=None, estimator_class=None):
        super().__init__(methods or {})
        self.estimator_class = estimator_class
