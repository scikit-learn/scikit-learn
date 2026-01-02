# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
import inspect
import re
from functools import lru_cache
from urllib.parse import quote

from sklearn.externals._numpydoc import docscrape


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


@lru_cache
def _scrape_estimator_docstring(docstring):
    return docscrape.NumpyDocString(docstring)


def _get_docstring(estimator_class, section_name, item):
    estimator_class_docs = inspect.getdoc(estimator_class)
    if estimator_class_docs and (
        structured_docstring := _scrape_estimator_docstring(estimator_class_docs)
    ):
        docstring_map = {
            item_docstring.name: item_docstring
            for item_docstring in structured_docstring[section_name]
        }
    else:
        docstring_map = {}
    if item_numpydoc := docstring_map.get(item, None):
        item_description = (
            f"{html.escape(item_numpydoc.name)}: "
            f"{html.escape(item_numpydoc.type)}<br><br>"
            f"{'<br>'.join(html.escape(line) for line in item_numpydoc.desc)}"
        )
    else:
        item_description = None
    return item_description
