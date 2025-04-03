# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
from collections import UserDict

from sklearn._config import get_config


class ParamsDict(UserDict):
    def __init__(self, data=None, non_default=tuple()):
        super().__init__(data or {})
        self.non_default = non_default

    @property
    def _repr_html_(self):
        # Taken from sklearn.base.BaseEstimator - (written 5 yrs ago)
        """HTML representation of estimator.

        This is redundant with the logic of `_repr_mimebundle_`. The latter
        should be favored in the long term, `_repr_html_` is only
        implemented for consumers who do not interpret `_repr_mimbundle_`.
        """
        if get_config()["display"] != "diagram":
            raise AttributeError(
                "_repr_html_ is only defined when the "
                "'display' configuration option is set to "
                "'diagram'"
            )
        return self._repr_html_inner

    def _repr_html_inner(self):
        return _html_template(self)

    def _repr_mimebundle_(self, **kwargs):
        """Mime bundle used by jupyter kernels to display estimator"""
        output = {"text/plain": repr(self)}
        if get_config()["display"] == "diagram":
            output["text/html"] = _html_template(self)


def _html_template(data):

    HTML_TEMPLATE = """
        <div class="estimator-table">
            <details>
                <summary>Parameters</summary>
                <table>
                  <tbody>
                    {rows}
                  </tbody>
                </table>
            </details>
        </div>
    """

    ROW_TEMPLATE = """
        <tr class="{row_class}">
            <td><i class="fa-regular fa-copy"
                 onclick="copyToClipboard('{param_name}',
                          this.parentElement.nextElementSibling)"
            ></i></td>
            <td class="param">{param_name}&nbsp;</td>
            <td class="value">{param_value}</td>
        </tr>
    """

    rows = []
    for x, y in data.items():
        if y != "deprecated" and isinstance(y, str):
            modified_y = f'"{y}"'
        else:
            modified_y = html.escape(str(y))

        row_class = "user-set" if x in data.non_default else "default"
        rows.append(
            ROW_TEMPLATE.format(
                row_class=row_class, param_name=x, param_value=modified_y
            )
        )

    return HTML_TEMPLATE.format(rows="\n".join(rows))
