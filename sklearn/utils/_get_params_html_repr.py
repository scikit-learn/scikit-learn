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

    out = ""
    html_start = """

        <div class="estimator-table">
            <details>
                <summary>Parameters</summary>
                <table>
                  <tbody>
        """

    for x, y in data.items():
        if y != "deprecated" and isinstance(y, str):
            modified_y = "".join(['"', y, '"'])
        else:
            modified_y = html.escape(str(y))

        if x in data.non_default:
            out += """
                    <tr class="user-set">
            """
        else:
            out += """
                    <tr class="default">
               """
        out += f"""
                        <td><i class="fa-regular fa-copy"
                         onclick="copyToClipboard('{x}',
                                  this.parentElement.nextElementSibling)"
                        </i></td>
                        <td class="param">{x}&nbsp;</td>
                        <td class="value">{modified_y}</td>
                    </tr>

                """
    html_end = """
        </tbody>
                </table>
            </details>
        </div>

    """
    return f"{html_start}{out}{html_end}"
