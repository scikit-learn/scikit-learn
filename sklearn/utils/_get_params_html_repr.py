# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import UserDict
from pathlib import Path

from sklearn._config import get_config


class ParamsDict(UserDict):
    @property
    def _repr_html_(self):
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


def _get_css_style():
    return Path(__file__).with_suffix(".css").read_text(encoding="utf-8")


def _html_template(data):

    style_template = _get_css_style()
    html_start = f"""
        <head><style>{style_template}</style>
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
        </head>
        <body>
        <div class="estimator-params">
            <details>
                <summary>Parameters</summary>
                <table>
                  <tbody>
        """
    out = ""
    for x, y in data.items():
        if x in data.non_default:
            out += f"""
                    <tr class="user-set">

                      <td><i class="fa-regular fa-copy"
                       onclick="copyToClipboard('{x}', this)"
                       style="color: gray; cursor: pointer;">
                      </i></td>
                      <td>{x}&nbsp;</td>
                      <td>{y}</td>
                    </tr>
            """
        else:
            out += f"""
                    <tr class="default">

                          <td><i class="fa-regular fa-copy"
                           onclick="copyToClipboard('{x}', this)"
                           style="color: gray; cursor: pointer; min-with: 5em">
                          </i></td>
                          <td>{x}&nbsp;</td>
                          <td>{y}</td>
                    </tr>
               """
    html_end = """
                  <tbody>
                </table>
            </details>

        </div>
        <script>
            function copyToClipboard(text, element) {{
                const parent = element.parentNode;
                const originalHTML = parent.innerHTML.replace('&nbsp;Copied!', '');

                navigator.clipboard.writeText(text)
                    .then(() => {
                        parent.innerHTML = originalHTML + "&nbsp;Copied!";
                        setTimeout(() => {
                            parent.innerHTML = originalHTML;
                        }, 1000);
                    })
                    .catch(err => console.error('Failed to copy:', err));
                return false;
            }}
        </script>
        </body>
    """
    html_template = f"{html_start}{out}{html_end}"
    # Remove the following:
    try:
        output_path = "get_params.html"
        with open(output_path, "w") as f:
            f.write(html_template)
    except Exception as e:
        print(f"Error saving HTML: {e}")
    return html_template
