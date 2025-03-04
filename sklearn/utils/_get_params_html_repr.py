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

        if y != "deprecated" and isinstance(y, str):
            y = "".join(['"', y, '"'])

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
                        <td class="value">{y}</td>
                    </tr>

                """
    html_end = """
                  <tbody>
                </table>
            </details>
        </div>
        <script>
            function copyToClipboard(text, element) {
                // Get the parameter prefix from the closest toggleable content
                const toggleableContent = element.closest('.sk-toggleable__content');
                const paramPrefix = toggleableContent ?
                                    toggleableContent.dataset.paramPrefix : '';
                const fullParamName = paramPrefix ? `${paramPrefix}${text}` : text;

                const originalStyle = element.style;
                const computedStyle = window.getComputedStyle(element);
                const originalWidth = computedStyle.width;
                const originalHTML = element.innerHTML.replace('Copied!', '');

                navigator.clipboard.writeText(fullParamName)
                    .then(() => {
                        element.style.width = originalWidth;
                        element.style.color = 'green';
                        element.innerHTML = "Copied!";

                        setTimeout(() => {
                            element.innerHTML = originalHTML;
                            element.style = originalStyle;
                        }, 2000);
                    })
                    .catch(err => console.error('Failed to copy:', err));
                return false;
            }
        </script>
        </body>
    """
    html_template = f"{html_start}{out}{html_end}"
    try:
        output_path = "get_params.html"
        with open(output_path, "w") as f:
            f.write(html_template)
    except Exception as e:
        print(f"Error saving HTML: {e}")
    return html_template
