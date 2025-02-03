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
        return html_template(self)

    def _repr_mimebundle_(self, **kwargs):
        """Mime bundle used by jupyter kernels to display estimator"""
        output = {"text/plain": repr(self)}
        if get_config()["display"] == "diagram":
            output["text/html"] = html_template(self)


def _get_css_style():
    return Path(__file__).with_suffix(".css").read_text(encoding="utf-8")


def html_template(data):

    num_parameters = len(data)
    style_template = _get_css_style()
    html_start = f"""
        <head><style>{style_template}</style>
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
        <script>
            function copyToClipboard(text) {{
                navigator.clipboard.writeText(text)
                    .then(() => console.log('Copied!'))
                    .catch(err => console.error('Failed to copy:', err));
                return false;
            }}
        </script>
        </head>
        <body>
        <div class="estimator-params">
             <ul>
            <details>
                <summary>Estimator x</summary>
                <ul>
        """
    out = ""
    for x, y in data.items():
        out += f"""
                      <li>{x}:{y}
                        <i class="fa-regular fa-copy"
                       onclick="copyToClipboard('{x}')"
                       style="color: #B197FC; cursor: pointer;">
                      </i>
                      </li>
           """
    html_end = """
            </details>
            </ul>
        </div>
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
