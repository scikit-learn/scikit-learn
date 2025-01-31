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
        # return f"<body><h1>Hello World aqui</h1></body>"
        return html_template(self)

    def _repr_mimebundle_(self, **kwargs):
        """Mime bundle used by jupyter kernels to display estimator"""
        output = {"text/plain": repr(self)}
        if get_config()["display"] == "diagram":
            # output["text/html"] = f"<body><h1>Hello World alla {self}</h1></body>"
            output["text/html"] = html_template(self)


def _get_css_style():
    return Path(__file__).with_suffix(".css").read_text(encoding="utf-8")


def html_template(data):
    style_template = _get_css_style()
    html_start = f"<head><style>{style_template}</style><body>"
    html_end = "</ul></body>"
    out = ""
    for x, y in data.items():
        out += f"""
                  <ul>
                      <li>{x}:{y}</li>
                  </ul>
           """
    return html_start + out + html_end
