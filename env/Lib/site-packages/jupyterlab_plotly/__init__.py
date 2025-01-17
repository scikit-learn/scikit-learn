def _jupyter_labextension_paths():
    return [{"src": "labextension", "dest": "jupyterlab-plotly"}]


def _jupyter_nbextension_paths():
    return [
        {
            "section": "notebook",
            "src": "nbextension",
            "dest": "jupyterlab-plotly",
            "require": "jupyterlab-plotly/extension",
        }
    ]
