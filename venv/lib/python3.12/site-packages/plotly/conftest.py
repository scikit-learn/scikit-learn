import os


def pytest_ignore_collect(path):
    # Ignored files, most of them are raising a chart studio error
    ignored_paths = [
        "exploding_module.py",
        "chunked_requests.py",
        "v2.py",
        "v1.py",
        "presentation_objs.py",
        "widgets.py",
        "dashboard_objs.py",
        "grid_objs.py",
        "config.py",
        "presentation_objs.py",
        "session.py",
    ]
    if (
        os.path.basename(str(path)) in ignored_paths
        or "plotly/plotly/plotly/__init__.py" in str(path)
        or "plotly/api/utils.py" in str(path)
    ):
        return True
