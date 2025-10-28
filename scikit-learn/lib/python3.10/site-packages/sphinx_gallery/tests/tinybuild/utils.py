"""Utility functions for doc building."""

import os.path as op

from sphinx_gallery.scrapers import matplotlib_scraper


def notebook_modification_function(notebook_content, notebook_filename):
    """Implement JupyterLite-specific modifications of notebooks."""
    source = f"JupyterLite-specific change for {notebook_filename}"
    markdown_cell = {"cell_type": "markdown", "metadata": {}, "source": source}
    notebook_content["cells"] = [markdown_cell] + notebook_content["cells"]


class MatplotlibFormatScraper:
    """Calls Matplotlib scraper, passing required `format` kwarg for testing."""

    def __repr__(self):
        return self.__class__.__name__

    def __call__(self, block, block_vars, gallery_conf):
        """Call Matplotlib scraper with required `format` kwarg for testing."""
        kwargs = dict()
        if (
            op.basename(block_vars["target_file"]) == "plot_svg.py"
            and gallery_conf["builder_name"] != "latex"
        ):
            kwargs["format"] = "svg"
        elif (
            op.basename(block_vars["target_file"]) == "plot_webp.py"
            and gallery_conf["builder_name"] != "latex"
        ):
            kwargs["format"] = "webp"
        return matplotlib_scraper(block, block_vars, gallery_conf, **kwargs)


class ResetArgv:
    """Provide `reset_argv` callable returning required `sys.argv` for test."""

    def __repr__(self):
        return "ResetArgv"

    def __call__(self, sphinx_gallery_conf, script_vars):
        """Return 'plot' arg if 'plot_command_line_args' example, for testing."""
        if "plot_command_line_args.py" in script_vars["src_file"]:
            return ["plot"]
        else:
            return []


def _raise(*args, **kwargs):
    import matplotlib.pyplot as plt

    plt.close("all")
    raise ValueError(
        "zero-size array to reduction operation minimum which has no identity"
    )


class MockScrapeProblem:
    """Used in 'reset_modules' to mock error during scraping."""

    def __init__(self):
        from matplotlib.colors import colorConverter

        self._orig = colorConverter.to_rgba

    def __repr__(self):
        return "MockScrapeProblem"

    def __call__(self, gallery_conf, fname):
        """Raise error for 'scraper_broken' example."""
        from matplotlib.colors import colorConverter

        if "scraper_broken" in fname:
            colorConverter.to_rgba = _raise
        else:
            colorConverter.to_rgba = self._orig


class MockSort:
    """Fake sort used to test that mini-gallery sort is correct."""

    def __repr__(self):
        return "MockSort"

    def __call__(self, f):
        """Sort plot_sub* for one test case."""
        if "subdir2" in f:
            return 0
        if "subdir1" in f:
            return 1
        return f


mock_scrape_problem = MockScrapeProblem()
matplotlib_format_scraper = MatplotlibFormatScraper()
reset_argv = ResetArgv()
mock_sort = MockSort()


def noop_key(x):
    """Sortkey that passes x."""
    return x
