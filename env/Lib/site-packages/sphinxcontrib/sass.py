from os import PathLike
from pathlib import Path
from typing import Dict, Optional, Union

import sass
from sphinx.application import Sphinx
from sphinx.util import logging


logger = logging.getLogger(__name__)


Targets = Dict[PathLike, PathLike]


def configure_path(conf_dir: str, src: Optional[Union[PathLike, Path]]) -> Path:
    if src is None:
        target = Path(conf_dir)
    elif isinstance(src, str):
        target = Path(src)
    else:
        target = src
    if not target.is_absolute():
        target = Path(conf_dir) / target
    return target


def build_sass_sources(app: Sphinx, env):
    logger.debug("Building stylesheet files")
    src_dir = configure_path(app.confdir, app.config.sass_src_dir)
    out_dir = configure_path(app.confdir, app.config.sass_out_dir)
    include_paths = [str(p) for p in app.config.sass_include_paths]
    targets: Targets = app.config.sass_targets
    output_style = app.config.sass_output_style
    # Create output directory
    out_dir.mkdir(exist_ok=True, parents=True)
    # Build css files
    for src, dst in targets.items():
        src_ = src_dir / src
        content = src_.read_text()
        css = sass.compile(
            string=content,
            output_style=output_style,
            include_paths=[str(src_.parent)] + include_paths,
        )
        out_path = out_dir / dst
        out_path.parent.mkdir(exist_ok=True, parents=True)
        (out_dir / dst).write_text(css)


def setup(app: Sphinx):
    """
    Setup function for this extension.
    """
    logger.debug(f"Using {__name__}")
    app.add_config_value("sass_include_paths", [], "html")
    app.add_config_value("sass_src_dir", None, "html")
    app.add_config_value("sass_out_dir", None, "html")
    app.add_config_value("sass_targets", {}, "html")
    app.add_config_value("sass_output_style", "nested", "html")
    app.connect("env-updated", build_sass_sources)
