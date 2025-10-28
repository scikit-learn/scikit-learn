"""A small sphinx extension to add "copy" buttons to code blocks."""
from pathlib import Path
from sphinx.util import logging

__version__ = "0.5.2"

logger = logging.getLogger(__name__)


def scb_static_path(app):
    app.config.html_static_path.append(
        str(Path(__file__).parent.joinpath("_static").absolute())
    )


def add_to_context(app, config):
    # Update the global context
    config.html_context.update(
        {"copybutton_prompt_text": config.copybutton_prompt_text}
    )
    config.html_context.update(
        {"copybutton_prompt_is_regexp": config.copybutton_prompt_is_regexp}
    )
    config.html_context.update(
        {"copybutton_only_copy_prompt_lines": config.copybutton_only_copy_prompt_lines}
    )
    config.html_context.update(
        {"copybutton_remove_prompts": config.copybutton_remove_prompts}
    )
    config.html_context.update(
        {"copybutton_copy_empty_lines": config.copybutton_copy_empty_lines}
    )
    config.html_context.update(
        {
            "copybutton_line_continuation_character": (
                config.copybutton_line_continuation_character
            )
        }
    )
    config.html_context.update(
        {"copybutton_here_doc_delimiter": config.copybutton_here_doc_delimiter}
    )

    # Old image path deprecation
    # REMOVE after next release
    if config.copybutton_image_path:
        path = Path(app.srcdir) / config.copybutton_image_path
        logger.warning("copybutton_image_path is deprecated, use copybutton_image_svg")
        if not path.exists():
            raise ValueError("copybutton_img_path does not exist")
        if not path.suffix == ".svg":
            raise ValueError("copybutton_img_path must be an SVG")
        config.copybutton_image_svg = path.read_text()

    config.html_context.update({"copybutton_image_svg": config.copybutton_image_svg})
    config.html_context.update({"copybutton_selector": config.copybutton_selector})
    config.html_context.update(
        {
            "copybutton_format_func": Path(__file__)
            .parent.joinpath("_static", "copybutton_funcs.js")
            .read_text()
            .replace("export function", "function")
        }
    )
    config.html_context.update({"copybutton_exclude": config.copybutton_exclude})


def setup(app):
    logger.verbose("Adding copy buttons to code blocks...")
    # Add our static path
    app.connect("builder-inited", scb_static_path)

    # configuration for this tool
    app.add_config_value("copybutton_prompt_text", "", "html")
    app.add_config_value("copybutton_prompt_is_regexp", False, "html")
    app.add_config_value("copybutton_only_copy_prompt_lines", True, "html")
    app.add_config_value("copybutton_remove_prompts", True, "html")
    app.add_config_value("copybutton_copy_empty_lines", True, "html")
    app.add_config_value("copybutton_line_continuation_character", "", "html")
    app.add_config_value("copybutton_here_doc_delimiter", "", "html")
    app.add_config_value("copybutton_image_svg", "", "html")
    app.add_config_value("copybutton_selector", "div.highlight pre", "html")
    app.add_config_value("copybutton_exclude", ".linenos", "html")

    # DEPRECATE THIS AFTER THE NEXT RELEASE
    app.add_config_value("copybutton_image_path", "", "html")

    # Add configuration value to the template
    app.connect("config-inited", add_to_context)

    # Add relevant code to headers
    app.add_css_file("copybutton.css")
    app.add_js_file("clipboard.min.js")
    app.add_js_file("copybutton.js")
    return {
        "version": __version__,
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
