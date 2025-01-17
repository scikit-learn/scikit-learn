"""The CLI for numpydoc."""

import argparse
import ast
from collections.abc import Sequence
from pathlib import Path
from typing import List, Union

from .docscrape_sphinx import get_doc_object
from .hooks import utils, validate_docstrings
from .validate import ERROR_MSGS, Validator, validate


def render_object(import_path: str, config: Union[List[str], None] = None) -> int:
    """Test numpydoc docstring generation for a given object."""
    # TODO: Move Validator._load_obj to a better place than validate
    print(get_doc_object(Validator._load_obj(import_path), config=dict(config or [])))
    return 0


def validate_object(import_path: str) -> int:
    """Run numpydoc docstring validation for a given object."""
    exit_status = 0
    results = validate(import_path)
    for err_code, err_desc in results["errors"]:
        exit_status += 1
        print(":".join([import_path, err_code, err_desc]))
    return exit_status


def get_parser() -> argparse.ArgumentParser:
    """
    Build an argument parser.

    Returns
    -------
    argparse.ArgumentParser
        The argument parser.
    """
    ap = argparse.ArgumentParser(prog="numpydoc", description=__doc__)
    subparsers = ap.add_subparsers(title="subcommands")

    def _parse_config(s):
        key, _, value = s.partition("=")
        value = ast.literal_eval(value)
        return key, value

    render = subparsers.add_parser(
        "render",
        description="Generate an expanded RST-version of the docstring.",
        help="generate the RST docstring with numpydoc",
    )
    render.add_argument("import_path", help="e.g. numpy.ndarray")
    render.add_argument(
        "-c",
        "--config",
        type=_parse_config,
        action="append",
        help="key=val where val will be parsed by literal_eval, "
        "e.g. -c use_plots=True. Multiple -c can be used.",
    )
    render.set_defaults(func=render_object)

    validate = subparsers.add_parser(
        "validate",
        description="Validate an object's docstring against the numpydoc standard.",
        help="validate the object's docstring and report errors",
    )
    validate.add_argument("import_path", help="e.g. numpy.ndarray")
    validate.set_defaults(func=validate_object)

    project_root_from_cwd, config_file = utils.find_project_root(["."])
    config_options = validate_docstrings.parse_config(project_root_from_cwd)
    ignored_checks = [
        f"- {check}: {ERROR_MSGS[check]}"
        for check in set(ERROR_MSGS.keys()) - config_options["checks"]
    ]
    ignored_checks_text = "\n  " + "\n  ".join(ignored_checks) + "\n"

    lint_parser = subparsers.add_parser(
        "lint",
        description="Run numpydoc validation on files with option to ignore individual checks.",
        help="validate all docstrings in file(s) using the abstract syntax tree",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    lint_parser.add_argument(
        "files", type=str, nargs="+", help="File(s) to run numpydoc validation on."
    )
    lint_parser.add_argument(
        "--config",
        type=str,
        help=(
            "Path to a directory containing a pyproject.toml or setup.cfg file.\n"
            "The hook will look for it in the root project directory.\n"
            "If both are present, only pyproject.toml will be used.\n"
            "Options must be placed under\n"
            "    - [tool:numpydoc_validation] for setup.cfg files and\n"
            "    - [tool.numpydoc_validation] for pyproject.toml files."
        ),
    )
    lint_parser.add_argument(
        "--ignore",
        type=str,
        nargs="*",
        help=(
            f"""Check codes to ignore.{
                ' Currently ignoring the following from '
                f'{Path(project_root_from_cwd) / config_file}: {ignored_checks_text}'
                'Values provided here will be in addition to the above, unless an alternate config is provided.'
                if ignored_checks else ''
            }"""
        ),
    )
    lint_parser.set_defaults(func=validate_docstrings.run_hook)

    return ap


def main(argv: Union[Sequence[str], None] = None) -> int:
    """CLI for numpydoc."""
    ap = get_parser()

    args = vars(ap.parse_args(argv))

    try:
        func = args.pop("func")
        return func(**args)
    except KeyError:
        ap.exit(status=2, message=ap.format_help())
