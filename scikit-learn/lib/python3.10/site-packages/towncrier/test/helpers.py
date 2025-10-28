from __future__ import annotations

import textwrap

from functools import wraps
from importlib import resources
from pathlib import Path
from subprocess import call
from typing import Any, Callable

from click.testing import CliRunner


def read(filename: str | Path) -> str:
    return Path(filename).read_text()


def write(path: str | Path, contents: str, dedent: bool = False) -> None:
    """
    Create a file with given contents including any missing parent directories
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if dedent:
        contents = textwrap.dedent(contents)
    p.write_text(contents)


def read_pkg_resource(path: str) -> str:
    """
    Read *path* from the towncrier package.
    """
    return (resources.files("towncrier") / path).read_text("utf-8")


def with_isolated_runner(fn: Callable[..., Any]) -> Callable[..., Any]:
    """
    Run *fn* within an isolated filesystem and add the kwarg *runner* to its
    arguments.
    """

    @wraps(fn)
    def test(*args: Any, **kw: Any) -> Any:
        runner = CliRunner()
        with runner.isolated_filesystem():
            return fn(*args, runner=runner, **kw)

    return test


def setup_simple_project(
    *,
    config: str | None = None,
    extra_config: str = "",
    pyproject_path: str = "pyproject.toml",
    mkdir_newsfragments: bool = True,
) -> None:
    if config is None:
        config = "[tool.towncrier]\n" 'package = "foo"\n' + extra_config
    else:
        config = textwrap.dedent(config)
    Path(pyproject_path).write_text(config)
    Path("foo").mkdir()
    Path("foo/__init__.py").write_text('__version__ = "1.2.3"\n')

    if mkdir_newsfragments:
        Path("foo/newsfragments").mkdir()


def with_project(
    *,
    config: str | None = None,
    pyproject_path: str = "pyproject.toml",
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Decorator to run a test with an isolated directory containing a simple
    project.

    The files are not managed by git.

    `config` is the content of the config file.
    It will be automatically dedented.

    `pyproject_path` is the path where to store the config file.
    """

    def decorator(fn: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(fn)
        def test(*args: Any, **kw: Any) -> Any:
            runner = CliRunner()
            with runner.isolated_filesystem():
                setup_simple_project(
                    config=config,
                    pyproject_path=pyproject_path,
                )

                return fn(*args, runner=runner, **kw)

        return test

    return decorator


def with_git_project(
    *,
    config: str | None = None,
    pyproject_path: str = "pyproject.toml",
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Decorator to run a test with an isolated directory containing a simple
    project checked into git.
    Use `config` to tweak the content of the config file.
    Use `pyproject_path` to tweak the location of the config file.
    """

    def decorator(fn: Callable[..., Any]) -> Callable[..., Any]:
        def _commit() -> None:
            call(["git", "add", "."])
            call(["git", "commit", "-m", "Second Commit"])

        @wraps(fn)
        def test(*args: Any, **kw: Any) -> Any:
            runner = CliRunner()
            with runner.isolated_filesystem():
                setup_simple_project(
                    config=config,
                    pyproject_path=pyproject_path,
                )

                call(["git", "init"])
                call(["git", "config", "user.name", "user"])
                call(["git", "config", "user.email", "user@example.com"])
                call(["git", "config", "commit.gpgSign", "false"])
                call(["git", "add", "."])
                call(["git", "commit", "-m", "Initial Commit"])

                return fn(*args, runner=runner, commit=_commit, **kw)

        return test

    return decorator
