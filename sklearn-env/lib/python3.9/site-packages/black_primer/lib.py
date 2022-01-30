import asyncio
import errno
import json
import logging
import os
import stat
import sys
from functools import partial
from pathlib import Path
from platform import system
from shutil import rmtree, which
from subprocess import CalledProcessError
from sys import version_info
from tempfile import TemporaryDirectory
from typing import Any, Callable, Dict, NamedTuple, Optional, Sequence, Tuple
from urllib.parse import urlparse

import click


WINDOWS = system() == "Windows"
BLACK_BINARY = "black.exe" if WINDOWS else "black"
GIT_BINARY = "git.exe" if WINDOWS else "git"
LOG = logging.getLogger(__name__)


# Windows needs a ProactorEventLoop if you want to exec subprocesses
# Starting with 3.8 this is the default - can remove when Black >= 3.8
# mypy only respects sys.platform if directly in the evaluation
# https://mypy.readthedocs.io/en/latest/common_issues.html#python-version-and-system-platform-checks  # noqa: B950
if sys.platform == "win32":
    asyncio.set_event_loop(asyncio.ProactorEventLoop())


class Results(NamedTuple):
    stats: Dict[str, int] = {}
    failed_projects: Dict[str, CalledProcessError] = {}


async def _gen_check_output(
    cmd: Sequence[str],
    timeout: float = 600,
    env: Optional[Dict[str, str]] = None,
    cwd: Optional[Path] = None,
) -> Tuple[bytes, bytes]:
    process = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
        env=env,
        cwd=cwd,
    )
    try:
        (stdout, stderr) = await asyncio.wait_for(process.communicate(), timeout)
    except asyncio.TimeoutError:
        process.kill()
        await process.wait()
        raise

    # A non-optional timeout was supplied to asyncio.wait_for, guaranteeing
    # a timeout or completed process.  A terminated Python process will have a
    # non-empty returncode value.
    assert process.returncode is not None

    if process.returncode != 0:
        cmd_str = " ".join(cmd)
        raise CalledProcessError(
            process.returncode, cmd_str, output=stdout, stderr=stderr
        )

    return (stdout, stderr)


def analyze_results(project_count: int, results: Results) -> int:
    failed_pct = round(((results.stats["failed"] / project_count) * 100), 2)
    success_pct = round(((results.stats["success"] / project_count) * 100), 2)

    click.secho("-- primer results 📊 --\n", bold=True)
    click.secho(
        f"{results.stats['success']} / {project_count} succeeded ({success_pct}%) ✅",
        bold=True,
        fg="green",
    )
    click.secho(
        f"{results.stats['failed']} / {project_count} FAILED ({failed_pct}%) 💩",
        bold=bool(results.stats["failed"]),
        fg="red",
    )
    s = "" if results.stats["disabled"] == 1 else "s"
    click.echo(f" - {results.stats['disabled']} project{s} disabled by config")
    s = "" if results.stats["wrong_py_ver"] == 1 else "s"
    click.echo(
        f" - {results.stats['wrong_py_ver']} project{s} skipped due to Python version"
    )
    click.echo(
        f" - {results.stats['skipped_long_checkout']} skipped due to long checkout"
    )

    if results.failed_projects:
        click.secho("\nFailed projects:\n", bold=True)

    for project_name, project_cpe in results.failed_projects.items():
        print(f"## {project_name}:")
        print(f" - Returned {project_cpe.returncode}")
        if project_cpe.stderr:
            print(f" - stderr:\n{project_cpe.stderr.decode('utf8')}")
        if project_cpe.stdout:
            print(f" - stdout:\n{project_cpe.stdout.decode('utf8')}")
        print("")

    return results.stats["failed"]


async def black_run(
    repo_path: Path,
    project_config: Dict[str, Any],
    results: Results,
    no_diff: bool = False,
) -> None:
    """Run Black and record failures"""
    cmd = [str(which(BLACK_BINARY))]
    if "cli_arguments" in project_config and project_config["cli_arguments"]:
        cmd.extend(project_config["cli_arguments"])
    cmd.append("--check")
    if no_diff:
        cmd.append(".")
    else:
        cmd.extend(["--diff", "."])

    with TemporaryDirectory() as tmp_path:
        # Prevent reading top-level user configs by manipulating environment variables
        env = {
            **os.environ,
            "XDG_CONFIG_HOME": tmp_path,  # Unix-like
            "USERPROFILE": tmp_path,  # Windows (changes `Path.home()` output)
        }

        try:
            _stdout, _stderr = await _gen_check_output(cmd, cwd=repo_path, env=env)
        except asyncio.TimeoutError:
            results.stats["failed"] += 1
            LOG.error(f"Running black for {repo_path} timed out ({cmd})")
        except CalledProcessError as cpe:
            # TODO: Tune for smarter for higher signal
            # If any other return value than 1 we raise - can disable project in config
            if cpe.returncode == 1:
                if not project_config["expect_formatting_changes"]:
                    results.stats["failed"] += 1
                    results.failed_projects[repo_path.name] = cpe
                else:
                    results.stats["success"] += 1
                return
            elif cpe.returncode > 1:
                results.stats["failed"] += 1
                results.failed_projects[repo_path.name] = cpe
                return

            LOG.error(f"Unknown error with {repo_path}")
            raise

    # If we get here and expect formatting changes something is up
    if project_config["expect_formatting_changes"]:
        results.stats["failed"] += 1
        results.failed_projects[repo_path.name] = CalledProcessError(
            0, cmd, b"Expected formatting changes but didn't get any!", b""
        )
        return

    results.stats["success"] += 1


async def git_checkout_or_rebase(
    work_path: Path,
    project_config: Dict[str, Any],
    rebase: bool = False,
    *,
    depth: int = 1,
) -> Optional[Path]:
    """git Clone project or rebase"""
    git_bin = str(which(GIT_BINARY))
    if not git_bin:
        LOG.error("No git binary found")
        return None

    repo_url_parts = urlparse(project_config["git_clone_url"])
    path_parts = repo_url_parts.path[1:].split("/", maxsplit=1)

    repo_path: Path = work_path / path_parts[1].replace(".git", "")
    cmd = [git_bin, "clone", "--depth", str(depth), project_config["git_clone_url"]]
    cwd = work_path
    if repo_path.exists() and rebase:
        cmd = [git_bin, "pull", "--rebase"]
        cwd = repo_path
    elif repo_path.exists():
        return repo_path

    try:
        _stdout, _stderr = await _gen_check_output(cmd, cwd=cwd)
    except (asyncio.TimeoutError, CalledProcessError) as e:
        LOG.error(f"Unable to git clone / pull {project_config['git_clone_url']}: {e}")
        return None

    return repo_path


def handle_PermissionError(
    func: Callable, path: Path, exc: Tuple[Any, Any, Any]
) -> None:
    """
    Handle PermissionError during shutil.rmtree.

    This checks if the erroring function is either 'os.rmdir' or 'os.unlink', and that
    the error was EACCES (i.e. Permission denied). If true, the path is set writable,
    readable, and executable by everyone. Finally, it tries the error causing delete
    operation again.

    If the check is false, then the original error will be reraised as this function
    can't handle it.
    """
    excvalue = exc[1]
    LOG.debug(f"Handling {excvalue} from {func.__name__}... ")
    if func in (os.rmdir, os.unlink) and excvalue.errno == errno.EACCES:
        LOG.debug(f"Setting {path} writable, readable, and executable by everyone... ")
        os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # chmod 0777
        func(path)  # Try the error causing delete operation again
    else:
        raise


async def load_projects_queue(
    config_path: Path,
) -> Tuple[Dict[str, Any], asyncio.Queue]:
    """Load project config and fill queue with all the project names"""
    with config_path.open("r") as cfp:
        config = json.load(cfp)

    # TODO: Offer more options here
    # e.g. Run on X random packages or specific sub list etc.
    project_names = sorted(config["projects"].keys())
    queue: asyncio.Queue = asyncio.Queue(maxsize=len(project_names))
    for project in project_names:
        await queue.put(project)

    return config, queue


async def project_runner(
    idx: int,
    config: Dict[str, Any],
    queue: asyncio.Queue,
    work_path: Path,
    results: Results,
    long_checkouts: bool = False,
    rebase: bool = False,
    keep: bool = False,
    no_diff: bool = False,
) -> None:
    """Check out project and run Black on it + record result"""
    loop = asyncio.get_event_loop()
    py_version = f"{version_info[0]}.{version_info[1]}"
    while True:
        try:
            project_name = queue.get_nowait()
        except asyncio.QueueEmpty:
            LOG.debug(f"project_runner {idx} exiting")
            return
        LOG.debug(f"worker {idx} working on {project_name}")

        project_config = config["projects"][project_name]

        # Check if disabled by config
        if "disabled" in project_config and project_config["disabled"]:
            results.stats["disabled"] += 1
            LOG.info(f"Skipping {project_name} as it's disabled via config")
            continue

        # Check if we should run on this version of Python
        if (
            "all" not in project_config["py_versions"]
            and py_version not in project_config["py_versions"]
        ):
            results.stats["wrong_py_ver"] += 1
            LOG.debug(f"Skipping {project_name} as it's not enabled for {py_version}")
            continue

        # Check if we're doing big projects / long checkouts
        if not long_checkouts and project_config["long_checkout"]:
            results.stats["skipped_long_checkout"] += 1
            LOG.debug(f"Skipping {project_name} as it's configured as a long checkout")
            continue

        repo_path = await git_checkout_or_rebase(work_path, project_config, rebase)
        if not repo_path:
            continue
        await black_run(repo_path, project_config, results, no_diff)

        if not keep:
            LOG.debug(f"Removing {repo_path}")
            rmtree_partial = partial(
                rmtree, path=repo_path, onerror=handle_PermissionError
            )
            await loop.run_in_executor(None, rmtree_partial)

        LOG.info(f"Finished {project_name}")


async def process_queue(
    config_file: str,
    work_path: Path,
    workers: int,
    keep: bool = False,
    long_checkouts: bool = False,
    rebase: bool = False,
    no_diff: bool = False,
) -> int:
    """
    Process the queue with X workers and evaluate results
    - Success is guaged via the config "expect_formatting_changes"

    Integer return equals the number of failed projects
    """
    results = Results()
    results.stats["disabled"] = 0
    results.stats["failed"] = 0
    results.stats["skipped_long_checkout"] = 0
    results.stats["success"] = 0
    results.stats["wrong_py_ver"] = 0

    config, queue = await load_projects_queue(Path(config_file))
    project_count = queue.qsize()
    s = "" if project_count == 1 else "s"
    LOG.info(f"{project_count} project{s} to run Black over")
    if project_count < 1:
        return -1

    s = "" if workers == 1 else "s"
    LOG.debug(f"Using {workers} parallel worker{s} to run Black")
    # Wait until we finish running all the projects before analyzing
    await asyncio.gather(
        *[
            project_runner(
                i,
                config,
                queue,
                work_path,
                results,
                long_checkouts,
                rebase,
                keep,
                no_diff,
            )
            for i in range(workers)
        ]
    )

    LOG.info("Analyzing results")
    return analyze_results(project_count, results)


if __name__ == "__main__":  # pragma: nocover
    raise NotImplementedError("lib is a library, funnily enough.")
