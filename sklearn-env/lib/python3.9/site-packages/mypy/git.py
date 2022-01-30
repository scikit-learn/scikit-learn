"""Git utilities."""

# Used also from setup.py, so don't pull in anything additional here (like mypy or typing):
import os
import subprocess


def is_git_repo(dir: str) -> bool:
    """Is the given directory version-controlled with git?"""
    return os.path.exists(os.path.join(dir, ".git"))


def have_git() -> bool:
    """Can we run the git executable?"""
    try:
        subprocess.check_output(["git", "--help"])
        return True
    except subprocess.CalledProcessError:
        return False
    except OSError:
        return False


def git_revision(dir: str) -> bytes:
    """Get the SHA-1 of the HEAD of a git repository."""
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=dir).strip()


def is_dirty(dir: str) -> bool:
    """Check whether a git repository has uncommitted changes."""
    output = subprocess.check_output(["git", "status", "-uno", "--porcelain"], cwd=dir)
    return output.strip() != b""


def has_extra_files(dir: str) -> bool:
    """Check whether a git repository has untracked files."""
    output = subprocess.check_output(["git", "clean", "--dry-run", "-d"], cwd=dir)
    return output.strip() != b""
