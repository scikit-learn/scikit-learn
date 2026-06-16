"""Synchronise the pixi minimum-dependency pins with ``_min_dependencies.py``.

The pixi configuration in ``pyproject.toml`` describes the CI environments. Two
of them (``min`` and ``doc-min``) pin a subset of the dependencies to the
minimum supported version. Those pins are the single source of truth in
``sklearn/_min_dependencies.py`` and are mirrored into ``pyproject.toml`` on the
lines that end with a ``# min`` marker, e.g.::

    numpy = "==1.24.1"  # min

This script rewrites the version of every ``# min`` line so that it matches
``sklearn/_min_dependencies.py``. Run it after bumping a minimum dependency::

    python build_tools/update_pixi_min_pins.py

Use ``--check`` to verify (without modifying the file) that the pins are in
sync, which is convenient as a CI guard::

    python build_tools/update_pixi_min_pins.py --check

Lines that pin a version but do *not* carry the ``# min`` marker (for instance
the ``doc-min`` ``pandas`` pin, which intentionally differs from the documented
minimum because conda-forge has no matching build) are left untouched.
"""

import argparse
import re
import sys
from pathlib import Path

# Import the single source of truth for the minimum supported versions.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sklearn._min_dependencies import dependent_packages

PYPROJECT_PATH = Path(__file__).resolve().parent.parent / "pyproject.toml"

# A line such as: ``numpy = "==1.24.1"  # min``. The trailing ``# min`` marker
# is what flags the pin as auto-generated from ``_min_dependencies.py``.
MIN_PIN_RE = re.compile(
    r"^(?P<prefix>\s*)"
    r"(?P<name>[A-Za-z0-9._-]+)"
    r'(?P<sep>\s*=\s*")=='
    r'(?P<version>[^"]+)'
    r'(?P<suffix>"\s*#\s*min\s*)$'
)


def _normalize(name):
    """Normalise a package name for case-insensitive comparison."""
    return name.lower().replace("_", "-")


def _min_versions():
    """Return a ``{normalized_name: min_version}`` mapping."""
    return {
        _normalize(name): min_version
        for name, (min_version, _extras) in dependent_packages.items()
    }


def update_min_pins(content, min_versions):
    """Return ``(new_content, changes)`` after rewriting ``# min`` pins.

    ``changes`` is a list of ``(name, old_version, new_version)`` tuples for the
    pins whose version changed.
    """
    changes = []
    new_lines = []
    for line in content.splitlines(keepends=True):
        match = MIN_PIN_RE.match(line.rstrip("\n"))
        if match is None:
            new_lines.append(line)
            continue

        name = match["name"]
        normalized = _normalize(name)
        if normalized not in min_versions:
            raise SystemExit(
                f"Package {name!r} is pinned with a '# min' marker in "
                f"pyproject.toml but is missing from "
                f"sklearn/_min_dependencies.py."
            )

        old_version = match["version"]
        new_version = min_versions[normalized]
        newline = "\n" if line.endswith("\n") else ""
        new_line = (
            f"{match['prefix']}{name}{match['sep']}=={new_version}"
            f"{match['suffix'].rstrip()}{newline}"
        )
        new_lines.append(new_line)
        if old_version != new_version:
            changes.append((name, old_version, new_version))

    return "".join(new_lines), changes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Only check that the pins are in sync; exit non-zero otherwise.",
    )
    args = parser.parse_args()

    content = PYPROJECT_PATH.read_text(encoding="utf-8")
    min_versions = _min_versions()
    new_content, changes = update_min_pins(content, min_versions)

    if args.check:
        if changes:
            print("pixi minimum-dependency pins are out of sync with")
            print("sklearn/_min_dependencies.py:")
            for name, old, new in changes:
                print(f"  - {name}: {old} -> {new}")
            print("\nRun `python build_tools/update_pixi_min_pins.py` to fix.")
            raise SystemExit(1)
        print("pixi minimum-dependency pins are in sync.")
        return

    if not changes:
        print("pixi minimum-dependency pins are already up to date.")
        return

    PYPROJECT_PATH.write_text(new_content, encoding="utf-8")
    print("Updated pixi minimum-dependency pins in pyproject.toml:")
    for name, old, new in changes:
        print(f"  - {name}: {old} -> {new}")
    print("\nRegenerate the lock file with `pixi lock` and commit the result.")


if __name__ == "__main__":
    main()
