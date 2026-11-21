"""Run Vale on the given files and report only alerts on changed lines.

Used by the pre-commit hook. Vale checks whole files, so running it on every
file a commit touches would also report problems in text the commit did not
change. This script keeps only the alerts on lines that are part of the change:

- the lines changed between ``PRE_COMMIT_FROM_REF`` and ``PRE_COMMIT_TO_REF``
  when pre-commit sets them (``pre-commit run --from-ref ... --to-ref ...``);
- otherwise the staged lines (a normal ``git commit``);
- otherwise the lines changed since the branch left ``main``, which is what
  ``pre-commit run --all-files`` and ``pre-commit run --files`` report.

Pass ``--all-lines`` to report every alert. The exit code is 1 when an
error-level alert remains, so warnings and suggestions never block a commit.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
HUNK = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", re.M)


def git(*args):
    return subprocess.run(
        ["git", *args], capture_output=True, text=True, check=False
    ).stdout


def diff_base():
    for ref in ("upstream/main", "origin/main", "main"):
        base = git("merge-base", "HEAD", ref).strip()
        if base:
            return [base]
    return ["HEAD"]


def changed_lines(files):
    """Map each file to the set of line numbers it adds or modifies."""
    from_ref = os.environ.get("PRE_COMMIT_FROM_REF")
    to_ref = os.environ.get("PRE_COMMIT_TO_REF")
    if from_ref and to_ref:
        refs = [f"{from_ref}...{to_ref}"]
    elif git("diff", "--cached", "--name-only", "--", *files).strip():
        refs = ["--cached"]
    else:
        refs = diff_base()

    lines = defaultdict(set)
    out = git("diff", "-U0", "--no-color", "--no-ext-diff", *refs, "--", *files)
    for chunk in re.split(r"^diff --git ", out, flags=re.M)[1:]:
        match = re.search(r"^\+\+\+ b/(.*)$", chunk, re.M)
        if not match:
            continue
        for start, count in HUNK.findall(chunk):
            start, count = int(start), int(count or 1)
            lines[match.group(1)].update(range(start, start + count))
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--all-lines", action="store_true")
    parser.add_argument("files", nargs="*")
    args = parser.parse_args()
    if not args.files:
        return 0

    vale = shutil.which("vale")
    if vale is None:
        sys.exit("vale not found: install it with 'pip install vale'")
    # Vale looks up rst2html on the PATH; ours must win, and run with the
    # Python that has docutils installed.
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join(
        [str(HERE), str(Path(sys.executable).parent), env["PATH"]]
    )
    result = subprocess.run(
        [vale, "--output=JSON", *args.files], capture_output=True, text=True, env=env
    )
    try:
        alerts = json.loads(result.stdout or "{}")
    except json.JSONDecodeError:
        sys.exit(result.stdout + result.stderr)
    if "Code" in alerts:  # Vale reports runtime errors as a single object
        sys.exit(alerts.get("Text", result.stdout))

    keep = None if args.all_lines else changed_lines(args.files)
    errors = 0
    for path in sorted(alerts):
        for alert in alerts[path]:
            if keep is not None and alert["Line"] not in keep.get(path, ()):
                continue
            errors += alert["Severity"] == "error"
            col = alert["Span"][0]
            link = f" ({alert['Link']})" if alert.get("Link") else ""
            print(
                f"{path}:{alert['Line']}:{col}: {alert['Severity']}: "
                f"{alert['Message']} [{alert['Check']}]{link}"
            )
    return int(errors > 0)


if __name__ == "__main__":
    sys.exit(main())
