"""
GitHub helper for scikit-learn contributions.

Creates issues and pull requests via the GitHub REST API.
Reads the fork/upstream remotes from git config automatically.

Usage
-----
export GITHUB_TOKEN=ghp_your_token_here

# Create an issue
python build_tools/github_helper.py issue \
    --title "PERF: avoid dia_array in _rescale_data" \
    --body-file ISSUE_perf_rescale_data.md \
    --labels "Enhancement,Needs Triage,module:linear_model"

# Create a pull request
python build_tools/github_helper.py pr \
    --title "PERF avoid n_samples x n_samples dia_array in _rescale_data" \
    --body-file PR_body.md \
    --base main \
    --head perf/avoid-dia-array-in-rescale-data \
    --issue 12345
"""

import argparse
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request


# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------

def _git(*args):
    result = subprocess.run(
        ["git"] + list(args),
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def _parse_github_remote(remote_url):
    """Extract 'owner/repo' from an https or ssh remote URL."""
    patterns = [
        r"github\.com[:/]([^/]+/[^/]+?)(?:\.git)?$",
    ]
    for pattern in patterns:
        match = re.search(pattern, remote_url)
        if match:
            return match.group(1)
    return None


def _detect_repos():
    """Return (fork_owner_repo, upstream_owner_repo) from git remotes."""
    remotes = {}
    for line in _git("remote", "-v").splitlines():
        parts = line.split()
        if len(parts) >= 2 and "(push)" in line:
            remotes[parts[0]] = _parse_github_remote(parts[1])

    upstream = remotes.get("upstream")
    origin = remotes.get("origin")

    if upstream is None:
        print(
            "WARNING: no 'upstream' remote found. "
            "PRs will target 'origin' as both base and head repo."
        )
        upstream = origin

    return origin, upstream


def _current_branch():
    return _git("rev-parse", "--abbrev-ref", "HEAD")


# ---------------------------------------------------------------------------
# GitHub API
# ---------------------------------------------------------------------------

def _api_request(method, path, token, payload=None):
    url = f"https://api.github.com{path}"
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "Content-Type": "application/json",
    }
    data = json.dumps(payload).encode() if payload else None
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as exc:
        body = exc.read().decode()
        print(f"GitHub API error {exc.code}: {body}", file=sys.stderr)
        sys.exit(1)


def create_issue(token, upstream_repo, title, body, labels):
    labels_list = [l.strip() for l in labels.split(",") if l.strip()]
    payload = {"title": title, "body": body, "labels": labels_list}
    result = _api_request("POST", f"/repos/{upstream_repo}/issues", token, payload)
    return result["number"], result["html_url"]


def create_pr(token, upstream_repo, fork_repo, title, body, base, head, issue_number):
    fork_owner = fork_repo.split("/")[0]
    full_head = f"{fork_owner}:{head}"
    payload = {
        "title": title,
        "body": body,
        "head": full_head,
        "base": base,
        "draft": False,
    }
    result = _api_request("POST", f"/repos/{upstream_repo}/pulls", token, payload)
    pr_number = result["number"]
    pr_url = result["html_url"]

    # Link issue to PR if provided
    if issue_number:
        _api_request(
            "POST",
            f"/repos/{upstream_repo}/issues/{pr_number}/comments",
            token,
            {"body": f"Closes #{issue_number}"},
        )

    return pr_number, pr_url


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _read_body(args):
    if args.body_file:
        with open(args.body_file, encoding="utf-8") as f:
            return f.read()
    if args.body:
        return args.body
    print("ERROR: provide --body or --body-file", file=sys.stderr)
    sys.exit(1)


def cmd_issue(args, token, origin, upstream):
    body = _read_body(args)
    number, url = create_issue(
        token,
        upstream,
        args.title,
        body,
        args.labels or "",
    )
    print(f"Issue created: #{number}")
    print(f"URL: {url}")
    print(f"\nNext step: open PR and use 'Closes #{number}' in the PR body.")
    return number


def cmd_pr(args, token, origin, upstream):
    body = _read_body(args)
    head = args.head or _current_branch()
    base = args.base or "main"

    pr_number, pr_url = create_pr(
        token,
        upstream,
        origin,
        args.title,
        body,
        base,
        head,
        args.issue,
    )
    print(f"PR created: #{pr_number}")
    print(f"URL: {pr_url}")

    if args.issue:
        print(f"Linked to issue #{args.issue}")

    print(
        f"\nNext step: rename changelog file:\n"
        f"  git mv doc/whats_new/upcoming_changes/**/PENDING.*.rst \\\n"
        f"         doc/whats_new/upcoming_changes/**/{pr_number}.*.rst\n"
        f"  git commit -m 'MNT update changelog filename with PR number (#{pr_number})'\n"
        f"  git push"
    )
    return pr_number


def main():
    token = os.environ.get("GITHUB_TOKEN")
    if not token:
        print(
            "ERROR: set GITHUB_TOKEN environment variable.\n"
            "  export GITHUB_TOKEN=ghp_your_token_here",
            file=sys.stderr,
        )
        sys.exit(1)

    origin, upstream = _detect_repos()

    parser = argparse.ArgumentParser(
        description="scikit-learn GitHub contribution helper"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- issue ---
    p_issue = sub.add_parser("issue", help="Create a GitHub issue")
    p_issue.add_argument("--title", required=True, help="Issue title")
    p_issue.add_argument("--body", help="Issue body (inline string)")
    p_issue.add_argument("--body-file", help="Path to a markdown file for the body")
    p_issue.add_argument(
        "--labels",
        default="Enhancement,Needs Triage",
        help="Comma-separated label names (default: 'Enhancement,Needs Triage')",
    )

    # --- pr ---
    p_pr = sub.add_parser("pr", help="Create a GitHub pull request")
    p_pr.add_argument("--title", required=True, help="PR title")
    p_pr.add_argument("--body", help="PR body (inline string)")
    p_pr.add_argument("--body-file", help="Path to a markdown file for the body")
    p_pr.add_argument(
        "--base", default="main", help="Target branch on upstream (default: main)"
    )
    p_pr.add_argument(
        "--head",
        help="Source branch on fork (default: current branch)",
    )
    p_pr.add_argument(
        "--issue",
        type=int,
        help="Issue number this PR closes (optional)",
    )

    args = parser.parse_args()

    if args.command == "issue":
        cmd_issue(args, token, origin, upstream)
    elif args.command == "pr":
        cmd_pr(args, token, origin, upstream)


if __name__ == "__main__":
    main()
