"""Gets the Pull Request Number from a workflow's run ID.

This script depends on
- `PyGithub` for interacting with GitHub
"""
import argparse

from github import Github

parser = argparse.ArgumentParser()
parser.add_argument("token", help="GitHub API token")
parser.add_argument("repo", help="Repo to query")
parser.add_argument("run_id", help="Workflow run ID", type=int)
args = parser.parse_args()

gh = Github(args.token)
sk_repo = gh.get_repo(args.repo)
run = sk_repo.get_workflow_run(args.run_id)

head = f"{run.head_repository.owner.login}:{run.head_branch}"
prs = list(sk_repo.get_pulls(state="all", sort="updated", direction="desc", head=head))
print(prs[0].number)
