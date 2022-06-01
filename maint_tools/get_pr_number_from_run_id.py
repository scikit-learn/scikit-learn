import argparse
from github import Github

parser = argparse.ArgumentParser()
parser.add_argument("repo")
parser.add_argument("run_id")

args = parser.parse_args()

gh = Github()

sk_repo = gh.get_repo(args.repo)
run = sk_repo.get_workflow_run(args.run_id)

head = f"{run.head_repository.owner.login}:{run.head_branch}"
prs = list(sk_repo.get_pulls(state="all", head=head))

pr = prs[0]

print(pr.number)
