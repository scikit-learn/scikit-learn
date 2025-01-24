# This script uses starlark for configuring when a cirrus CI job runs:
# https://cirrus-ci.org/guide/programming-tasks/

load("cirrus", "env", "fs", "http")

def main(ctx):
    # Only run for scikit-learn/scikit-learn. For debugging on a fork, you can
    # comment out the following condition.
    if env.get("CIRRUS_REPO_FULL_NAME") != "scikit-learn/scikit-learn":
        return []

    arm_tests_yaml = "build_tools/cirrus/arm_tests.yml"

    # Nightly jobs always run
    if env.get("CIRRUS_CRON", "") == "nightly":
        return fs.read(arm_tests_yaml)

    # Get commit message for event. We can not use `git` here because there is
    # no command line access in starlark. Thus we need to query the GitHub API
    # for the commit message. Note that `CIRRUS_CHANGE_MESSAGE` can not be used
    # because it is set to the PR's title and not the latest commit message.
    SHA = env.get("CIRRUS_CHANGE_IN_REPO")
    REPO = env.get("CIRRUS_REPO_FULL_NAME")
    url = "https://api.github.com/repos/" + REPO + "/git/commits/" + SHA
    response = http.get(url).json()
    commit_msg = response["message"]

    jobs_to_run = ""

    if "[cirrus arm]" in commit_msg:
        jobs_to_run += fs.read(arm_tests_yaml)

    return jobs_to_run
