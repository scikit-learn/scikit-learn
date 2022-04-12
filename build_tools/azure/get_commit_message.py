import os
import subprocess


def get_commit_message():
    """Retrieve the commit message."""
    build_source_version_message = os.environ["BUILD_SOURCEVERSIONMESSAGE"]

    if os.environ["BUILD_REASON"] == "PullRequest":
        # By default pull requests use refs/pull/PULL_ID/merge as the source branch
        # which has a "Merge ID into ID" as a commit message. The latest commit
        # message is the second to last commit
        commit_id = build_source_version_message.split()[1]
        git_cmd = "git log {} -1 --pretty=%B ".format(commit_id)
        commit_message = subprocess.run(
            git_cmd, shell=True, capture_output=True, text=True
        ).stdout.strip()
    else:
        commit_message = build_source_version_message

    return commit_message


if __name__ == "__main__":
    # set the environment variable to be propagated to other steps
    message = get_commit_message()
    print("##vso[task.setvariable variable=message;isOutput=true]{}".format(message))
