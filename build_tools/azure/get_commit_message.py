import argparse
import os
import subprocess
import warnings


def get_commit_message():
    """Retrieve the commit message."""
    build_source_version_message = os.environ.get("BUILD_SOURCEVERSIONMESSAGE")
    if build_source_version_message is None:
        # We are not on Azure: behaviour based on commit-message is not
        # supported for now.
        # TODO: this should be implemented at one point for GHA.
        warnings.warn(
            "get_commit_message not supported outside Azure for now, "
            "returning empty commit message"
        )
        return ""

    if os.environ["BUILD_REASON"] == "PullRequest":
        # By default pull requests use refs/pull/PULL_ID/merge as the source branch
        # which has a "Merge ID into ID" as a commit message. The latest commit
        # message is the second to last commit
        commit_id = build_source_version_message.split()[1]
        git_cmd = ["git", "log", commit_id, "-1", "--pretty=%B"]
        commit_message = subprocess.run(
            git_cmd, capture_output=True, text=True
        ).stdout.strip()
    else:
        commit_message = build_source_version_message

    # Sanitize the commit message to avoid introducing a vulnerability: a PR
    # submitter could include the "##vso" special marker in their commit
    # message to attempt to obfuscate the injection of arbitrary commands in
    # the Azure pipeline.
    #
    # This can be a problem if the PR reviewers do not pay close enough
    # attention to the full commit message prior to clicking the merge button
    # and as a result make the inject code run in a protected branch with
    # elevated access to CI secrets. On a protected branch, Azure
    # already sanitizes `BUILD_SOURCEVERSIONMESSAGE`, but the message
    # will still be sanitized here out of precaution.
    commit_message = commit_message.replace("##vso", "..vso")

    return commit_message


def parsed_args():
    parser = argparse.ArgumentParser(
        description=(
            "Show commit message that triggered the build in Azure DevOps pipeline"
        )
    )
    parser.add_argument(
        "--only-show-message",
        action="store_true",
        default=False,
        help=(
            "Only print commit message. Useful for direct use in scripts rather than"
            " setting output variable of the Azure job"
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parsed_args()
    commit_message = get_commit_message()

    if args.only_show_message:
        print(commit_message)
    else:
        # set the environment variable to be propagated to other steps
        print(f"##vso[task.setvariable variable=message;isOutput=true]{commit_message}")
        print(f"commit message: {commit_message}")  # helps debugging
