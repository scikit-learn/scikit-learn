import os
import subprocess


def get_commit_message():
    
    # By default pull requests use refs/pull/PULL_ID/merge as the source branch
    # which has a "Merge ID into ID" as a commit message. The latest commit
    # message is the second to last commit
    build_sourceversionmessage = os.environ["BUILD_SOURCEVERSIONMESSAGE"]
    commit_id = build_sourceversionmessage.split()[1]

    git_cmd = f"git log -1 --pretty=%B {commit_id}"
    commit_message = subprocess.run(
        git_cmd, shell=True, capture_output=True, text=True
    ).stdout.strip()

    print(commit_message)
    return commit_message
    # By default pull requests use refs/pull/PULL_ID/merge as the source branch
    # which has a "Merge ID into ID" as a commit message. The latest commit
    # message is the second to last commit
    # COMMIT_ID=$(echo $BUILD_SOURCEVERSIONMESSAGE | awk '{print $2}')
    # message=$(git log $COMMIT_ID -1 --pretty=%B)
    # message=$(echo $message)

    # echo "##vso[task.setvariable variable=message;isOutput=true]$message"


def get_selected_tests():
    """Parse the commit message to check if pytest should run only specific tests.

    If so, selected tests will be run with SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all".

    The commit message must take the form:
        <title> [all random seeds]
        <test_name_1>
        <test_name_2>
        ...
    """
    commit_message = get_commit_message()

    if "[all random seeds]" in commit_message:
        selected_tests = commit_message.split("[all random seeds]")[1].strip()
        selected_tests = selected_tests.replace("\n", " or ")
    else:
        selected_tests = ""

    # set the environment variable to be propagated to other steps
    print(f"##vso[task.setvariable variable=SELECTED_TESTS]'{selected_tests}'")


if __name__ == "__main__":
    get_selected_tests()
