import os

from get_commit_message import get_commit_message


def get_selected_tests():
    """Parse the commit message to check if pytest should run only specific tests.

    If so, selected tests will be run with SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all".

    The commit message must take the form:
        <title> [all random seeds]
        <test_name_1>
        <test_name_2>
        ...
    """
    if "SELECTED_TESTS" in os.environ or os.environ.get("GITHUB_ACTIONS", False):
        raise RuntimeError(
            f"It seems that the `SELECTED_TESTS` environment variable is set or you are"
            f"in a GitHub Actions context "
            f"(GITHUB_ACTIONS={os.environ.get('GITHUB_ACTIONS', False)})."
            f"Instead, please use directly `SELECTED_TESTS`."
        )

    commit_message = get_commit_message()

    if "[all random seeds]" in commit_message:
        selected_tests = commit_message.split("[all random seeds]")[1].strip()
        selected_tests = selected_tests.replace("\n", " or ")
    else:
        selected_tests = ""

    return selected_tests


if __name__ == "__main__":
    selected_tests = get_selected_tests()

    # set the environment variable to be propagated to other steps
    print(f"##vso[task.setvariable variable=SELECTED_TESTS]'{selected_tests}'")
    print(f"selected tests: {selected_tests}")  # helps debugging
