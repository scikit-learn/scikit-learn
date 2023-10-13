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
    commit_message = get_commit_message()

    if "[all random seeds]" in commit_message:
        selected_tests = commit_message.split("[all random seeds]")[1].strip()
        selected_tests = selected_tests.replace("\n", " or ")
    else:
        selected_tests = ""

    return selected_tests


if __name__ == "__main__":
    # set the environment variable to be propagated to other steps
    selected_tests = get_selected_tests()

    if selected_tests:
        print(f"##vso[task.setvariable variable=SELECTED_TESTS]'{selected_tests}'")
        print(f"selected tests: {selected_tests}")  # helps debugging
    else:
        print("no selected tests")
