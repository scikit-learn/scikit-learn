import os


def get_selected_tests():
    """Parse the commit message to check if pytest should run only specific tests.

    If so, selected tests will be run with SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all".

    The commit message must take the form:
        <title> [all random seeds]
        <test_name_1>
        <test_name_2>
        ...
    """
    commit_message = os.environ["COMMIT_MESSAGE"]

    if "[all random seeds]" in commit_message:
        selected_tests = commit_message.split("[all random seeds]")[1].strip()
        selected_tests = selected_tests.replace(" ", " or ")
    else:
        selected_tests = ""

    # set the environment variable to be propagated to other steps
    print("##vso[task.setvariable variable=SELECTED_TESTS]'{}'".format(selected_tests))


if __name__ == "__main__":
    get_selected_tests()
