import os


def get_selected_tests():
    """Parse the commit message to check if pytest should run only specific tests.

    The commit message must take the form:
        <title> [all random seeds]
        <test_name_1>
        <test_name_2>
        ...
    """
    commit_message = os.environ["COMMIT_MESSAGE"]
    print(f"commit message: {commit_message}")

    if "[all random seeds]" in commit_message:
        selected_tests = commit_message.split("[all random seeds]")[1].strip()
        selected_tests = selected_tests.replace(" ", " or ")
    else:
        selected_tests = ""

    print(f"##vso[task.setvariable variable=SELECTED_TESTS]'{selected_tests}'")


if __name__ == "__main__":
    get_selected_tests()
