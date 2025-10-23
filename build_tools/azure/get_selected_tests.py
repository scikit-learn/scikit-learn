import argparse
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
    if "SELECTED_TESTS" in os.environ:
        raise RuntimeError(
            "This legacy script should only be used on Azure. "
            "On GitHub actions, use the 'SELECTED_TESTS' environment variable"
        )

    commit_message = get_commit_message()

    if "[all random seeds]" in commit_message:
        selected_tests = commit_message.split("[all random seeds]")[1].strip()
        selected_tests = selected_tests.replace("\n", " or ")
    else:
        selected_tests = ""

    return selected_tests


def parsed_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--only-show-selected-tests",
        action="store_true",
        default=False,
        help=(
            "Only print selected tests. Useful for direct use in scripts rather than"
            " setting output variable of the Azure job"
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parsed_args()
    selected_tests = get_selected_tests()

    if args.only_show_selected_tests:
        print(selected_tests)
    else:
        # set the environment variable to be propagated to other steps
        print(f"##vso[task.setvariable variable=SELECTED_TESTS]'{selected_tests}'")
        print(f"selected tests: {selected_tests}")  # helps debugging
