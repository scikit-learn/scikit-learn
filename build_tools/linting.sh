#!/bin/bash

# Note that any change in this file, adding or removing steps or changing the
# printed messages, should be also reflected in the `get_comment.py` file.

# This script shouldn't exit if a command / pipeline fails
set +e
# pipefail is necessary to propagate exit codes
set -o pipefail

global_status=0

echo -e "### Running black ###\n"
black --check --diff .
status=$?

if [[ $status -eq 0 ]]
then
    echo -e "No problem detected by black\n"
else
    echo -e "Problems detected by black, please run black and commit the result\n"
    global_status=1
fi

echo -e "### Running ruff ###\n"
ruff check --output-format=full .
status=$?
if [[ $status -eq 0 ]]
then
    echo -e "No problem detected by ruff\n"
else
    echo -e "Problems detected by ruff, please fix them\n"
    global_status=1
fi

echo -e "### Running mypy ###\n"
mypy sklearn/
status=$?
if [[ $status -eq 0 ]]
then
    echo -e "No problem detected by mypy\n"
else
    echo -e "Problems detected by mypy, please fix them\n"
    global_status=1
fi

echo -e "### Running cython-lint ###\n"
cython-lint sklearn/
status=$?
if [[ $status -eq 0 ]]
then
    echo -e "No problem detected by cython-lint\n"
else
    echo -e "Problems detected by cython-lint, please fix them\n"
    global_status=1
fi

# For docstrings and warnings of deprecated attributes to be rendered
# properly, the `deprecated` decorator must come before the `property` decorator
# (else they are treated as functions)

echo -e "### Checking for bad deprecation order ###\n"
bad_deprecation_property_order=`git grep -A 10 "@property"  -- "*.py" | awk '/@property/,/def /' | grep -B1 "@deprecated"`

if [ ! -z "$bad_deprecation_property_order" ]
then
    echo "deprecated decorator should come before property decorator"
    echo "found the following occurrences:"
    echo $bad_deprecation_property_order
    echo -e "\nProblems detected by deprecation order check\n"
    global_status=1
else
    echo -e "No problems detected related to deprecation order\n"
fi

# Check for default doctest directives ELLIPSIS and NORMALIZE_WHITESPACE

echo -e "### Checking for default doctest directives ###\n"
doctest_directive="$(git grep -nw -E "# doctest\: \+(ELLIPSIS|NORMALIZE_WHITESPACE)")"

if [ ! -z "$doctest_directive" ]
then
    echo "ELLIPSIS and NORMALIZE_WHITESPACE doctest directives are enabled by default, but were found in:"
    echo "$doctest_directive"
    echo -e "\nProblems detected by doctest directive check\n"
    global_status=1
else
    echo -e "No problems detected related to doctest directives\n"
fi

# Check for joblib.delayed and joblib.Parallel imports
# TODO(1.7): remove ":!sklearn/utils/_joblib.py"
echo -e "### Checking for joblib imports ###\n"
joblib_status=0
joblib_delayed_import="$(git grep -l -A 10 -E "joblib import.+delayed" -- "*.py" ":!sklearn/utils/_joblib.py" ":!sklearn/utils/parallel.py")"
if [ ! -z "$joblib_delayed_import" ]; then
    echo "Use from sklearn.utils.parallel import delayed instead of joblib delayed. The following files contains imports to joblib.delayed:"
    echo "$joblib_delayed_import"
    joblib_status=1
fi
joblib_Parallel_import="$(git grep -l -A 10 -E "joblib import.+Parallel" -- "*.py" ":!sklearn/utils/_joblib.py" ":!sklearn/utils/parallel.py")"
if [ ! -z "$joblib_Parallel_import" ]; then
    echo "Use from sklearn.utils.parallel import Parallel instead of joblib Parallel. The following files contains imports to joblib.Parallel:"
    echo "$joblib_Parallel_import"
    joblib_status=1
fi

if [[ $joblib_status -eq 0 ]]
then
    echo -e "No problems detected related to joblib imports\n"
else
    echo -e "\nProblems detected by joblib import check\n"
    global_status=1
fi

echo -e "### Linting completed ###\n"

if [[ $global_status -eq 1 ]]
then
    echo -e "Linting failed\n"
    exit 1
else
    echo -e "Linting passed\n"
    exit 0
fi
