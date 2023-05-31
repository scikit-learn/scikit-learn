#!/bin/bash

set -e
# pipefail is necessary to propagate exit codes
set -o pipefail

black --check --diff .
echo -e "No problem detected by black\n"

flake8 --show-source .
echo -e "No problem detected by flake8\n"

mypy sklearn/
echo -e "No problem detected by mypy\n"

cython-lint sklearn/
echo -e "No problem detected by cython-lint\n"

# For docstrings and warnings of deprecated attributes to be rendered
# properly, the property decorator must come before the deprecated decorator
# (else they are treated as functions)

# do not error when grep -B1 "@property" finds nothing
set +e
bad_deprecation_property_order=`git grep -A 10 "@property"  -- "*.py" | awk '/@property/,/def /' | grep -B1 "@deprecated"`

if [ ! -z "$bad_deprecation_property_order" ]
then
    echo "property decorator should come before deprecated decorator"
    echo "found the following occurrences:"
    echo $bad_deprecation_property_order
    exit 1
fi

# Check for default doctest directives ELLIPSIS and NORMALIZE_WHITESPACE

doctest_directive="$(git grep -nw -E "# doctest\: \+(ELLIPSIS|NORMALIZE_WHITESPACE)")"

if [ ! -z "$doctest_directive" ]
then
    echo "ELLIPSIS and NORMALIZE_WHITESPACE doctest directives are enabled by default, but were found in:"
    echo "$doctest_directive"
    exit 1
fi

joblib_delayed_import="$(git grep -l -A 10 -E "joblib import.+delayed" -- "*.py" ":!sklearn/utils/_joblib.py" ":!sklearn/utils/parallel.py")"
if [ ! -z "$joblib_delayed_import" ]; then
    echo "Use from sklearn.utils.parallel import delayed instead of joblib delayed. The following files contains imports to joblib.delayed:"
    echo "$joblib_delayed_import"
    exit 1
fi
joblib_Parallel_import="$(git grep -l -A 10 -E "joblib import.+Parallel" -- "*.py" ":!sklearn/utils/_joblib.py" ":!sklearn/utils/parallel.py")"
if [ ! -z "$joblib_Parallel_import" ]; then
    echo "Use from sklearn.utils.parallel import Parallel instead of joblib Parallel. The following files contains imports to joblib.Parallel:"
    echo "$joblib_Parallel_import"
    exit 1
fi
