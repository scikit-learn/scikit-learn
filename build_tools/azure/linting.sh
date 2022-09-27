#!/bin/bash

set -e
# pipefail is necessary to propagate exit codes
set -o pipefail

flake8 --show-source .
echo -e "No problem detected by flake8\n"

# For docstrings and warnings of deprecated attributes to be rendered
# properly, the property decorator must come before the deprecated decorator
# (else they are treated as functions)

# do not error when grep -B1 "@property" finds nothing
set +e
bad_deprecation_property_order=`git grep -A 10 "@property"  -- "*.py" | awk '/@property/,/def /' | grep -B1 "@deprecated"`

if [ ! -z "$bad_deprecation_property_order" ]
then
    echo "property decorator should come before deprecated decorator"
    echo "found the following occurrencies:"
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

joblib_import="$(git grep -l -A 10 -E "joblib import.+delayed" -- "*.py" ":!sklearn/utils/_joblib.py" ":!sklearn/utils/fixes.py")"

if [ ! -z "$joblib_import" ]; then
    echo "Use from sklearn.utils.fixes import delayed instead of joblib delayed. The following files contains imports to joblib.delayed:"
    echo "$joblib_import"
    exit 1
fi
