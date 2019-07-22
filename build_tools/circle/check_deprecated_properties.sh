#!/bin/bash

# For docstrings and warnings of deprecated attributes to be rendered
# properly, the property decorator must come before the deprecated decorator
# (else they are treated as functions)
bad_deprecation_property_order=`git grep -A 10 "@property" | awk '/@property/,/def /' | grep -B1 "@deprecated"`
# exclude this file from the matches
bad_deprecation_property_order=`echo $bad_deprecation_property_order | grep -v check_deprecated_properties`

if [ ! -z "$bad_deprecation_property_order" ]
then
    echo "property decorator should come before deprecated decorator"
    echo "found the following occurrencies:"
    echo $bad_deprecation_property_order
    exit 1
fi
