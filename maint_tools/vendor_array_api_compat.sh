#!/bin/bash

# Vendors https://github.com/data-apis/array-api-compat/ into sklearn/externals

ARRAY_API_COMPAT_SHA="b32a5b32892f5f4b5052ef54a04b8ed51936b008"
URL="https://github.com/data-apis/array-api-compat/archive/$ARRAY_API_COMPAT_SHA.tar.gz"

rm -rf sklearn/externals/_array_api_compat

curl -s -L $URL |
    tar xvz --strip-components=1 -C sklearn/externals array-api-compat-$ARRAY_API_COMPAT_SHA/array_api_compat

mv sklearn/externals/array_api_compat sklearn/externals/_array_api_compat
