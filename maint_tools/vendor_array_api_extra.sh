#!/bin/bash

# Vendors https://github.com/data-apis/array-api-extra/ into sklearn/externals

VERSION="0.3.3"
URL="https://github.com/data-apis/array-api-extra/archive/refs/tags/v0.3.3.tar.gz"
ROOT_DIR=sklearn/externals/_array_api_extra

rm -rf $ROOT_DIR/*

curl -s -L $URL |
    tar xvz --strip-components=1 -C sklearn/externals/_array_api_extra \
        array-api-extra-$VERSION/src/array_api_extra/__init__.py \
        array-api-extra-$VERSION/src/array_api_extra/_funcs.py \
        array-api-extra-$VERSION/src/array_api_extra/py.typed \
        array-api-extra-$VERSION/src/array_api_extra/_lib/_compat.py \
        array-api-extra-$VERSION/src/array_api_extra/_lib/_typing.py \
        array-api-extra-$VERSION/src/array_api_extra/_lib/_utils.py \
        array-api-extra-$VERSION/LICENSE

mv $ROOT_DIR/src/array_api_extra/* $ROOT_DIR
rm -rf $ROOT_DIR/src 

echo "Update this directory using maint_tools/vendor_array_api_extra.sh" >$ROOT_DIR/README.md
