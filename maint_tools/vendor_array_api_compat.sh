#!/bin/bash

# Vendors https://github.com/data-apis/array-api-compat/ into sklearn/externals

set -o nounset
set -o errexit

URL="https://github.com/data-apis/array-api-compat.git"
VERSION="1.12"

ROOT_DIR=sklearn/externals/array_api_compat

rm -rf $ROOT_DIR
mkdir $ROOT_DIR
mkdir $ROOT_DIR/.tmp
git clone $URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $VERSION
popd
mv -v $ROOT_DIR/.tmp/array_api_compat/* $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/LICENSE $ROOT_DIR/
rm -rf $ROOT_DIR/.tmp

echo "Update this directory using maint_tools/vendor_array_api_compat.sh" >$ROOT_DIR/README.md
