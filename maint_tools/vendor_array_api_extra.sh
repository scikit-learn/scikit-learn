#!/bin/bash

# Vendors https://github.com/data-apis/array-api-extra/ into sklearn/externals

set -o nounset
set -o errexit

URL="https://github.com/bruAristimunha/array-api-extra.git"
VERSION="e34d415ec22eaa3f069fa0d3a35b62ed1d4e181d"

ROOT_DIR=sklearn/externals/array_api_extra

rm -rf $ROOT_DIR
mkdir $ROOT_DIR
mkdir $ROOT_DIR/.tmp
git clone $URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $VERSION
popd
mv -v $ROOT_DIR/.tmp/src/array_api_extra/* $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/LICENSE $ROOT_DIR/
rm -rf $ROOT_DIR/.tmp

echo "Update this directory using maint_tools/vendor_array_api_extra.sh" >$ROOT_DIR/README.md
