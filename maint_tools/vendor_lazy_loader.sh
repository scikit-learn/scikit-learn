#!/bin/bash

# Vendors https://github.com/scientific-python/lazy_loader/ into sklearn/externals

SHA="2334bd279d40e0dadd3af48fe4ec494d3ce7f47d"
URL="https://github.com/scientific-python/lazy_loader/archive/2334bd279d40e0dadd3af48fe4ec494d3ce7f47d.tar.gz"
ROOT_DIR=sklearn/externals/_lazy_loader

rm -rf $ROOT_DIR/*

curl -s -L $URL |
    tar xvz --strip-components=1 -C sklearn/externals/_lazy_loader \
        lazy_loader-$SHA/lazy_loader/__init__.py \
        lazy_loader-$SHA/LICENSE.md

mv $ROOT_DIR/lazy_loader/__init__.py $ROOT_DIR/__init__.py
rmdir $ROOT_DIR/lazy_loader

echo "Update this directory using maint_tools/vendor_lazy_loader.sh" >$ROOT_DIR/README.md
