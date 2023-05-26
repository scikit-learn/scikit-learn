#!/bin/bash

# Vendors https://github.com/scientific-python/lazy_loader/ into sklearn/externals

SHA="2334bd279d40e0dadd3af48fe4ec494d3ce7f47d"
URL="https://github.com/scientific-python/lazy_loader/archive/$SHA.tar.gz"

rm -rf sklearn/externals/_lazy_loader

curl -s -L $URL |
    tar xvz --strip-components=1 -C sklearn/externals lazy_loader-$SHA/lazy_loader/__init__.py

mv sklearn/externals/lazy_loader sklearn/externals/_lazy_loader
