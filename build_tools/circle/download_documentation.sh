#!/bin/bash

set -e
set -x

wget $GHA_META
unzip doc.zip -d doc
