#!/bin/sh
# Script to do a local install of sphinx-gallery
rm -rf tmp sphinx_gallery

easy_install -Zeab tmp sphinx-gallery

cp -vru tmp/sphinx-gallery/sphinx_gallery/ .

echo "Remember to add sphinx_gallery to your version control"
echo "Use in case of git:"
echo "$ git add sphinx_gallery"
