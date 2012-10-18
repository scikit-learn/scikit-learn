Maintainer / core-developer information
========================================

Making a release
------------------

1. Update docs:

    - edit the doc/whats_new.rst file to add release title and commit
      statistics. You can retrieve commit statistics with::

        $ git shortlog -ns 0.998..

    - edit the doc/conf.py to increase the version number

    - edit the doc/themes/scikit-learn/layout.html to change the 'News'
      entry of the front page.

2. Update the version number in sklearn/__init__.py, the __version__
   variable

3. Create the tag and push it::

    $ git tag 0.999

    $ git push origin --tags

4. create tarballs:

   - Wipe clean your repo::

       $ git clean -xfd

   - Register and upload on PyPI::

       $ python setup.py sdist register upload

   - Upload manually the tarbal on sourceforge:
     https://sourceforge.net/projects/scikit-learn/files/

5. Push the documentation to the website (see README in doc folder)


6. Build binaries for windows and push them to PyPI::

    $ python setup.py bdist_wininst upload

   And upload them also to sourceforge
