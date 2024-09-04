Install conda using the
`miniforge installers <https://github.com/conda-forge/miniforge#miniforge>`__ (no
administrator permission required). Then run:

.. prompt:: bash

  conda create -n sklearn-env -c conda-forge scikit-learn
  conda activate sklearn-env

In order to check your installation, you can use:

.. prompt:: bash

  conda list scikit-learn  # show scikit-learn version and location
  conda list               # show all installed packages in the environment
  python -c "import sklearn; sklearn.show_versions()"
