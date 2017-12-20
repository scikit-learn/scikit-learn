# Requirements files

Dependencies for scikit-learn can be installed using either `pip` or `conda`.

* [default.txt](default.txt) &mdash; Default requirements
* [docs.txt](docs.txt) &mdash; Requirements for building documentation
* [test.txt](test.txt) &mdash; Requirements for running tests

## Example usage with `pip`

Testing dependencies can be installed with

```bash
pip install -r requirements/test.txt
```

and the required dependencies for building the documentation can be installed with

```bash
pip install -r requirements/docs.txt
```

## Example usage with `conda`

Testing dependencies can be installed with

```bash
conda install --file requirements/test.txt
```

However, one of the dependencies for building the documentation, `sphinx-gallery`, is missing from the default conda channel and should be installed from conda-forge.

```bash
conda config --append channels conda-forge
conda install --file requirements/docs.txt
```
