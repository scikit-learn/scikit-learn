# Requirements files

Dependencies can be installed using either `pip` or `conda`.

## `pip`

* [pip/default.txt](pip/default.txt) &mdash; Default requirements
* [pip/docs.txt](pip/docs.txt) &mdash; Documentation requirements
* [pip/test.txt](pip/test.txt) &mdash; Requirements for running tests

#### Example usage

Install the required dependencies for building the documentation.

```bash
pip install -r requirements/pip/docs.txt
```

## `conda`

* [conda/default.yml](conda/default.yml) &mdash; Default requirements
* [conda/docs.yml](conda/docs.yml) &mdash; Documentation requirements
* [conda/test.yml](conda/test.yml) &mdash; Requirements for running tests

#### Example usage

Install the required dependencies for building the documentation.

```bash
conda env update --file requirements/conda/docs.yml
```
