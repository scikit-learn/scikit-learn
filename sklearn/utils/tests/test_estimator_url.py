from urllib.parse import urlparse

from sklearn.base import BaseEstimator


class CustomValidEstimator(BaseEstimator):
    def fit(self):
        return self


class CustomEstimatorTemplateOverride(CustomValidEstimator):
    """
    With this custom estimator,
    we want to check that just changing the template yields the correct URL.
    """

    # Private values, used only for the test
    _domain = "example.com"

    _doc_link = f"https://{_domain}/{{major}}.{{minor}}/docs/{{estimator_name}}.html"


class CustomEstimatorMethodOverride(CustomValidEstimator):
    """
    With this custom estimator,
    we want to check that overriding `_get_url_link` yield the correct URL.
    """

    # Private values, used only for the test
    _domain = "example.com"

    def _get_url_link(self):
        estimator_name = self.__class__.__name__
        return f"https://{self._domain}/docs/{estimator_name}.html"


def test_defaults():
    """
    Tests that the default behavior works as expected.

    We want to assert that:
    - The returned URL is valid
    - The domain is `scikit-learn.org`
    - The URL contains the name of the estimator
    """
    estimator = CustomValidEstimator().fit()
    url = estimator._get_url_link()
    parsed_url = urlparse(url)
    # Perform assertions to check the expected behavior
    assert parsed_url.scheme in {"http", "https"}
    assert parsed_url.netloc == "scikit-learn.org"
    assert CustomValidEstimator.__name__ in parsed_url.path


def test_template_override():
    """
    Tests that overriding the URL template works as expected.

    We want to assert that:
    - The returned URL is valid
    - The domain is `example.org`
    - The URL contains the name of the estimator
    """
    estimator = CustomEstimatorTemplateOverride()
    url = estimator._get_url_link()
    parsed_url = urlparse(url)
    # Perform assertions to check the expected behavior
    assert parsed_url.scheme in {"http", "https"}
    assert parsed_url.netloc == CustomEstimatorTemplateOverride._domain
    assert CustomEstimatorTemplateOverride.__name__ in parsed_url.path


def test_method_override():
    """
    Tests that overriding the `_get_url_link` method works as expected.

    We want to assert that:
    - The returned URL is valid
    - The domain is the same as the one specified in the class
    - The URL contains the name of the estimator
    """
    estimator = CustomEstimatorMethodOverride()
    url = estimator._get_url_link()
    parsed_url = urlparse(url)
    # Perform assertions to check the expected behavior
    assert parsed_url.scheme in {"http", "https"}
    assert parsed_url.netloc == CustomEstimatorMethodOverride._domain
    assert CustomEstimatorTemplateOverride.__name__ in parsed_url.path
