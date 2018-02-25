# Author: Jonas Eschle 'Mayou36@jonas.eschle.com'
# License: BSD 3 clause
# Testing the absolute import of attributes like sklearn.module or sklearn.attribute


from sklearn.utils.testing import assert_equal


def test_abs_import_skl():
    # Test if the modules/attributes have been imported correctly
    import sklearn

    failed_attributes = [att for att in sklearn.__all__ if not hasattr(sklearn, att)]
    assert_equal(failed_attributes, [])
