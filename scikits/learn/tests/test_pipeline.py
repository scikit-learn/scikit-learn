"""
Test the pipeline module.
"""

from nose.tools import assert_raises, assert_equal

from ..pipeline import Pipeline
from ..svm import SVC
from ..feature_selection import SelectKBest, f_classif

def test_pipeline_init():
    """ Test the various init parameters of the pipeline.
    """
    assert_raises(AssertionError, Pipeline)
    clf = SVC()
    # Smoke test with only an estimator
    pipe = Pipeline(estimator=clf)
    assert_equal(pipe._get_param_names(), clf._get_param_names())
    assert_equal(pipe._get_params(), clf._get_params())

    # Check that params are set
    pipe._set_params(C=0.1)
    assert_equal(clf.C, 0.1)
    # Smoke test the repr:
    repr(pipe)
    
    # Test with two objects
    clf = SVC()
    filter1 = SelectKBest(f_classif)
    pipe = Pipeline(transforms=[filter1], estimator=clf)
    assert_equal(set(pipe._get_param_names()
                        ).intersection(clf._get_param_names()), 
                 set(clf._get_param_names()))

    # Check that params are set
    pipe._set_params(C=0.1)
    assert_equal(clf.C, 0.1)
    # Smoke test the repr:
    repr(pipe)

    # Test with two named objects
    clf = SVC()
    filter1 = SelectKBest(f_classif)
    pipe = Pipeline(transforms=[filter1], estimator=clf, 
                    names=['anova', 'svm'])

    # Check that params are set
    assert_raises(AssertionError, pipe._set_params, anova_C=0.1)
    pipe._set_params(svm_C=0.1)
    assert_equal(clf.C, 0.1)

    # Check that we can't instanciate a Pipeline with overlapping
    # parameters
    assert_raises(AssertionError, Pipeline, 
                transforms=[filter1, filter1], estimator=clf)
    # But that it's OK when specifying names
    pipe = Pipeline(transforms=[filter1, filter1], estimator=clf, 
                    names=['filter1', 'filter2'])
    # Smoke test the repr:
    repr(pipe)

    # Test reinit
    pipe2 = pipe._reinit()
    assert_equal(pipe._get_params(), pipe2._get_params())
