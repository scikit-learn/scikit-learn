import pytest
import pandas as pd
import pandas.util.testing as tm
from pandas.core.dtypes.missing import array_equivalent


# Fixtures
# ========
@pytest.fixture
def df():
    """DataFrame with columns 'L1', 'L2', and 'L3' """
    return pd.DataFrame({'L1': [1, 2, 3],
                         'L2': [11, 12, 13],
                         'L3': ['A', 'B', 'C']})


@pytest.fixture(params=[[], ['L1'], ['L1', 'L2'], ['L1', 'L2', 'L3']])
def df_levels(request, df):
    """DataFrame with columns or index levels 'L1', 'L2', and 'L3' """
    levels = request.param

    if levels:
        df = df.set_index(levels)

    return df


@pytest.fixture
def df_ambig(df):
    """DataFrame with levels 'L1' and 'L2' and labels 'L1' and 'L3' """
    df = df.set_index(['L1', 'L2'])

    df['L1'] = df['L3']

    return df


@pytest.fixture
def df_duplabels(df):
    """DataFrame with level 'L1' and labels 'L2', 'L3', and 'L2' """
    df = df.set_index(['L1'])
    df = pd.concat([df, df['L2']], axis=1)

    return df


@pytest.fixture
def panel():
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        return pd.Panel()


# Test is label/level reference
# =============================
def get_labels_levels(df_levels):
    expected_labels = list(df_levels.columns)
    expected_levels = [name for name in df_levels.index.names
                       if name is not None]
    return expected_labels, expected_levels


def assert_label_reference(frame, labels, axis):
    for label in labels:
        assert frame._is_label_reference(label, axis=axis)
        assert not frame._is_level_reference(label, axis=axis)
        assert frame._is_label_or_level_reference(label, axis=axis)


def assert_level_reference(frame, levels, axis):
    for level in levels:
        assert frame._is_level_reference(level, axis=axis)
        assert not frame._is_label_reference(level, axis=axis)
        assert frame._is_label_or_level_reference(level, axis=axis)


# DataFrame
# ---------
@pytest.mark.parametrize('axis', [0, 1])
def test_is_level_or_label_reference_df_simple(df_levels, axis):

    # Compute expected labels and levels
    expected_labels, expected_levels = get_labels_levels(df_levels)

    # Transpose frame if axis == 1
    if axis == 1:
        df_levels = df_levels.T

    # Perform checks
    assert_level_reference(df_levels, expected_levels, axis=axis)
    assert_label_reference(df_levels, expected_labels, axis=axis)


@pytest.mark.parametrize('axis', [0, 1])
def test_is_level_reference_df_ambig(df_ambig, axis):

    # Transpose frame if axis == 1
    if axis == 1:
        df_ambig = df_ambig.T

    # df has both an on-axis level and off-axis label named L1
    # Therefore L1 should reference the label, not the level
    assert_label_reference(df_ambig, ['L1'], axis=axis)

    # df has an on-axis level named L2 and it is not ambiguous
    # Therefore L2 is an level reference
    assert_level_reference(df_ambig, ['L2'], axis=axis)

    # df has a column named L3 and it not an level reference
    assert_label_reference(df_ambig, ['L3'], axis=axis)


# Series
# ------
def test_is_level_reference_series_simple_axis0(df):

    # Make series with L1 as index
    s = df.set_index('L1').L2
    assert_level_reference(s, ['L1'], axis=0)
    assert not s._is_level_reference('L2')

    # Make series with L1 and L2 as index
    s = df.set_index(['L1', 'L2']).L3
    assert_level_reference(s, ['L1', 'L2'], axis=0)
    assert not s._is_level_reference('L3')


def test_is_level_reference_series_axis1_error(df):

    # Make series with L1 as index
    s = df.set_index('L1').L2

    with tm.assert_raises_regex(ValueError, "No axis named 1"):
        s._is_level_reference('L1', axis=1)


# Panel
# -----
def test_is_level_reference_panel_error(panel):
    msg = ("_is_level_reference is not implemented for {type}"
           .format(type=type(panel)))

    with tm.assert_raises_regex(NotImplementedError, msg):
        panel._is_level_reference('L1', axis=0)


def test_is_label_reference_panel_error(panel):
    msg = ("_is_label_reference is not implemented for {type}"
           .format(type=type(panel)))

    with tm.assert_raises_regex(NotImplementedError, msg):
        panel._is_label_reference('L1', axis=0)


def test_is_label_or_level_reference_panel_error(panel):
    msg = ("_is_label_or_level_reference is not implemented for {type}"
           .format(type=type(panel)))

    with tm.assert_raises_regex(NotImplementedError, msg):
        panel._is_label_or_level_reference('L1', axis=0)


# Test _check_label_or_level_ambiguity_df
# =======================================

# DataFrame
# ---------
@pytest.mark.parametrize('axis', [0, 1])
def test_check_label_or_level_ambiguity_df(df_ambig, axis):

    # Transpose frame if axis == 1
    if axis == 1:
        df_ambig = df_ambig.T

    # df_ambig has both an on-axis level and off-axis label named L1
    # Therefore L1 is ambiguous
    with tm.assert_produces_warning(FutureWarning,
                                    clear=True) as w:

        assert df_ambig._check_label_or_level_ambiguity('L1', axis=axis)
        warning_msg = w[0].message.args[0]
        if axis == 0:
            assert warning_msg.startswith("'L1' is both an index level "
                                          "and a column label")
        else:
            assert warning_msg.startswith("'L1' is both a column level "
                                          "and an index label")

    # df_ambig has an on-axis level named L2 and it is not ambiguous
    # No warning should be raised
    with tm.assert_produces_warning(None):
        assert not df_ambig._check_label_or_level_ambiguity('L2', axis=axis)

    # df_ambig has an off-axis label named L3 and it is not ambiguous
    with tm.assert_produces_warning(None):
        assert not df_ambig._is_level_reference('L3', axis=axis)


# Series
# ------
def test_check_label_or_level_ambiguity_series(df):

    # A series has no columns and therefore references are never ambiguous

    # Make series with L1 as index
    s = df.set_index('L1').L2
    with tm.assert_produces_warning(None):
        assert not s._check_label_or_level_ambiguity('L1', axis=0)
        assert not s._check_label_or_level_ambiguity('L2', axis=0)

    # Make series with L1 and L2 as index
    s = df.set_index(['L1', 'L2']).L3
    with tm.assert_produces_warning(None):
        assert not s._check_label_or_level_ambiguity('L1', axis=0)
        assert not s._check_label_or_level_ambiguity('L2', axis=0)
        assert not s._check_label_or_level_ambiguity('L3', axis=0)


def test_check_label_or_level_ambiguity_series_axis1_error(df):

    # Make series with L1 as index
    s = df.set_index('L1').L2

    with tm.assert_raises_regex(ValueError, "No axis named 1"):
        s._check_label_or_level_ambiguity('L1', axis=1)


# Panel
# -----
def test_check_label_or_level_ambiguity_panel_error(panel):
    msg = ("_check_label_or_level_ambiguity is not implemented for {type}"
           .format(type=type(panel)))

    with tm.assert_raises_regex(NotImplementedError, msg):
        panel._check_label_or_level_ambiguity('L1', axis=0)


# Test _get_label_or_level_values
# ===============================
def assert_label_values(frame, labels, axis):
    for label in labels:
        if axis == 0:
            expected = frame[label]._values
        else:
            expected = frame.loc[label]._values

        result = frame._get_label_or_level_values(label, axis=axis,
                                                  stacklevel=2)
        assert array_equivalent(expected, result)


def assert_level_values(frame, levels, axis):
    for level in levels:
        if axis == 0:
            expected = frame.index.get_level_values(level=level)._values
        else:
            expected = (frame.columns
                        .get_level_values(level=level)
                        ._values)

        result = frame._get_label_or_level_values(level, axis=axis)
        assert array_equivalent(expected, result)


# DataFrame
# ---------
@pytest.mark.parametrize('axis', [0, 1])
def test_get_label_or_level_values_df_simple(df_levels, axis):

    # Compute expected labels and levels
    expected_labels, expected_levels = get_labels_levels(df_levels)

    # Transpose frame if axis == 1
    if axis == 1:
        df_levels = df_levels.T

    # Perform checks
    assert_label_values(df_levels, expected_labels, axis=axis)
    assert_level_values(df_levels, expected_levels, axis=axis)


@pytest.mark.parametrize('axis', [0, 1])
def test_get_label_or_level_values_df_ambig(df_ambig, axis):

    # Transpose frame if axis == 1
    if axis == 1:
        df_ambig = df_ambig.T

    # df has both an on-axis level and off-axis label named L1
    # Therefore L1 is ambiguous but will default to label
    with tm.assert_produces_warning(FutureWarning):
        assert_label_values(df_ambig, ['L1'], axis=axis)

    # df has an on-axis level named L2 and it is not ambiguous
    with tm.assert_produces_warning(None):
        assert_level_values(df_ambig, ['L2'], axis=axis)

    # df has an off-axis label named L3 and it is not ambiguous
    with tm.assert_produces_warning(None):
        assert_label_values(df_ambig, ['L3'], axis=axis)


@pytest.mark.parametrize('axis', [0, 1])
def test_get_label_or_level_values_df_duplabels(df_duplabels, axis):

    # Transpose frame if axis == 1
    if axis == 1:
        df_duplabels = df_duplabels.T

    # df has unambiguous level 'L1'
    assert_level_values(df_duplabels, ['L1'], axis=axis)

    # df has unique label 'L3'
    assert_label_values(df_duplabels, ['L3'], axis=axis)

    # df has duplicate labels 'L2'
    if axis == 0:
        expected_msg = "The column label 'L2' is not unique"
    else:
        expected_msg = "The index label 'L2' is not unique"

    with tm.assert_raises_regex(ValueError, expected_msg):
        assert_label_values(df_duplabels, ['L2'], axis=axis)


# Series
# ------
def test_get_label_or_level_values_series_axis0(df):

    # Make series with L1 as index
    s = df.set_index('L1').L2
    assert_level_values(s, ['L1'], axis=0)

    # Make series with L1 and L2 as index
    s = df.set_index(['L1', 'L2']).L3
    assert_level_values(s, ['L1', 'L2'], axis=0)


def test_get_label_or_level_values_series_axis1_error(df):

    # Make series with L1 as index
    s = df.set_index('L1').L2

    with tm.assert_raises_regex(ValueError, "No axis named 1"):
        s._get_label_or_level_values('L1', axis=1)


# Panel
# -----
def test_get_label_or_level_values_panel_error(panel):
    msg = ("_get_label_or_level_values is not implemented for {type}"
           .format(type=type(panel)))

    with tm.assert_raises_regex(NotImplementedError, msg):
        panel._get_label_or_level_values('L1', axis=0)


# Test _drop_labels_or_levels
# ===========================
def assert_labels_dropped(frame, labels, axis):
    for label in labels:
        df_dropped = frame._drop_labels_or_levels(label, axis=axis)

        if axis == 0:
            assert label in frame.columns
            assert label not in df_dropped.columns
        else:
            assert label in frame.index
            assert label not in df_dropped.index


def assert_levels_dropped(frame, levels, axis):
    for level in levels:
        df_dropped = frame._drop_labels_or_levels(level, axis=axis)

        if axis == 0:
            assert level in frame.index.names
            assert level not in df_dropped.index.names
        else:
            assert level in frame.columns.names
            assert level not in df_dropped.columns.names


# DataFrame
# ---------
@pytest.mark.parametrize('axis', [0, 1])
def test_drop_labels_or_levels_df(df_levels, axis):

    # Compute expected labels and levels
    expected_labels, expected_levels = get_labels_levels(df_levels)

    # Transpose frame if axis == 1
    if axis == 1:
        df_levels = df_levels.T

    # Perform checks
    assert_labels_dropped(df_levels, expected_labels, axis=axis)
    assert_levels_dropped(df_levels, expected_levels, axis=axis)

    with tm.assert_raises_regex(ValueError, "not valid labels or levels"):
        df_levels._drop_labels_or_levels('L4', axis=axis)


# Series
# ------
def test_drop_labels_or_levels_series(df):

    # Make series with L1 as index
    s = df.set_index('L1').L2
    assert_levels_dropped(s, ['L1'], axis=0)

    with tm.assert_raises_regex(ValueError, "not valid labels or levels"):
        s._drop_labels_or_levels('L4', axis=0)

    # Make series with L1 and L2 as index
    s = df.set_index(['L1', 'L2']).L3
    assert_levels_dropped(s, ['L1', 'L2'], axis=0)

    with tm.assert_raises_regex(ValueError, "not valid labels or levels"):
        s._drop_labels_or_levels('L4', axis=0)


# Panel
# -----
def test_drop_labels_or_levels_panel_error(panel):
    msg = ("_drop_labels_or_levels is not implemented for {type}"
           .format(type=type(panel)))

    with tm.assert_raises_regex(NotImplementedError, msg):
        panel._drop_labels_or_levels('L1', axis=0)
