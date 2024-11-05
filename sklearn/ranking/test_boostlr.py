import pytest
from sklearn.ranking.BoostingLRWrapper import BoostingLRWrapper
from sklearn.ranking.utils import start_jvm, stop_jvm, load_dataset_as_Instances, kendalls_tau

############################################
# Fixtures
############################################

# Fixture to start and stop the JVM
@pytest.fixture(scope="session", autouse=True)
def jvm_manager():
    """Start the JVM before tests and stop it after."""
    start_jvm()
    yield  # Let the tests run
    stop_jvm()  # Stop the JVM after all tests are done

# Fixtures to load the small dataset
@pytest.fixture
def train_data():
    """Fixture to load training data as Weka Instances."""
    return load_dataset_as_Instances("./datasets/sushi_mini_train.xarff")

@pytest.fixture
def test_data():
    """Fixture to load test data as Weka Instances."""
    return load_dataset_as_Instances("./datasets/sushi_mini_test.xarff")



# Fixtures to load the same ranked dataset
@pytest.fixture
def train_data_same_rank():
    """Fixture to load training data as Weka Instances."""
    return load_dataset_as_Instances("./datasets/sushi_mini_same_rank_train.xarff")

@pytest.fixture
def test_data_same_rank():
    """Fixture to load test data as Weka Instances."""
    return load_dataset_as_Instances("./datasets/sushi_mini_same_rank_test.xarff")


############################################
# Tests
############################################

# Initialization Test
def test_initialization():
    """Test the initialization of BoostingLRWrapper."""
    model = BoostingLRWrapper(max_iterations=100, seed=7)
    assert model.max_iterations == 100
    assert model.seed == 7

# Fit Method Test
def test_fit(train_data):
    """Test the fit method of BoostingLRWrapper."""
    model = BoostingLRWrapper(max_iterations=50, seed=7)
    try:
        model.fit(train_data)
    except Exception as e:
        pytest.fail(f"Fit method failed with exception: {e}")

# Predict Method Test
def test_predict(train_data, test_data):
    """Test the predict method of BoostingLRWrapper."""

    model = BoostingLRWrapper(max_iterations=50, seed=42)
    model.fit(train_data)
    predictions = model.predict(test_data)
    
    for pred in predictions:
        # Ensure that the number of predictions matches the number of instances
        assert len(pred) == 10

def test_predict_same_rank_dataset(train_data_same_rank, test_data_same_rank):
    """Test the predict method return the same rank like the instances as test_data"""
    model = BoostingLRWrapper(max_iterations=50, seed=7)
    model.fit(train_data_same_rank)
    predictions = model.predict(test_data_same_rank)
    
    for pred in predictions:
        # Ensure that the number of predictions matches the number of instances
        assert pred == pytest.approx([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], rel=1e-9)

# Score Method Test
def test_score(train_data, test_data):
    """Test the score method of BoostingLRWrapper."""
    model = BoostingLRWrapper(max_iterations=50, seed=42)
    model.fit(train_data)
    score = model.score(test_data)
    
    # Check if the score is a float and within a reasonable range
    assert isinstance(score, float)
    print("score: ", score)
    assert 0.0 <= score <= 1.0, "Score should be between 0 and 1"

def test_score_same_rank_dataset(train_data_same_rank, test_data_same_rank):
    """Test the score method of BoostingLRWrapper."""
    model = BoostingLRWrapper(max_iterations=50, seed=42)
    model.fit(train_data_same_rank)
    score = model.score(test_data_same_rank)
    
    # Check if the score is a float and be approximately 1.0
    assert isinstance(score, float)
    print("score: ", score)
    assert score == pytest.approx(1.0, rel=1e-9), "Score should be approximately 1.0"

# Additional Tests
def test_kendalls_tau():
    """Test the kendalls_tau utility function"""
    
    # Test 1: Perfect disagreement
    prefs = [1, 2, 3, 4, 5]
    preds = [5, 4, 3, 2, 1]
    tau = kendalls_tau(prefs, preds)
    assert tau == pytest.approx(-1.0, rel=1e-9)

    # Test 2: Perfect agreement
    prefs = [1, 2, 3, 4, 5]
    preds = [1, 2, 3, 4, 5]
    tau = kendalls_tau(prefs, preds)
    assert tau == pytest.approx(1.0, rel=1e-9)

if __name__ == "__main__":
    pytest.main([__file__])

