import numpy as np
import pytest
from sklearn.model_selection import ParameterSampler
from sklearn.utils import check_random_state


def test_parameter_sampler_distribution_weights():
    """Test the distribution_weights parameter of ParameterSampler."""
    # Test with uniform weights (default behavior)
    param_distributions = [
        {"a": [1, 2, 3]},
        {"b": [4, 5, 6]},
        {"c": [7, 8, 9]}
    ]
    
    # Uniform sampling
    sampler = ParameterSampler(param_distributions, n_iter=100, random_state=42)
    samples = list(sampler)
    
    # Check that we get samples from all distributions
    dist_indices = []
    for sample in samples:
        if "a" in sample:
            dist_indices.append(0)
        elif "b" in sample:
            dist_indices.append(1)
        elif "c" in sample:
            dist_indices.append(2)
    
    # With uniform sampling, we should see samples from all distributions
    assert len(set(dist_indices)) == 3
    
    # Test with custom weights
    weights = [0.1, 0.1, 0.8]  # Prefer the third distribution
    sampler = ParameterSampler(
        param_distributions, 
        n_iter=100, 
        random_state=42,
        distribution_weights=weights
    )
    samples = list(sampler)
    
    # Check that we get samples from all distributions
    dist_indices = []
    for sample in samples:
        if "a" in sample:
            dist_indices.append(0)
        elif "b" in sample:
            dist_indices.append(1)
        elif "c" in sample:
            dist_indices.append(2)
    
    # With weighted sampling, we should see more samples from the third distribution
    # Count occurrences
    counts = [dist_indices.count(i) for i in range(3)]
    # The third distribution should be sampled most frequently
    assert counts[2] > counts[0]
    assert counts[2] > counts[1]
    
    # Test with zero weights
    weights = [0.5, 0.5, 0.0]  # Never sample from the third distribution
    sampler = ParameterSampler(
        param_distributions, 
        n_iter=100, 
        random_state=42,
        distribution_weights=weights
    )
    samples = list(sampler)
    
    # Check that we get samples only from first two distributions
    dist_indices = []
    for sample in samples:
        if "a" in sample:
            dist_indices.append(0)
        elif "b" in sample:
            dist_indices.append(1)
        elif "c" in sample:
            dist_indices.append(2)
    
    # We should never sample from the third distribution
    assert 2 not in dist_indices
    # But we should sample from first two
    assert 0 in dist_indices
    assert 1 in dist_indices


def test_parameter_sampler_distribution_weights_validation():
    """Test validation of distribution_weights parameter."""
    param_distributions = [
        {"a": [1, 2, 3]},
        {"b": [4, 5, 6]},
        {"c": [7, 8, 9]}
    ]
    
    # Test with wrong length
    with pytest.raises(ValueError, match="distribution_weights must have the same length"):
        ParameterSampler(
            param_distributions, 
            n_iter=10, 
            distribution_weights=[0.5, 0.5]  # Wrong length
        )
    
    # Test with negative weights
    with pytest.raises(ValueError, match="distribution_weights must contain only non-negative values"):
        ParameterSampler(
            param_distributions, 
            n_iter=10, 
            distribution_weights=[0.5, 0.5, -0.1]
        )
    
    # Test with all zero weights
    with pytest.raises(ValueError, match="distribution_weights must contain at least one positive value"):
        ParameterSampler(
            param_distributions, 
            n_iter=10, 
            distribution_weights=[0.0, 0.0, 0.0]
        )