def test_fastfood_performance_comparison_between_methods():
    """compares the performance of Fastfood and RKS"""
    # generate data
    X = rng.random_sample(size=(1000, 4096))
    Y = rng.random_sample(size=(10000, 4096))
    X /= X.sum(axis=1)[:, np.newaxis]
    Y /= Y.sum(axis=1)[:, np.newaxis]

    # calculate feature maps
    gamma = 10.
    sigma = np.sqrt(1 / (2 * gamma))
    number_of_features_to_generate = 4096*4

    exact_start = datetime.datetime.utcnow()
    # original rbf kernel method:
    # rbf_kernel(X, X, gamma=gamma)
    # rbf_kernel(X, Y, gamma=gamma)
    exact_end = datetime.datetime.utcnow()
    exact_spent_time = exact_end - exact_start
    print "Timimg exact rbf: \t\t", exact_spent_time

    rbf_transform = Fastfood(sigma=sigma,
                             n_components=number_of_features_to_generate,
                             tradeoff_mem_accuracy='mem',
                             random_state=42)
    _ = rbf_transform.fit(X)
    fastfood_fast_vec_start = datetime.datetime.utcnow()
    # Fastfood: approximate kernel mapping
    _ = rbf_transform.transform(X)
    _ = rbf_transform.transform(Y)
    fastfood_fast_vec_end = datetime.datetime.utcnow()
    fastfood_fast_vec_spent_time = fastfood_fast_vec_end - \
        fastfood_fast_vec_start
    print "Timimg fastfood fast vectorized: \t\t", fastfood_fast_vec_spent_time

    rks_rbf_transform = RBFSampler(gamma=gamma,
                                   n_components=number_of_features_to_generate,
                                   random_state=42)
    _ = rks_rbf_transform.fit(X)
    rks_start = datetime.datetime.utcnow()
    # Random Kitchens Sinks: approximate kernel mapping
    _ = rks_rbf_transform.transform(X)
    _ = rks_rbf_transform.transform(Y)
    rks_end = datetime.datetime.utcnow()
    rks_spent_time = rks_end - rks_start
    print "Timimg rks: \t\t\t", rks_spent_time

    assert_greater(rks_spent_time, fastfood_fast_vec_spent_time)


