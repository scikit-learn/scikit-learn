BoostLR Algorithm
=================

.. currentmodule:: sklearn.ranking

Overview
--------
The BoostLR algorithm is a supervised learning method designed for label ranking tasks.
It leverages boosting techniques to improve the predictive performance of label ranking models.

Features
--------
- Implements boosting-based label ranking using Weka and JPype.
- Integrates seamlessly with the `scikit-learn` API.
- Supports custom parameters, such as the maximum number of boosting iterations and seed.

Requirements
------------
- **Java Version**: Ensure that you have Java 8 installed and configured on your system, as Weka requires Java 8 for proper functionality.

Dependencies
------------
- `JPype1`: Used to interface with Weka’s Java implementation.

Note: Unlike other Weka-based Python integrations, this implementation does not use `python-weka-wrapper3`.

Usage Example
-------------
Here’s a simple example of how to use the `BoostingLRWrapper` class:

.. code-block:: python

    from sklearn.ranking import BoostingLRWrapper
    import jpype

    # Load your Weka Instances data
    train_data = load_dataset_as_Instances("path/to/train_data.arff")
    test_data = load_dataset_as_Instances("path/to/test_data.arff")

    # Initialize and train the model
    model = BoostingLRWrapper(max_iterations=100, seed=42)
    model.fit(train_data)

    # Make predictions
    predictions = model.predict(test_data)
    print("Predicted Rankings:", predictions)

    # Evaluate the model (implementing your own scoring metric)
    score = model.score(test_data)
    print("Model Score:", score)

API Reference
-------------
.. autoclass:: BoostingLRWrapper
    :members:
    :undoc-members:
    :show-inheritance:

Parameters
----------
- **max_iterations**: `int`, default=50
  The maximum number of boosting iterations.

- **seed**: `int`
  The seed for the random number generator. Used for reproducibility.

Methods
-------
- **fit(train_data)**: Train the model on Weka Instances data.
- **predict(test_data)**: Predict rankings for Weka Instances data.
- **score(test_data)**: Compute the score of the model using a ranking evaluation metric.

Notes
-----
- Make sure that Weka is properly installed and configured on your system.
- **Java Requirement**: The BoostLR algorithm requires Java 8. Ensure your `JAVA_HOME` environment variable points to the Java 8 installation directory.
- This implementation is specifically designed for label ranking tasks and may require specific data preprocessing.


