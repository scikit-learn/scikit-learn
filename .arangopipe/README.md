# Machine Learning Metadata Management with ArangoDB

Machine learning model development requires extensive experimentation. Experiments are required for:

1. Understanding the characteristics of the data

2. Testing hypothesis about the data.
3.  Testing candidate models for a machine learning task.

4. Determining the best set of hyper-parameters for a machine learning task.

5. Determining the best optimization algorithm for a machine learning task.

6. Determining the best parameters for an optimization algorithm.

7. Determine the effects of various pre-procrossesing tasks on algorithm performance.

8. Feature selection and feature extraction.

9. Best loss function for a machine learning task.

10. Best model metric to match a business metric.

The above list is by no means exhaustive. They are just some of the motivations for performing experiments for machine learning model development. When these experiments are performed, there is a need to capture the results of the experiments and key results in the form of standard machine learning artifacts, for example a notebook used for performing this experiment. There is also a need to organize the data captured from the above experiments in a cohesive manner. Related project activities must be connected. Personel working of machine learning projects must be able to query this data in an adhoc manner to obtain information that they or their co-workers performed. This is the need for which Arangopipe was created.

Arangopipe is a tool for machine learning metadata management using ArangoDB. ArangoDB is a multi-model database and supports a graph and document oriented data model. This makes it a very convinient tool to store meta-data from machine learning experiments. Arangopipe is an API for caturing and analyzing meta-data from machine learning experiments. In the notebook we provide a brief overview of tracking machine learning metadata using Arangopipe. For more information about Arangopipe, please see the project page.

**NOTE**:

You should be able to run this notebook without installing ArangoDB. The notebook connects to a managed service (cloud) instance of ArangoDB.