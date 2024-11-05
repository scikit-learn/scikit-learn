from jpype import JClass
from sklearn.ranking.BoostingLR import BoostingLR
from .utils import *

class BoostingLRWrapper:
    def __init__(self, max_iterations=50, seed=None):
        self.max_iterations = max_iterations
        self.seed = seed
        self.boosting_lr = BoostingLR(max_iterations=self.max_iterations, seed=self.seed)

    def fit(self, train_data):
        """
        Train the BoostingLR model on the provided training data.
        :param train_data: Weka Instances object containing the training data
        """
        # Use the provided training data to build the classifier
        lrt = JClass("weka.classifiers.labelranking.LRT")()
        self.boosting_lr.build_classifier(train_data, lrt)

    def predict(self, test_data):
        """
        Predict rankings for the provided test data.
        :param test_data: Weka Instances object containing the test data
        :return: List of predicted rankings for each instance in the test data
        """
        predictions = []
        for i in range(test_data.numInstances()):
            instance = test_data.instance(i)
            preds = self.boosting_lr.predict(instance)
            predictions.append(preds)
        return predictions
    
    def score(self, test_instances):
        """
        Compute the score of the model on the test data.
        :param test_data: Weka Instances object containing the test data
        :param true_rankings: Weka Instances object containing the true rankings
        :return: Normalized Discounted Cumulative Gain (NDCG) score
        """
        # Predict the rankings using the model
        predicted_rankings = self.predict(test_instances)

        total_kt = 0.0
        for i in range(test_instances.numInstances()):
            instance = test_instances.instance(i)

            prefs = self.boosting_lr.preferences(instance)
            preds = predicted_rankings[i]

            total_kt += kendalls_tau(preds, prefs)

        # Calculate the average Kendall's Tau for this fold
        return total_kt / len(predicted_rankings)
