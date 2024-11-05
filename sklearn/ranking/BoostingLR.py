from jpype import JClass, JObject
from jpype.types import JDouble, JArray
import numpy as np
from .utils import kendalls_tau


# Implement the BoostLR class in Python
class BoostingLR:
    def __init__(self, max_iterations=50, seed=None):
        self.max_iterations = max_iterations
        self.m_Seed = seed
        self.m_Classifiers = None
        self.iters = 0
        print("seed: ", self.m_Seed)

    def build_classifier(self, data, lrt):

        # Dynamically load the Java classes
        Instances = JClass('weka.core.Instances')
        AbstractClassifier = JClass('weka.classifiers.AbstractClassifier')

        # Initialize
        self.m_Classifiers = AbstractClassifier.makeCopies(lrt, self.max_iterations)
        self.m_Alphas = np.zeros(self.max_iterations)
        
        # Use Java's Random class
        random_instance = JClass("java.util.Random")(self.m_Seed)

        
        training = Instances(data)
        num_instances = training.numInstances()

        t = 0

        # Initialize weights
        for i in range(num_instances):
            instance = training.instance(i)
            instance.setWeight(1)

        # Normalize weights
        sum_of_weights = training.sumOfWeights()
        for i in range(num_instances):
            instance = training.instance(i)
            instance.setWeight(instance.weight() / sum_of_weights)

        L_avg = 0

        while t < self.max_iterations:
            print("#", end="", flush=True)
            self.iters += 1

            # Obtain current weights
            weights = np.array([training.instance(i).weight() for i in range(num_instances)])

            # Sample and train
            sample = training.resampleWithWeights(random_instance, JArray(JDouble)(weights.tolist()))
            self.m_Classifiers[t].buildClassifier(sample)

            # Calculate losses and sum of losses
            Ls = np.zeros(num_instances)
            sum_of_losses = 0
            for i in range(num_instances):
                instance = training.instance(i)

                prefs = self.preferences(instance)
                preds = self.m_Classifiers[t].distributionForInstance(instance)

                l = 1 - kendalls_tau(prefs, preds)
                Ls[i] = l
                sum_of_losses += l

            if sum_of_losses > 0:
                # Normalize losses
                Ls /= sum_of_losses

                # Calculate L_avg
                L_avg = np.dot(Ls, weights)

                # Calculate Beta, and Alpha
                Beta = L_avg / (1.0 - L_avg)
                Alpha = np.log(1.0 / Beta)

                self.m_Alphas[t] = Alpha

                # Update new weights
                for i in range(num_instances):
                    instance = training.instance(i)
                    instance.setWeight(weights[i] * (Beta ** (1 - Ls[i])))

                # Normalize weights
                sum_of_weights = training.sumOfWeights()
                for i in range(num_instances):
                    instance = training.instance(i)
                    instance.setWeight(instance.weight() / sum_of_weights)
            else:
                Alpha = 1000
                self.m_Alphas[t] = Alpha
                break

            t += 1

        for i in range(t + 1, self.max_iterations):
            self.m_Classifiers[i] = self.m_Classifiers[t]
            self.m_Alphas[i] = 0

        print()

    def preferences(self, instance):
        # Cast the instance to PreferenceDenseInstance
        pdi = JClass("weka.core.labelranking.PreferenceDenseInstance")(instance)

        # Print all available keys in the HashMap for debugging
        all_keys = pdi.getHashMap().keySet()
        # print(f"Available keys in the HashMap: {all_keys}")

        # Try to find the first key with a non-null value in the hashmap
        for key in all_keys:
            prefs = pdi.getHashMap().get(key)
            if prefs is not None:
                break
        else:
            raise ValueError("No preference matrix found in the instance.")

        # Sum the rows and store the result in a list
        temp = []
        for row in prefs:
            sum_row = sum(row)
            temp.append(len(prefs) - 1 - sum_row)

        # Convert the result to a numpy array (or a standard Python list)
        result = np.array(temp, dtype=float)

        return result

    def calc_predss(self, instance):
        predss = []
        for i in range(self.iters):
            preds = self.m_Classifiers[i].distributionForInstance(instance)
            predss.append(preds)
        return predss


    def weighted_borda(self, rankings, weights):
        
        num_labels = len(rankings[0])
        
        # Initialize weightedPreds with zeros
        weighted_preds = np.zeros(num_labels)
        
        # Accumulate weighted predictions
        for i in range(len(rankings)):
            preds = rankings[i]
            weight = weights[i]
            
            for j in range(num_labels):
                weighted_preds[j] += weight * preds[j]
        
        # Determine final predictions based on weighted Borda count
        final_preds = np.zeros(num_labels)
        for j in range(num_labels):
            min_index = np.argmin(weighted_preds)
            final_preds[min_index] = float(j)
            weighted_preds[min_index] = np.inf  # Set to infinity to exclude this index in the next iteration
        
        return final_preds
    
    def predict(self, instance):
        predss = self.calc_predss(instance)
        return self.weighted_borda(predss, self.m_Alphas)
        
   