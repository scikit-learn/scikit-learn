# Hi everyone!!!
from random import randint

class BisectingKMeans():
    def __init__(self, max_n_clusters):
       self.max_n_clusters = max_n_clusters # num of clusters
       self.labels = None  # np array which stores the cluster info of each data point
       self.centroids = np.zeros(max_n_clusters) # np array which stores the center location of each cluster
       self.scores = np.zeros(max_n_clusters) # np array which stores the cost of each cluster

    # #Jun and Vishwa
    # def fit(self, X):
    #    """
    #    :param X:
    #    :return:
    #    """
    #    X = np.ndarray(X)
    #    labels = np.zeros(X.shape[0])
    #    largestLabel = 0

    #    for i in range(self.max_n_clusters-1):
    #        kmeans = KMeans(n_clusters=2, init='k-means++')
    #        target_label = self._next_cluster_to_split()

    #        # find corresponding points in the target_cluster call sub_X
    #        # use labels to create a mask and filter X
           
    #        target_label_indices = np.where(labels == target_label)
    #        print(target_label_indices, "are the target label indices")
           
    #        # using the target_label_indices get the corresponding X values into sub_X
    #        sub_X = X[target_label_indices]

    #        kmeans.fit(sub_X)
    #        sub_labels = kmeans.labels_  # np.ndarray with lengh equal to sub_X, labels 0 or 1
    #        score = kmeans.score()
    #        sub_centroids = kmeans.cluster_centers_

    #        # end of this loop[]
           
    #        # update X with new labels
    #        # [0, 1] -> map elements with 0 label to targetLabel and 
    #        # elements with 1 label to largestLabel+1
    #        # update labels using sub_labels
    #        self._update_labels(sub_labels, X, target_label_indices, target_label)
    #        # update self.centroids self.scores with new labels and sub_scores

    #        # if early_stopping:
    #        #     # check stopping criteria
    #        #     break

    #    return self

    def _update_labels(self, sub_labels, X, target_label_indices, target_label):
        """
        Update the labels in X based on sub_labels outputted by the KMeans cluster split
        
        sub_labels: contains 0 or 1 for sub_X inputted into the KMeans call
        target_label_indices: contains indices within X that have the label value as target_label
        
        """
        # sub_labels=[0,0,1]-> 1 -> target_label_indices 
        # map the sub_labels to actual indices in the target_label_indices
        one_label_indices = target_label_indices[np.where(sub_labels == 1)]
        print(one_label_indices, "one_label_indices")

        # set value target_label + 1 for indices where sub_labels == 1
        X[one_label_indices] = target_label + 1
        print(X, "updated labels!")
    

    def _next_cluster_to_split(self):
        """
        check self.clusters
        the cluster with the lowest score will the next cluster to split
        return the next cluster (label) to split
        :return:
        """
        
        minScore = float('inf')
        cluster = 0
        print(self.scores)
        for i in range(0, len(self.scores)):
            if minScore > self.scores[i]:
                minScore = self.scores[i]
                cluster = i
                
        return cluster
         
       
    #Shruti and Vishwaa
    def _set_cluster_location(self, labelIndex, X, isFirst=False):
        """
            Set the cluster location for the cluster with label value
            
            Parameters:
            
            labelIndex: int (3)
            
            DataSet
            X: [[a1,,,.....,a5,a6,...,an], [b1,b2,b3...,bn], .... [n1,n2,n3,...,nn]] 
            
            isFirst: default False
            
            use findRandomRow in this function

            Return:
            updated self.centroids for the label X

        """
        #Select a random point from the data points as a cluster location
        #add to self.centroids

        cordinate = []
        #NOTE: Assuming each labelIndex comes in order
        if isFirst:
            value = randint(0, len(X)-1)
            cordinate = X[value]
        # randompoint = select random from DataSet
            #generate integer from 0 to len(dataset)
            #get the coordinate at the integer generated
            #this is the randompoint
        #set the first index of the clusters to randompoint
        else:
            dimension_size = len(X[0])
            for i in range(0, dimension_size):
                value = randint(0, 1000000) #todo: what is the coordinate boundary?
                cordinate.append(value)
            
            #labeindex is the start of a new cluster
            # any random cordinate
            #append to clusters
        self.centroids[labelIndex] = cordinate
        #current state
            #Centroids: [[3,2],[5,6]] ---- this function is going to set this (we assume its not done yet)
            #Lables: [0, 0, 1, 1] --- 2 clusters             
            #dataset: [[1,2],[3,4],[5,6],[7,8]]

            #_set_cluster_location()
    
    #Ryland and Aliza
    def inertia_cost(self, x1, x2):
        '''
        >>> inertia_cost([1, 0, 1], [0, 1, 1])
            2
        '''
        inertia = 0
        print("x1", x1)
        print("x2", x2)
        for a, b in zip(x1, x2):
            inertia += (a - b)**2
        return inertia
   
    #Shruti and Vishwaa
    def set_cluster_cost(self, cluster, X):
        """
        Set the cluster cost for the cluster with label value
        """
        
        # check what they did in kmeans, the cost function should have been implemented
        # you only need to update the self.scores here
        print(cluster)
        center = self.centroids[cluster]
        cost = 0
        count = 0.0
        
        for i in range(0, len(self.labels)):
            if self.labels[i] == cluster:
                
                cost += self.inertia_cost(X[i], center)
                count+=1.0
        print("cost", cost)
        print(cost/count)
        self.scores[cluster] = float(cost/count)

        
    
        
    #Ryland and Aliza
    def predict(self, x):
        """Predict the closest cluster each sample in X belongs to.

            Parameters
            ----------
            X : {array-like, sparse matrix} of shape (n_samples, n_features)
                New data to predict.

            Returns
            -------
            labels : ndarray of shape (n_samples,)
                Ind   
            """
        predictions = []

        # for each point in X, check which centroid is closest
        for x in X:
            # inertia = sum of squared distance to cluster centroid
            min_inertia = float('inf')
            label = 0

            for i in range(len(self.centroids)):
                centroid = self.centroids[i]
                inertia = self.inertia_cost(x, centroid)

                if inertia < min_inertia:
                    min_inertia = inertia
                    label = i
            
            predictions.append(label)

        return predictions

# test
if __name__ == "__main__":
    #from sklearn.cluster import BisectingKMeans
    import numpy as np

    X = np.array([[1],[2],[3],
               [4],[5],[6], [7],[8],[10]])

                
    clf = BisectingKMeans(max_n_clusters=2)
    
    clf.labels = np.array([0,0,0,1,1,1, 2,2,2])
    clf.centroids = np.array([[2], [5], [8.33]])
    clf.scores = np.array([2, 2, 4],dtype=np.float64)

    y = clf._next_cluster_to_split()
    z = clf.set_cluster_cost(0, X)
    a = clf._set_cluster_location(0, X)

    print(clf.centroids)