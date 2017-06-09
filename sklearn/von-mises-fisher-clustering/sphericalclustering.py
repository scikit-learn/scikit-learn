# -*- coding: utf-8 -*-
"""
This module clusters the data using the VMF distribution
    
GUIDE : 
    Prof. Anand A Joshi - ajoshi@sipi.usc.edu
TEAM:
    Bhavana Ganesh - bhavanag@usc.edu
    Mahathi Vatsal Salopanthula - salopant@usc.edu
    Sayali Ghume - ghume@usc.edu
    
Contact any of the members for queries and bugs.
"""

import numpy as np
import mpmath

#getting the size of the data
def sphericalknn (data,no_clusters):
    
    H=len(data)
    W=len(data[0])
    
#calculating the mean to remove the dc component
    sum_sample=sum(data)
    norm_sample=np.linalg.norm(sum_sample)
    mean_sample=(sum_sample)/(norm_sample)

#calculating global mean to get centroids of the clusters
    deviation=0.01
    mean_global = np.zeros([no_clusters,W])
    for i in range (0,no_clusters):
        random_sample=np.random.rand(1,W)-0.5
        random_norm=deviation*(np.random.rand())
        random_sample2=(random_norm*random_sample)/np.linalg.norm(random_sample)
        temp = mean_sample+random_sample2
        mean_global[i,:] = temp/np.linalg.norm(temp)

  
#calculating mean from spherical kmeans

    sum_sample3 = np.zeros([1,W])
    difference=1
    epsilon=0.01
    number=100
    iteration=0

    while (difference>epsilon):
        iteration=iteration+1
        number2=number
        #computing the nearest neighbour and assigning the points
        #E Step in EM algorithm
        mean_global2 = np.transpose(mean_global)
        value=np.dot(data, mean_global2)
        value_max=value.max(1)
        clusters=np.argmax(value,axis=1)
    
    #computing value of the function
        number=sum(value_max) 
    #print(number)

    #computing centroids for the clusters
    #M step in EM algorithm
        for i in range(0,no_clusters):
            sum_sample3=sum(data[np.where(clusters==i)])
            if(mpmath.norm(sum_sample3) != 0):
                temp2=sum_sample3/mpmath.norm(sum_sample3)#np.linalg.norm(sum_sample3) #Check this
                mean_global[i,:] = temp2

        
        difference=abs(number-number2)
    
    return clusters