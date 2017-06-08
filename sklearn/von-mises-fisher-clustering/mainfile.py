# -*- coding: utf-8 -*-
"""
PROJECT : Clustering data using von Mises-Fisher Distribution
Reference: http://jmlr.csail.mit.edu/papers/v6/banerjee05a.html
    
GUIDE : 
    Prof. Anand A Joshi - ajoshi@sipi.usc.edu
TEAM:
    Bhavana Ganesh - bhavanag@usc.edu
    Mahathi Vatsal Salopanthula - salopant@usc.edu
    Sayali Ghume - ghume@usc.edu

NOTE: Threading does not work on IPython console so this code would run only
      on a python console
      
Contact any of the members for queries and bugs.
"""

import vonmisesGenerate as vmg
import sphericalclustering as snn
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import itertools
import time
import multiprocessing as mp
from multiprocessing import Pool

def clust (params):
    """ Returns the cluster labels and exracted raw data
    
    This function extracts data from file whose path is the parameter to the 
    function.Samples from von mises-fisher distribution with each row in the 
    data set to be the mean vector and clusters the data using 
    Spherical Kmeans clustering
    
    The given function ''clust'' must have the following signature::
        
        output = clust(params) where, 
        
        fpath = ['.../data.npz']
        par2 = [no_clusters]
        params = zip(fpath, par2)
        
        where params is file address and no_of clusters zipped
        
    Dependencies::
        The function requires the below mentioned libraries to be installed
        1. Numpy
        2. Mpmath
        
        The function requires the module ''vonmisesGenerate.py'' and 
        ''sphericalclustering.py'' to be in Path
 
    """
    fileadd, no_clusters = params
    dat = np.load(fileadd)
    #brain = dat.f.data  #Comment while debugging
    brain2 = dat.f.data #Uncomment while debugging
    brain = brain2[1:3,:] #Uncomment while debugging
    del dat

    H = len(brain)
    W = len(brain[0])
    data = np.zeros([H,W])
    #generating random samples using the VMF distribution
    for i in range (0,H):
        mu = np.ndarray.tolist(brain[i,:])
        val = np.linalg.norm(mu)
        data[i,:] = vmg.randVMF(1,mu/val,1)

    #Clustering using the sphericalclustering function
    clusters = snn.sphericalknn(data,no_clusters)
    return (clusters,brain)



def clustplot (clusters,brain,no_clusters):
    """ Plots the clusters in 3D
    
    This function creates two plots one displaying the clusters and another 
    a 3D plot displaying the first 3 dimensions of the extarcted data with
    markers colored based on clustering
    
    The given function ''clustplot'' must have the following signature::
        
        clustplot(clusters,data)
        
    Dependencies::
        The function requires the below mentioned libraries to be installed
        1. Matplotlib
        2. mpl_toolkits
        3. itertools    
    
    """
    H = len(brain)
    W = len(brain[0])
    #displaying the clusters
    plt.plot(clusters, marker = '.', linestyle = '')
    #Plot first 3 colums of data
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors = itertools.cycle(["r", "b", "g", "k", "m", "y"])
    marker = itertools.cycle(["o", "*", "+", "x", "D", "s", "^"])

    for i in range (0,no_clusters):
        a = np.zeros([1,W])
        j=0
        for p in range (0,H):
            if clusters[p] == i:
                a = np.vstack((a,brain[p,:]))
                j=j+1
        a = np.delete(a,(0),axis=0)
        ax.scatter(a[:,147],a[:,148],a[:,149], c=next(colors), marker=next(marker))
    plt.show()
    
    
    
################################ MAIN FUNCTION#############################    

if __name__ == '__main__':
    """ Processes multiple data in parallel using multiprocessing library
    
    The main function can process multiple data at the same time by threading 
    the process and then plots the first dataset.
    
    """
    start_time = time.clock() #Start clock
    #Define Parameters
    fpath = ['D:/Codes/vonmises/data.npz']
    no_clusters = 200
    par2 = [no_clusters]
    params = zip(fpath, par2) #Zip the file address and number of clusters to be formed
    pool = Pool()
    '''
    Pool is a convenient means of parallelizing the execution of a function 
    across multiple input values, distributing the input data across processes 
    (data parallelism)
    Source: https://docs.python.org/2/library/multiprocessing.html
    '''
    #print(mp.current_process())
    out = pool.map(clust,params)
    
    #plotting the first dataset
    label = (out[0])[0]
    brain = (out[0])[1]
    clustplot(label,brain,no_clusters)
    
    print((time.clock() - start_time)/60,"minutes") #Print Elapsed time
    