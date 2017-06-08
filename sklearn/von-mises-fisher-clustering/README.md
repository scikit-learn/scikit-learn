# sphere-kit
Computational tools for processing and analysis on hypersphere

von-mises-fischer-clustering functions :

vonmisesGenerate.py
   This module contains all the functions that are used to sample from the VMF distribution.

sphericalclustering.py
   This module clusters the data using the VMF distribution.

- clust:
    Returns the cluster labels and exracted raw data
    
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

- clustplot:
    Plots the clusters in 3D
    
    This function creates two plots one displaying the clusters and another 
    a 3D plot displaying the first 3 dimensions of the extarcted data with
    markers colored based on clustering
    
    The given function ''clustplot'' must have the following signature::
        
        clustplot(clusters,data)
        
    Dependencies:
        The function requires the below mentioned libraries to be installed
        1. Matplotlib
        2. mpl_toolkits
        3. itertools  

- randVMF:
  
    This function generates samples from the Von Mises Fisher Distribution
    The inputs are:
        N=Number of samples
        mu=Mean of the VMF distrivution. It should always be 1
        k=Kappa parameter
    The output is N x dimension(mu) matrix of samples generated from VMF distribution.

- randVMFMeanDir:
  
    This function generates random samples from VMF distribution using rejection sampling.
    The inputs are:
        N=Number of samples
        k=Kappa parameter  
        p=Dimension of VMF distribution
    The output is N x 1 vector of random samples from VMF tangent distribution

- VMFMeanDirDensity: 
  
    This function is the tangent direction density of VMF distribution.
    The inputs are:
        x=Tangent value between [-1 1]
        k=kappa parameter
        p=Dimension of the VMF distribution
    The output is y=density of the VMF tangent density

- randUniformSphere:
    
    This function generate the random samples uniformly on a d-dimensional sphere
    The inputs are :
        N=Number of samples
        p=Dimension of the VMF distribution
    The output is N x d dim matrix which are the N random samples generated 
    on the unit p-dimensional sphere
    
- nullspace:

    This function generates the null space of a matrix
    The source from where it is taken is:  
        http://scipy-cookbook.readthedocs.io/items/RankNullspace.html


  
