"""
Polygonal Coordinate System (PCS)
"""

# author: Caio Flexa <caio.rodrigues@icen.ufpa.br>
# author: Walisson Gomes <walisson.gomes@icen.ufpa.br>
# author: Isaac Elgrably <isaac.elgrably@icen.ufpa.br>
# License: BSD 3 clause (C) 2019

from ..base import BaseEstimator

import numpy as np

class PCS( BaseEstimator ):
    """PCS Embedding
    
    This code presents a geometric approach for data embedding called Polygonal 
    Coordinate System (PCS), which is able to efficiently represent 
    multidimensional data into a 2D plane by preserving the global structure of 
    the data. For this purpose, data are represented across a regular polygon 
    or interface between the high dimensionality and the two-dimensional data. 
    PCS can properly handle massive amounts of data (Big Data) by adopting an 
    incremental and linear-time complexity dimensionality reduction.
    
    Parameters
    ----------
    u : float, optional (default = None)
        Size of each side in the regular polygon used to project the
        M-dimensional data into 2-dimensional data. If no u is specified, it is
        calculated based on the range of values in which the data is presented.
        
    origin : array of float of size 2, optional (default = [ 0, 0])
        represents where the polygon is located. If not specified, the first
        vertex of the polygon is in [0, 0].

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.manifold import PCS
    >>> target = load_iris().target
    >>> X = load_iris().data
    >>> X.shape
    (150, 4)
    >>>
    >>> embedding = PCS()
    >>> Y = embedding.fit_transform( X )
    >>> Y.shape
    (150, 2)
    """

    def __init__( self, u = None, origin = [ 0, 0 ] ):
        self.u = u
        self.origin = origin
    
    def _fit_transform( self, X ):
        """Compute the embedding vectors for data X

        Parameters
        ----------
        X : array, shape = (n_samples, n_features), in the form of a
            numpy array.

        Returns
        -------
        X_new : array, shape ( n_samples, 2 )
            Embedding of the training data in a low-dimensional space.
        """
        
        # (Part 1) - Data normalization
        # The NxS input vector X can be in three states:
        #   1. Normalized on [0,1]
        #   2. Normalized, but on a different range of values
        #   3. Non-normalized on [0,1]
        # Since all columns must be in the same interval, we need to normalize
        # the data, if it is not. u is defined accordingly.
        
        N, S = X.shape
        
        self.N = N
        self.S = S
        
        # Here, u is determined, if not specified
        if self.u is None:
            flag = True
            Xmin = X.min( axis = 0 )
            Xmax = X.max( axis = 0 )
            
            for j in range( S - 1 ):
                if Xmin[ j ] != Xmin[ j + 1 ] or Xmax[ j ] != Xmax[ j + 1 ]:
                    flag = False
                    break
            if flag:
                self.u = Xmax[ 0 ] - Xmin[ 0 ]
            else:
                X = ( X - X.min( axis = 0 ) ) /\
                                      ( X.max( axis = 0 ) - X.min( axis = 0 ) )
                self.u = 1
        else:
            X = ( X - X.min( axis = 0 ) ) /\
                                      ( X.max( axis = 0 ) - X.min( axis = 0 ) )
            X = self.u * X
        
        u = self.u
        auxpair = self.origin
        
        # (Part 2) - Data projection
        # Since data are in the same interval, it is necessary to rotate and to
        # translate the metavectors to their final position in the polygons sides.
        
        # Interior angle of the polygon, auxiiar variable, and omega.
        alpha = 180 * ( S - 2 ) / S
        tau = ( 180 - alpha ) / ( 2 - S % 2 )
        omega = X.min()
		
        # Meta-vector angles and metavector coordinates.
        angs = np.zeros( S )
        pairs = np.zeros( [ S + 1, 2 ] )
        
        # Iteratively, find the starting point of each metavector.
        for j in range( 1, S + 1 ):
            beta = 180 * ( tau + j * ( S - 2 ) ) / S % 180
			
            if beta == 0:
                beta = 180
							
            angs[ j - 1 ] = np.deg2rad( beta )
			
            pairs[ j - 1, : ] = auxpair;
			
            # Coordinate pairs update.
            if j < ( S / 2 + 1 ):
                auxpair = [ u * np.cos( angs[ j - 1 ] ) + auxpair[ 0 ],
						           u * np.sin( angs[ j - 1 ] ) + auxpair[ 1 ] ]
            else:
                auxpair = [ u * -np.cos( angs[ j - 1 ] ) + auxpair[ 0 ],
						          u * -np.sin( angs[ j - 1 ] ) + auxpair[ 1 ] ]
		
        pairs[ -1, : ] = pairs[ 0, : ];
		
        # Every column is projected by a rotation + translation procedure.
        W1 = np.zeros( [ N, S ] );
        W2 = np.zeros( [ N, S ] );
        for i in range( S ):
            Mr = np.array( [ [ np.cos( angs[ i ] ), -np.sin( angs[ i ] ) ],
                              [ np.sin( angs[ i ] ), np.cos( angs[ i ] ) ] ] );
			
            bidimX = np.column_stack( ( X[ :, i ] - omega,
                                                   np.zeros( X.shape[ 0 ] ) ) )
					
            if i < ( S / 2 ):
                auxW = np.matmul( bidimX, Mr.T ) + pairs[ i, : ];
            else:
                auxW = np.matmul( bidimX, Mr.T ) + pairs[ i + 1, : ];
			
            W1[ :, i ] = auxW[ :, 0 ]
            W2[ :, i ] = auxW[ :, 1 ]
        
        # (Part 3) - 2D mapping
        # Y represents the final mapping for each sample in the 2D plane.
        # self.projection stores the data position in the interspace.
        Y = np.array( [ W1.mean( axis = 1 ), W2.mean( axis = 1 ) ] )
        
        self.projection_ = ( W1, W2 )
        self.embedding_ = Y.T

    def fit( self, X, y = None ):
        """Compute the embedding vectors for data X

        Parameters
        ----------
        X : array, shape = (n_samples, n_features), in the form of a
            numpy array.

        y : Ignored

        Returns
        -------
        self : returns an instance of self.
        """
        
        self._fit_transform( X )
        return self

    def fit_transform( self, X, y = None ):
        """Compute the embedding vectors for data X

        Parameters
        ----------
        X : array, shape = (n_samples, n_features), in the form of a
            numpy array.

        y : Ignored

        Returns
        -------
        X_new : array, shape (n_samples, 2)
            Embedding of the training data in a low-dimensional space.
        """
        self._fit_transform( X )
        return self.embedding_
    