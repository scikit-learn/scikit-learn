import numpy as np
from exceptions import ValueError
from scipy import stats

class PLS(object):
    """
    Projection to Latent Structures (aka Partial Least Squares)
    """
    
    def __init__(self, n_comp=2, tol=1e-6, scale=True):
        """
        Parameters
        ----------
        n_comp: int, optional (default=2)
            number of components to use
                
        tol: float, optional (default=1e-6)
            Convergence tolerance for NIPALS
        
        scale: bool, optional (default=True)
            If true, scale all variables to unit variance
        
        Attributes
        ----------
        Xmean: numpy array, shape = [1, n_variables]
            The means of each variable in the training data
        
        Ymean: numpy array, shape = [1, n_responses]
            The means of each response in the training data
        
        Xstd: numpy array, shape = [1, n_variables]
            The standard deviation of each of the variables in the training data
            
        Ystd: numpy array, shape = [1, n_responses]
            The standard deviation of each of the responses in the training data
        
        P: numpy array, shape = [n_variables, n_comp]
            The X space loadings
            
        Q: numpy array, shape = [n_responses, n_comp]
            The Y space loadings (aka C)
        
        W: numpy array, shape = [n_variables, n_comp]
            The deflated weight vector
            
        Ws: numpy array, shape = [n_variables, n_comp]
            The undeflated weight vector (aka W* or R).  Used for interpretation of weightings
        
        Examples
        --------
        >>> model = PLS()
        >>> model.fit()
        []
        
        """
        assert n_comp>0, ValueError('cannot less than 1 component')
        self.n_comp = n_comp
        self.tol = tol
        self.scale = scale
    
    def __calc_new_u(self,u,X,Y):
        """Helper function to calculate updated u vector for NIPALS"""    
        w = X.T*u/(u.T*u)
        w = w/np.linalg.norm(w)
        t = X*w/(w.T*w)
        q = Y.T*t/(t.T*t)
        unew =  Y*q/(q.T*q)
        return w,t,q,unew  

    def fit(self, x, y): 
        """
        Fit the PLS model according to the given training data and parameters
        
        Parameters
        ----------
        x: numpy array, shape=[n_samples, n_variables]
            Training vector, where n_samples is the number of samples and 
            n_variables is the number of variables.
        y: numpy array, shape=[n_samples, n_responses]
            Target values (classification uses binary variables, real numbers
            in regression)
        
        """

        self.n_samples = x.shape[0]
        self.n_variables = x.shape[1]
        self.n_responses = y.shape[1]
        
        #assert self.n_comp<=self.n_variables, ValueError('Number of components > variables')
        
        A = self.n_comp
        #Centre Data
        Xmean = x.mean(axis=0)
        Xc = x - Xmean
        Ymean = y.mean(axis=0)
        Yc = y - Ymean
        
        #Scale to unit variance
        #TODO: Custom scaling
        if self.scale:
            Xstd = x.std(axis=0)
            Ystd = y.std(axis=0)
        else:
            Xstd = 1.
            Ystd = 1.
            
        Xs = Xc/Xstd
        Ys = Yc/Ystd   
        
        self.Xmean = Xmean
        self.Ymean = Ymean
        self.Xstd = Xstd
        self.Ystd = Ystd
        self.intercept_ = self.Ymean
        
        #Convert data to matrices to simplify matrix operations
        X = np.matrix(Xs)
        Y = np.matrix(Ys)
    
        #Calculate metrics
        SSX = np.zeros((A+1,))
        SSY = np.zeros((A+1,))
        SSX[0] = np.array(X.T*X)[0][0]
        SSY[0] = np.array(Y.T*Y)[0][0]
    
        #Pick an initial guess for u.
        #Currently selects first column in Y.
        u = np.matrix(Y[:,0])
        
        for i in range(A):
            
            #Initialize convergence
            conv = 1 
            
            while conv > self.tol:
                w,t,q,unew=self.__calc_new_u(u,X,Y)
                conv = np.linalg.norm(u-unew)
                u = unew
    
            #Calculate loadings for current component
            p = X.T*t/(t.T*t)
            
            if i == 0:
                W=w
                T=t
                Q=q
                P=p
                U=u
                
            else:
                P = np.concatenate((P,p),axis=1)
                T = np.concatenate((T,t),axis=1)
                W = np.concatenate((W,w),axis=1)
                Q = np.concatenate((Q,q),axis=1)
                U = np.concatenate((U,u),axis=1)
            
            #Deflation of X and Y
            X_hat = t*p.T
                
            E = X-X_hat
            X = E
            Y_hat = t*q.T
            F = Y-Y_hat
            Y = F
    
            ##Calculation of the sum of squares
            SSX[i+1]= np.array(X.T*X)[0][0]
            SSY[i+1]= np.array(Y.T*Y)[0][0]
            
            #Reset u for next component
            u = np.matrix(Y[:,0])
        
        #Calculate coefficients
        B = W*np.linalg.inv(P.T*W)*Q.T
        Ws = W*np.linalg.inv(P.T*W) 
        
        self.coef_ = np.array(B)
        self.P = np.array(P)
        self.Q = np.array(Q)
        self.W = np.array(W)
        self.Ws = np.array(Ws)
        self.SSX = SSX
        self.SSY = SSY
        
    def predict(self, x):
        """
        The function predicts the responses using set of test vectors, x.
          
        Parameters
        ----------
        
        x: numpy array, shape = [>=1, n_variables]
            The test vector
            
        Note
        ----
        Alternatively, a prediction may be made using:
            
            Ynew =  np.dot(X,self.coef_)*self.Ystd+self.intercept_
            
        """
        T = self._scores(x)
        return np.dot(T,self.Q.T)
    
    def transform(self,x):
        """
        Transforms the data to match the transformations on the training data.
        
        Parameters
        ----------
        x: numpy array, shape = [>=1, n_variables]
            The raw input matrix
        
        Returns
        -------
        X: numpy array, shape = [same as x, n_variables]
            The transformed input matrix
        
        """
        return (x - self.Xmean)/self.Xstd 
    
    def _scores(self,x):
        """
        Function to return the X scores
        
        Parameters
        ----------
        
        x: numpy array, shape = [>=1, n_variables]
            The test vector
        
        Returns
        -------
        
        T: numpy array, shape = [n_samples in x, n_comp]
            The X scores
        """
        
        X = self.transform(x)#Scale X values
        T = np.dot(X,self.Ws) #New Scores
        return T
    
    def _statistics(self,x):
        """
        Function to return the Squared Prediction Error in X (SPEx)
        and Hotelling's T-squared values.
        
        Parameters
        ----------
        
        x: numpy array, shape = [>=1, n_variables]
            The test vector
        
        Returns
        -------
        
        SPEx: numpy 1D array, shape = [n_samples in x]
            The Squared Prediction Error in X
        
        HotT2: numpy 1D array, shape = [n_samples in x]
            Hotelling's T-squared values
        """
        
        X = self.transform(x)
        Xnew = self.predict(X)
        Ex = X - Xnew
        T = self._scores(x)
        SPEx = np.sum(np.array(Ex)**2,axis=1)
        eig = np.var(T,axis=0)
        HotT2 = np.sum(np.array(T)**2/np.array(eig),axis=1)
        return SPEx, HotT2
    
    def _VIP(self):
        """
        Determine the Variables Important to Projection (VIP)
        
        Returns
        -------
        VIP: numpy 1D array, [n_variables]
        
        """
        #number of variables
        k = self.W.shape[0]
        #number of components
        a = self.n_comp
        W = self.W
        SSY = self.SSY
        return np.sqrt(np.sum(W**2*(SSY[:a]-SSY[1:a+1]),axis=1)*k/(SSY[0]-SSY[a])) 
    
def statistics_limits(model, x):
    """
    Calculate 95% and 99% limits on the Squared Prediction Error
    and the T squared values.
    
    Parameters:
    
    model: PLS.PLS
    
    x: numpy.array, shape = [>=1, model.n_variables]
    
    """
    out = model._statistics(x)
    SPEx = out[0]
    n = SPEx.shape[0] 
    n_comp = model.n_comp
    
    #Limits for Squared Prediction Error in  X
    v = SPEx.var()
    m = SPEx.mean()
    SPE95 = v/(2*m)*stats.chi2.ppf(0.95, 2*m**2/v)
    SPE99 = v/(2*m)*stats.chi2.ppf(0.99, 2*m**2/v)
    
    #Limits for Hotelling's T Squared
    T2_mult = (n-1)*(n+1)*n_comp/((n-n_comp)*n)
    df1 = n_comp    
    df2 = n-n_comp
    T295 = T2_mult*stats.f.ppf(0.95,df1,df2)
    T299 = T2_mult*stats.f.ppf(0.99,df1,df2)
    
    return SPE95, SPE99, T295, T299

def plot_SPE_T2(plt,model,X,obs):
    """
    Plot Squared prediction Error (SPE) and T squared (T2) values.
    
    Parameters
    ----------
    plt: matplotlib.pyplot object
    
    model: PLS.PLS object
    
    X: numpy array, shape = [>=1, model.n_variables]
        Vector for which to calculate the score values in the model
    
    obs: list with IDs for each observation in X

    Returns
    -------
    
    fig: matplotlib.figure.Figure   
    
    """
    SPE95, SPE99, T295, T299 = statistics_limits(model, X)
    SPEx, HotT2 = model._statistics(X)
    
    n =  model.n_samples
    
    fig = plt.figure()
    axSPE = fig.add_subplot(211)
    axSPE.plot(SPEx,marker='o',markerfacecolor='None',c='k', picker=5)
    for x,y,s in zip(range(n),SPEx,obs):
        if y>SPE95:
            axSPE.text(x,y,s, fontsize='x-small')
        
    axSPE.axhline(y=SPE95,color='g',ls='dashed',lw=1, label='95%')
    axSPE.axhline(y=SPE99,color='r',ls='solid',lw=1, label='99%')
    axSPE.set_xlabel('Observation')
    axSPE.set_ylabel('SPEx')
    leg = axSPE.legend(loc='best')
    leg.draggable()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='small')
    axHotT2 = fig.add_subplot(212)
    axHotT2.plot(HotT2,marker='o',markerfacecolor='None',c='k', picker=5)
    for x,y,s in zip(range(n),HotT2,obs):
        if y>T295:
            axHotT2.text(x,y,s, fontsize='x-small')

    axHotT2.axhline(y=T295,color='g',ls='dashed',lw=1, label='95%')
    axHotT2.axhline(y=T299,color='r',ls='solid',lw=1, label='99%')
    axHotT2.set_xlabel('Observation')
    axHotT2.set_ylabel(r'$Hotelling\'s$ $T^2$')
    leg = axHotT2.legend(loc='best')
    leg.draggable()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='small')
    return fig

def plot_VIP(plt, model, var):
    """
    Plot Variables Important to Projection (VIP).
    
    Parameters
    ----------
    plt: matplotlib.pyplot object
    
    model: PLS.PLS object
    
    var: list with variable names for each observation in X
    
    Returns
    -------
    
    fig: matplotlib.figure.Figure
    """
        
    #Plotting VIPs
    vip = model._VIP()
    k = model.n_variables
    fig = plt.figure(figsize=(10,4))
    fig.subplots_adjust(bottom=0.3)
    ivip = np.argsort(vip, axis=0)
    VIP_labels = var[ivip[::-1]] 
    ind = np.arange(1,k+1)
    width = 0.35
    ax = fig.add_subplot(111)
    ax.bar(ind, vip[ivip[::-1]], color='b')
    ax.set_xticks(ind+width)
    ax.set_xticklabels(VIP_labels,size='small',rotation='90')
    
    return fig

def plot_WQbar(plt, model, var, var_ind=[], index=0):
    """
    Bar Chart of W*,Q
    
    Parameters
    ----------
    plt: matplotlib.pyplot object
    
    model: PLS.PLS object
    
    var: list with variable names for each observation in X
    
    var_ind: list, optional (default=[])
        integers representing the indices of the variables to use
        if empty, use all variables 
        
    index: int, optional (default=0)
        The index of the component you wish to plot
        
    Returns
    -------
    
    fig: matplotlib.figure.Figure
    """
    
    n_res = model.n_responses
    
    if var_ind == []:
        ivar = np.arange(model.n_variables)
    else:
        ivar = np.array(var_ind, dtype='uint16')
    
    Ws = model.Ws[ivar,:]
    Q = model.Q
    var = var[ivar]
    n_var = var.shape[0]
    
    fig = plt.figure()
    fig.subplots_adjust(bottom=0.4)
    ax = fig.add_subplot(111)
    ind = np.arange(n_var+n_res)
    width = 0.35
    Y_label = np.array(["Y"])
    labels = np.concatenate((Y_label,var))
    ax.bar(ind, np.concatenate((Q[:,index],Ws[:,index])), color=['r']*n_res+['b']*n_var)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(labels,size='small',rotation='90')
    return fig

def plot_WQ(plt, model, var, var_ind=[],ind1=0, ind2=1):
    """
    Scatter plot of Ws*Q for any two components in W and Q
    
    Parameters
    ----------
    plt: matplotlib.pyplot
    
    model: PLS.PLS object
    
    var: list
        names of the variables in the model
    
    var_ind: list, optional(default=[])
        integers representing the indices of the variables to use
        if empty, use all variables 
        
    ind1: int, optional (default=0)
        The first component index.
    
    ind2: int, optional (default=1)
        The second component index.
    
    Returns
    -------
    
    fig: matplotlib.figure.Figure
    """
    

    n_var = model.n_variables
    
    if var_ind == []:
        ivar = np.arange(n_var)
    else:
        ivar = np.array(var_ind, dtype='uint16')
    
    Ws = model.Ws[ivar,:]
    Q = model.Q
    var = var[ivar]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(Ws[:,ind1],Ws[:,ind2],ls='None',marker='o',c='k',markerfacecolor='k',markeredgecolor='k')
    for j,var in enumerate(var):
        ax.text(Ws[j,ind1],Ws[j,ind2],var,fontsize='x-small',clip_on=True)
    ax.plot(Q[:,ind1],Q[:,ind2],ls='None',marker='s',c='r',markerfacecolor='r',markeredgecolor='r')
    ax.text(Q[:,ind1],Q[:,ind2],"Response",fontsize='x-small',clip_on=True)
    ax.set_xlabel("Component %d"%ind1)
    ax.set_ylabel("Component %d"%ind2)
    ax.grid()
    ax.axhline(y=0)
    ax.axvline(x=0)   
    
    return fig 