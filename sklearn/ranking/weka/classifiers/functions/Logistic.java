/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    Logistic.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Optimization;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.RemoveUseless;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class for building and using a multinomial logistic regression model with a ridge estimator.<br/>
 * <br/>
 * There are some modifications, however, compared to the paper of leCessie and van Houwelingen(1992): <br/>
 * <br/>
 * If there are k classes for n instances with m attributes, the parameter matrix B to be calculated will be an m*(k-1) matrix.<br/>
 * <br/>
 * The probability for class j with the exception of the last class is<br/>
 * <br/>
 * Pj(Xi) = exp(XiBj)/((sum[j=1..(k-1)]exp(Xi*Bj))+1) <br/>
 * <br/>
 * The last class has probability<br/>
 * <br/>
 * 1-(sum[j=1..(k-1)]Pj(Xi)) <br/>
 * 	= 1/((sum[j=1..(k-1)]exp(Xi*Bj))+1)<br/>
 * <br/>
 * The (negative) multinomial log-likelihood is thus: <br/>
 * <br/>
 * L = -sum[i=1..n]{<br/>
 * 	sum[j=1..(k-1)](Yij * ln(Pj(Xi)))<br/>
 * 	+(1 - (sum[j=1..(k-1)]Yij)) <br/>
 * 	* ln(1 - sum[j=1..(k-1)]Pj(Xi))<br/>
 * 	} + ridge * (B^2)<br/>
 * <br/>
 * In order to find the matrix B for which L is minimised, a Quasi-Newton Method is used to search for the optimized values of the m*(k-1) variables.  Note that before we use the optimization procedure, we 'squeeze' the matrix B into a m*(k-1) vector.  For details of the optimization procedure, please check weka.core.Optimization class.<br/>
 * <br/>
 * Although original Logistic Regression does not deal with instance weights, we modify the algorithm a little bit to handle the instance weights.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * le Cessie, S., van Houwelingen, J.C. (1992). Ridge Estimators in Logistic Regression. Applied Statistics. 41(1):191-201.<br/>
 * <br/>
 * Note: Missing values are replaced using a ReplaceMissingValuesFilter, and nominal attributes are transformed into numeric attributes using a NominalToBinaryFilter.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{leCessie1992,
 *    author = {le Cessie, S. and van Houwelingen, J.C.},
 *    journal = {Applied Statistics},
 *    number = {1},
 *    pages = {191-201},
 *    title = {Ridge Estimators in Logistic Regression},
 *    volume = {41},
 *    year = {1992}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turn on debugging output.</pre>
 * 
 * <pre> -R &lt;ridge&gt;
 *  Set the ridge in the log-likelihood.</pre>
 * 
 * <pre> -M &lt;number&gt;
 *  Set the maximum number of iterations (default -1, until convergence).</pre>
 * 
 <!-- options-end -->
 *
 * @author Xin Xu (xx5@cs.waikato.ac.nz)
 * @version $Revision: 5928 $
 */
public class Logistic extends AbstractClassifier 
  implements OptionHandler, WeightedInstancesHandler, TechnicalInformationHandler {
  
  /** for serialization */
  static final long serialVersionUID = 3932117032546553727L;
  
  /** The coefficients (optimized parameters) of the model */
  protected double [][] m_Par;
    
  /** The data saved as a matrix */
  protected double [][] m_Data;
    
  /** The number of attributes in the model */
  protected int m_NumPredictors;
    
  /** The index of the class attribute */
  protected int m_ClassIndex;
    
  /** The number of the class labels */
  protected int m_NumClasses;
    
  /** The ridge parameter. */
  protected double m_Ridge = 1e-8;
    
  /** An attribute filter */
  private RemoveUseless m_AttFilter;
    
  /** The filter used to make attributes numeric. */
  private NominalToBinary m_NominalToBinary;
    
  /** The filter used to get rid of missing values. */
  private ReplaceMissingValues m_ReplaceMissingValues;
    
  /** Debugging output */
  protected boolean m_Debug;

  /** Log-likelihood of the searched model */
  protected double m_LL;
    
  /** The maximum number of iterations. */
  private int m_MaxIts = -1;

  private Instances m_structure;
    
  /**
   * Returns a string describing this classifier
   * @return a description of the classifier suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Class for building and using a multinomial logistic "
      +"regression model with a ridge estimator.\n\n"
      +"There are some modifications, however, compared to the paper of "
      +"leCessie and van Houwelingen(1992): \n\n" 
      +"If there are k classes for n instances with m attributes, the "
      +"parameter matrix B to be calculated will be an m*(k-1) matrix.\n\n"
      +"The probability for class j with the exception of the last class is\n\n"
      +"Pj(Xi) = exp(XiBj)/((sum[j=1..(k-1)]exp(Xi*Bj))+1) \n\n"
      +"The last class has probability\n\n"
      +"1-(sum[j=1..(k-1)]Pj(Xi)) \n\t= 1/((sum[j=1..(k-1)]exp(Xi*Bj))+1)\n\n"
      +"The (negative) multinomial log-likelihood is thus: \n\n"
      +"L = -sum[i=1..n]{\n\tsum[j=1..(k-1)](Yij * ln(Pj(Xi)))"
      +"\n\t+(1 - (sum[j=1..(k-1)]Yij)) \n\t* ln(1 - sum[j=1..(k-1)]Pj(Xi))"
      +"\n\t} + ridge * (B^2)\n\n"
      +"In order to find the matrix B for which L is minimised, a "
      +"Quasi-Newton Method is used to search for the optimized values of "
      +"the m*(k-1) variables.  Note that before we use the optimization "
      +"procedure, we 'squeeze' the matrix B into a m*(k-1) vector.  For "
      +"details of the optimization procedure, please check "
      +"weka.core.Optimization class.\n\n"
      +"Although original Logistic Regression does not deal with instance "
      +"weights, we modify the algorithm a little bit to handle the "
      +"instance weights.\n\n"
      +"For more information see:\n\n"
      + getTechnicalInformation().toString() + "\n\n"
      +"Note: Missing values are replaced using a ReplaceMissingValuesFilter, and "
      +"nominal attributes are transformed into numeric attributes using a "
      +"NominalToBinaryFilter.";
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "le Cessie, S. and van Houwelingen, J.C.");
    result.setValue(Field.YEAR, "1992");
    result.setValue(Field.TITLE, "Ridge Estimators in Logistic Regression");
    result.setValue(Field.JOURNAL, "Applied Statistics");
    result.setValue(Field.VOLUME, "41");
    result.setValue(Field.NUMBER, "1");
    result.setValue(Field.PAGES, "191-201");
    
    return result;
  }

  /**
   * Returns an enumeration describing the available options
   *
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {
    Vector newVector = new Vector(3);
    newVector.addElement(new Option("\tTurn on debugging output.",
				    "D", 0, "-D"));
    newVector.addElement(new Option("\tSet the ridge in the log-likelihood.",
				    "R", 1, "-R <ridge>"));
    newVector.addElement(new Option("\tSet the maximum number of iterations"+
				    " (default -1, until convergence).",
				    "M", 1, "-M <number>"));
    return newVector.elements();
  }
    
  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Turn on debugging output.</pre>
   * 
   * <pre> -R &lt;ridge&gt;
   *  Set the ridge in the log-likelihood.</pre>
   * 
   * <pre> -M &lt;number&gt;
   *  Set the maximum number of iterations (default -1, until convergence).</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    setDebug(Utils.getFlag('D', options));

    String ridgeString = Utils.getOption('R', options);
    if (ridgeString.length() != 0) 
      m_Ridge = Double.parseDouble(ridgeString);
    else 
      m_Ridge = 1.0e-8;
	
    String maxItsString = Utils.getOption('M', options);
    if (maxItsString.length() != 0) 
      m_MaxIts = Integer.parseInt(maxItsString);
    else 
      m_MaxIts = -1;
  }
    
  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {
	
    String [] options = new String [5];
    int current = 0;
	
    if (getDebug()) 
      options[current++] = "-D";
    options[current++] = "-R";
    options[current++] = ""+m_Ridge;	
    options[current++] = "-M";
    options[current++] = ""+m_MaxIts;
    while (current < options.length) 
      options[current++] = "";
    return options;
  }
   
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String debugTipText() {
    return "Output debug information to the console.";
  }

  /**
   * Sets whether debugging output will be printed.
   *
   * @param debug true if debugging output should be printed
   */
  public void setDebug(boolean debug) {
    m_Debug = debug;
  }
    
  /**
   * Gets whether debugging output will be printed.
   *
   * @return true if debugging output will be printed
   */
  public boolean getDebug() {
    return m_Debug;
  }      

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String ridgeTipText() {
    return "Set the Ridge value in the log-likelihood.";
  }

  /**
   * Sets the ridge in the log-likelihood.
   *
   * @param ridge the ridge
   */
  public void setRidge(double ridge) {
    m_Ridge = ridge;
  }
    
  /**
   * Gets the ridge in the log-likelihood.
   *
   * @return the ridge
   */
  public double getRidge() {
    return m_Ridge;
  }
   
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String maxItsTipText() {
    return "Maximum number of iterations to perform.";
  }

  /**
   * Get the value of MaxIts.
   *
   * @return Value of MaxIts.
   */
  public int getMaxIts() {
	
    return m_MaxIts;
  }
    
  /**
   * Set the value of MaxIts.
   *
   * @param newMaxIts Value to assign to MaxIts.
   */
  public void setMaxIts(int newMaxIts) {
	
    m_MaxIts = newMaxIts;
  }    
    
  private class OptEng extends Optimization{
    /** Weights of instances in the data */
    private double[] weights;

    /** Class labels of instances */
    private int[] cls;
	
    /** 
     * Set the weights of instances
     * @param w the weights to be set
     */ 
    public void setWeights(double[] w) {
      weights = w;
    }
	
    /** 
     * Set the class labels of instances
     * @param c the class labels to be set
     */ 
    public void setClassLabels(int[] c) {
      cls = c;
    }
	
    /** 
     * Evaluate objective function
     * @param x the current values of variables
     * @return the value of the objective function 
     */
    protected double objectiveFunction(double[] x){
      double nll = 0; // -LogLikelihood
      int dim = m_NumPredictors+1; // Number of variables per class
	    
      for(int i=0; i<cls.length; i++){ // ith instance

	double[] exp = new double[m_NumClasses-1];
	int index;
	for(int offset=0; offset<m_NumClasses-1; offset++){ 
	  index = offset * dim;
	  for(int j=0; j<dim; j++)
	    exp[offset] += m_Data[i][j]*x[index + j];
	}
	double max = exp[Utils.maxIndex(exp)];
	double denom = Math.exp(-max);
	double num;
	if (cls[i] == m_NumClasses - 1) { // Class of this instance
	  num = -max;
	} else {
	  num = exp[cls[i]] - max;
	}
	for(int offset=0; offset<m_NumClasses-1; offset++){
	  denom += Math.exp(exp[offset] - max);
	}
		
	nll -= weights[i]*(num - Math.log(denom)); // Weighted NLL
      }
	    
      // Ridge: note that intercepts NOT included
      for(int offset=0; offset<m_NumClasses-1; offset++){
	for(int r=1; r<dim; r++)
	  nll += m_Ridge*x[offset*dim+r]*x[offset*dim+r];
      }
	    
      return nll;
    }

    /** 
     * Evaluate Jacobian vector
     * @param x the current values of variables
     * @return the gradient vector 
     */
    protected double[] evaluateGradient(double[] x){
      double[] grad = new double[x.length];
      int dim = m_NumPredictors+1; // Number of variables per class
	    
      for(int i=0; i<cls.length; i++){ // ith instance
	double[] num=new double[m_NumClasses-1]; // numerator of [-log(1+sum(exp))]'
	int index;
	for(int offset=0; offset<m_NumClasses-1; offset++){ // Which part of x
	  double exp=0.0;
	  index = offset * dim;
	  for(int j=0; j<dim; j++)
	    exp += m_Data[i][j]*x[index + j];
	  num[offset] = exp;
	}

	double max = num[Utils.maxIndex(num)];
	double denom = Math.exp(-max); // Denominator of [-log(1+sum(exp))]'
	for(int offset=0; offset<m_NumClasses-1; offset++){
	  num[offset] = Math.exp(num[offset] - max);
	  denom += num[offset];
	}
	Utils.normalize(num, denom);
		
	// Update denominator of the gradient of -log(Posterior)
	double firstTerm;
	for(int offset=0; offset<m_NumClasses-1; offset++){ // Which part of x
	  index = offset * dim;
	  firstTerm = weights[i] * num[offset];
	  for(int q=0; q<dim; q++){
	    grad[index + q] += firstTerm * m_Data[i][q];
	  }
	}
		
	if(cls[i] != m_NumClasses-1){ // Not the last class
	  for(int p=0; p<dim; p++){
	    grad[cls[i]*dim+p] -= weights[i]*m_Data[i][p]; 
	  }
	}
      }
	    
      // Ridge: note that intercepts NOT included
      for(int offset=0; offset<m_NumClasses-1; offset++){
	for(int r=1; r<dim; r++)
	  grad[offset*dim+r] += 2*m_Ridge*x[offset*dim+r];
      }
	    
      return grad;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 5928 $");
    }
  }

  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.DATE_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }
    
  /**
   * Builds the classifier
   *
   * @param train the training data to be used for generating the
   * boosted classifier.
   * @throws Exception if the classifier could not be built successfully
   */
  public void buildClassifier(Instances train) throws Exception {
    // can classifier handle the data?
    getCapabilities().testWithFail(train);

    // remove instances with missing class
    train = new Instances(train);
    train.deleteWithMissingClass();
    
    // Replace missing values	
    m_ReplaceMissingValues = new ReplaceMissingValues();
    m_ReplaceMissingValues.setInputFormat(train);
    train = Filter.useFilter(train, m_ReplaceMissingValues);

    // Remove useless attributes
    m_AttFilter = new RemoveUseless();
    m_AttFilter.setInputFormat(train);
    train = Filter.useFilter(train, m_AttFilter);
	
    // Transform attributes
    m_NominalToBinary = new NominalToBinary();
    m_NominalToBinary.setInputFormat(train);
    train = Filter.useFilter(train, m_NominalToBinary);
    
    // Save the structure for printing the model
    m_structure = new Instances(train, 0);
	
    // Extract data
    m_ClassIndex = train.classIndex();
    m_NumClasses = train.numClasses();

    int nK = m_NumClasses - 1;                     // Only K-1 class labels needed 
    int nR = m_NumPredictors = train.numAttributes() - 1;
    int nC = train.numInstances();
	
    m_Data = new double[nC][nR + 1];               // Data values
    int [] Y  = new int[nC];                       // Class labels
    double [] xMean= new double[nR + 1];           // Attribute means
    double [] xSD  = new double[nR + 1];           // Attribute stddev's
    double [] sY = new double[nK + 1];             // Number of classes
    double [] weights = new double[nC];            // Weights of instances
    double totWeights = 0;                         // Total weights of the instances
    m_Par = new double[nR + 1][nK];                // Optimized parameter values
	
    if (m_Debug) {
      System.out.println("Extracting data...");
    }
	
    for (int i = 0; i < nC; i++) {
      // initialize X[][]
      Instance current = train.instance(i);
      Y[i] = (int)current.classValue();  // Class value starts from 0
      weights[i] = current.weight();     // Dealing with weights
      totWeights += weights[i];
	    
      m_Data[i][0] = 1;
      int j = 1;
      for (int k = 0; k <= nR; k++) {
	if (k != m_ClassIndex) {
	  double x = current.value(k);
	  m_Data[i][j] = x;
	  xMean[j] += weights[i]*x;
	  xSD[j] += weights[i]*x*x;
	  j++;
	}
      }
	    
      // Class count
      sY[Y[i]]++;	
    }
	
    if((totWeights <= 1) && (nC > 1))
      throw new Exception("Sum of weights of instances less than 1, please reweight!");

    xMean[0] = 0; xSD[0] = 1;
    for (int j = 1; j <= nR; j++) {
      xMean[j] = xMean[j] / totWeights;
      if(totWeights > 1)
	xSD[j] = Math.sqrt(Math.abs(xSD[j] - totWeights*xMean[j]*xMean[j])/(totWeights-1));
      else
	xSD[j] = 0;
    }

    if (m_Debug) {	    
      // Output stats about input data
      System.out.println("Descriptives...");
      for (int m = 0; m <= nK; m++)
	System.out.println(sY[m] + " cases have class " + m);
      System.out.println("\n Variable     Avg       SD    ");
      for (int j = 1; j <= nR; j++) 
	System.out.println(Utils.doubleToString(j,8,4) 
			   + Utils.doubleToString(xMean[j], 10, 4) 
			   + Utils.doubleToString(xSD[j], 10, 4)
			   );
    }
	
    // Normalise input data 
    for (int i = 0; i < nC; i++) {
      for (int j = 0; j <= nR; j++) {
	if (xSD[j] != 0) {
	  m_Data[i][j] = (m_Data[i][j] - xMean[j]) / xSD[j];
	}
      }
    }
	
    if (m_Debug) {
      System.out.println("\nIteration History..." );
    }
	
    double x[] = new double[(nR+1)*nK];
    double[][] b = new double[2][x.length]; // Boundary constraints, N/A here

    // Initialize
    for(int p=0; p<nK; p++){
      int offset=p*(nR+1);	 
      x[offset] =  Math.log(sY[p]+1.0) - Math.log(sY[nK]+1.0); // Null model
      b[0][offset] = Double.NaN;
      b[1][offset] = Double.NaN;   
      for (int q=1; q <= nR; q++){
	x[offset+q] = 0.0;		
	b[0][offset+q] = Double.NaN;
	b[1][offset+q] = Double.NaN;
      }	
    }
	
    OptEng opt = new OptEng();	
    opt.setDebug(m_Debug);
    opt.setWeights(weights);
    opt.setClassLabels(Y);

    if(m_MaxIts == -1){  // Search until convergence
      x = opt.findArgmin(x, b);
      while(x==null){
	x = opt.getVarbValues();
	if (m_Debug)
	  System.out.println("200 iterations finished, not enough!");
	x = opt.findArgmin(x, b);
      }
      if (m_Debug)
	System.out.println(" -------------<Converged>--------------");
    }
    else{
      opt.setMaxIteration(m_MaxIts);
      x = opt.findArgmin(x, b);
      if(x==null) // Not enough, but use the current value
	x = opt.getVarbValues();
    }
	
    m_LL = -opt.getMinFunction(); // Log-likelihood

    // Don't need data matrix anymore
    m_Data = null;
	    
    // Convert coefficients back to non-normalized attribute units
    for(int i=0; i < nK; i++){
      m_Par[0][i] = x[i*(nR+1)];
      for(int j = 1; j <= nR; j++) {
	m_Par[j][i] = x[i*(nR+1)+j];
	if (xSD[j] != 0) {
	  m_Par[j][i] /= xSD[j];
	  m_Par[0][i] -= m_Par[j][i] * xMean[j];
	}
      }
    }
  }		
    
  /**
   * Computes the distribution for a given instance
   *
   * @param instance the instance for which distribution is computed
   * @return the distribution
   * @throws Exception if the distribution can't be computed successfully
   */
  public double [] distributionForInstance(Instance instance) 
    throws Exception {
	
    m_ReplaceMissingValues.input(instance);
    instance = m_ReplaceMissingValues.output();
    m_AttFilter.input(instance);
    instance = m_AttFilter.output();
    m_NominalToBinary.input(instance);
    instance = m_NominalToBinary.output();
	
    // Extract the predictor columns into an array
    double [] instDat = new double [m_NumPredictors + 1];
    int j = 1;
    instDat[0] = 1;
    for (int k = 0; k <= m_NumPredictors; k++) {
      if (k != m_ClassIndex) {
	instDat[j++] = instance.value(k);
      }
    }
	
    double [] distribution = evaluateProbability(instDat);
    return distribution;
  }

  /**
   * Compute the posterior distribution using optimized parameter values
   * and the testing instance.
   * @param data the testing instance
   * @return the posterior probability distribution
   */ 
  private double[] evaluateProbability(double[] data){
    double[] prob = new double[m_NumClasses],
      v = new double[m_NumClasses];

    // Log-posterior before normalizing
    for(int j = 0; j < m_NumClasses-1; j++){
      for(int k = 0; k <= m_NumPredictors; k++){
	v[j] += m_Par[k][j] * data[k];
      }
    }
    v[m_NumClasses-1] = 0;
	
    // Do so to avoid scaling problems
    for(int m=0; m < m_NumClasses; m++){
      double sum = 0;
      for(int n=0; n < m_NumClasses-1; n++)
	sum += Math.exp(v[n] - v[m]);
      prob[m] = 1 / (sum + Math.exp(-v[m]));
    }
	
    return prob;
  } 

  /**
   * Returns the coefficients for this logistic model.
   * The first dimension indexes the attributes, and
   * the second the classes.
   * 
   * @return the coefficients for this logistic model
   */
  public double [][] coefficients() {
    return m_Par;
  }
    
  /**
   * Gets a string describing the classifier.
   *
   * @return a string describing the classifer built.
   */
  public String toString() {
    StringBuffer temp = new StringBuffer();

    String result = "";
    temp.append("Logistic Regression with ridge parameter of " + m_Ridge);
    if (m_Par == null) {
      return result + ": No model built yet.";
    }

    // find longest attribute name
    int attLength = 0;
    for (int i = 0; i < m_structure.numAttributes(); i++) {
      if (i != m_structure.classIndex() && 
          m_structure.attribute(i).name().length() > attLength) {
        attLength = m_structure.attribute(i).name().length();
      }
    }

    if ("Intercept".length() > attLength) {
      attLength = "Intercept".length();
    }

    if ("Variable".length() > attLength) {
      attLength = "Variable".length();
    }
    attLength += 2;

    int colWidth = 0;
    // check length of class names
    for (int i = 0; i < m_structure.classAttribute().numValues() - 1; i++) {
      if (m_structure.classAttribute().value(i).length() > colWidth) {
        colWidth = m_structure.classAttribute().value(i).length();
      }
    }

    // check against coefficients and odds ratios
    for (int j = 1; j <= m_NumPredictors; j++) {
      for (int k = 0; k < m_NumClasses - 1; k++) {
        if (Utils.doubleToString(m_Par[j][k], 12, 4).trim().length() > colWidth) {
          colWidth = Utils.doubleToString(m_Par[j][k], 12, 4).trim().length();
        }
        double ORc = Math.exp(m_Par[j][k]);
	String t = " " + ((ORc > 1e10) ?  "" + ORc : Utils.doubleToString(ORc, 12, 4));
        if (t.trim().length() > colWidth) {
          colWidth = t.trim().length();
        }
      }
    }

    if ("Class".length() > colWidth) {
      colWidth = "Class".length();
    }
    colWidth += 2;
    
    
    temp.append("\nCoefficients...\n");
    temp.append(Utils.padLeft(" ", attLength) + Utils.padLeft("Class", colWidth) + "\n");
    temp.append(Utils.padRight("Variable", attLength));

    for (int i = 0; i < m_NumClasses - 1; i++) {
      String className = m_structure.classAttribute().value(i);
      temp.append(Utils.padLeft(className, colWidth));
    }
    temp.append("\n");
    int separatorL = attLength + ((m_NumClasses - 1) * colWidth);
    for (int i = 0; i < separatorL; i++) {
      temp.append("=");
    }
    temp.append("\n");
                
    int j = 1;
    for (int i = 0; i < m_structure.numAttributes(); i++) {
      if (i != m_structure.classIndex()) {
        temp.append(Utils.padRight(m_structure.attribute(i).name(), attLength));
        for (int k = 0; k < m_NumClasses-1; k++) {
          temp.append(Utils.padLeft(Utils.doubleToString(m_Par[j][k], 12, 4).trim(), colWidth));
        }
        temp.append("\n");
        j++;
      }
    }
	
    temp.append(Utils.padRight("Intercept", attLength));
    for (int k = 0; k < m_NumClasses-1; k++) {
      temp.append(Utils.padLeft(Utils.doubleToString(m_Par[0][k], 10, 4).trim(), colWidth)); 
    }
    temp.append("\n");
	
    temp.append("\n\nOdds Ratios...\n");
    temp.append(Utils.padLeft(" ", attLength) + Utils.padLeft("Class", colWidth) + "\n");
    temp.append(Utils.padRight("Variable", attLength));

    for (int i = 0; i < m_NumClasses - 1; i++) {
      String className = m_structure.classAttribute().value(i);
      temp.append(Utils.padLeft(className, colWidth));
    }
    temp.append("\n");
    for (int i = 0; i < separatorL; i++) {
      temp.append("=");
    }
    temp.append("\n");

    j = 1;
    for (int i = 0; i < m_structure.numAttributes(); i++) {
      if (i != m_structure.classIndex()) {
        temp.append(Utils.padRight(m_structure.attribute(i).name(), attLength));
        for (int k = 0; k < m_NumClasses-1; k++) {
          double ORc = Math.exp(m_Par[j][k]);
          String ORs = " " + ((ORc > 1e10) ?  "" + ORc : Utils.doubleToString(ORc, 12, 4));
          temp.append(Utils.padLeft(ORs.trim(), colWidth));
        }
        temp.append("\n");
        j++;
      }
    }

    return temp.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5928 $");
  }
    
  /**
   * Main method for testing this class.
   *
   * @param argv should contain the command line arguments to the
   * scheme (see Evaluation)
   */
  public static void main(String [] argv) {
    runClassifier(new Logistic(), argv);
  }
}
