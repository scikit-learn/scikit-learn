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
 *    SMO.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.functions.supportVector.Kernel;
import weka.classifiers.functions.supportVector.PolyKernel;
import weka.classifiers.functions.supportVector.SMOset;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.SerializedObject;
import weka.core.Tag;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.Normalize;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;
import weka.filters.unsupervised.attribute.Standardize;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Implements John Platt's sequential minimal optimization algorithm for training a support vector classifier.<br/>
 * <br/>
 * This implementation globally replaces all missing values and transforms nominal attributes into binary ones. It also normalizes all attributes by default. (In that case the coefficients in the output are based on the normalized data, not the original data --- this is important for interpreting the classifier.)<br/>
 * <br/>
 * Multi-class problems are solved using pairwise classification (1-vs-1 and if logistic models are built pairwise coupling according to Hastie and Tibshirani, 1998).<br/>
 * <br/>
 * To obtain proper probability estimates, use the option that fits logistic regression models to the outputs of the support vector machine. In the multi-class case the predicted probabilities are coupled using Hastie and Tibshirani's pairwise coupling method.<br/>
 * <br/>
 * Note: for improved speed normalization should be turned off when operating on SparseInstances.<br/>
 * <br/>
 * For more information on the SMO algorithm, see<br/>
 * <br/>
 * J. Platt: Fast Training of Support Vector Machines using Sequential Minimal Optimization. In B. Schoelkopf and C. Burges and A. Smola, editors, Advances in Kernel Methods - Support Vector Learning, 1998.<br/>
 * <br/>
 * S.S. Keerthi, S.K. Shevade, C. Bhattacharyya, K.R.K. Murthy (2001). Improvements to Platt's SMO Algorithm for SVM Classifier Design. Neural Computation. 13(3):637-649.<br/>
 * <br/>
 * Trevor Hastie, Robert Tibshirani: Classification by Pairwise Coupling. In: Advances in Neural Information Processing Systems, 1998.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;incollection{Platt1998,
 *    author = {J. Platt},
 *    booktitle = {Advances in Kernel Methods - Support Vector Learning},
 *    editor = {B. Schoelkopf and C. Burges and A. Smola},
 *    publisher = {MIT Press},
 *    title = {Fast Training of Support Vector Machines using Sequential Minimal Optimization},
 *    year = {1998},
 *    URL = {http://research.microsoft.com/\~jplatt/smo.html},
 *    PS = {http://research.microsoft.com/\~jplatt/smo-book.ps.gz},
 *    PDF = {http://research.microsoft.com/\~jplatt/smo-book.pdf}
 * }
 * 
 * &#64;article{Keerthi2001,
 *    author = {S.S. Keerthi and S.K. Shevade and C. Bhattacharyya and K.R.K. Murthy},
 *    journal = {Neural Computation},
 *    number = {3},
 *    pages = {637-649},
 *    title = {Improvements to Platt's SMO Algorithm for SVM Classifier Design},
 *    volume = {13},
 *    year = {2001},
 *    PS = {http://guppy.mpe.nus.edu.sg/\~mpessk/svm/smo_mod_nc.ps.gz}
 * }
 * 
 * &#64;inproceedings{Hastie1998,
 *    author = {Trevor Hastie and Robert Tibshirani},
 *    booktitle = {Advances in Neural Information Processing Systems},
 *    editor = {Michael I. Jordan and Michael J. Kearns and Sara A. Solla},
 *    publisher = {MIT Press},
 *    title = {Classification by Pairwise Coupling},
 *    volume = {10},
 *    year = {1998},
 *    PS = {http://www-stat.stanford.edu/\~hastie/Papers/2class.ps}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 * <pre> -no-checks
 *  Turns off all checks - use with caution!
 *  Turning them off assumes that data is purely numeric, doesn't
 *  contain any missing values, and has a nominal class. Turning them
 *  off also means that no header information will be stored if the
 *  machine is linear. Finally, it also assumes that no instance has
 *  a weight equal to 0.
 *  (default: checks on)</pre>
 * 
 * <pre> -C &lt;double&gt;
 *  The complexity constant C. (default 1)</pre>
 * 
 * <pre> -N
 *  Whether to 0=normalize/1=standardize/2=neither. (default 0=normalize)</pre>
 * 
 * <pre> -L &lt;double&gt;
 *  The tolerance parameter. (default 1.0e-3)</pre>
 * 
 * <pre> -P &lt;double&gt;
 *  The epsilon for round-off error. (default 1.0e-12)</pre>
 * 
 * <pre> -M
 *  Fit logistic models to SVM outputs. </pre>
 * 
 * <pre> -V &lt;double&gt;
 *  The number of folds for the internal
 *  cross-validation. (default -1, use training data)</pre>
 * 
 * <pre> -W &lt;double&gt;
 *  The random number seed. (default 1)</pre>
 * 
 * <pre> -K &lt;classname and parameters&gt;
 *  The Kernel to use.
 *  (default: weka.classifiers.functions.supportVector.PolyKernel)</pre>
 * 
 * <pre> 
 * Options specific to kernel weka.classifiers.functions.supportVector.PolyKernel:
 * </pre>
 * 
 * <pre> -D
 *  Enables debugging output (if available) to be printed.
 *  (default: off)</pre>
 * 
 * <pre> -no-checks
 *  Turns off all checks - use with caution!
 *  (default: checks on)</pre>
 * 
 * <pre> -C &lt;num&gt;
 *  The size of the cache (a prime number), 0 for full cache and 
 *  -1 to turn it off.
 *  (default: 250007)</pre>
 * 
 * <pre> -E &lt;num&gt;
 *  The Exponent to use.
 *  (default: 1.0)</pre>
 * 
 * <pre> -L
 *  Use lower-order terms.
 *  (default: no)</pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Shane Legg (shane@intelligenesis.net) (sparse vector code)
 * @author Stuart Inglis (stuart@reeltwo.com) (sparse vector code)
 * @version $Revision: 6024 $
 */
public class SMO 
  extends AbstractClassifier 
  implements WeightedInstancesHandler, TechnicalInformationHandler {

  /** for serialization */
  static final long serialVersionUID = -6585883636378691736L;
  
  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return  "Implements John Platt's sequential minimal optimization "
      + "algorithm for training a support vector classifier.\n\n"
      + "This implementation globally replaces all missing values and "
      + "transforms nominal attributes into binary ones. It also "
      + "normalizes all attributes by default. (In that case the coefficients "
      + "in the output are based on the normalized data, not the "
      + "original data --- this is important for interpreting the classifier.)\n\n"
      + "Multi-class problems are solved using pairwise classification "
      + "(1-vs-1 and if logistic models are built pairwise coupling "
      + "according to Hastie and Tibshirani, 1998).\n\n"
      + "To obtain proper probability estimates, use the option that fits "
      + "logistic regression models to the outputs of the support vector "
      + "machine. In the multi-class case the predicted probabilities "
      + "are coupled using Hastie and Tibshirani's pairwise coupling "
      + "method.\n\n"
      + "Note: for improved speed normalization should be turned off when "
      + "operating on SparseInstances.\n\n"
      + "For more information on the SMO algorithm, see\n\n"
      + getTechnicalInformation().toString();
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
    TechnicalInformation 	additional;
    
    result = new TechnicalInformation(Type.INCOLLECTION);
    result.setValue(Field.AUTHOR, "J. Platt");
    result.setValue(Field.YEAR, "1998");
    result.setValue(Field.TITLE, "Fast Training of Support Vector Machines using Sequential Minimal Optimization");
    result.setValue(Field.BOOKTITLE, "Advances in Kernel Methods - Support Vector Learning");
    result.setValue(Field.EDITOR, "B. Schoelkopf and C. Burges and A. Smola");
    result.setValue(Field.PUBLISHER, "MIT Press");
    result.setValue(Field.URL, "http://research.microsoft.com/~jplatt/smo.html");
    result.setValue(Field.PDF, "http://research.microsoft.com/~jplatt/smo-book.pdf");
    result.setValue(Field.PS, "http://research.microsoft.com/~jplatt/smo-book.ps.gz");
    
    additional = result.add(Type.ARTICLE);
    additional.setValue(Field.AUTHOR, "S.S. Keerthi and S.K. Shevade and C. Bhattacharyya and K.R.K. Murthy");
    additional.setValue(Field.YEAR, "2001");
    additional.setValue(Field.TITLE, "Improvements to Platt's SMO Algorithm for SVM Classifier Design");
    additional.setValue(Field.JOURNAL, "Neural Computation");
    additional.setValue(Field.VOLUME, "13");
    additional.setValue(Field.NUMBER, "3");
    additional.setValue(Field.PAGES, "637-649");
    additional.setValue(Field.PS, "http://guppy.mpe.nus.edu.sg/~mpessk/svm/smo_mod_nc.ps.gz");
    
    additional = result.add(Type.INPROCEEDINGS);
    additional.setValue(Field.AUTHOR, "Trevor Hastie and Robert Tibshirani");
    additional.setValue(Field.YEAR, "1998");
    additional.setValue(Field.TITLE, "Classification by Pairwise Coupling");
    additional.setValue(Field.BOOKTITLE, "Advances in Neural Information Processing Systems");
    additional.setValue(Field.VOLUME, "10");
    additional.setValue(Field.PUBLISHER, "MIT Press");
    additional.setValue(Field.EDITOR, "Michael I. Jordan and Michael J. Kearns and Sara A. Solla");
    additional.setValue(Field.PS, "http://www-stat.stanford.edu/~hastie/Papers/2class.ps");
    
    return result;
  }

  /**
   * Class for building a binary support vector machine.
   */
  public class BinarySMO 
    implements Serializable {
    
    /** for serialization */
    static final long serialVersionUID = -8246163625699362456L;
    
    /** The Lagrange multipliers. */
    protected double[] m_alpha;

    /** The thresholds. */
    protected double m_b, m_bLow, m_bUp;

    /** The indices for m_bLow and m_bUp */
    protected int m_iLow, m_iUp;

    /** The training data. */
    protected Instances m_data;

    /** Weight vector for linear machine. */
    protected double[] m_weights;

    /** Variables to hold weight vector in sparse form.
	(To reduce storage requirements.) */
    protected double[] m_sparseWeights;
    protected int[] m_sparseIndices;

    /** Kernel to use **/
    protected Kernel m_kernel;

    /** The transformed class values. */
    protected double[] m_class;

    /** The current set of errors for all non-bound examples. */
    protected double[] m_errors;

    /* The five different sets used by the algorithm. */
    /** {i: 0 < m_alpha[i] < C} */
    protected SMOset m_I0;
    /**  {i: m_class[i] = 1, m_alpha[i] = 0} */
    protected SMOset m_I1; 
    /**  {i: m_class[i] = -1, m_alpha[i] =C} */
    protected SMOset m_I2; 
    /** {i: m_class[i] = 1, m_alpha[i] = C} */
    protected SMOset m_I3;
    /**  {i: m_class[i] = -1, m_alpha[i] = 0} */
    protected SMOset m_I4; 

    /** The set of support vectors */
    protected SMOset m_supportVectors; // {i: 0 < m_alpha[i]}

    /** Stores logistic regression model for probability estimate */
    protected Logistic m_logistic = null;

    /** Stores the weight of the training instances */
    protected double m_sumOfWeights = 0;

    /**
     * Fits logistic regression model to SVM outputs analogue
     * to John Platt's method.  
     *
     * @param insts the set of training instances
     * @param cl1 the first class' index
     * @param cl2 the second class' index
     * @param numFolds the number of folds for cross-validation
     * @param random for randomizing the data
     * @throws Exception if the sigmoid can't be fit successfully
     */
    protected void fitLogistic(Instances insts, int cl1, int cl2,
			     int numFolds, Random random) 
      throws Exception {

      // Create header of instances object
      FastVector atts = new FastVector(2);
      atts.addElement(new Attribute("pred"));
      FastVector attVals = new FastVector(2);
      attVals.addElement(insts.classAttribute().value(cl1));
      attVals.addElement(insts.classAttribute().value(cl2));
      atts.addElement(new Attribute("class", attVals));
      Instances data = new Instances("data", atts, insts.numInstances());
      data.setClassIndex(1);

      // Collect data for fitting the logistic model
      if (numFolds <= 0) {

	// Use training data
	for (int j = 0; j < insts.numInstances(); j++) {
	  Instance inst = insts.instance(j);
	  double[] vals = new double[2];
	  vals[0] = SVMOutput(-1, inst);
	  if (inst.classValue() == cl2) {
	    vals[1] = 1;
	  }
	  data.add(new DenseInstance(inst.weight(), vals));
	}
      } else {

	// Check whether number of folds too large
	if (numFolds > insts.numInstances()) {
	  numFolds = insts.numInstances();
	}

	// Make copy of instances because we will shuffle them around
	insts = new Instances(insts);
	
	// Perform three-fold cross-validation to collect
	// unbiased predictions
	insts.randomize(random);
	insts.stratify(numFolds);
	for (int i = 0; i < numFolds; i++) {
	  Instances train = insts.trainCV(numFolds, i, random);
          /*	  SerializedObject so = new SerializedObject(this);
                  BinarySMO smo = (BinarySMO)so.getObject(); */
          BinarySMO smo = new BinarySMO();
          smo.setKernel(Kernel.makeCopy(SMO.this.m_kernel));
          smo.buildClassifier(train, cl1, cl2, false, -1, -1);
	  Instances test = insts.testCV(numFolds, i);
	  for (int j = 0; j < test.numInstances(); j++) {
	    double[] vals = new double[2];
	    vals[0] = smo.SVMOutput(-1, test.instance(j));
	    if (test.instance(j).classValue() == cl2) {
	      vals[1] = 1;
	    }
	    data.add(new DenseInstance(test.instance(j).weight(), vals));
	  }
	}
      }

      // Build logistic regression model
      m_logistic = new Logistic();
      m_logistic.buildClassifier(data);
    }
    
    /**
     * sets the kernel to use
     * 
     * @param value	the kernel to use
     */
    public void setKernel(Kernel value) {
      m_kernel = value;
    }
    
    /**
     * Returns the kernel to use
     * 
     * @return 		the current kernel
     */
    public Kernel getKernel() {
      return m_kernel;
    }

    /**
     * Method for building the binary classifier.
     *
     * @param insts the set of training instances
     * @param cl1 the first class' index
     * @param cl2 the second class' index
     * @param fitLogistic true if logistic model is to be fit
     * @param numFolds number of folds for internal cross-validation
     * @param randomSeed random number generator for cross-validation
     * @throws Exception if the classifier can't be built successfully
     */
    protected void buildClassifier(Instances insts, int cl1, int cl2,
				 boolean fitLogistic, int numFolds,
				 int randomSeed) throws Exception {
      
      // Initialize some variables
      m_bUp = -1; m_bLow = 1; m_b = 0; 
      m_alpha = null; m_data = null; m_weights = null; m_errors = null;
      m_logistic = null; m_I0 = null; m_I1 = null; m_I2 = null;
      m_I3 = null; m_I4 = null;	m_sparseWeights = null; m_sparseIndices = null;

      // Store the sum of weights
      m_sumOfWeights = insts.sumOfWeights();
      
      // Set class values
      m_class = new double[insts.numInstances()];
      m_iUp = -1; m_iLow = -1;
      for (int i = 0; i < m_class.length; i++) {
	if ((int) insts.instance(i).classValue() == cl1) {
	  m_class[i] = -1; m_iLow = i;
	} else if ((int) insts.instance(i).classValue() == cl2) {
	  m_class[i] = 1; m_iUp = i;
	} else {
	  throw new Exception ("This should never happen!");
	}
      }

      // Check whether one or both classes are missing
      if ((m_iUp == -1) || (m_iLow == -1)) {
	if (m_iUp != -1) {
	  m_b = -1;
	} else if (m_iLow != -1) {
	  m_b = 1;
	} else {
	  m_class = null;
	  return;
	}
	if (m_KernelIsLinear) {
	  m_sparseWeights = new double[0];
	  m_sparseIndices = new int[0];
	  m_class = null;
	} else {
	  m_supportVectors = new SMOset(0);
	  m_alpha = new double[0];
	  m_class = new double[0];
	}

	// Fit sigmoid if requested
	if (fitLogistic) {
	  fitLogistic(insts, cl1, cl2, numFolds, new Random(randomSeed));
	}
	return;
      }
      
      // Set the reference to the data
      m_data = insts;

      // If machine is linear, reserve space for weights
      if (m_KernelIsLinear) {
	m_weights = new double[m_data.numAttributes()];
      } else {
	m_weights = null;
      }
      
      // Initialize alpha array to zero
      m_alpha = new double[m_data.numInstances()];
      
      // Initialize sets
      m_supportVectors = new SMOset(m_data.numInstances());
      m_I0 = new SMOset(m_data.numInstances());
      m_I1 = new SMOset(m_data.numInstances());
      m_I2 = new SMOset(m_data.numInstances());
      m_I3 = new SMOset(m_data.numInstances());
      m_I4 = new SMOset(m_data.numInstances());

      // Clean out some instance variables
      m_sparseWeights = null;
      m_sparseIndices = null;
      
      // init kernel
      m_kernel.buildKernel(m_data);
      
      // Initialize error cache
      m_errors = new double[m_data.numInstances()];
      m_errors[m_iLow] = 1; m_errors[m_iUp] = -1;
     
      // Build up I1 and I4
      for (int i = 0; i < m_class.length; i++ ) {
	if (m_class[i] == 1) {
	  m_I1.insert(i);
	} else {
	  m_I4.insert(i);
	}
      }
      
      // Loop to find all the support vectors
      int numChanged = 0;
      boolean examineAll = true;
      while ((numChanged > 0) || examineAll) {
	numChanged = 0;
	if (examineAll) {
	  for (int i = 0; i < m_alpha.length; i++) {
	    if (examineExample(i)) {
	      numChanged++;
	    }
	  }
	} else {
	  
	  // This code implements Modification 1 from Keerthi et al.'s paper
	  for (int i = 0; i < m_alpha.length; i++) {
	    if ((m_alpha[i] > 0) &&  
		(m_alpha[i] < m_C * m_data.instance(i).weight())) {
	      if (examineExample(i)) {
		numChanged++;
	      }
	      
	      // Is optimality on unbound vectors obtained?
	      if (m_bUp > m_bLow - 2 * m_tol) {
		numChanged = 0;
		break;
	      }
	    }
	  }
	  
	  //This is the code for Modification 2 from Keerthi et al.'s paper
	  /*boolean innerLoopSuccess = true; 
	    numChanged = 0;
	    while ((m_bUp < m_bLow - 2 * m_tol) && (innerLoopSuccess == true)) {
	    innerLoopSuccess = takeStep(m_iUp, m_iLow, m_errors[m_iLow]);
	    }*/
	}
	
	if (examineAll) {
	  examineAll = false;
	} else if (numChanged == 0) {
	  examineAll = true;
	}
      }
      
      // Set threshold
      m_b = (m_bLow + m_bUp) / 2.0;
      
      // Save memory
      m_kernel.clean(); 
      
      m_errors = null;
      m_I0 = m_I1 = m_I2 = m_I3 = m_I4 = null;
      
      // If machine is linear, delete training data
      // and store weight vector in sparse format
      if (m_KernelIsLinear) {
	
	// We don't need to store the set of support vectors
	m_supportVectors = null;

	// We don't need to store the class values either
	m_class = null;
	
	// Clean out training data
	if (!m_checksTurnedOff) {
	  m_data = new Instances(m_data, 0);
	} else {
	  m_data = null;
	}
	
	// Convert weight vector
	double[] sparseWeights = new double[m_weights.length];
	int[] sparseIndices = new int[m_weights.length];
	int counter = 0;
	for (int i = 0; i < m_weights.length; i++) {
	  if (m_weights[i] != 0.0) {
	    sparseWeights[counter] = m_weights[i];
	    sparseIndices[counter] = i;
	    counter++;
	  }
	}
	m_sparseWeights = new double[counter];
	m_sparseIndices = new int[counter];
	System.arraycopy(sparseWeights, 0, m_sparseWeights, 0, counter);
	System.arraycopy(sparseIndices, 0, m_sparseIndices, 0, counter);
	
	// Clean out weight vector
	m_weights = null;
	
	// We don't need the alphas in the linear case
	m_alpha = null;
      }
      
      // Fit sigmoid if requested
      if (fitLogistic) {
	fitLogistic(insts, cl1, cl2, numFolds, new Random(randomSeed));
      }

    }
    
    /**
     * Computes SVM output for given instance.
     *
     * @param index the instance for which output is to be computed
     * @param inst the instance 
     * @return the output of the SVM for the given instance
     * @throws Exception in case of an error
     */
    public double SVMOutput(int index, Instance inst) throws Exception {
      
      double result = 0;
      
      // Is the machine linear?
      if (m_KernelIsLinear) {
	
	// Is weight vector stored in sparse format?
	if (m_sparseWeights == null) {
	  int n1 = inst.numValues(); 
	  for (int p = 0; p < n1; p++) {
	    if (inst.index(p) != m_classIndex) {
	      result += m_weights[inst.index(p)] * inst.valueSparse(p);
	    }
	  }
	} else {
	  int n1 = inst.numValues(); int n2 = m_sparseWeights.length;
	  for (int p1 = 0, p2 = 0; p1 < n1 && p2 < n2;) {
	    int ind1 = inst.index(p1); 
	    int ind2 = m_sparseIndices[p2];
	    if (ind1 == ind2) {
	      if (ind1 != m_classIndex) {
		result += inst.valueSparse(p1) * m_sparseWeights[p2];
	      }
	      p1++; p2++;
	    } else if (ind1 > ind2) {
	      p2++;
	    } else { 
	      p1++;
	    }
	  }
	}
      } else {
	for (int i = m_supportVectors.getNext(-1); i != -1; 
	     i = m_supportVectors.getNext(i)) {
	  result += m_class[i] * m_alpha[i] * m_kernel.eval(index, i, inst);
	}
      }
      result -= m_b;
      
      return result;
    }

    /**
     * Prints out the classifier.
     *
     * @return a description of the classifier as a string
     */
    public String toString() {

      StringBuffer text = new StringBuffer();
      int printed = 0;

      if ((m_alpha == null) && (m_sparseWeights == null)) {
	return "BinarySMO: No model built yet.\n";
      }
      try {
	text.append("BinarySMO\n\n");

	// If machine linear, print weight vector
	if (m_KernelIsLinear) {
	  text.append("Machine linear: showing attribute weights, ");
	  text.append("not support vectors.\n\n");

	  // We can assume that the weight vector is stored in sparse
	  // format because the classifier has been built
	  for (int i = 0; i < m_sparseWeights.length; i++) {
	    if (m_sparseIndices[i] != (int)m_classIndex) {
	      if (printed > 0) {
		text.append(" + ");
	      } else {
		text.append("   ");
	      }
	      text.append(Utils.doubleToString(m_sparseWeights[i], 12, 4) +
			  " * ");
	      if (m_filterType == FILTER_STANDARDIZE) {
		text.append("(standardized) ");
	      } else if (m_filterType == FILTER_NORMALIZE) {
		text.append("(normalized) ");
	      }
	      if (!m_checksTurnedOff) {
		text.append(m_data.attribute(m_sparseIndices[i]).name()+"\n");
	      } else {
		text.append("attribute with index " + 
			    m_sparseIndices[i] +"\n");
	      }
	      printed++;
	    }
	  }
	} else {
	  for (int i = 0; i < m_alpha.length; i++) {
	    if (m_supportVectors.contains(i)) {
	      double val = m_alpha[i];
	      if (m_class[i] == 1) {
		if (printed > 0) {
		  text.append(" + ");
		}
	      } else {
		text.append(" - ");
	      }
	      text.append(Utils.doubleToString(val, 12, 4) 
			  + " * <");
	      for (int j = 0; j < m_data.numAttributes(); j++) {
		if (j != m_data.classIndex()) {
		  text.append(m_data.instance(i).toString(j));
		}
		if (j != m_data.numAttributes() - 1) {
		  text.append(" ");
		}
	      }
	      text.append("> * X]\n");
	      printed++;
	    }
	  }
	}
	if (m_b > 0) {
	  text.append(" - " + Utils.doubleToString(m_b, 12, 4));
	} else {
	  text.append(" + " + Utils.doubleToString(-m_b, 12, 4));
	}

	if (!m_KernelIsLinear) {
	  text.append("\n\nNumber of support vectors: " + 
		      m_supportVectors.numElements());
	}
	int numEval = 0;
	int numCacheHits = -1;
	if (m_kernel != null) {
	  numEval = m_kernel.numEvals();
	  numCacheHits = m_kernel.numCacheHits();
	}
	text.append("\n\nNumber of kernel evaluations: " + numEval);
	if (numCacheHits >= 0 && numEval > 0) {
	  double hitRatio = 1 - numEval*1.0/(numCacheHits+numEval);
	  text.append(" (" + Utils.doubleToString(hitRatio*100, 7, 3).trim() + "% cached)");
	}

      } catch (Exception e) {
	e.printStackTrace();

	return "Can't print BinarySMO classifier.";
      }
    
      return text.toString();
    }

    /**
     * Examines instance.
     *
     * @param i2 index of instance to examine
     * @return true if examination was successfull
     * @throws Exception if something goes wrong
     */
    protected boolean examineExample(int i2) throws Exception {
    
      double y2, F2;
      int i1 = -1;
    
      y2 = m_class[i2];
      if (m_I0.contains(i2)) {
	F2 = m_errors[i2];
      } else {
	F2 = SVMOutput(i2, m_data.instance(i2)) + m_b - y2;
	m_errors[i2] = F2;
      
	// Update thresholds
	if ((m_I1.contains(i2) || m_I2.contains(i2)) && (F2 < m_bUp)) {
	  m_bUp = F2; m_iUp = i2;
	} else if ((m_I3.contains(i2) || m_I4.contains(i2)) && (F2 > m_bLow)) {
	  m_bLow = F2; m_iLow = i2;
	}
      }

      // Check optimality using current bLow and bUp and, if
      // violated, find an index i1 to do joint optimization
      // with i2...
      boolean optimal = true;
      if (m_I0.contains(i2) || m_I1.contains(i2) || m_I2.contains(i2)) {
	if (m_bLow - F2 > 2 * m_tol) {
	  optimal = false; i1 = m_iLow;
	}
      }
      if (m_I0.contains(i2) || m_I3.contains(i2) || m_I4.contains(i2)) {
	if (F2 - m_bUp > 2 * m_tol) {
	  optimal = false; i1 = m_iUp;
	}
      }
      if (optimal) {
	return false;
      }

      // For i2 unbound choose the better i1...
      if (m_I0.contains(i2)) {
	if (m_bLow - F2 > F2 - m_bUp) {
	  i1 = m_iLow;
	} else {
	  i1 = m_iUp;
	}
      }
      if (i1 == -1) {
	throw new Exception("This should never happen!");
      }
      return takeStep(i1, i2, F2);
    }

    /**
     * Method solving for the Lagrange multipliers for
     * two instances.
     *
     * @param i1 index of the first instance
     * @param i2 index of the second instance
     * @param F2
     * @return true if multipliers could be found
     * @throws Exception if something goes wrong
     */
    protected boolean takeStep(int i1, int i2, double F2) throws Exception {

      double alph1, alph2, y1, y2, F1, s, L, H, k11, k12, k22, eta,
	a1, a2, f1, f2, v1, v2, Lobj, Hobj;
      double C1 = m_C * m_data.instance(i1).weight();
      double C2 = m_C * m_data.instance(i2).weight();

      // Don't do anything if the two instances are the same
      if (i1 == i2) {
	return false;
      }

      // Initialize variables
      alph1 = m_alpha[i1]; alph2 = m_alpha[i2];
      y1 = m_class[i1]; y2 = m_class[i2];
      F1 = m_errors[i1];
      s = y1 * y2;

      // Find the constraints on a2
      if (y1 != y2) {
	L = Math.max(0, alph2 - alph1); 
	H = Math.min(C2, C1 + alph2 - alph1);
      } else {
	L = Math.max(0, alph1 + alph2 - C1);
	H = Math.min(C2, alph1 + alph2);
      }
      if (L >= H) {
	return false;
      }

      // Compute second derivative of objective function
      k11 = m_kernel.eval(i1, i1, m_data.instance(i1));
      k12 = m_kernel.eval(i1, i2, m_data.instance(i1));
      k22 = m_kernel.eval(i2, i2, m_data.instance(i2));
      eta = 2 * k12 - k11 - k22;

      // Check if second derivative is negative
      if (eta < 0) {

	// Compute unconstrained maximum
	a2 = alph2 - y2 * (F1 - F2) / eta;

	// Compute constrained maximum
	if (a2 < L) {
	  a2 = L;
	} else if (a2 > H) {
	  a2 = H;
	}
      } else {

	// Look at endpoints of diagonal
	f1 = SVMOutput(i1, m_data.instance(i1));
	f2 = SVMOutput(i2, m_data.instance(i2));
	v1 = f1 + m_b - y1 * alph1 * k11 - y2 * alph2 * k12; 
	v2 = f2 + m_b - y1 * alph1 * k12 - y2 * alph2 * k22; 
	double gamma = alph1 + s * alph2;
	Lobj = (gamma - s * L) + L - 0.5 * k11 * (gamma - s * L) * (gamma - s * L) - 
	  0.5 * k22 * L * L - s * k12 * (gamma - s * L) * L - 
	  y1 * (gamma - s * L) * v1 - y2 * L * v2;
	Hobj = (gamma - s * H) + H - 0.5 * k11 * (gamma - s * H) * (gamma - s * H) - 
	  0.5 * k22 * H * H - s * k12 * (gamma - s * H) * H - 
	  y1 * (gamma - s * H) * v1 - y2 * H * v2;
	if (Lobj > Hobj + m_eps) {
	  a2 = L;
	} else if (Lobj < Hobj - m_eps) {
	  a2 = H;
	} else {
	  a2 = alph2;
	}
      }
      if (Math.abs(a2 - alph2) < m_eps * (a2 + alph2 + m_eps)) {
	return false;
      }
      
      // To prevent precision problems
      if (a2 > C2 - m_Del * C2) {
	a2 = C2;
      } else if (a2 <= m_Del * C2) {
	a2 = 0;
      }
      
      // Recompute a1
      a1 = alph1 + s * (alph2 - a2);
      
      // To prevent precision problems
      if (a1 > C1 - m_Del * C1) {
	a1 = C1;
      } else if (a1 <= m_Del * C1) {
	a1 = 0;
      }
      
      // Update sets
      if (a1 > 0) {
	m_supportVectors.insert(i1);
      } else {
	m_supportVectors.delete(i1);
      }
      if ((a1 > 0) && (a1 < C1)) {
	m_I0.insert(i1);
      } else {
	m_I0.delete(i1);
      }
      if ((y1 == 1) && (a1 == 0)) {
	m_I1.insert(i1);
      } else {
	m_I1.delete(i1);
      }
      if ((y1 == -1) && (a1 == C1)) {
	m_I2.insert(i1);
      } else {
	m_I2.delete(i1);
      }
      if ((y1 == 1) && (a1 == C1)) {
	m_I3.insert(i1);
      } else {
	m_I3.delete(i1);
      }
      if ((y1 == -1) && (a1 == 0)) {
	m_I4.insert(i1);
      } else {
	m_I4.delete(i1);
      }
      if (a2 > 0) {
	m_supportVectors.insert(i2);
      } else {
	m_supportVectors.delete(i2);
      }
      if ((a2 > 0) && (a2 < C2)) {
	m_I0.insert(i2);
      } else {
	m_I0.delete(i2);
      }
      if ((y2 == 1) && (a2 == 0)) {
	m_I1.insert(i2);
      } else {
	m_I1.delete(i2);
      }
      if ((y2 == -1) && (a2 == C2)) {
	m_I2.insert(i2);
      } else {
	m_I2.delete(i2);
      }
      if ((y2 == 1) && (a2 == C2)) {
	m_I3.insert(i2);
      } else {
	m_I3.delete(i2);
      }
      if ((y2 == -1) && (a2 == 0)) {
	m_I4.insert(i2);
      } else {
	m_I4.delete(i2);
      }
      
      // Update weight vector to reflect change a1 and a2, if linear SVM
      if (m_KernelIsLinear) {
	Instance inst1 = m_data.instance(i1);
	for (int p1 = 0; p1 < inst1.numValues(); p1++) {
	  if (inst1.index(p1) != m_data.classIndex()) {
	    m_weights[inst1.index(p1)] += 
	      y1 * (a1 - alph1) * inst1.valueSparse(p1);
	  }
	}
	Instance inst2 = m_data.instance(i2);
	for (int p2 = 0; p2 < inst2.numValues(); p2++) {
	  if (inst2.index(p2) != m_data.classIndex()) {
	    m_weights[inst2.index(p2)] += 
	      y2 * (a2 - alph2) * inst2.valueSparse(p2);
	  }
	}
      }
      
      // Update error cache using new Lagrange multipliers
      for (int j = m_I0.getNext(-1); j != -1; j = m_I0.getNext(j)) {
	if ((j != i1) && (j != i2)) {
	  m_errors[j] += 
	    y1 * (a1 - alph1) * m_kernel.eval(i1, j, m_data.instance(i1)) + 
	    y2 * (a2 - alph2) * m_kernel.eval(i2, j, m_data.instance(i2));
	}
      }
      
      // Update error cache for i1 and i2
      m_errors[i1] += y1 * (a1 - alph1) * k11 + y2 * (a2 - alph2) * k12;
      m_errors[i2] += y1 * (a1 - alph1) * k12 + y2 * (a2 - alph2) * k22;
      
      // Update array with Lagrange multipliers
      m_alpha[i1] = a1;
      m_alpha[i2] = a2;
      
      // Update thresholds
      m_bLow = -Double.MAX_VALUE; m_bUp = Double.MAX_VALUE;
      m_iLow = -1; m_iUp = -1;
      for (int j = m_I0.getNext(-1); j != -1; j = m_I0.getNext(j)) {
	if (m_errors[j] < m_bUp) {
	  m_bUp = m_errors[j]; m_iUp = j;
	}
	if (m_errors[j] > m_bLow) {
	  m_bLow = m_errors[j]; m_iLow = j;
	}
      }
      if (!m_I0.contains(i1)) {
	if (m_I3.contains(i1) || m_I4.contains(i1)) {
	  if (m_errors[i1] > m_bLow) {
	    m_bLow = m_errors[i1]; m_iLow = i1;
	  } 
	} else {
	  if (m_errors[i1] < m_bUp) {
	    m_bUp = m_errors[i1]; m_iUp = i1;
	  }
	}
      }
      if (!m_I0.contains(i2)) {
	if (m_I3.contains(i2) || m_I4.contains(i2)) {
	  if (m_errors[i2] > m_bLow) {
	    m_bLow = m_errors[i2]; m_iLow = i2;
	  }
	} else {
	  if (m_errors[i2] < m_bUp) {
	    m_bUp = m_errors[i2]; m_iUp = i2;
	  }
	}
      }
      if ((m_iLow == -1) || (m_iUp == -1)) {
	throw new Exception("This should never happen!");
      }

      // Made some progress.
      return true;
    }
  
    /**
     * Quick and dirty check whether the quadratic programming problem is solved.
     * 
     * @throws Exception if checking fails
     */
    protected void checkClassifier() throws Exception {

      double sum = 0;
      for (int i = 0; i < m_alpha.length; i++) {
	if (m_alpha[i] > 0) {
	  sum += m_class[i] * m_alpha[i];
	}
      }
      System.err.println("Sum of y(i) * alpha(i): " + sum);

      for (int i = 0; i < m_alpha.length; i++) {
	double output = SVMOutput(i, m_data.instance(i));
	if (Utils.eq(m_alpha[i], 0)) {
	  if (Utils.sm(m_class[i] * output, 1)) {
	    System.err.println("KKT condition 1 violated: " + m_class[i] * output);
	  }
	} 
	if (Utils.gr(m_alpha[i], 0) && 
	    Utils.sm(m_alpha[i], m_C * m_data.instance(i).weight())) {
	  if (!Utils.eq(m_class[i] * output, 1)) {
	    System.err.println("KKT condition 2 violated: " + m_class[i] * output);
	  }
	} 
	if (Utils.eq(m_alpha[i], m_C * m_data.instance(i).weight())) {
	  if (Utils.gr(m_class[i] * output, 1)) {
	    System.err.println("KKT condition 3 violated: " + m_class[i] * output);
	  }
	} 
      }
    }  
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 6024 $");
    }
  }

  /** filter: Normalize training data */
  public static final int FILTER_NORMALIZE = 0;
  /** filter: Standardize training data */
  public static final int FILTER_STANDARDIZE = 1;
  /** filter: No normalization/standardization */
  public static final int FILTER_NONE = 2;
  /** The filter to apply to the training data */
  public static final Tag [] TAGS_FILTER = {
    new Tag(FILTER_NORMALIZE, "Normalize training data"),
    new Tag(FILTER_STANDARDIZE, "Standardize training data"),
    new Tag(FILTER_NONE, "No normalization/standardization"),
  };

  /** The binary classifier(s) */
  protected BinarySMO[][] m_classifiers = null;
  
  /** The complexity parameter. */
  protected double m_C = 1.0;
  
  /** Epsilon for rounding. */
  protected double m_eps = 1.0e-12;
  
  /** Tolerance for accuracy of result. */
  protected double m_tol = 1.0e-3;

  /** Whether to normalize/standardize/neither */
  protected int m_filterType = FILTER_NORMALIZE;

  /** The filter used to make attributes numeric. */
  protected NominalToBinary m_NominalToBinary;

  /** The filter used to standardize/normalize all values. */
  protected Filter m_Filter = null;

  /** The filter used to get rid of missing values. */
  protected ReplaceMissingValues m_Missing;

  /** The class index from the training data */
  protected int m_classIndex = -1;

  /** The class attribute */
  protected Attribute m_classAttribute;
  
  /** whether the kernel is a linear one */
  protected boolean m_KernelIsLinear = false;

  /** Turn off all checks and conversions? Turning them off assumes
      that data is purely numeric, doesn't contain any missing values,
      and has a nominal class. Turning them off also means that
      no header information will be stored if the machine is linear. 
      Finally, it also assumes that no instance has a weight equal to 0.*/
  protected boolean m_checksTurnedOff;

  /** Precision constant for updating sets */
  protected static double m_Del = 1000 * Double.MIN_VALUE;

  /** Whether logistic models are to be fit */
  protected boolean m_fitLogisticModels = false;

  /** The number of folds for the internal cross-validation */
  protected int m_numFolds = -1;

  /** The random number seed  */
  protected int m_randomSeed = 1;

  /** the kernel to use */
  protected Kernel m_kernel = new PolyKernel();
  
  /**
   * Turns off checks for missing values, etc. Use with caution.
   */
  public void turnChecksOff() {

    m_checksTurnedOff = true;
  }

  /**
   * Turns on checks for missing values, etc.
   */
  public void turnChecksOn() {

    m_checksTurnedOff = false;
  }

  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = getKernel().getCapabilities();
    result.setOwner(this);
    
    // attribute
    result.enableAllAttributeDependencies();
    // with NominalToBinary we can also handle nominal attributes, but only
    // if the kernel can handle numeric attributes
    if (result.handles(Capability.NUMERIC_ATTRIBUTES))
      result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.disableAllClasses();
    result.disableAllClassDependencies();
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }

  /**
   * Method for building the classifier. Implements a one-against-one
   * wrapper for multi-class problems.
   *
   * @param insts the set of training instances
   * @throws Exception if the classifier can't be built successfully
   */
  public void buildClassifier(Instances insts) throws Exception {

    if (!m_checksTurnedOff) {
      // can classifier handle the data?
      getCapabilities().testWithFail(insts);

      // remove instances with missing class
      insts = new Instances(insts);
      insts.deleteWithMissingClass();
      
      /* Removes all the instances with weight equal to 0.
       MUST be done since condition (8) of Keerthi's paper 
       is made with the assertion Ci > 0 (See equation (3a). */
      Instances data = new Instances(insts, insts.numInstances());
      for(int i = 0; i < insts.numInstances(); i++){
        if(insts.instance(i).weight() > 0)
          data.add(insts.instance(i));
      }
      if (data.numInstances() == 0) {
        throw new Exception("No training instances left after removing " + 
        "instances with weight 0!");
      }
      insts = data;
    }

    if (!m_checksTurnedOff) {
      m_Missing = new ReplaceMissingValues();
      m_Missing.setInputFormat(insts);
      insts = Filter.useFilter(insts, m_Missing); 
    } else {
      m_Missing = null;
    }

    if (getCapabilities().handles(Capability.NUMERIC_ATTRIBUTES)) {
      boolean onlyNumeric = true;
      if (!m_checksTurnedOff) {
	for (int i = 0; i < insts.numAttributes(); i++) {
	  if (i != insts.classIndex()) {
	    if (!insts.attribute(i).isNumeric()) {
	      onlyNumeric = false;
	      break;
	    }
	  }
	}
      }
      
      if (!onlyNumeric) {
	m_NominalToBinary = new NominalToBinary();
	m_NominalToBinary.setInputFormat(insts);
	insts = Filter.useFilter(insts, m_NominalToBinary);
      } 
      else {
	m_NominalToBinary = null;
      }
    }
    else {
      m_NominalToBinary = null;
    }

    if (m_filterType == FILTER_STANDARDIZE) {
      m_Filter = new Standardize();
      m_Filter.setInputFormat(insts);
      insts = Filter.useFilter(insts, m_Filter); 
    } else if (m_filterType == FILTER_NORMALIZE) {
      m_Filter = new Normalize();
      m_Filter.setInputFormat(insts);
      insts = Filter.useFilter(insts, m_Filter); 
    } else {
      m_Filter = null;
    }

    m_classIndex = insts.classIndex();
    m_classAttribute = insts.classAttribute();
    m_KernelIsLinear = (m_kernel instanceof PolyKernel) && (((PolyKernel) m_kernel).getExponent() == 1.0);
    
    // Generate subsets representing each class
    Instances[] subsets = new Instances[insts.numClasses()];
    for (int i = 0; i < insts.numClasses(); i++) {
      subsets[i] = new Instances(insts, insts.numInstances());
    }
    for (int j = 0; j < insts.numInstances(); j++) {
      Instance inst = insts.instance(j);
      subsets[(int)inst.classValue()].add(inst);
    }
    for (int i = 0; i < insts.numClasses(); i++) {
      subsets[i].compactify();
    }

    // Build the binary classifiers
    Random rand = new Random(m_randomSeed);
    m_classifiers = new BinarySMO[insts.numClasses()][insts.numClasses()];
    for (int i = 0; i < insts.numClasses(); i++) {
      for (int j = i + 1; j < insts.numClasses(); j++) {
	m_classifiers[i][j] = new BinarySMO();
	m_classifiers[i][j].setKernel(Kernel.makeCopy(getKernel()));
	Instances data = new Instances(insts, insts.numInstances());
	for (int k = 0; k < subsets[i].numInstances(); k++) {
	  data.add(subsets[i].instance(k));
	}
	for (int k = 0; k < subsets[j].numInstances(); k++) {
	  data.add(subsets[j].instance(k));
	}
	data.compactify();
	data.randomize(rand);
	m_classifiers[i][j].buildClassifier(data, i, j, 
					    m_fitLogisticModels,
					    m_numFolds, m_randomSeed);
      }
    }
  }

  /**
   * Estimates class probabilities for given instance.
   * 
   * @param inst the instance to compute the probabilities for
   * @throws Exception in case of an error
   */
  public double[] distributionForInstance(Instance inst) throws Exception {

    // Filter instance
    if (!m_checksTurnedOff) {
      m_Missing.input(inst);
      m_Missing.batchFinished();
      inst = m_Missing.output();
    }

    if (m_NominalToBinary != null) {
      m_NominalToBinary.input(inst);
      m_NominalToBinary.batchFinished();
      inst = m_NominalToBinary.output();
    }
    
    if (m_Filter != null) {
      m_Filter.input(inst);
      m_Filter.batchFinished();
      inst = m_Filter.output();
    }
    
    if (!m_fitLogisticModels) {
      double[] result = new double[inst.numClasses()];
      for (int i = 0; i < inst.numClasses(); i++) {
	for (int j = i + 1; j < inst.numClasses(); j++) {
	  if ((m_classifiers[i][j].m_alpha != null) || 
	      (m_classifiers[i][j].m_sparseWeights != null)) {
	    double output = m_classifiers[i][j].SVMOutput(-1, inst);
	    if (output > 0) {
	      result[j] += 1;
	    } else {
	      result[i] += 1;
	    }
	  }
	} 
      }
      Utils.normalize(result);
      return result;
    } else {

      // We only need to do pairwise coupling if there are more
      // then two classes.
      if (inst.numClasses() == 2) {
	double[] newInst = new double[2];
	newInst[0] = m_classifiers[0][1].SVMOutput(-1, inst);
	newInst[1] = Utils.missingValue();
	return m_classifiers[0][1].m_logistic.
	  distributionForInstance(new DenseInstance(1, newInst));
      }
      double[][] r = new double[inst.numClasses()][inst.numClasses()];
      double[][] n = new double[inst.numClasses()][inst.numClasses()];
      for (int i = 0; i < inst.numClasses(); i++) {
	for (int j = i + 1; j < inst.numClasses(); j++) {
	  if ((m_classifiers[i][j].m_alpha != null) || 
	      (m_classifiers[i][j].m_sparseWeights != null)) {
	    double[] newInst = new double[2];
	    newInst[0] = m_classifiers[i][j].SVMOutput(-1, inst);
	    newInst[1] = Utils.missingValue();
	    r[i][j] = m_classifiers[i][j].m_logistic.
	      distributionForInstance(new DenseInstance(1, newInst))[0];
	    n[i][j] = m_classifiers[i][j].m_sumOfWeights;
	  }
	}
      }
      return weka.classifiers.meta.MultiClassClassifier.pairwiseCoupling(n, r);
    }
  }

  /**
   * Returns an array of votes for the given instance.
   * @param inst the instance
   * @return array of votex
   * @throws Exception if something goes wrong
   */
  public int[] obtainVotes(Instance inst) throws Exception {

    // Filter instance
    if (!m_checksTurnedOff) {
      m_Missing.input(inst);
      m_Missing.batchFinished();
      inst = m_Missing.output();
    }

    if (m_NominalToBinary != null) {
      m_NominalToBinary.input(inst);
      m_NominalToBinary.batchFinished();
      inst = m_NominalToBinary.output();
    }
    
    if (m_Filter != null) {
      m_Filter.input(inst);
      m_Filter.batchFinished();
      inst = m_Filter.output();
    }

    int[] votes = new int[inst.numClasses()];
    for (int i = 0; i < inst.numClasses(); i++) {
      for (int j = i + 1; j < inst.numClasses(); j++) {
	double output = m_classifiers[i][j].SVMOutput(-1, inst);
	if (output > 0) {
	  votes[j] += 1;
	} else {
	  votes[i] += 1;
	}
      }
    }
    return votes;
  }

  /**
   * Returns the weights in sparse format.
   */
  public double [][][] sparseWeights() {
    
    int numValues = m_classAttribute.numValues();
    double [][][] sparseWeights = new double[numValues][numValues][];
    
    for (int i = 0; i < numValues; i++) {
      for (int j = i + 1; j < numValues; j++) {
	sparseWeights[i][j] = m_classifiers[i][j].m_sparseWeights;
      }
    }
    
    return sparseWeights;
  }
  
  /**
   * Returns the indices in sparse format.
   */
  public int [][][] sparseIndices() {
    
    int numValues = m_classAttribute.numValues();
    int [][][] sparseIndices = new int[numValues][numValues][];

    for (int i = 0; i < numValues; i++) {
      for (int j = i + 1; j < numValues; j++) {
	sparseIndices[i][j] = m_classifiers[i][j].m_sparseIndices;
      }
    }
    
    return sparseIndices;
  }
  
  /**
   * Returns the bias of each binary SMO.
   */
  public double [][] bias() {
    
    int numValues = m_classAttribute.numValues();
    double [][] bias = new double[numValues][numValues];

    for (int i = 0; i < numValues; i++) {
      for (int j = i + 1; j < numValues; j++) {
	bias[i][j] = m_classifiers[i][j].m_b;
      }
    }
    
    return bias;
  }
  
  /*
   * Returns the number of values of the class attribute.
   */
  public int numClassAttributeValues() {

    return m_classAttribute.numValues();
  }
  
  /*
   * Returns the names of the class attributes.
   */
  public String [] classAttributeNames() {

    int numValues = m_classAttribute.numValues();
    
    String [] classAttributeNames = new String[numValues];
    
    for (int i = 0; i < numValues; i++) {
      classAttributeNames[i] = m_classAttribute.value(i);
    }
    
    return classAttributeNames;
  }
  
  /**
   * Returns the attribute names.
   */
  public String [][][] attributeNames() {
    
    int numValues = m_classAttribute.numValues();
    String [][][] attributeNames = new String[numValues][numValues][];
    
    for (int i = 0; i < numValues; i++) {
      for (int j = i + 1; j < numValues; j++) {
        //	int numAttributes = m_classifiers[i][j].m_data.numAttributes();
	int numAttributes = m_classifiers[i][j].m_sparseIndices.length;
	String [] attrNames = new String[numAttributes];
	for (int k = 0; k < numAttributes; k++) {
	  attrNames[k] = m_classifiers[i][j].
            m_data.attribute(m_classifiers[i][j].m_sparseIndices[k]).name();
	}
	attributeNames[i][j] = attrNames;          
      }
    }
    return attributeNames;
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector result = new Vector();

    Enumeration enm = super.listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    result.addElement(new Option(
	"\tTurns off all checks - use with caution!\n"
	+ "\tTurning them off assumes that data is purely numeric, doesn't\n"
	+ "\tcontain any missing values, and has a nominal class. Turning them\n"
	+ "\toff also means that no header information will be stored if the\n"
	+ "\tmachine is linear. Finally, it also assumes that no instance has\n"
	+ "\ta weight equal to 0.\n"
	+ "\t(default: checks on)",
	"no-checks", 0, "-no-checks"));

    result.addElement(new Option(
	"\tThe complexity constant C. (default 1)",
	"C", 1, "-C <double>"));
    
    result.addElement(new Option(
	"\tWhether to 0=normalize/1=standardize/2=neither. " +
	"(default 0=normalize)",
	"N", 1, "-N"));
    
    result.addElement(new Option(
	"\tThe tolerance parameter. " +
	"(default 1.0e-3)",
	"L", 1, "-L <double>"));
    
    result.addElement(new Option(
	"\tThe epsilon for round-off error. " +
	"(default 1.0e-12)",
	"P", 1, "-P <double>"));
    
    result.addElement(new Option(
	"\tFit logistic models to SVM outputs. ",
	"M", 0, "-M"));
    
    result.addElement(new Option(
	"\tThe number of folds for the internal\n" +
	"\tcross-validation. " +
	"(default -1, use training data)",
	"V", 1, "-V <double>"));
    
    result.addElement(new Option(
	"\tThe random number seed. " +
	"(default 1)",
	"W", 1, "-W <double>"));
    
    result.addElement(new Option(
	"\tThe Kernel to use.\n"
	+ "\t(default: weka.classifiers.functions.supportVector.PolyKernel)",
	"K", 1, "-K <classname and parameters>"));

    result.addElement(new Option(
	"",
	"", 0, "\nOptions specific to kernel "
	+ getKernel().getClass().getName() + ":"));
    
    enm = ((OptionHandler) getKernel()).listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    return result.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  If set, classifier is run in debug mode and
   *  may output additional info to the console</pre>
   * 
   * <pre> -no-checks
   *  Turns off all checks - use with caution!
   *  Turning them off assumes that data is purely numeric, doesn't
   *  contain any missing values, and has a nominal class. Turning them
   *  off also means that no header information will be stored if the
   *  machine is linear. Finally, it also assumes that no instance has
   *  a weight equal to 0.
   *  (default: checks on)</pre>
   * 
   * <pre> -C &lt;double&gt;
   *  The complexity constant C. (default 1)</pre>
   * 
   * <pre> -N
   *  Whether to 0=normalize/1=standardize/2=neither. (default 0=normalize)</pre>
   * 
   * <pre> -L &lt;double&gt;
   *  The tolerance parameter. (default 1.0e-3)</pre>
   * 
   * <pre> -P &lt;double&gt;
   *  The epsilon for round-off error. (default 1.0e-12)</pre>
   * 
   * <pre> -M
   *  Fit logistic models to SVM outputs. </pre>
   * 
   * <pre> -V &lt;double&gt;
   *  The number of folds for the internal
   *  cross-validation. (default -1, use training data)</pre>
   * 
   * <pre> -W &lt;double&gt;
   *  The random number seed. (default 1)</pre>
   * 
   * <pre> -K &lt;classname and parameters&gt;
   *  The Kernel to use.
   *  (default: weka.classifiers.functions.supportVector.PolyKernel)</pre>
   * 
   * <pre> 
   * Options specific to kernel weka.classifiers.functions.supportVector.PolyKernel:
   * </pre>
   * 
   * <pre> -D
   *  Enables debugging output (if available) to be printed.
   *  (default: off)</pre>
   * 
   * <pre> -no-checks
   *  Turns off all checks - use with caution!
   *  (default: checks on)</pre>
   * 
   * <pre> -C &lt;num&gt;
   *  The size of the cache (a prime number), 0 for full cache and 
   *  -1 to turn it off.
   *  (default: 250007)</pre>
   * 
   * <pre> -E &lt;num&gt;
   *  The Exponent to use.
   *  (default: 1.0)</pre>
   * 
   * <pre> -L
   *  Use lower-order terms.
   *  (default: no)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported 
   */
  public void setOptions(String[] options) throws Exception {
    String 	tmpStr;
    String[]	tmpOptions;
    
    setChecksTurnedOff(Utils.getFlag("no-checks", options));

    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0)
      setC(Double.parseDouble(tmpStr));
    else
      setC(1.0);

    tmpStr = Utils.getOption('L', options);
    if (tmpStr.length() != 0)
      setToleranceParameter(Double.parseDouble(tmpStr));
    else
      setToleranceParameter(1.0e-3);
    
    tmpStr = Utils.getOption('P', options);
    if (tmpStr.length() != 0)
      setEpsilon(Double.parseDouble(tmpStr));
    else
      setEpsilon(1.0e-12);
    
    tmpStr = Utils.getOption('N', options);
    if (tmpStr.length() != 0)
      setFilterType(new SelectedTag(Integer.parseInt(tmpStr), TAGS_FILTER));
    else
      setFilterType(new SelectedTag(FILTER_NORMALIZE, TAGS_FILTER));
    
    setBuildLogisticModels(Utils.getFlag('M', options));
    
    tmpStr = Utils.getOption('V', options);
    if (tmpStr.length() != 0)
      setNumFolds(Integer.parseInt(tmpStr));
    else
      setNumFolds(-1);
    
    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() != 0)
      setRandomSeed(Integer.parseInt(tmpStr));
    else
      setRandomSeed(1);

    tmpStr     = Utils.getOption('K', options);
    tmpOptions = Utils.splitOptions(tmpStr);
    if (tmpOptions.length != 0) {
      tmpStr        = tmpOptions[0];
      tmpOptions[0] = "";
      setKernel(Kernel.forName(tmpStr, tmpOptions));
    }
    
    super.setOptions(options);
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int       i;
    Vector    result;
    String[]  options;

    result = new Vector();
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    if (getChecksTurnedOff())
      result.add("-no-checks");

    result.add("-C");
    result.add("" + getC());
    
    result.add("-L");
    result.add("" + getToleranceParameter());
    
    result.add("-P");
    result.add("" + getEpsilon());
    
    result.add("-N");
    result.add("" + m_filterType);
    
    if (getBuildLogisticModels())
      result.add("-M");
    
    result.add("-V");
    result.add("" + getNumFolds());
    
    result.add("-W");
    result.add("" + getRandomSeed());

    result.add("-K");
    result.add("" + getKernel().getClass().getName() + " " + Utils.joinOptions(getKernel().getOptions()));
    
    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * Disables or enables the checks (which could be time-consuming). Use with
   * caution!
   * 
   * @param value	if true turns off all checks
   */
  public void setChecksTurnedOff(boolean value) {
    if (value)
      turnChecksOff();
    else
      turnChecksOn();
  }
  
  /**
   * Returns whether the checks are turned off or not.
   * 
   * @return		true if the checks are turned off
   */
  public boolean getChecksTurnedOff() {
    return m_checksTurnedOff;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String checksTurnedOffTipText() {
    return "Turns time-consuming checks off - use with caution.";
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String kernelTipText() {
    return "The kernel to use.";
  }
  
  /**
   * sets the kernel to use
   * 
   * @param value	the kernel to use
   */
  public void setKernel(Kernel value) {
    m_kernel = value;
  }
  
  /**
   * Returns the kernel to use
   * 
   * @return 		the current kernel
   */
  public Kernel getKernel() {
    return m_kernel;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String cTipText() {
    return "The complexity parameter C.";
  }
  
  /**
   * Get the value of C.
   *
   * @return Value of C.
   */
  public double getC() {
    
    return m_C;
  }
  
  /**
   * Set the value of C.
   *
   * @param v  Value to assign to C.
   */
  public void setC(double v) {
    
    m_C = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String toleranceParameterTipText() {
    return "The tolerance parameter (shouldn't be changed).";
  }
  
  /**
   * Get the value of tolerance parameter.
   * @return Value of tolerance parameter.
   */
  public double getToleranceParameter() {
    
    return m_tol;
  }
  
  /**
   * Set the value of tolerance parameter.
   * @param v  Value to assign to tolerance parameter.
   */
  public void setToleranceParameter(double v) {
    
    m_tol = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String epsilonTipText() {
    return "The epsilon for round-off error (shouldn't be changed).";
  }
  
  /**
   * Get the value of epsilon.
   * @return Value of epsilon.
   */
  public double getEpsilon() {
    
    return m_eps;
  }
  
  /**
   * Set the value of epsilon.
   * @param v  Value to assign to epsilon.
   */
  public void setEpsilon(double v) {
    
    m_eps = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String filterTypeTipText() {
    return "Determines how/if the data will be transformed.";
  }
  
  /**
   * Gets how the training data will be transformed. Will be one of
   * FILTER_NORMALIZE, FILTER_STANDARDIZE, FILTER_NONE.
   *
   * @return the filtering mode
   */
  public SelectedTag getFilterType() {

    return new SelectedTag(m_filterType, TAGS_FILTER);
  }
  
  /**
   * Sets how the training data will be transformed. Should be one of
   * FILTER_NORMALIZE, FILTER_STANDARDIZE, FILTER_NONE.
   *
   * @param newType the new filtering mode
   */
  public void setFilterType(SelectedTag newType) {
    
    if (newType.getTags() == TAGS_FILTER) {
      m_filterType = newType.getSelectedTag().getID();
    }
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String buildLogisticModelsTipText() {
    return "Whether to fit logistic models to the outputs (for proper "
      + "probability estimates).";
  }

  /**
   * Get the value of buildLogisticModels.
   *
   * @return Value of buildLogisticModels.
   */
  public boolean getBuildLogisticModels() {
    
    return m_fitLogisticModels;
  }
  
  /**
   * Set the value of buildLogisticModels.
   *
   * @param newbuildLogisticModels Value to assign to buildLogisticModels.
   */
  public void setBuildLogisticModels(boolean newbuildLogisticModels) {
    
    m_fitLogisticModels = newbuildLogisticModels;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numFoldsTipText() {
    return "The number of folds for cross-validation used to generate "
      + "training data for logistic models (-1 means use training data).";
  }
  
  /**
   * Get the value of numFolds.
   *
   * @return Value of numFolds.
   */
  public int getNumFolds() {
    
    return m_numFolds;
  }
  
  /**
   * Set the value of numFolds.
   *
   * @param newnumFolds Value to assign to numFolds.
   */
  public void setNumFolds(int newnumFolds) {
    
    m_numFolds = newnumFolds;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String randomSeedTipText() {
    return "Random number seed for the cross-validation.";
  }
  
  /**
   * Get the value of randomSeed.
   *
   * @return Value of randomSeed.
   */
  public int getRandomSeed() {
    
    return m_randomSeed;
  }
  
  /**
   * Set the value of randomSeed.
   *
   * @param newrandomSeed Value to assign to randomSeed.
   */
  public void setRandomSeed(int newrandomSeed) {
    
    m_randomSeed = newrandomSeed;
  }
  
  /**
   * Prints out the classifier.
   *
   * @return a description of the classifier as a string
   */
  public String toString() {
    
    StringBuffer text = new StringBuffer();
    
    if ((m_classAttribute == null)) {
      return "SMO: No model built yet.";
    }
    try {
      text.append("SMO\n\n");
      text.append("Kernel used:\n  " + m_kernel.toString() + "\n\n");
      
      for (int i = 0; i < m_classAttribute.numValues(); i++) {
	for (int j = i + 1; j < m_classAttribute.numValues(); j++) {
	  text.append("Classifier for classes: " + 
		      m_classAttribute.value(i) + ", " +
		      m_classAttribute.value(j) + "\n\n");
	  text.append(m_classifiers[i][j]);
	  if (m_fitLogisticModels) {
	    text.append("\n\n");
	    if ( m_classifiers[i][j].m_logistic == null) {
	      text.append("No logistic model has been fit.\n");
	    } else {
	      text.append(m_classifiers[i][j].m_logistic);
	    }
	  }
	  text.append("\n\n");
	}
      }
    } catch (Exception e) {
      return "Can't print SMO classifier.";
    }
    
    return text.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6024 $");
  }
  
  /**
   * Main method for testing this class.
   */
  public static void main(String[] argv) {
    runClassifier(new SMO(), argv);
  }
}

