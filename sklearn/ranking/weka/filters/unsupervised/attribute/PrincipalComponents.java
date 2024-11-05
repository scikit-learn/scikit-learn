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
 * PrincipalComponents.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.unsupervised.attribute;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.core.matrix.EigenvalueDecomposition;
import weka.core.matrix.Matrix;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Performs a principal components analysis and transformation of the data.<br/>
 * Dimensionality reduction is accomplished by choosing enough eigenvectors to account for some percentage of the variance in the original data -- default 0.95 (95%).<br/>
 * Based on code of the attribute selection scheme 'PrincipalComponents' by Mark Hall and Gabi Schmidberger.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C
 *  Center (rather than standardize) the
 *  data and compute PCA using the covariance (rather
 *   than the correlation) matrix.</pre>
 * 
 * <pre> -R &lt;num&gt;
 *  Retain enough PC attributes to account
 *  for this proportion of variance in the original data.
 *  (default: 0.95)</pre>
 * 
 * <pre> -A &lt;num&gt;
 *  Maximum number of attributes to include in 
 *  transformed attribute names.
 *  (-1 = include all, default: 5)</pre>
 * 
 * <pre> -M &lt;num&gt;
 *  Maximum number of PC attributes to retain.
 *  (-1 = include all, default: -1)</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz) -- attribute selection code
 * @author Gabi Schmidberger (gabi@cs.waikato.ac.nz) -- attribute selection code
 * @author fracpete (fracpete at waikato dot ac dot nz) -- filter code
 * @version $Revision: 6714 $
 */
public class PrincipalComponents
  extends Filter
  implements OptionHandler, UnsupervisedFilter {

  /** for serialization. */
  private static final long serialVersionUID = -5649876869480249303L;

  /** The data to transform analyse/transform. */
  protected Instances m_TrainInstances;

  /** Keep a copy for the class attribute (if set). */
  protected Instances m_TrainCopy;

  /** The header for the transformed data format. */
  protected Instances m_TransformedFormat;

  /** Data has a class set. */
  protected boolean m_HasClass;

  /** Class index. */
  protected int m_ClassIndex;

  /** Number of attributes. */
  protected int m_NumAttribs;

  /** Number of instances. */
  protected int m_NumInstances;

  /** Correlation matrix for the original data. */
  protected double[][] m_Correlation;
  
  /** 
   * If true, center (rather than standardize) the data and
   * compute PCA from covariance (rather than correlation)
   * matrix.
   */
  private boolean m_center = false;

  /** Will hold the unordered linear transformations of the (normalized)
      original data. */
  protected double[][] m_Eigenvectors;

  /** Eigenvalues for the corresponding eigenvectors. */
  protected double[] m_Eigenvalues = null;

  /** Sorted eigenvalues. */
  protected int[] m_SortedEigens;

  /** sum of the eigenvalues. */
  protected double m_SumOfEigenValues = 0.0;

  /** Filters for replacing missing values. */
  protected ReplaceMissingValues m_ReplaceMissingFilter;
  
  /** Filter for turning nominal values into numeric ones. */
  protected NominalToBinary m_NominalToBinaryFilter;
  
  /** Filter for removing class attribute, nominal attributes with 0 or 1 value. */
  protected Remove m_AttributeFilter;
  
  /** Filter for standardizing the data */
  protected Standardize m_standardizeFilter;
  
  /** Filter for centering the data */
  protected Center m_centerFilter;

  /** The number of attributes in the pc transformed data. */
  protected int m_OutputNumAtts = -1;  

  /** the amount of varaince to cover in the original data when
      retaining the best n PC's. */
  protected double m_CoverVariance = 0.95;

  /** maximum number of attributes in the transformed attribute name. */
  protected int m_MaxAttrsInName = 5;

  /** maximum number of attributes in the transformed data (-1 for all). */
  protected int m_MaxAttributes = -1;

  /**
   * Returns a string describing this filter.
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Performs a principal components analysis and transformation of "
      + "the data.\n"
      + "Dimensionality reduction is accomplished by choosing enough eigenvectors "
      + "to account for some percentage of the variance in the original data -- "
      + "default 0.95 (95%).\n"
      + "Based on code of the attribute selection scheme 'PrincipalComponents' "
      + "by Mark Hall and Gabi Schmidberger.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
    
    result.addElement(new Option("\tCenter (rather than standardize) the" +
        "\n\tdata and compute PCA using the covariance (rather" +
        "\n\t than the correlation) matrix.",
        "C", 0, "-C"));

    result.addElement(new Option(
	"\tRetain enough PC attributes to account\n"
	+"\tfor this proportion of variance in the original data.\n"
	+ "\t(default: 0.95)",
	"R", 1, "-R <num>"));

    result.addElement(new Option(
	"\tMaximum number of attributes to include in \n"
	+ "\ttransformed attribute names.\n"
	+ "\t(-1 = include all, default: 5)", 
	"A", 1, "-A <num>"));

    result.addElement(new Option(
	"\tMaximum number of PC attributes to retain.\n"
	+ "\t(-1 = include all, default: -1)", 
	"M", 1, "-M <num>"));

    return result.elements();
  }

  /**
   * Parses a list of options for this object. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -C
   *  Center (rather than standardize) the
   *  data and compute PCA using the covariance (rather
   *   than the correlation) matrix.</pre>
   * 
   * <pre> -R &lt;num&gt;
   *  Retain enough PC attributes to account
   *  for this proportion of variance in the original data.
   *  (default: 0.95)</pre>
   * 
   * <pre> -A &lt;num&gt;
   *  Maximum number of attributes to include in 
   *  transformed attribute names.
   *  (-1 = include all, default: 5)</pre>
   * 
   * <pre> -M &lt;num&gt;
   *  Maximum number of PC attributes to retain.
   *  (-1 = include all, default: -1)</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String        tmpStr;

    tmpStr = Utils.getOption('R', options);
    if (tmpStr.length() != 0)
      setVarianceCovered(Double.parseDouble(tmpStr));
    else
      setVarianceCovered(0.95);

    tmpStr = Utils.getOption('A', options);
    if (tmpStr.length() != 0)
      setMaximumAttributeNames(Integer.parseInt(tmpStr));
    else
      setMaximumAttributeNames(5);

    tmpStr = Utils.getOption('M', options);
    if (tmpStr.length() != 0)
      setMaximumAttributes(Integer.parseInt(tmpStr));
    else
      setMaximumAttributes(-1);

    setCenterData(Utils.getFlag('C', options));    
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;

    result = new Vector<String>();

    result.add("-R");
    result.add("" + getVarianceCovered());

    result.add("-A");
    result.add("" + getMaximumAttributeNames());

    result.add("-M");
    result.add("" + getMaximumAttributes());

    if (getCenterData())
      result.add("-C");

    return result.toArray(new String[result.size()]);
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String centerDataTipText() {
    return "Center (rather than standardize) the data. PCA will "
      + "be computed from the covariance (rather than correlation) "
      + "matrix";
  }
  
  /**
   * Set whether to center (rather than standardize)
   * the data. If set to true then PCA is computed
   * from the covariance rather than correlation matrix.
   * 
   * @param center true if the data is to be
   * centered rather than standardized
   */
  public void setCenterData(boolean center) {
    m_center = center;
  }
  
  /**
   * Get whether to center (rather than standardize)
   * the data. If true then PCA is computed
   * from the covariance rather than correlation matrix. 
   * 
   * @return true if the data is to be centered rather
   * than standardized.
   */
  public boolean getCenterData() {
    return m_center;
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String varianceCoveredTipText() {
    return "Retain enough PC attributes to account for this proportion of variance.";
  }

  /**
   * Sets the amount of variance to account for when retaining
   * principal components.
   * 
   * @param value 	the proportion of total variance to account for
   */
  public void setVarianceCovered(double value) {
    m_CoverVariance = value;
  }

  /**
   * Gets the proportion of total variance to account for when
   * retaining principal components.
   * 
   * @return 		the proportion of variance to account for
   */
  public double getVarianceCovered() {
    return m_CoverVariance;
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String maximumAttributeNamesTipText() {
    return "The maximum number of attributes to include in transformed attribute names.";
  }

  /**
   * Sets maximum number of attributes to include in
   * transformed attribute names.
   * 
   * @param value 	the maximum number of attributes
   */
  public void setMaximumAttributeNames(int value) {
    m_MaxAttrsInName = value;
  }

  /**
   * Gets maximum number of attributes to include in
   * transformed attribute names.
   * 
   * @return 		the maximum number of attributes
   */
  public int getMaximumAttributeNames() {
    return m_MaxAttrsInName;
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String maximumAttributesTipText() {
    return "The maximum number of PC attributes to retain.";
  }

  /**
   * Sets maximum number of PC attributes to retain.
   * 
   * @param value 	the maximum number of attributes
   */
  public void setMaximumAttributes(int value) {
    m_MaxAttributes = value;
  }

  /**
   * Gets maximum number of PC attributes to retain.
   * 
   * @return 		the maximum number of attributes
   */
  public int getMaximumAttributes() {
    return m_MaxAttributes;
  }

  /**
   * Returns the capabilities of this evaluator.
   *
   * @return            the capabilities of this evaluator
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enable(Capability.PREFERENCE_ATTRIBUTE);
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.DATE_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enable(Capability.RANKING);
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    result.enable(Capability.NO_CLASS);

    return result;
  }

  /**
   * Determines the output format based on the input format and returns 
   * this. In case the output format cannot be returned immediately, i.e.,
   * immediateOutputFormat() returns false, then this method will be called
   * from batchFinished().
   *
   * @param inputFormat     the input format to base the output format on
   * @return                the output format
   * @throws Exception      in case the determination goes wrong
   * @see   #hasImmediateOutputFormat()
   * @see   #batchFinished()
   */
  protected Instances determineOutputFormat(Instances inputFormat) throws Exception {
    double 		cumulative;
    FastVector 		attributes;
    int 		i;
    int 		j;
    StringBuffer 	attName;
    double[] 		coeff_mags;
    int 		num_attrs;
    int[] 		coeff_inds;
    double 		coeff_value;
    int			numAttsLowerBound;
    
    if (m_Eigenvalues == null)
      return inputFormat;

    if (m_MaxAttributes > 0)
      numAttsLowerBound = m_NumAttribs - m_MaxAttributes;
    else
      numAttsLowerBound = 0;
    if (numAttsLowerBound < 0)
      numAttsLowerBound = 0;
    
    cumulative = 0.0;
    attributes = new FastVector();
    for (i = m_NumAttribs - 1; i >= numAttsLowerBound; i--) {
      attName = new StringBuffer();
      // build array of coefficients
      coeff_mags = new double[m_NumAttribs];
      for (j = 0; j < m_NumAttribs; j++)
	coeff_mags[j] = -Math.abs(m_Eigenvectors[j][m_SortedEigens[i]]);
      num_attrs = (m_MaxAttrsInName > 0) ? Math.min(m_NumAttribs, m_MaxAttrsInName) : m_NumAttribs;

      // this array contains the sorted indices of the coefficients
      if (m_NumAttribs > 0) {
	// if m_maxAttrsInName > 0, sort coefficients by decreasing magnitude
	coeff_inds = Utils.sort(coeff_mags);
      }
      else {
	// if  m_maxAttrsInName <= 0, use all coeffs in original order
	coeff_inds = new int[m_NumAttribs];
	for (j = 0; j < m_NumAttribs; j++)
	  coeff_inds[j] = j;
      }
      // build final attName string
      for (j = 0; j < num_attrs; j++) {
	coeff_value = m_Eigenvectors[coeff_inds[j]][m_SortedEigens[i]];
	if (j > 0 && coeff_value >= 0)
	  attName.append("+");
	attName.append(
	    Utils.doubleToString(coeff_value,5,3) 
	    + inputFormat.attribute(coeff_inds[j]).name());
      }
      if (num_attrs < m_NumAttribs)
	attName.append("...");

      attributes.addElement(new Attribute(attName.toString()));
      cumulative += m_Eigenvalues[m_SortedEigens[i]];

      if ((cumulative / m_SumOfEigenValues) >= m_CoverVariance)
	break;
    }

    if (m_HasClass)
      attributes.addElement(m_TrainCopy.classAttribute().copy());

    Instances outputFormat = 
      new Instances(
	  m_TrainCopy.relationName() + "_principal components", attributes, 0);

    // set the class to be the last attribute if necessary
    if (m_HasClass)
      outputFormat.setClassIndex(outputFormat.numAttributes() - 1);

    m_OutputNumAtts = outputFormat.numAttributes();
    
    return outputFormat;
  }
  
  protected void fillCovariance() throws Exception {    
    
    if (!m_center) {
      fillCorrelation();
      return;
    }
    
    double[] att = new double[m_TrainInstances.numInstances()];
    
    // now center the data by subtracting the mean
    m_centerFilter = new Center();
    m_centerFilter.setInputFormat(m_TrainInstances);
    m_TrainInstances = Filter.useFilter(m_TrainInstances, m_centerFilter);
    
    // now compute the covariance matrix
    m_Correlation = new double[m_NumAttribs][m_NumAttribs];
    
    for (int i = 0; i < m_NumAttribs; i++) {
      for (int j = 0; j < m_NumAttribs; j++) {
        
        double cov = 0;
        for (int k = 0; k < m_NumInstances; k++) {
       
          if (i == j) {
            cov += (m_TrainInstances.instance(k).value(i) *
                m_TrainInstances.instance(k).value(i));
          } else {
          cov += (m_TrainInstances.instance(k).value(i) *
              m_TrainInstances.instance(k).value(j));
          }
        }
        
        cov /= (double)(m_TrainInstances.numInstances() - 1);
        m_Correlation[i][j] = cov;
        m_Correlation[j][i] = cov;                
      }
    }
  }

  /**
   * Fill the correlation matrix.
   */
  protected void fillCorrelation() throws Exception {
    int		i;
    int		j;
    int		k;
    double[] 	att1;
    double[] 	att2;
    double 	corr;
    
    m_Correlation = new double[m_NumAttribs][m_NumAttribs];
    att1          = new double [m_NumInstances];
    att2          = new double [m_NumInstances];

    for (i = 0; i < m_NumAttribs; i++) {
      for (j = 0; j < m_NumAttribs; j++) {
        for (k = 0; k < m_NumInstances; k++) {
          att1[k] = m_TrainInstances.instance(k).value(i);
          att2[k] = m_TrainInstances.instance(k).value(j);
        }
	if (i == j) {
	  m_Correlation[i][j] = 1.0;
	}
	else {	  
	  corr = Utils.correlation(att1,att2,m_NumInstances);
	  m_Correlation[i][j] = corr;
	  m_Correlation[j][i] = corr;
	}
      }
    }
    
    // now standardize the input data
    m_standardizeFilter = new Standardize();
    m_standardizeFilter.setInputFormat(m_TrainInstances);
    m_TrainInstances = Filter.useFilter(m_TrainInstances, m_standardizeFilter);
  }

  /**
   * Transform an instance in original (unormalized) format.
   * 
   * @param instance 	an instance in the original (unormalized) format
   * @return 		a transformed instance
   * @throws Exception 	if instance can't be transformed
   */
  protected Instance convertInstance(Instance instance) throws Exception {
    Instance	result;
    double[] 	newVals;
    Instance 	tempInst;
    double 	cumulative;
    int		i;
    int		j;
    double 	tempval;
    int		numAttsLowerBound;
    
    newVals  = new double[m_OutputNumAtts];
    tempInst = (Instance) instance.copy();

    m_ReplaceMissingFilter.input(tempInst);
    m_ReplaceMissingFilter.batchFinished();
    tempInst = m_ReplaceMissingFilter.output();    

    m_NominalToBinaryFilter.input(tempInst);
    m_NominalToBinaryFilter.batchFinished();
    tempInst = m_NominalToBinaryFilter.output();

    if (m_AttributeFilter != null) {
      m_AttributeFilter.input(tempInst);
      m_AttributeFilter.batchFinished();
      tempInst = m_AttributeFilter.output();
    }
    
    if (!m_center) {
      m_standardizeFilter.input(tempInst);
      m_standardizeFilter.batchFinished();
      tempInst = m_standardizeFilter.output();
    } else {
      m_centerFilter.input(tempInst);
      m_centerFilter.batchFinished();
      tempInst = m_centerFilter.output();
    }

    if (m_HasClass)
      newVals[m_OutputNumAtts - 1] = instance.value(instance.classIndex());

    if (m_MaxAttributes > 0)
      numAttsLowerBound = m_NumAttribs - m_MaxAttributes;
    else
      numAttsLowerBound = 0;
    if (numAttsLowerBound < 0)
      numAttsLowerBound = 0;
    
    cumulative = 0;
    for (i = m_NumAttribs - 1; i >= numAttsLowerBound; i--) {
      tempval = 0.0;
      for (j = 0; j < m_NumAttribs; j++)
	tempval += m_Eigenvectors[j][m_SortedEigens[i]] * tempInst.value(j);

      newVals[m_NumAttribs - i - 1] = tempval;
      cumulative += m_Eigenvalues[m_SortedEigens[i]];
      if ((cumulative / m_SumOfEigenValues) >= m_CoverVariance)
	break;
    }

    // create instance
    if (instance instanceof SparseInstance)
      result = new SparseInstance(instance.weight(), newVals);
    //RANKING BEGIN
    if(instance instanceof PreferenceDenseInstance){
    	PreferenceDenseInstance instanceCopy = (PreferenceDenseInstance)instance;
    	result = new PreferenceDenseInstance(instance.weight(), newVals, instanceCopy.getHashMap());
    }
    //RANKING END
    else{
      result = new DenseInstance(instance.weight(), newVals);
    }
    
    return result;
  }

  /**
   * Initializes the filter with the given input data.
   *
   * @param instances   the data to process
   * @throws Exception  in case the processing goes wrong
   * @see               #batchFinished()
   */
  protected void setup(Instances instances) throws Exception {
    int				i;
    int				j;
    Vector<Integer> 		deleteCols;
    int[] 			todelete;
    double[][] 			v;
    Matrix 			corr;
    EigenvalueDecomposition 	eig;
    Matrix 			V;
    
    m_TrainInstances = new Instances(instances);

    // make a copy of the training data so that we can get the class
    // column to append to the transformed data (if necessary)
    m_TrainCopy = new Instances(m_TrainInstances, 0);

    m_ReplaceMissingFilter = new ReplaceMissingValues();
    m_ReplaceMissingFilter.setInputFormat(m_TrainInstances);
    m_TrainInstances = Filter.useFilter(m_TrainInstances, m_ReplaceMissingFilter);

    m_NominalToBinaryFilter = new NominalToBinary();
    m_NominalToBinaryFilter.setInputFormat(m_TrainInstances);
    m_TrainInstances = Filter.useFilter(m_TrainInstances, m_NominalToBinaryFilter);

    // delete any attributes with only one distinct value or are all missing
    deleteCols = new Vector<Integer>();
    for (i = 0; i < m_TrainInstances.numAttributes(); i++) {
      if (m_TrainInstances.numDistinctValues(i) <= 1)
	deleteCols.addElement(i);
    }

    if (m_TrainInstances.classIndex() >=0) {
      // get rid of the class column
      m_HasClass = true;
      m_ClassIndex = m_TrainInstances.classIndex();
      deleteCols.addElement(new Integer(m_ClassIndex));
    }

    // remove columns from the data if necessary
    if (deleteCols.size() > 0) {
      m_AttributeFilter = new Remove();
      todelete = new int [deleteCols.size()];
      for (i = 0; i < deleteCols.size(); i++)
	todelete[i] = ((Integer)(deleteCols.elementAt(i))).intValue();
      m_AttributeFilter.setAttributeIndicesArray(todelete);
      m_AttributeFilter.setInvertSelection(false);
      m_AttributeFilter.setInputFormat(m_TrainInstances);
      m_TrainInstances = Filter.useFilter(m_TrainInstances, m_AttributeFilter);
    }

    // can evaluator handle the processed data ? e.g., enough attributes?
    getCapabilities().testWithFail(m_TrainInstances);

    m_NumInstances = m_TrainInstances.numInstances();
    m_NumAttribs   = m_TrainInstances.numAttributes();

    //fillCorrelation();
    fillCovariance();

    // get eigen vectors/values
    corr = new Matrix(m_Correlation);
    eig  = corr.eig();
    V    = eig.getV();
    v    = new double[m_NumAttribs][m_NumAttribs];
    for (i = 0; i < v.length; i++) {
      for (j = 0; j < v[0].length; j++)
        v[i][j] = V.get(i, j);
    }
    m_Eigenvectors = (double[][]) v.clone();
    m_Eigenvalues  = (double[]) eig.getRealEigenvalues().clone();

    // any eigenvalues less than 0 are not worth anything --- change to 0
    for (i = 0; i < m_Eigenvalues.length; i++) {
      if (m_Eigenvalues[i] < 0)
	m_Eigenvalues[i] = 0.0;
    }
    m_SortedEigens     = Utils.sort(m_Eigenvalues);
    m_SumOfEigenValues = Utils.sum(m_Eigenvalues);

    m_TransformedFormat = determineOutputFormat(m_TrainInstances);
    setOutputFormat(m_TransformedFormat);
    
    m_TrainInstances = null;
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo 	an Instances object containing the input 
   * 				instance structure (any instances contained 
   * 				in the object are ignored - only the structure 
   * 				is required).
   * @return 			true if the outputFormat may be collected 
   * 				immediately
   * @throws Exception 		if the input format can't be set successfully
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {
    super.setInputFormat(instanceInfo);

    m_Eigenvalues           = null;
    m_OutputNumAtts         = -1;
    m_AttributeFilter       = null;
    m_NominalToBinaryFilter = null;
    m_SumOfEigenValues      = 0.0;
    
    return false;
  }

  /**
   * Input an instance for filtering. Filter requires all
   * training instances be read before producing output.
   *
   * @param instance 			the input instance
   * @return 				true if the filtered instance may now be
   * 					collected with output().
   * @throws IllegalStateException 	if no input format has been set
   * @throws Exception 			if conversion fails
   */
  public boolean input(Instance instance) throws Exception {
    Instance 	inst;
    
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    if (isNewBatch()) {
      resetQueue();
      m_NewBatch = false;
    }
    
    if (isFirstBatchDone()) {
      inst = convertInstance(instance);
      inst.setDataset(getOutputFormat());
      push(inst);
      return true;
    }
    else {
      bufferInput(instance);
      return false;
    }
  }

  /**
   * Signify that this batch of input to the filter is finished.
   *
   * @return true 			if there are instances pending output
   * @throws NullPointerException 	if no input structure has been defined,
   * @throws Exception 			if there was a problem finishing the batch.
   */
  public boolean batchFinished() throws Exception {
    int		i;
    Instances	insts;
    Instance	inst;
    
    if (getInputFormat() == null)
      throw new NullPointerException("No input instance format defined");

    insts = getInputFormat();

    if (!isFirstBatchDone())
      setup(insts);
    
    for (i = 0; i < insts.numInstances(); i++) {
      inst = convertInstance(insts.instance(i));
      inst.setDataset(getOutputFormat());
      push(inst);
    }
    
    flushInput();
    m_NewBatch       = true;
    m_FirstBatchDone = true;
    
    return (numPendingOutput() != 0);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6714 $");
  }

  /**
   * Main method for running this filter.
   *
   * @param args 	should contain arguments to the filter: use -h for help
   */
  public static void main(String[] args) {
    runFilter(new PrincipalComponents(), args);
  }
}

