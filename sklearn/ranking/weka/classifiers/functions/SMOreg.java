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
 *    SMOreg.java
 *    Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.functions.supportVector.Kernel;
import weka.classifiers.functions.supportVector.PolyKernel;
import weka.classifiers.functions.supportVector.RegOptimizer;
import weka.classifiers.functions.supportVector.RegSMOImproved;
import weka.core.AdditionalMeasureProducer;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
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

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * SMOreg implements the support vector machine for regression. The parameters can be learned using various algorithms. The algorithm is selected by setting the RegOptimizer. The most popular algorithm (RegSMOImproved) is due to Shevade, Keerthi et al and this is the default RegOptimizer.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * S.K. Shevade, S.S. Keerthi, C. Bhattacharyya, K.R.K. Murthy: Improvements to the SMO Algorithm for SVM Regression. In: IEEE Transactions on Neural Networks, 1999.<br/>
 * <br/>
 * A.J. Smola, B. Schoelkopf (1998). A tutorial on support vector regression.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Shevade1999,
 *    author = {S.K. Shevade and S.S. Keerthi and C. Bhattacharyya and K.R.K. Murthy},
 *    booktitle = {IEEE Transactions on Neural Networks},
 *    title = {Improvements to the SMO Algorithm for SVM Regression},
 *    year = {1999},
 *    PS = {http://guppy.mpe.nus.edu.sg/\~mpessk/svm/ieee_smo_reg.ps.gz}
 * }
 * 
 * &#64;techreport{Smola1998,
 *    author = {A.J. Smola and B. Schoelkopf},
 *    note = {NeuroCOLT2 Technical Report NC2-TR-1998-030},
 *    title = {A tutorial on support vector regression},
 *    year = {1998}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C &lt;double&gt;
 *  The complexity constant C.
 *  (default 1)</pre>
 * 
 * <pre> -N
 *  Whether to 0=normalize/1=standardize/2=neither.
 *  (default 0=normalize)</pre>
 * 
 * <pre> -I &lt;classname and parameters&gt;
 *  Optimizer class used for solving quadratic optimization problem
 *  (default weka.classifiers.functions.supportVector.RegSMOImproved)</pre>
 * 
 * <pre> -K &lt;classname and parameters&gt;
 *  The Kernel to use.
 *  (default: weka.classifiers.functions.supportVector.PolyKernel)</pre>
 * 
 * <pre> 
 * Options specific to optimizer ('-I') weka.classifiers.functions.supportVector.RegSMOImproved:
 * </pre>
 * 
 * <pre> -T &lt;double&gt;
 *  The tolerance parameter for checking the stopping criterion.
 *  (default 0.001)</pre>
 * 
 * <pre> -V
 *  Use variant 1 of the algorithm when true, otherwise use variant 2.
 *  (default true)</pre>
 * 
 * <pre> -P &lt;double&gt;
 *  The epsilon for round-off error.
 *  (default 1.0e-12)</pre>
 * 
 * <pre> -L &lt;double&gt;
 *  The epsilon parameter in epsilon-insensitive loss function.
 *  (default 1.0e-3)</pre>
 * 
 * <pre> -W &lt;double&gt;
 *  The random number seed.
 *  (default 1)</pre>
 * 
 * <pre> 
 * Options specific to kernel ('-K') weka.classifiers.functions.supportVector.PolyKernel:
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
 * @author  Remco Bouckaert (remco@cs.waikato.ac.nz,rrb@xm.co.nz)
 * @version $Revision: 5928 $
 */
public class SMOreg 
  extends AbstractClassifier 
  implements WeightedInstancesHandler, AdditionalMeasureProducer, 
             TechnicalInformationHandler {
  
  /** for serialization */
  private static final long serialVersionUID = -7149606251113102827L;
  
  /** The filter to apply to the training data: Normalzie */
  public static final int FILTER_NORMALIZE = 0;
  /** The filter to apply to the training data: Standardize */
  public static final int FILTER_STANDARDIZE = 1;
  /** The filter to apply to the training data: None */
  public static final int FILTER_NONE = 2;
  /** The filter to apply to the training data */
  public static final Tag[] TAGS_FILTER =
  {
    new Tag(FILTER_NORMALIZE, "Normalize training data"),
    new Tag(FILTER_STANDARDIZE, "Standardize training data"),
    new Tag(FILTER_NONE, "No normalization/standardization"),
  };
  
  /** Whether to normalize/standardize/neither */
  protected int m_filterType = FILTER_NORMALIZE;
  
  /** The filter used to make attributes numeric. */
  protected NominalToBinary m_NominalToBinary;
  
  /** The filter used to standardize/normalize all values. */
  protected Filter m_Filter = null;
  
  /** The filter used to get rid of missing values. */
  protected ReplaceMissingValues m_Missing;
  
  /** Only numeric attributes in the dataset? If so, less need to filter */
  protected boolean m_onlyNumeric;
  
  /** capacity parameter **/
  protected double m_C = 1.0;
  
  /** coefficients used by normalization filter for doing its linear transformation 
   * so that result = svmoutput * m_x1 + m_x0 **/
  protected double m_x1 = 1.0;
  protected double m_x0 = 0.0;
  
  /** contains the algorithm used for learning **/
  protected RegOptimizer m_optimizer = new RegSMOImproved();

  /** the configured kernel */
  protected Kernel m_kernel = new PolyKernel();
  
  /**
   * Returns a string describing classifier
   * 
   * @return 		a description suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return
        "SMOreg implements the support vector machine for regression. "
      + "The parameters can be learned using various algorithms. The "
      + "algorithm is selected by setting the RegOptimizer. The most "
      + "popular algorithm (" 
      + RegSMOImproved.class.getName().replaceAll(".*\\.", "") 
      + ") is due to Shevade, Keerthi " 
      + "et al and this is the default RegOptimizer.\n\n"
      + "For more information see:\n\n"
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
    TechnicalInformation	additional;
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "S.K. Shevade and S.S. Keerthi and C. Bhattacharyya and K.R.K. Murthy");
    result.setValue(Field.TITLE, "Improvements to the SMO Algorithm for SVM Regression");
    result.setValue(Field.BOOKTITLE, "IEEE Transactions on Neural Networks");
    result.setValue(Field.YEAR, "1999");
    result.setValue(Field.PS, "http://guppy.mpe.nus.edu.sg/~mpessk/svm/ieee_smo_reg.ps.gz");

    additional = result.add(Type.TECHREPORT);
    additional.setValue(Field.AUTHOR, "A.J. Smola and B. Schoelkopf");
    additional.setValue(Field.TITLE, "A tutorial on support vector regression");
    additional.setValue(Field.NOTE, "NeuroCOLT2 Technical Report NC2-TR-1998-030");
    additional.setValue(Field.YEAR, "1998");

    return result;
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Enumeration enm;
    Vector result = new Vector();
    
    result.addElement(new Option(
	"\tThe complexity constant C.\n"
	+ "\t(default 1)", 
	"C", 1, "-C <double>"));
    
    result.addElement(new Option(
	"\tWhether to 0=normalize/1=standardize/2=neither.\n" 
	+ "\t(default 0=normalize)", 
	"N", 1, "-N"));
    
    result.addElement(new Option(
	"\tOptimizer class used for solving quadratic optimization problem\n" 
	+ "\t(default " + RegSMOImproved.class.getName() + ")",
	"I", 1, "-I <classname and parameters>"));
    
    result.addElement(new Option(
	"\tThe Kernel to use.\n"
	+ "\t(default: weka.classifiers.functions.supportVector.PolyKernel)",
	"K", 1, "-K <classname and parameters>"));

    result.addElement(new Option(
	"",
	"", 0, "\nOptions specific to optimizer ('-I') "
	+ getRegOptimizer().getClass().getName() + ":"));

    enm = ((OptionHandler) getRegOptimizer()).listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    result.addElement(new Option(
	"",
	"", 0, "\nOptions specific to kernel ('-K') "
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
   * <pre> -C &lt;double&gt;
   *  The complexity constant C.
   *  (default 1)</pre>
   * 
   * <pre> -N
   *  Whether to 0=normalize/1=standardize/2=neither.
   *  (default 0=normalize)</pre>
   * 
   * <pre> -I &lt;classname and parameters&gt;
   *  Optimizer class used for solving quadratic optimization problem
   *  (default weka.classifiers.functions.supportVector.RegSMOImproved)</pre>
   * 
   * <pre> -K &lt;classname and parameters&gt;
   *  The Kernel to use.
   *  (default: weka.classifiers.functions.supportVector.PolyKernel)</pre>
   * 
   * <pre> 
   * Options specific to optimizer ('-I') weka.classifiers.functions.supportVector.RegSMOImproved:
   * </pre>
   * 
   * <pre> -T &lt;double&gt;
   *  The tolerance parameter for checking the stopping criterion.
   *  (default 0.001)</pre>
   * 
   * <pre> -V
   *  Use variant 1 of the algorithm when true, otherwise use variant 2.
   *  (default true)</pre>
   * 
   * <pre> -P &lt;double&gt;
   *  The epsilon for round-off error.
   *  (default 1.0e-12)</pre>
   * 
   * <pre> -L &lt;double&gt;
   *  The epsilon parameter in epsilon-insensitive loss function.
   *  (default 1.0e-3)</pre>
   * 
   * <pre> -W &lt;double&gt;
   *  The random number seed.
   *  (default 1)</pre>
   * 
   * <pre> 
   * Options specific to kernel ('-K') weka.classifiers.functions.supportVector.PolyKernel:
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
    String	tmpStr;
    String[]	tmpOptions;
    
    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0) {
      setC(Double.parseDouble(tmpStr));
    } else {
      setC(1.0);
    }
    
    String nString = Utils.getOption('N', options);
    if (nString.length() != 0) {
      setFilterType(new SelectedTag(Integer.parseInt(nString), TAGS_FILTER));
    } else {
      setFilterType(new SelectedTag(FILTER_NORMALIZE, TAGS_FILTER));
    }
    
    tmpStr = Utils.getOption('I', options);
    tmpOptions = Utils.splitOptions(tmpStr);
    if (tmpOptions.length != 0) {
      tmpStr        = tmpOptions[0];
      tmpOptions[0] = "";
      setRegOptimizer(
	  (RegOptimizer) Utils.forName(RegOptimizer.class, tmpStr, tmpOptions));
    }
    else {
      setRegOptimizer(new RegSMOImproved());
    }

    tmpStr     = Utils.getOption('K', options);
    tmpOptions = Utils.splitOptions(tmpStr);
    if (tmpOptions.length != 0) {
      tmpStr        = tmpOptions[0];
      tmpOptions[0] = "";
      setKernel(Kernel.forName(tmpStr, tmpOptions));
    }
    else {
      setKernel(new PolyKernel());
    }
  }
  
  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int       	i;
    Vector    	result;
    String[]  	options;

    result = new Vector();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    result.add("-C");
    result.add("" + getC());
    
    result.add("-N");
    result.add("" + m_filterType);
    
    result.add("-I");
    result.add("" + getRegOptimizer().getClass().getName() + " " + Utils.joinOptions(getRegOptimizer().getOptions()));

    result.add("-K");
    result.add("" + getKernel().getClass().getName() + " " + Utils.joinOptions(getKernel().getOptions()));

    return (String[]) result.toArray(new String[result.size()]);	  
  }
  
  /**
   * Returns default capabilities of the classifier.
   *
   * @return		the capabilities of this classifier
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
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }
  
  /**
   * Method for building the classifier.
   *
   * @param instances the set of training instances
   * @throws Exception if the classifier can't be built successfully
   */
  public void buildClassifier(Instances instances) throws Exception {
    // can classifier handle the data?
    getCapabilities().testWithFail(instances);

    // remove instances with missing class
    instances = new Instances(instances);
    instances.deleteWithMissingClass();
    
    // Removes all the instances with weight equal to 0.
    // MUST be done since condition (8) of Keerthi's paper 
    // is made with the assertion Ci > 0 (See equation (3a).
    Instances data = new Instances(instances, 0);
    for (int i = 0; i < instances.numInstances(); i++) {
      if (instances.instance(i).weight() > 0) {
	data.add(instances.instance(i));
      }
    }
    
    if (data.numInstances() == 0) {
      throw new Exception("No training instances left after removing " + "instance with either a weight null or a missing class!");
    }
    instances = data;
    
    m_onlyNumeric = true;
    for (int i = 0; i < instances.numAttributes(); i++) {
      if (i != instances.classIndex()) {
	if (!instances.attribute(i).isNumeric()) {
	  m_onlyNumeric = false;
	  break;
	}
      }
    }
    m_Missing = new ReplaceMissingValues();
    m_Missing.setInputFormat(instances);
    instances = Filter.useFilter(instances, m_Missing);
    
    if (!m_onlyNumeric) {
      m_NominalToBinary = new NominalToBinary();
      m_NominalToBinary.setInputFormat(instances);
      instances = Filter.useFilter(instances, m_NominalToBinary);
    } else {
      m_NominalToBinary = null;
    }
    
    // retrieve two different class values used to determine filter transformation
    double y0 = instances.instance(0).classValue();
    int index = 1;
    while (index < instances.numInstances() && instances.instance(index).classValue() == y0) {
      index++;
    }
    if (index == instances.numInstances()) {
      // degenerate case, all class values are equal
      // we don't want to deal with this, too much hassle
      throw new Exception("All class values are the same. At least two class values should be different");
    }
    double y1 = instances.instance(index).classValue();
    
    // apply filters
    if (m_filterType == FILTER_STANDARDIZE) {
      m_Filter = new Standardize();
      ((Standardize)m_Filter).setIgnoreClass(true);
      m_Filter.setInputFormat(instances);
      instances = Filter.useFilter(instances, m_Filter);			
    } else if (m_filterType == FILTER_NORMALIZE) {
      m_Filter = new Normalize();
      ((Normalize)m_Filter).setIgnoreClass(true);
      m_Filter.setInputFormat(instances);
      instances = Filter.useFilter(instances, m_Filter);
    } else {
      m_Filter = null;
    }
    if (m_Filter != null) {
      double z0 = instances.instance(0).classValue();
      double z1 = instances.instance(index).classValue();
      m_x1 = (y0-y1) / (z0 - z1); // no division by zero, since y0 != y1 guaranteed => z0 != z1 ???
      m_x0 = (y0 - m_x1 * z0); // = y1 - m_x1 * z1
    } else {
      m_x1 = 1.0;
      m_x0 = 0.0; 			
    }
    
    m_optimizer.setSMOReg(this);
    m_optimizer.buildClassifier(instances);
  }
  
  /**
   * Classifies the given instance using the linear regression function.
   *
   * @param instance the test instance
   * @return the classification
   * @throws Exception if classification can't be done successfully
   */
  public double classifyInstance(Instance instance) throws Exception {
    // Filter instance
    m_Missing.input(instance);
    m_Missing.batchFinished();
    instance = m_Missing.output();
    
    if (!m_onlyNumeric) {
      m_NominalToBinary.input(instance);
      m_NominalToBinary.batchFinished();
      instance = m_NominalToBinary.output();
    }
    
    if (m_Filter != null) {
      m_Filter.input(instance);
      m_Filter.batchFinished();
      instance = m_Filter.output();
    }
    
    double result = m_optimizer.SVMOutput(instance);
    return result * m_x1 + m_x0;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String regOptimizerTipText() {
    return "The learning algorithm.";
  }
  
  /**
   * sets the learning algorithm
   * 
   * @param regOptimizer	the learning algorithm
   */
  public void setRegOptimizer(RegOptimizer regOptimizer) {
    m_optimizer = regOptimizer;
  }
  
  /**
   * returns the learning algorithm
   * 
   * @return		the learning algorithm
   */
  public RegOptimizer getRegOptimizer() {
    return m_optimizer;
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
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String cTipText() {
    return "The complexity parameter C.";
  }
  
  /**
   * Get the value of C.
   *
   * @return 		Value of C.
   */
  public double getC() {
    return m_C;
  }
  
  /**
   * Set the value of C.
   *
   * @param v  		Value to assign to C.
   */
  public void setC(double v) {
    m_C = v;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String filterTypeTipText() {
    return "Determines how/if the data will be transformed.";
  }
  
  /**
   * Gets how the training data will be transformed. Will be one of
   * FILTER_NORMALIZE, FILTER_STANDARDIZE, FILTER_NONE.
   *
   * @return 		the filtering mode
   */
  public SelectedTag getFilterType() {
    return new SelectedTag(m_filterType, TAGS_FILTER);
  }
  
  /**
   * Sets how the training data will be transformed. Should be one of
   * FILTER_NORMALIZE, FILTER_STANDARDIZE, FILTER_NONE.
   *
   * @param newType 	the new filtering mode
   */
  public void setFilterType(SelectedTag newType) {
    if (newType.getTags() == TAGS_FILTER) {
      m_filterType = newType.getSelectedTag().getID();
    }
  }
  
  /**
   * Prints out the classifier.
   *
   * @return a description of the classifier as a string
   */
  public String toString() {
    StringBuffer text = new StringBuffer();
    
    if (m_optimizer == null || !m_optimizer.modelBuilt()) {
      return "SMOreg: No model built yet.";
    }
    
    try {
      text.append(m_optimizer.toString());
    } 
    catch (Exception e) {
      return "Can't print SMVreg classifier.";
    }
    
    return text.toString();
  }
  
  /**
   * Returns an enumeration of the measure names. Additional measures
   * must follow the naming convention of starting with "measure", eg.
   * double measureBlah()
   * 
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector result = new Vector();
    
    result.addElement("measureKernelEvaluations");
    result.addElement("measureCacheHits");
    
    return result.elements();
  }
  
  /**
   * Returns the value of the named measure
   * @param measureName the name of the measure to query for its value
   * @return the value of the named measure
   * @throws IllegalArgumentException if the named measure is not supported
   */
  public double getMeasure(String measureName) {
    if (measureName.equals("measureKernelEvaluations"))
      return measureKernelEvaluations();
    else if (measureName.equals("measureCacheHits"))
      return measureCacheHits();
    else
      throw new IllegalArgumentException("Measure '" +  measureName + "' is not supported!");
  }
  
  /** 
   * number of kernel evaluations used in learing
   * 
   * @return		the number of kernel evaluations
   */
  protected double measureKernelEvaluations() {
    if (m_optimizer != null) {
      return m_optimizer.getKernelEvaluations();
    } else {
      return 0; 
    }
  }
  
  /** 
   * number of kernel cache hits used during learing
   * 
   * @return		the number of kernel cache hits
   */
  protected double measureCacheHits() {
    if (m_optimizer != null) {
      return m_optimizer.getCacheHits();
    } else {
      return 0; 
    }
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
   * Main method for running this classifier.
   * 
   * @param args	the commandline options
   */
  public static void main(String[] args) {
    runClassifier(new SMOreg(), args);
  }
}
