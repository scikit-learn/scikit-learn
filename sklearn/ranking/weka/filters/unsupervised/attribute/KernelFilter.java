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
 * KernelFilter.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import weka.classifiers.functions.supportVector.Kernel;
import weka.classifiers.functions.supportVector.PolyKernel;
import weka.classifiers.functions.supportVector.RBFKernel;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.MathematicalExpression;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SingleIndex;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.converters.ConverterUtils.DataSource;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.AllFilter;
import weka.filters.Filter;
import weka.filters.SimpleBatchFilter;
import weka.filters.UnsupervisedFilter;

import java.io.File;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Converts the given set of predictor variables into a kernel matrix. The class value remains unchangedm, as long as the preprocessing filter doesn't change it.<br/>
 * By default, the data is preprocessed with the Center filter, but the user can choose any filter (NB: one must be careful that the filter does not alter the class attribute unintentionally). With weka.filters.AllFilter the preprocessing gets disabled.<br/>
 * <br/>
 * For more information regarding preprocessing the data, see:<br/>
 * <br/>
 * K.P. Bennett, M.J. Embrechts: An Optimization Perspective on Kernel Partial Least Squares Regression. In: Advances in Learning Theory: Methods, Models and Applications, 227-249, 2003.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Bennett2003,
 *    author = {K.P. Bennett and M.J. Embrechts},
 *    booktitle = {Advances in Learning Theory: Methods, Models and Applications},
 *    editor = {J. Suykens et al.},
 *    pages = {227-249},
 *    publisher = {IOS Press, Amsterdam, The Netherlands},
 *    series = {NATO Science Series, Series III: Computer and System Sciences},
 *    title = {An Optimization Perspective on Kernel Partial Least Squares Regression},
 *    volume = {190},
 *    year = {2003}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns on output of debugging information.</pre>
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
 * <pre> -F &lt;filename&gt;
 *  The file to initialize the filter with (optional).</pre>
 * 
 * <pre> -C &lt;num&gt;
 *  The class index for the file to initialize with,
 *  First and last are valid (optional, default: last).</pre>
 * 
 * <pre> -K &lt;classname and parameters&gt;
 *  The Kernel to use.
 *  (default: weka.classifiers.functions.supportVector.PolyKernel)</pre>
 * 
 * <pre> -kernel-factor
 *  Defines a factor for the kernel.
 *   - RBFKernel: a factor for gamma
 *    Standardize: 1/(2*N)
 *    Normalize..: 6/N
 *  Available parameters are:
 *   N for # of instances, A for # of attributes
 *  (default: 1)</pre>
 * 
 * <pre> -P &lt;classname and parameters&gt;
 *  The Filter used for preprocessing (use weka.filters.AllFilter
 *  to disable preprocessing).
 *  (default: weka.filters.unsupervised.attribute.Center)</pre>
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
 * <pre> 
 * Options specific to preprocessing filter weka.filters.unsupervised.attribute.Center:
 * </pre>
 * 
 * <pre> -unset-class-temporarily
 *  Unsets the class index temporarily before the filter is
 *  applied to the data.
 *  (default: no)</pre>
 * 
 <!-- options-end -->
 *
 * @author Jonathan Miles (jdm18@cs.waikato.ac.nz) 
 * @author FracPete (fracpete at waikato dot ac dot nz) 
 * @version $Revision: 5987 $
 */
public class KernelFilter
  extends SimpleBatchFilter 
  implements UnsupervisedFilter, TechnicalInformationHandler {

  /** for serialization */
  static final long serialVersionUID = 213800899640387499L;

  /** The number of instances in the training data. */
  protected int m_NumTrainInstances;

  /** Kernel to use **/
  protected Kernel m_Kernel = new PolyKernel();

  /** the Kernel which is actually used for computation */
  protected Kernel m_ActualKernel = null;

  /** Turn off all checks and conversions? Turning them off assumes
      that data is purely numeric, doesn't contain any missing values,
      and has a nominal class. Turning them off also means that
      no header information will be stored if the machine is linear. 
      Finally, it also assumes that no instance has a weight equal to 0.*/
  protected boolean m_checksTurnedOff;

  /** The filter used to make attributes numeric. */
  protected NominalToBinary m_NominalToBinary;

  /** The filter used to get rid of missing values. */
  protected ReplaceMissingValues m_Missing;

  /** The dataset to initialize the filter with */
  protected File m_InitFile = new File(System.getProperty("user.dir"));

  /** the class index for the file to initialized with 
   * @see #m_InitFile */
  protected SingleIndex m_InitFileClassIndex = new SingleIndex("last");
  
  /** whether the filter was initialized */
  protected boolean m_Initialized = false;

  /** optimizes the kernel with this formula 
   * (A = # of attributes, N = # of instances)*/
  protected String m_KernelFactorExpression = "1";

  /** the calculated kernel factor
   * @see #m_KernelFactorExpression */
  protected double m_KernelFactor = 1.0;
  
  /** for centering/standardizing the data */
  protected Filter m_Filter = new Center();
  
  /** for centering/standardizing the data (the actual filter to use) */
  protected Filter m_ActualFilter = null;
  
  /**
   * Returns a string describing this filter.
   *
   * @return      a description of the filter suitable for
   *              displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Converts the given set of predictor variables into a kernel matrix. "
      + "The class value remains unchangedm, as long as the preprocessing "
      + "filter doesn't change it.\n"
      + "By default, the data is preprocessed with the Center filter, but the "
      + "user can choose any filter (NB: one must be careful that the filter "
      + "does not alter the class attribute unintentionally). With "
      + "weka.filters.AllFilter the preprocessing gets disabled.\n\n"
      + "For more information regarding preprocessing the data, see:\n\n"
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
    TechnicalInformation	result;
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "K.P. Bennett and M.J. Embrechts");
    result.setValue(Field.TITLE, "An Optimization Perspective on Kernel Partial Least Squares Regression");
    result.setValue(Field.YEAR, "2003");
    result.setValue(Field.EDITOR, "J. Suykens et al.");
    result.setValue(Field.BOOKTITLE, "Advances in Learning Theory: Methods, Models and Applications");
    result.setValue(Field.PAGES, "227-249");
    result.setValue(Field.PUBLISHER, "IOS Press, Amsterdam, The Netherlands");
    result.setValue(Field.SERIES, "NATO Science Series, Series III: Computer and System Sciences");
    result.setValue(Field.VOLUME, "190");
    
    return result;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector        result;
    Enumeration   enm;

    result = new Vector();

    enm = super.listOptions();
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
	"\tThe file to initialize the filter with (optional).",
	"F", 1, "-F <filename>"));

    result.addElement(new Option(
	"\tThe class index for the file to initialize with,\n"
	+ "\tFirst and last are valid (optional, default: last).",
	"C", 1, "-C <num>"));

    result.addElement(new Option(
	"\tThe Kernel to use.\n"
	+ "\t(default: weka.classifiers.functions.supportVector.PolyKernel)",
	"K", 1, "-K <classname and parameters>"));

    result.addElement(new Option(
	"\tDefines a factor for the kernel.\n"
	+ "\t\t- RBFKernel: a factor for gamma\n"
	+ "\t\t\tStandardize: 1/(2*N)\n"
	+ "\t\t\tNormalize..: 6/N\n"
	+ "\tAvailable parameters are:\n"
	+ "\t\tN for # of instances, A for # of attributes\n"
	+ "\t(default: 1)",
	"kernel-factor", 0, "-kernel-factor"));

    result.addElement(new Option(
	"\tThe Filter used for preprocessing (use weka.filters.AllFilter\n"
	+ "\tto disable preprocessing).\n"
	+ "\t(default: " + Center.class.getName() + ")",
	"P", 1, "-P <classname and parameters>"));

    // kernel options
    result.addElement(new Option(
	"",
	"", 0, "\nOptions specific to kernel "
	+ getKernel().getClass().getName() + ":"));
    
    enm = ((OptionHandler) getKernel()).listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    // filter options
    if (getPreprocessing() instanceof OptionHandler) {
      result.addElement(new Option(
	  "",
	  "", 0, "\nOptions specific to preprocessing filter "
	  + getPreprocessing().getClass().getName() + ":"));

      enm = ((OptionHandler) getPreprocessing()).listOptions();
      while (enm.hasMoreElements())
	result.addElement(enm.nextElement());
    }
    
    return result.elements();
  }	  

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int		i;
    Vector	result;
    String[]	options;
    String	tmpStr;

    result = new Vector();
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    if (getChecksTurnedOff())
      result.add("-no-checks");

    if ((getInitFile() != null) && getInitFile().isFile()) {
      result.add("-F");
      result.add("" + getInitFile().getAbsolutePath());

      result.add("-C");
      result.add("" + getInitFileClassIndex());
    }

    result.add("-K");
    result.add("" + getKernel().getClass().getName() + " " + Utils.joinOptions(getKernel().getOptions()));

    result.add("-kernel-factor");
    result.add("" + getKernelFactorExpression());

    result.add("-P");
    tmpStr = getPreprocessing().getClass().getName();
    if (getPreprocessing() instanceof OptionHandler)
      tmpStr += " " + Utils.joinOptions(((OptionHandler) getPreprocessing()).getOptions());
    result.add("" + tmpStr);

    return (String[]) result.toArray(new String[result.size()]);	  
  }	  

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Turns on output of debugging information.</pre>
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
   * <pre> -F &lt;filename&gt;
   *  The file to initialize the filter with (optional).</pre>
   * 
   * <pre> -C &lt;num&gt;
   *  The class index for the file to initialize with,
   *  First and last are valid (optional, default: last).</pre>
   * 
   * <pre> -K &lt;classname and parameters&gt;
   *  The Kernel to use.
   *  (default: weka.classifiers.functions.supportVector.PolyKernel)</pre>
   * 
   * <pre> -kernel-factor
   *  Defines a factor for the kernel.
   *   - RBFKernel: a factor for gamma
   *    Standardize: 1/(2*N)
   *    Normalize..: 6/N
   *  Available parameters are:
   *   N for # of instances, A for # of attributes
   *  (default: 1)</pre>
   * 
   * <pre> -P &lt;classname and parameters&gt;
   *  The Filter used for preprocessing (use weka.filters.AllFilter
   *  to disable preprocessing).
   *  (default: weka.filters.unsupervised.attribute.Center)</pre>
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
   * <pre> 
   * Options specific to preprocessing filter weka.filters.unsupervised.attribute.Center:
   * </pre>
   * 
   * <pre> -unset-class-temporarily
   *  Unsets the class index temporarily before the filter is
   *  applied to the data.
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
    
    setChecksTurnedOff(Utils.getFlag("no-checks", options));

    tmpStr = Utils.getOption('F', options);
    if (tmpStr.length() != 0)
      setInitFile(new File(tmpStr));
    else 
      setInitFile(null);

    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0)
      setInitFileClassIndex(tmpStr);
    else 
      setInitFileClassIndex("last");

    tmpStr     = Utils.getOption('K', options);
    tmpOptions = Utils.splitOptions(tmpStr);
    if (tmpOptions.length != 0) {
      tmpStr        = tmpOptions[0];
      tmpOptions[0] = "";
      setKernel(Kernel.forName(tmpStr, tmpOptions));
    }
    
    tmpStr = Utils.getOption("kernel-factor", options);
    if (tmpStr.length() != 0)
      setKernelFactorExpression(tmpStr);
    else 
      setKernelFactorExpression("1");
    
    tmpStr = Utils.getOption("P", options);
    tmpOptions = Utils.splitOptions(tmpStr);
    if (tmpOptions.length != 0) {
      tmpStr        = tmpOptions[0];
      tmpOptions[0] = "";
      setPreprocessing((Filter) Utils.forName(Filter.class, tmpStr, tmpOptions));
    }
    else {
      setPreprocessing(new Center());
    }

    super.setOptions(options);
  }	  
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String initFileTipText() {
    return "The dataset to initialize the filter with.";
  }

  /**
   * Gets the file to initialize the filter with, can be null.
   *
   * @return 		the file
   */
  public File getInitFile() {
    return m_InitFile;
  }
    
  /**
   * Sets the file to initialize the filter with, can be null.
   *
   * @param value	the file
   */
  public void setInitFile(File value) {
    m_InitFile = value;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String initFileClassIndexTipText() {
    return "The class index of the dataset to initialize the filter with (first and last are valid).";
  }

  /**
   * Gets the class index of the file to initialize the filter with.
   *
   * @return 		the class index
   */
  public String getInitFileClassIndex() {
    return m_InitFileClassIndex.getSingleIndex();
  }
    
  /**
   * Sets class index of the file to initialize the filter with.
   *
   * @param value	the class index
   */
  public void setInitFileClassIndex(String value) {
    m_InitFileClassIndex.setSingleIndex(value);
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
   * Gets the kernel to use.
   *
   * @return 		the kernel
   */
  public Kernel getKernel() {
    return m_Kernel;
  }
    
  /**
   * Sets the kernel to use.
   *
   * @param value	the kernel
   */
  public void setKernel(Kernel value) {
    m_Kernel = value;
  }

  /**
   * Disables or enables the checks (which could be time-consuming). Use with
   * caution!
   * 
   * @param value	if true turns off all checks
   */
  public void setChecksTurnedOff(boolean value) {
    m_checksTurnedOff = value;
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
  public String kernelFactorExpressionTipText() {
    return "The factor for the kernel, with A = # of attributes and N = # of instances.";
  }

  /**
   * Gets the expression for the kernel.
   *
   * @return 		the expression
   */
  public String getKernelFactorExpression() {
    return m_KernelFactorExpression;
  }
    
  /**
   * Sets the expression for the kernel.
   *
   * @param value	the file
   */
  public void setKernelFactorExpression(String value) {
    m_KernelFactorExpression = value;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String preprocessingTipText() {
    return "Sets the filter to use for preprocessing (use the AllFilter for no preprocessing).";
  }

  /**
   * Sets the filter to use for preprocessing (use the AllFilter for no 
   * preprocessing)
   *
   * @param value 	the preprocessing filter
   */
  public void setPreprocessing(Filter value) {
    m_Filter       = value;
    m_ActualFilter = null;
  }

  /**
   * Gets the filter used for preprocessing
   *
   * @return 		the current preprocessing filter.
   */
  public Filter getPreprocessing() {
    return m_Filter;
  }

  /**
   * resets the filter, i.e., m_NewBatch to true and m_FirstBatchDone to
   * false.
   */
  protected void reset() {
    super.reset();
    
    m_Initialized = false;
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
    return new Instances(inputFormat);
  }
  
  /**
   * initializes the filter with the given dataset, i.e., the kernel gets
   * built. Needs to be called before the first call of Filter.useFilter or
   * batchFinished(), if not the -F option (or setInitFile(File) is used).
   * 
   * @param instances	the data to initialize with
   * @throws Exception	if building of kernel fails
   */
  public void initFilter(Instances instances) throws Exception {
    HashMap	symbols;
    
    // determine kernel factor
    symbols = new HashMap();
    symbols.put("A", new Double(instances.numAttributes()));
    symbols.put("N", new Double(instances.numInstances()));
    m_KernelFactor = MathematicalExpression.evaluate(getKernelFactorExpression(), symbols);
    
    // init filters
    if (!m_checksTurnedOff) {
      m_Missing = new ReplaceMissingValues();
      m_Missing.setInputFormat(instances);
      instances = Filter.useFilter(instances, m_Missing); 
    } 
    else {
      m_Missing = null;
    }

    if (getKernel().getCapabilities().handles(Capability.NUMERIC_ATTRIBUTES)) {
	boolean onlyNumeric = true;
	if (!m_checksTurnedOff) {
	  for (int i = 0; i < instances.numAttributes(); i++) {
	    if (i != instances.classIndex()) {
	      if (!instances.attribute(i).isNumeric()) {
		onlyNumeric = false;
		break;
	      }
	    }
	  }
	}
	
	if (!onlyNumeric) {
	  m_NominalToBinary = new NominalToBinary();
	  m_NominalToBinary.setInputFormat(instances);
	  instances = Filter.useFilter(instances, m_NominalToBinary);
	} 
	else {
	  m_NominalToBinary = null;
	}
    }
    else {
      m_NominalToBinary = null;
    }

    if ((m_Filter != null) && (m_Filter.getClass() != AllFilter.class)) {
      m_ActualFilter = Filter.makeCopy(m_Filter);
      m_ActualFilter.setInputFormat(instances);
      instances = Filter.useFilter(instances, m_ActualFilter);
    }
    else {
      m_ActualFilter = null;
    }

    m_NumTrainInstances = instances.numInstances();

    // set factor for kernel
    m_ActualKernel = Kernel.makeCopy(m_Kernel);
    if (m_ActualKernel instanceof RBFKernel)
      ((RBFKernel) m_ActualKernel).setGamma(
	  m_KernelFactor * ((RBFKernel) m_ActualKernel).getGamma());
    // build kernel
    m_ActualKernel.buildKernel(instances);

    m_Initialized = true;
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities 	result;
    
    if (getKernel() == null) {
      result = super.getCapabilities();
      result.disableAll();
    } else {
      result = getKernel().getCapabilities();
    }

    result.setMinimumNumberInstances(0);
    
    result.disable(Capability.RANKING);
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    
    return result;
  }

  /**
   * Processes the given data (may change the provided dataset) and returns
   * the modified version. This method is called in batchFinished().
   *
   * @param instances   the data to process
   * @return            the modified data
   * @throws Exception  in case the processing goes wrong
   * @see               #batchFinished()
   */
  protected Instances process(Instances instances) throws Exception {
    // initializing necessary?
    if (!m_Initialized) {
      // do we have a file to initialize with?
      if ((getInitFile() != null) && getInitFile().isFile()) {
	DataSource source = new DataSource(getInitFile().getAbsolutePath());
	Instances data = source.getDataSet();
	m_InitFileClassIndex.setUpper(data.numAttributes() - 1);
	data.setClassIndex(m_InitFileClassIndex.getIndex());
	initFilter(data);
      }
      else {
	initFilter(instances);
      }
    }

    // apply filters
    if (m_Missing != null)
      instances = Filter.useFilter(instances, m_Missing); 
    if (m_NominalToBinary != null)
      instances = Filter.useFilter(instances, m_NominalToBinary); 
    if (m_ActualFilter != null)
      instances = Filter.useFilter(instances, m_ActualFilter);

    // backup class attribute and remove it
    double[] classes = instances.attributeToDoubleArray(instances.classIndex());
    int classIndex = instances.classIndex();
    instances.setClassIndex(-1);
    instances.deleteAttributeAt(classIndex);

    // generate new header
    FastVector atts = new FastVector();
    for (int j = 0; j < m_NumTrainInstances; j++)
      atts.addElement(new Attribute("Kernel " + j));
    atts.addElement(new Attribute("Class"));
    Instances result = new Instances("Kernel", atts, 0);
    result.setClassIndex(result.numAttributes() - 1);

    // compute matrix
    for (int i = 0; i < instances.numInstances(); i++) {
      double[] k = new double[m_NumTrainInstances + 1];
      
      for (int j = 0; j < m_NumTrainInstances; j++) {
	double v = m_ActualKernel.eval(-1, j, instances.instance(i));
	k[j] = v;
      }
      k[k.length - 1] = classes[i];

      // create new instance
      //RANKING BEGIN
      Instance in;
      if(instances.instance(i) instanceof PreferenceDenseInstance){
    	  PreferenceDenseInstance pdi = (PreferenceDenseInstance)instances.instance(i);
    	  in = new PreferenceDenseInstance(1.0,k,pdi.getHashMap());
      }
      else
    	  in = new DenseInstance(1.0, k);
      //RANKING END
      result.add(in);    
    }

    if (!isFirstBatchDone())
      setOutputFormat(result);
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }

  /**
   * runs the filter with the given arguments
   *
   * @param args      the commandline arguments
   */
  public static void main(String[] args) {
    runFilter(new KernelFilter(), args);
  }
}
