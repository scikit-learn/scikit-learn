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
 * AddClassification.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.supervised.attribute;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Capabilities.Capability;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.WekaException;
import weka.filters.SimpleBatchFilter;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.ObjectInputStream;
import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * A filter for adding the classification, the class distribution and an error flag to a dataset with a classifier. The classifier is either trained on the data itself or provided as serialized model.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns on output of debugging information.</pre>
 * 
 * <pre> -W &lt;classifier specification&gt;
 *  Full class name of classifier to use, followed
 *  by scheme options. eg:
 *   "weka.classifiers.bayes.NaiveBayes -D"
 *  (default: weka.classifiers.rules.ZeroR)</pre>
 * 
 * <pre> -serialized &lt;file&gt;
 *  Instead of training a classifier on the data, one can also provide
 *  a serialized model and use that for tagging the data.</pre>
 * 
 * <pre> -classification
 *  Adds an attribute with the actual classification.
 *  (default: off)</pre>
 * 
 * <pre> -remove-old-class
 *  Removes the old class attribute.
 *  (default: off)</pre>
 * 
 * <pre> -distribution
 *  Adds attributes with the distribution for all classes 
 *  (for numeric classes this will be identical to the attribute 
 *  output with '-classification').
 *  (default: off)</pre>
 * 
 * <pre> -error
 *  Adds an attribute indicating whether the classifier output 
 *  a wrong classification (for numeric classes this is the numeric 
 *  difference).
 *  (default: off)</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class AddClassification
  extends SimpleBatchFilter {

  /** for serialization */
  private static final long serialVersionUID = -1931467132568441909L;

  /** The classifier template used to do the classification */
  protected Classifier m_Classifier = new weka.classifiers.rules.ZeroR();

  /** The file from which to load a serialized classifier */
  protected File m_SerializedClassifierFile = new File(System.getProperty("user.dir"));
  
  /** The actual classifier used to do the classification */
  protected Classifier m_ActualClassifier = null;

  /** whether to output the classification */
  protected boolean m_OutputClassification = false;

  /** whether to remove the old class attribute */
  protected boolean m_RemoveOldClass = false;
  
  /** whether to output the class distribution */
  protected boolean m_OutputDistribution = false;
  
  /** whether to output the error flag */
  protected boolean m_OutputErrorFlag = false;

  /**
   * Returns a string describing this filter
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "A filter for adding the classification, the class distribution and "
      + "an error flag to a dataset with a classifier. The classifier is "
      + "either trained on the data itself or provided as serialized model.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector        	result;
    Enumeration   	en;

    result = new Vector();

    en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());

    result.addElement(new Option(
	"\tFull class name of classifier to use, followed\n"
	+ "\tby scheme options. eg:\n"
	+ "\t\t\"weka.classifiers.bayes.NaiveBayes -D\"\n"
	+ "\t(default: weka.classifiers.rules.ZeroR)",
	"W", 1, "-W <classifier specification>"));

    result.addElement(new Option(
	"\tInstead of training a classifier on the data, one can also provide\n"
	+ "\ta serialized model and use that for tagging the data.",
	"serialized", 1, "-serialized <file>"));

    result.addElement(new Option(
	"\tAdds an attribute with the actual classification.\n"
	+ "\t(default: off)",
	"classification", 0, "-classification"));

    result.addElement(new Option(
	"\tRemoves the old class attribute.\n"
	+ "\t(default: off)",
	"remove-old-class", 0, "-remove-old-class"));

    result.addElement(new Option(
	"\tAdds attributes with the distribution for all classes \n"
        + "\t(for numeric classes this will be identical to the attribute \n"
        + "\toutput with '-classification').\n"
	+ "\t(default: off)",
	"distribution", 0, "-distribution"));

    result.addElement(new Option(
	"\tAdds an attribute indicating whether the classifier output \n"
        + "\ta wrong classification (for numeric classes this is the numeric \n"
        + "\tdifference).\n"
	+ "\t(default: off)",
	"error", 0, "-error"));

    return result.elements();
  }

  /**
   * Parses the options for this object. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Turns on output of debugging information.</pre>
   * 
   * <pre> -W &lt;classifier specification&gt;
   *  Full class name of classifier to use, followed
   *  by scheme options. eg:
   *   "weka.classifiers.bayes.NaiveBayes -D"
   *  (default: weka.classifiers.rules.ZeroR)</pre>
   * 
   * <pre> -serialized &lt;file&gt;
   *  Instead of training a classifier on the data, one can also provide
   *  a serialized model and use that for tagging the data.</pre>
   * 
   * <pre> -classification
   *  Adds an attribute with the actual classification.
   *  (default: off)</pre>
   * 
   * <pre> -remove-old-class
   *  Removes the old class attribute.
   *  (default: off)</pre>
   * 
   * <pre> -distribution
   *  Adds attributes with the distribution for all classes 
   *  (for numeric classes this will be identical to the attribute 
   *  output with '-classification').
   *  (default: off)</pre>
   * 
   * <pre> -error
   *  Adds an attribute indicating whether the classifier output 
   *  a wrong classification (for numeric classes this is the numeric 
   *  difference).
   *  (default: off)</pre>
   * 
   <!-- options-end -->
   *
   * @param options	the options to use
   * @throws Exception	if setting of options fails
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    String[] 	tmpOptions;
    File	file;
    boolean 	serializedModel;

    setOutputClassification(Utils.getFlag("classification", options));
    
    setRemoveOldClass(Utils.getFlag("remove-old-class", options));
    
    setOutputDistribution(Utils.getFlag("distribution", options));

    setOutputErrorFlag(Utils.getFlag("error", options));
    
    serializedModel = false;
    tmpStr = Utils.getOption("serialized", options);
    if (tmpStr.length() != 0) {
      file = new File(tmpStr);
      if (!file.exists())
	throw new FileNotFoundException(
	    "File '" + file.getAbsolutePath() + "' not found!");
      if (file.isDirectory())
	throw new FileNotFoundException(
	    "'" + file.getAbsolutePath() + "' points to a directory not a file!");
      setSerializedClassifierFile(file);
      serializedModel = true;
    }
    else {
      setSerializedClassifierFile(null);
    }
    
    if (!serializedModel) {
      tmpStr = Utils.getOption('W', options);
      if (tmpStr.length() == 0)
	tmpStr = weka.classifiers.rules.ZeroR.class.getName();
      tmpOptions = Utils.splitOptions(tmpStr);
      if (tmpOptions.length == 0)
	throw new Exception("Invalid classifier specification string");
      tmpStr = tmpOptions[0];
      tmpOptions[0] = "";
      setClassifier(AbstractClassifier.forName(tmpStr, tmpOptions));
    }

    super.setOptions(options);
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int		i;
    Vector	result;
    String[]	options;
    File	file;

    result = new Vector();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    if (getOutputClassification())
      result.add("-classification");

    if (getRemoveOldClass())
      result.add("-remove-old-class");

    if (getOutputDistribution())
      result.add("-distribution");

    if (getOutputErrorFlag())
      result.add("-error");

    file = getSerializedClassifierFile();
    if ((file != null) && (!file.isDirectory())) {
      result.add("-serialized");
      result.add(file.getAbsolutePath());
    }
    else {
      result.add("-W");
      result.add(getClassifierSpec());
    }
    
    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities 	result;
    
    if (getClassifier() == null) {
      result = super.getCapabilities();
      result.disableAll();
    } else {
      result = getClassifier().getCapabilities();
    }
    
    result.setMinimumNumberInstances(0);
    
    result.disable(Capability.RANKING);
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    
    return result;
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String classifierTipText() {
    return "The classifier to use for classification.";
  }

  /**
   * Sets the classifier to classify instances with.
   *
   * @param value 	The classifier to be used (with its options set).
   */
  public void setClassifier(Classifier value) {
    m_Classifier = value;
  }
  
  /**
   * Gets the classifier used by the filter.
   *
   * @return 		The classifier to be used.
   */
  public Classifier getClassifier() {
    return m_Classifier;
  }

  /**
   * Gets the classifier specification string, which contains the class name of
   * the classifier and any options to the classifier.
   *
   * @return 		the classifier string.
   */
  protected String getClassifierSpec() {
    String	result;
    Classifier 	c;
    
    c      = getClassifier();
    result = c.getClass().getName();
    if (c instanceof OptionHandler)
      result += " " + Utils.joinOptions(((OptionHandler) c).getOptions());
    
    return result;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String serializedClassifierFileTipText() {
    return "A file containing the serialized model of a trained classifier.";
  }

  /**
   * Gets the file pointing to a serialized, trained classifier. If it is
   * null or pointing to a directory it will not be used.
   * 
   * @return		the file the serialized, trained classifier is located 
   * 			in
   */
  public File getSerializedClassifierFile() {
    return m_SerializedClassifierFile;
  }

  /**
   * Sets the file pointing to a serialized, trained classifier. If the
   * argument is null, doesn't exist or pointing to a directory, then the 
   * value is ignored.
   * 
   * @param value	the file pointing to the serialized, trained classifier
   */
  public void setSerializedClassifierFile(File value) {
    if ((value == null) || (!value.exists()))
      value = new File(System.getProperty("user.dir"));

    m_SerializedClassifierFile = value;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String outputClassificationTipText() {
    return "Whether to add an attribute with the actual classification.";
  }

  /**
   * Get whether the classifiction of the classifier is output.
   *
   * @return 		true if the classification of the classifier is output.
   */
  public boolean getOutputClassification() {
    return m_OutputClassification;
  }
  
  /**
   * Set whether the classification of the classifier is output.
   *
   * @param value 	whether the classification of the classifier is output.
   */
  public void setOutputClassification(boolean value) {
    m_OutputClassification = value;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String removeOldClassTipText() {
    return "Whether to remove the old class attribute.";
  }

  /**
   * Get whether the old class attribute is removed.
   *
   * @return 		true if the old class attribute is removed.
   */
  public boolean getRemoveOldClass() {
    return m_RemoveOldClass;
  }
  
  /**
   * Set whether the old class attribute is removed.
   *
   * @param value 	whether the old class attribute is removed.
   */
  public void setRemoveOldClass(boolean value) {
    m_RemoveOldClass = value;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String outputDistributionTipText() {
    return 
        "Whether to add attributes with the distribution for all classes "
      + "(for numeric classes this will be identical to the attribute output "
      + "with 'outputClassification').";
  }

  /**
   * Get whether the classifiction of the classifier is output.
   *
   * @return 		true if the distribution of the classifier is output.
   */
  public boolean getOutputDistribution() {
    return m_OutputDistribution;
  }
  
  /**
   * Set whether the Distribution of the classifier is output.
   *
   * @param value 	whether the distribution of the classifier is output.
   */
  public void setOutputDistribution(boolean value) {
    m_OutputDistribution = value;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String outputErrorFlagTipText() {
    return 
        "Whether to add an attribute indicating whether the classifier output "
      + "a wrong classification (for numeric classes this is the numeric "
      + "difference).";
  }

  /**
   * Get whether the classifiction of the classifier is output.
   *
   * @return 		true if the classification of the classifier is output.
   */
  public boolean getOutputErrorFlag() {
    return m_OutputErrorFlag;
  }
  
  /**
   * Set whether the classification of the classifier is output.
   *
   * @param value 	whether the classification of the classifier is output.
   */
  public void setOutputErrorFlag(boolean value) {
    m_OutputErrorFlag = value;
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
  protected Instances determineOutputFormat(Instances inputFormat)
      throws Exception {
    
    Instances	result;
    FastVector	atts;
    int		i;
    FastVector	values;
    int		classindex;
    
    classindex = -1;
    
    // copy old attributes
    atts = new FastVector();
    for (i = 0; i < inputFormat.numAttributes(); i++) {
      // remove class?
      if ((i == inputFormat.classIndex()) && (getRemoveOldClass()) )
	continue;
      // record class index
      if (i == inputFormat.classIndex())
	classindex = i;
      atts.addElement(inputFormat.attribute(i).copy());
    }
    
    // add new attributes
    // 1. classification?
    if (getOutputClassification()) {
      // if old class got removed, use this one
      if (classindex == -1)
	classindex = atts.size();
      atts.addElement(inputFormat.classAttribute().copy("classification"));
    }
    
    // 2. distribution?
    if (getOutputDistribution()) {
      if (inputFormat.classAttribute().isNominal()) {
	for (i = 0; i < inputFormat.classAttribute().numValues(); i++) {
	  atts.addElement(new Attribute("distribution_" + inputFormat.classAttribute().value(i)));
	}
      }
      else {
	atts.addElement(new Attribute("distribution"));
      }
    }
    
    // 2. error flag?
    if (getOutputErrorFlag()) {
      if (inputFormat.classAttribute().isNominal()) {
	values = new FastVector();
	values.addElement("no");
	values.addElement("yes");
	atts.addElement(new Attribute("error", values));
      }
      else {
	atts.addElement(new Attribute("error"));
      }
    }
    
    // generate new header
    result = new Instances(inputFormat.relationName(), atts, 0);
    result.setClassIndex(classindex);
    
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
    Instances		result;
    double[]		newValues;
    double[]		oldValues;
    int			i;
    int			start;
    int			n;
    Instance		newInstance;
    Instance		oldInstance;
    Instances		header;
    double[]		distribution;
    File		file;
    ObjectInputStream 	ois;
    
    // load or train classifier
    if (!isFirstBatchDone()) {
      file = getSerializedClassifierFile();
      if (!file.isDirectory()) {
	ois = new ObjectInputStream(new FileInputStream(file));
	m_ActualClassifier = (Classifier) ois.readObject();
	header = null;
	// let's see whether there's an Instances header stored as well
	try {
	  header = (Instances) ois.readObject();
	}
	catch (Exception e) {
	  // ignored
	}
	ois.close();
	// same dataset format?
	if ((header != null) && (!header.equalHeaders(instances)))
	  throw new WekaException(
	      "Training header of classifier and filter dataset don't match:\n"
	      + header.equalHeadersMsg(instances));
      }
      else {
	m_ActualClassifier = AbstractClassifier.makeCopy(m_Classifier);
	m_ActualClassifier.buildClassifier(instances);
      }
    }
    
    result = getOutputFormat();
    
    // traverse all instances
    for (i = 0; i < instances.numInstances(); i++) {
      oldInstance = instances.instance(i);
      oldValues   = oldInstance.toDoubleArray();
      newValues   = new double[result.numAttributes()];
      
      start = oldValues.length;
      if (getRemoveOldClass())
	start--;

      // copy old values
      System.arraycopy(oldValues, 0, newValues, 0, start);
      
      // add new values:
      // 1. classification?
      if (getOutputClassification()) {
	newValues[start] = m_ActualClassifier.classifyInstance(oldInstance);
	start++;
      }
      
      // 2. distribution?
      if (getOutputDistribution()) {
	distribution = m_ActualClassifier.distributionForInstance(oldInstance);
	for (n = 0; n < distribution.length; n++) {
	  newValues[start] = distribution[n];
	  start++;
	}
      }
      
      // 3. error flag?
      if (getOutputErrorFlag()) {
	if (result.classAttribute().isNominal()) {
	  if (oldInstance.classValue() == m_ActualClassifier.classifyInstance(oldInstance))
	    newValues[start] = 0;
	  else
	    newValues[start] = 1;
	}
	else {
	  newValues[start] = m_ActualClassifier.classifyInstance(oldInstance) - oldInstance.classValue();
	}
	start++;
      }
      
      // create new instance
      if (oldInstance instanceof SparseInstance)
	newInstance = new SparseInstance(oldInstance.weight(), newValues);
      else
	newInstance = new DenseInstance(oldInstance.weight(), newValues);

      // copy string/relational values from input to output
      copyValues(newInstance, false, oldInstance.dataset(), getOutputFormat());

      result.add(newInstance);
    }
    
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
    runFilter(new AddClassification(), args);
  }
}
