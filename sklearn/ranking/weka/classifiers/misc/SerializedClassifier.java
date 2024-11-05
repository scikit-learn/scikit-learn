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
 * SerializedClassifier.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.misc;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.SerializationHelper;
import weka.core.Utils;
import weka.core.Capabilities.Capability;

import java.io.File;
import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * A wrapper around a serialized classifier model. This classifier loads a serialized models and uses it to make predictions.<br/>
 * <br/>
 * Warning: since the serialized model doesn't get changed, cross-validation cannot bet used with this classifier.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 * <pre> -model &lt;filename&gt;
 *  The file containing the serialized model.
 *  (required)</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5928 $
 */
public class SerializedClassifier
  extends AbstractClassifier {

  /** for serialization */
  private static final long serialVersionUID = 4599593909947628642L;

  /** the serialized classifier model used for making predictions */
  protected transient Classifier m_Model = null;
  
  /** the file where the serialized model is stored */
  protected File m_ModelFile = new File(System.getProperty("user.dir"));
  
  /**
   * Returns a string describing classifier
   * 
   * @return 		a description suitable for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "A wrapper around a serialized classifier model. This classifier loads "
      + "a serialized models and uses it to make predictions.\n\n"
      + "Warning: since the serialized model doesn't get changed, cross-validation "
      + "cannot bet used with this classifier.";
  }

  /**
   * Gets an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions(){
    Vector        	result;
    Enumeration   	enm;

    result = new Vector();

    enm = super.listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    result.addElement(new Option(
	"\tThe file containing the serialized model.\n"
	+ "\t(required)",
	"model", 1, "-model <filename>"));

    return result.elements();
  }
  
  /**
   * returns the options of the current setup
   *
   * @return		the current options
   */
  public String[] getOptions(){
    int       	i;
    Vector    	result;
    String[]  	options;

    result = new Vector();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-model");
    result.add("" + getModelFile());

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * Parses the options for this object. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  If set, classifier is run in debug mode and
   *  may output additional info to the console</pre>
   * 
   * <pre> -model &lt;filename&gt;
   *  The file containing the serialized model.
   *  (required)</pre>
   * 
   <!-- options-end -->
   *
   * @param options	the options to use
   * @throws Exception	if setting of options fails
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    super.setOptions(options);
    
    tmpStr = Utils.getOption("model", options);
    if (tmpStr.length() != 0)
      setModelFile(new File(tmpStr));
    else
      setModelFile(new File(System.getProperty("user.dir")));
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String modelFileTipText() {
    return "The serialized classifier model to use for predictions.";
  }

  /**
   * Gets the file containing the serialized model.
   *
   * @return 		the file.
   */
  public File getModelFile() {
    return m_ModelFile;
  }
  
  /**
   * Sets the file containing the serialized model.
   *
   * @param value 	the file.
   */
  public void setModelFile(File value) {
    m_ModelFile = value;
    
    if (value.exists() && value.isFile()) {
      try {
	initModel();
      }
      catch (Exception e) {
	throw new IllegalArgumentException("Cannot load model from file '" + value + "': " + e);
      }
    }
  }

  /**
   * Sets the fully built model to use, if one doesn't want to load a model
   * from a file or already deserialized a model from somewhere else.
   * 
   * @param value	the built model
   * @see		#getCurrentModel()
   */
  public void setModel(Classifier value) {
    m_Model = value;
  }
  
  /**
   * Gets the currently loaded model (can be null). Call buildClassifier method
   * to load model from file.
   * 
   * @return		the current model
   * @see		#setModel(Classifier)
   */
  public Classifier getCurrentModel() {
    return m_Model;
  }
  
  /**
   * loads the serialized model if necessary, throws an Exception if the
   * derserialization fails.
   * 
   * @throws Exception	if deserialization fails
   */
  protected void initModel() throws Exception {
    if (m_Model == null)
      m_Model = (Classifier) SerializationHelper.read(m_ModelFile.getAbsolutePath());
  }

  /**
   * Returns default capabilities of the base classifier.
   *
   * @return      the capabilities of the base classifier
   */
  public Capabilities getCapabilities() {
    Capabilities        result;

    // init model if necessary
    try {
      initModel();
    }
    catch (Exception e) {
      System.err.println(e);
    }

    if (m_Model != null) {
      result = m_Model.getCapabilities();
    } else {
      result = new Capabilities(this);
      result.disableAll();
    }
    
    // set dependencies
    for (Capability cap: Capability.values())
      result.enableDependency(cap);
    
    //RANKING BEGIN
    result.enable(Capability.PREFERENCE_ATTRIBUTE);
    result.enable(Capability.RANKING);
    //RANKING END
    
    result.setOwner(this);
    
    return result;
  }
  
  /**
   * Calculates the class membership probabilities for the given test
   * instance.
   *
   * @param instance the instance to be classified
   * @return preedicted class probability distribution
   * @throws Exception if distribution can't be computed successfully 
   */
  public double[] distributionForInstance(Instance instance) throws Exception {
    double[]	result;

    // init model if necessary
    initModel();
    
    result = m_Model.distributionForInstance(instance);
    
    return result;
  }
  
  /**
   * loads only the serialized classifier
   * 
   * @param data        the training instances
   * @throws Exception  if something goes wrong
   */
  public void buildClassifier(Instances data) throws Exception {
    // init model if necessary
    initModel();

    // can classifier handle the data?
    getCapabilities().testWithFail(data);
  }

  /**
   * Returns a string representation of the classifier
   * 
   * @return		the string representation of the classifier
   */
  public String toString() {
    StringBuffer	result;
    
    if (m_Model == null) {
      result = new StringBuffer("No model loaded yet.");
    }
    else {
      result = new StringBuffer();
      result.append("SerializedClassifier\n");
      result.append("====================\n\n");
      result.append("File: " + getModelFile() + "\n\n");
      result.append(m_Model.toString());
    }
    
    return result.toString();
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
   * Runs the classifier with the given options
   * 
   * @param args	the commandline options
   */
  public static void main(String[] args) {
    runClassifier(new SerializedClassifier(), args);
  }
}
