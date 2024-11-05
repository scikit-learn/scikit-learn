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
 * NominalToString.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 *
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
import weka.core.SingleIndex;
import weka.core.Range;
import weka.core.UnsupportedAttributeTypeException;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Converts a nominal attribute (i.e. set number of values) to string (i.e. unspecified number of values).
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C &lt;col&gt;
 *  Sets the range of attributes to convert (default last).</pre>
 * 
 <!-- options-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6246 $
 */
public class NominalToString
  extends Filter 
  implements UnsupervisedFilter, OptionHandler {

  /** for serialization */
  static final long serialVersionUID = 8655492378380068939L;
  
  /** The attribute's index setting. */
  private Range m_AttIndex = new Range("last");

  /**
   * Returns a string describing this filter
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Converts a nominal attribute (that is, a set number of values) to string "
      + "(that is, an unspecified number of values).";
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enableAllAttributes();
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
    result.enable(Capability.NO_CLASS);
    
    return result;
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo 	an Instances object containing the input 
   * 				instance structure (any instances contained 
   * 				in the object are ignored - only the 
   * 				structure is required).
   * @return 			true if the outputFormat may be collected 
   * 				immediately.
   * @throws UnsupportedAttributeTypeException 	if the selected attribute
   * 						a string attribute.
   * @throws Exception 		if the input format can't be set 
   * 				successfully.
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {
    super.setInputFormat(instanceInfo);

    m_AttIndex.setUpper(instanceInfo.numAttributes() - 1);
    
    /*    if (!instanceInfo.attribute(m_AttIndex.getIndex()).isNominal())
      throw new UnsupportedAttributeTypeException("Chosen attribute is not of "
      + "type nominal."); */

    return false;
  }

  /**
   * Input an instance for filtering. The instance is processed
   * and made available for output immediately.
   *
   * @param instance 		the input instance.
   * @return 			true if the filtered instance may now be
   * 				collected with output().
   * @throws IllegalStateException 	if no input structure has been defined.
   */
  public boolean input(Instance instance) {
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }

    if (isOutputFormatDefined()) {
      Instance newInstance = (Instance) instance.copy();
      push(newInstance);
      return true;
    }

    bufferInput(instance);
    return false;
  }


  /**
   * Signifies that this batch of input to the filter is finished. If the 
   * filter requires all instances prior to filtering, output() may now 
   * be called to retrieve the filtered instances.
   *
   * @return 		true if there are instances pending output.
   * @throws IllegalStateException 	if no input structure has been defined.
   */
  public boolean batchFinished() {
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    if (!isOutputFormatDefined()) {
      setOutputFormat();

      // Convert pending input instances
      for(int i = 0; i < getInputFormat().numInstances(); i++)
	push((Instance) getInputFormat().instance(i).copy());
    } 

    flushInput();
    m_NewBatch = true;

    return (numPendingOutput() != 0);
  }


  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    result.addElement(new Option(
	"\tSets the range of attributes to convert (default last).",
	"C", 1, "-C <col>"));

    return result.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -C &lt;col&gt;
   *  Sets the range of attributes to convert (default last).</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0)
      setAttributeIndexes(tmpStr);
    else
      setAttributeIndexes("last");
       
    if (getInputFormat() != null)
      setInputFormat(getInputFormat());
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector    result;

    result = new Vector();

    result.add("-C");
    result.add("" + (getAttributeIndexes()));

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String attributeIndexesTipText() {
    return "Sets a range attributes to process. Any non-nominal "
      + "attributes in the range are left untouched (\"first\" and \"last\" are valid values)";
  }

  /**
   * Get the index of the attribute used.
   *
   * @return 		the index of the attribute
   */
  public String getAttributeIndexes() {
    //    return m_AttIndex.getSingleIndex();
    return m_AttIndex.getRanges();
  }

  /**
   * Sets index of the attribute used.
   *
   * @param attIndex 	the index of the attribute
   */
  public void setAttributeIndexes(String attIndex) {
    //    m_AttIndex.setSingleIndex(attIndex);
    m_AttIndex.setRanges(attIndex);
  }

  /**
   * Set the output format. Takes the current average class values
   * and m_InputFormat and calls setOutputFormat(Instances) 
   * appropriately.
   */
  private void setOutputFormat() {
    Instances 	newData;
    FastVector 	newAtts;
      
    // Compute new attributes
    newAtts = new FastVector(getInputFormat().numAttributes());
    for (int j = 0; j < getInputFormat().numAttributes(); j++) {
      Attribute att = getInputFormat().attribute(j);

      if (!att.isNominal() && !att.isRanking() || !m_AttIndex.isInRange(j))
	newAtts.addElement(att); 
      else
	newAtts.addElement(new Attribute(att.name(), (FastVector) null));
    }
      
    // Construct new header
    newData = new Instances(getInputFormat().relationName(), newAtts, 0);
    newData.setClassIndex(getInputFormat().classIndex());

    setOutputFormat(newData);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6246 $");
  }
  
  /**
   * Main method for testing this class.
   *
   * @param args 	should contain arguments to the filter: use -h for help
   */
  public static void main(String [] args) {
    runFilter(new NominalToString(), args);
  }
}

