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
 * AddID.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.unsupervised.attribute;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SingleIndex;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * An instance filter that adds an ID attribute to the dataset. The new attribute contains a unique ID for each instance.<br/>
 * Note: The ID is not reset for the second batch of files (using -b and -r and -s).
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C &lt;index&gt;
 *  Specify where to insert the ID. First and last
 *  are valid indexes.(default first)</pre>
 * 
 * <pre> -N &lt;name&gt;
 *  Name of the new attribute.
 *  (default = 'ID')</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class AddID
  extends Filter
  implements UnsupervisedFilter, OptionHandler {
  
  /** for serialization */
  static final long serialVersionUID = 4734383199819293390L;

  /** the index of the attribute */
  protected SingleIndex m_Index = new SingleIndex("first");

  /** the name of the attribute */
  protected String m_Name = "ID";
  
  /** the counter for the ID */
  protected int m_Counter = -1;
  
  /**
   * Returns a string describing this filter
   *
   * @return            a description of the filter suitable for
   *                    displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "An instance filter that adds an ID attribute to the dataset. "
      + "The new attribute contains a unique ID for each instance.\n"
      + "Note: The ID is not reset for the second batch of files (using -b "
      + "and -r and -s).";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    result.addElement(new Option(
              "\tSpecify where to insert the ID. First and last\n"
              +"\tare valid indexes.(default first)",
              "C", 1, "-C <index>"));

    result.addElement(new Option(
              "\tName of the new attribute.\n"
              +"\t(default = 'ID')",
              "N", 1,"-N <name>"));

    return result.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -C &lt;index&gt;
   *  Specify where to insert the ID. First and last
   *  are valid indexes.(default first)</pre>
   * 
   * <pre> -N &lt;name&gt;
   *  Name of the new attribute.
   *  (default = 'ID')</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      tmpStr;

    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0)
      m_Index.setSingleIndex(tmpStr);
    else
      m_Index.setSingleIndex("first");
    
    tmpStr = Utils.getOption('N', options);
    if (tmpStr.length() != 0)
      m_Name = tmpStr;
    else
      m_Name = "ID";

    if (getInputFormat() != null)
      setInputFormat(getInputFormat());
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector      result;
    
    result = new Vector();

    result.add("-C");
    result.add(getIDIndex());

    result.add("-N");
    result.add(getAttributeName());
    
    return (String[]) result.toArray(new String[result.size()]);
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String attributeNameTipText() {
    return "Set the new attribute's name.";
  }

  /**
   * Get the name of the attribute to be created
   *
   * @return the current attribute name
   */
  public String getAttributeName() {
    return m_Name;
  }

  /** 
   * Set the new attribute's name
   *
   * @param value the new name
   */
  public void setAttributeName(String value) {
    m_Name = value;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String IDIndexTipText() {
    return 
        "The position (starting from 1) where the attribute will be inserted "
      + "(first and last are valid indices).";
  }

  /**
   * Get the index of the attribute used.
   *
   * @return the index of the attribute
   */
  public String getIDIndex() {
    return m_Index.getSingleIndex();
  }

  /**
   * Sets index of the attribute used.
   *
   * @param value the index of the attribute
   */
  public void setIDIndex(String value) {
    m_Index.setSingleIndex(value);
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
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if the format couldn't be set successfully
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {
    Instances           outputFormat;
    Attribute           newAttribute;

    super.setInputFormat(instanceInfo);

    m_Counter = -1;
    m_Index.setUpper(instanceInfo.numAttributes());
    outputFormat = new Instances(instanceInfo, 0);
    newAttribute = new Attribute(m_Name);

    if ((m_Index.getIndex() < 0) || 
        (m_Index.getIndex() > getInputFormat().numAttributes()))
      throw new IllegalArgumentException("Index out of range");
    
    outputFormat.insertAttributeAt(newAttribute, m_Index.getIndex());
    setOutputFormat(outputFormat);
    
    return true;
  }

  /**
   * Input an instance for filtering. Filter requires all
   * training instances be read before producing output.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been set.
   */
  public boolean input(Instance instance) {
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    
    if (!isFirstBatchDone()) {
      bufferInput(instance);
      return false;
    } 
    else {
      convertInstance(instance);
      return true;
    }
  }

  /**
   * Signify that this batch of input to the filter is finished. 
   * If the filter requires all instances prior to filtering,
   * output() may now be called to retrieve the filtered instances.
   *
   * @return true if there are instances pending output
   * @throws IllegalStateException if no input structure has been defined
   */
  public boolean batchFinished() {
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    if (!isFirstBatchDone()) {
      m_Counter = 0;
      
      // Convert pending input instances
      for (int i = 0; i < getInputFormat().numInstances(); i++)
        convertInstance(getInputFormat().instance(i));
    } 
    
    // Free memory
    flushInput();

    m_NewBatch = true;
    m_FirstBatchDone = true;
    
    return (numPendingOutput() != 0);
  }

  /**
   * Convert a single instance over. The converted instance is 
   * added to the end of the output queue.
   *
   * @param instance the instance to convert
   */
  protected void convertInstance(Instance instance) {
    Instance            inst;
    
    m_Counter++;

    // build instance
    try {
      inst = (Instance)instance.copy();

      // First copy string values from input to output
      copyValues(inst, true, inst.dataset(), getOutputFormat());

      // Insert the new attribute and reassign to output
      inst.setDataset(null);
      inst.insertAttributeAt(m_Index.getIndex());
      inst.setValue(m_Index.getIndex(), m_Counter);
      inst.setDataset(getOutputFormat());

      push(inst);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
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
   * Main method for testing this class.
   *
   * @param args should contain arguments to the filter: use -h for help
   */
  public static void main(String[] args) {
    runFilter(new AddID(), args);
  }
}
