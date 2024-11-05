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
 * AddValues.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
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
import weka.core.UnsupportedAttributeTypeException;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.Filter;
import weka.filters.StreamableFilter;
import weka.filters.UnsupervisedFilter;

import java.util.Collections;
import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Adds the labels from the given list to an attribute if they are missing. The labels can also be sorted in an ascending manner. If no labels are provided then only the (optional) sorting applies.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C &lt;col&gt;
 *  Sets the attribute index
 *  (default last).</pre>
 * 
 * <pre> -L &lt;label1,label2,...&gt;
 *  Comma-separated list of labels to add.
 *  (default: none)</pre>
 * 
 * <pre> -S
 *  Turns on the sorting of the labels.</pre>
 * 
 <!-- options-end -->
 *
 * Based on code from AddValues.
 *
 * @author  FracPete (fracpete at waikato dot ac dot nz) 
 * @version $Revision: 5987 $
 * @see     AddValues
 */
public class AddValues 
  extends Filter
  implements UnsupervisedFilter, StreamableFilter, OptionHandler {

  /** for serialization */
  private static final long serialVersionUID = -8100622241742393656L;

  /** The attribute's index setting. */
  protected SingleIndex m_AttIndex = new SingleIndex("last"); 

  /** The values to add. */
  protected Vector m_Labels = new Vector();

  /** Whether to sort the values. */
  protected boolean m_Sort = false;

  /** the array with the sorted label indices */
  protected int[] m_SortedIndices;
  
  /**
   * Returns a string describing this filter
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Adds the labels from the given list to an attribute if they are "
      + "missing. The labels can also be sorted in an ascending manner. "
      + "If no labels are provided then only the (optional) sorting applies.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector	result;
    
    result = new Vector();

    result.addElement(new Option(
	"\tSets the attribute index\n"
	+ "\t(default last).",
	"C", 1, "-C <col>"));

    result.addElement(new Option(
	"\tComma-separated list of labels to add.\n"
	+ "\t(default: none)",
	"L", 1, "-L <label1,label2,...>"));

    result.addElement(new Option(
	"\tTurns on the sorting of the labels.",
	"S", 0, "-S"));

    return result.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -C &lt;col&gt;
   *  Sets the attribute index
   *  (default last).</pre>
   * 
   * <pre> -L &lt;label1,label2,...&gt;
   *  Comma-separated list of labels to add.
   *  (default: none)</pre>
   * 
   * <pre> -S
   *  Turns on the sorting of the labels.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0)
      setAttributeIndex(tmpStr);
    else
      setAttributeIndex("last");

    tmpStr = Utils.getOption('L', options);
    if (tmpStr.length() != 0)
      setLabels(tmpStr);
    else
      setLabels("");

    setSort(Utils.getFlag('S', options));
   
    if (getInputFormat() != null)
      setInputFormat(getInputFormat());
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector	result;
    
    result = new Vector();

    result.add("-C");
    result.add("" + getAttributeIndex());
    
    result.add("-L");
    result.add("" + getLabels());
    
    if (getSort())
      result.add("-S");

    return (String[]) result.toArray(new String[result.size()]);
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
   * 				immediately
   * @throws Exception 		if the input format can't be set successfully
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {
    Attribute 	att;
    Attribute 	attNew;
    Vector 	allLabels;
    Enumeration enm;
    int 	i;
    FastVector	values;
    FastVector	atts;
    Instances	instNew;

    super.setInputFormat(instanceInfo);
    
    m_AttIndex.setUpper(instanceInfo.numAttributes() - 1);
    att = instanceInfo.attribute(m_AttIndex.getIndex());
    if (!att.isNominal() && !att.isRanking())
      throw new UnsupportedAttributeTypeException("Chosen attribute not nominal.");
    
    // merge labels
    allLabels = new Vector();
    enm = att.enumerateValues();
    while (enm.hasMoreElements())
      allLabels.add(enm.nextElement());
    for (i = 0; i < m_Labels.size(); i++) {
      if (!allLabels.contains(m_Labels.get(i)))
	allLabels.add(m_Labels.get(i));
    }
    
    // generate index array
    if (getSort())
      Collections.sort(allLabels);
    m_SortedIndices = new int[att.numValues()];
    enm             = att.enumerateValues();
    i               = 0;
    while (enm.hasMoreElements()) {
      m_SortedIndices[i] = allLabels.indexOf(enm.nextElement());
      i++;
    }
    
    // generate new header
    values = new FastVector();
    for (i = 0; i < allLabels.size(); i++)
      values.addElement(allLabels.get(i));
    attNew = new Attribute(att.name(), values);

    atts = new FastVector();
    for (i = 0; i < instanceInfo.numAttributes(); i++) {
      if (i == m_AttIndex.getIndex())
	atts.addElement(attNew);
      else
	atts.addElement(instanceInfo.attribute(i));
    }
    
    instNew = new Instances(instanceInfo.relationName(), atts, 0);
    instNew.setClassIndex(instanceInfo.classIndex());
    
    // set new format
    setOutputFormat(instNew);
    
    return true;
  }

  /**
   * Input an instance for filtering. The instance is processed
   * and made available for output immediately.
   *
   * @param instance 	the input instance
   * @return 		true if the filtered instance may now be
   * 			collected with output().
   * @throws IllegalStateException if no input format has been set.
   */
  public boolean input(Instance instance) {
    Instance 	newInstance;
    double[]	values;

    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    
    // generate new Instance
    values = instance.toDoubleArray();
    values[m_AttIndex.getIndex()] = m_SortedIndices[(int) values[m_AttIndex.getIndex()]];
    //RANKING BEGIN
    if(instance instanceof PreferenceDenseInstance){
    	PreferenceDenseInstance pdi = (PreferenceDenseInstance) instance;
    	newInstance = new PreferenceDenseInstance(instance.weight(), values, pdi.getHashMap());
    }
    else
    	newInstance = new DenseInstance(instance.weight(), values);    
    //RANKING END

    // copy string values etc. from input to output
    copyValues(instance, false, instance.dataset(), getOutputFormat());
    
    push(newInstance);
    
    return true;
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String attributeIndexTipText() {
    return "Sets which attribute to process. This "
      + "attribute must be nominal (\"first\" and \"last\" are valid values)";
  }

  /**
   * Get the index of the attribute used.
   *
   * @return 		the index of the attribute
   */
  public String getAttributeIndex() {
    return m_AttIndex.getSingleIndex();
  }

  /**
   * Sets index of the attribute used.
   *
   * @param attIndex 	the index of the attribute
   */
  public void setAttributeIndex(String attIndex) {
    m_AttIndex.setSingleIndex(attIndex);
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String labelsTipText() {
    return "Comma-separated list of lables to add.";
  }

  /**
   * Get the comma-separated list of labels that are added.
   *
   * @return 		the list of labels
   */
  public String getLabels() {
    String	result;
    int		i;

    result = "";
    for (i = 0; i < m_Labels.size(); i++) {
      if (i > 0)
	result += ",";
      result += Utils.quote((String) m_Labels.get(i));
    }
    
    return result;
  }

  /**
   * Sets the comma-separated list of labels.
   *
   * @param value	the list
   */
  public void setLabels(String value) {
    int		i;
    String	label;
    boolean	quoted;
    boolean	add;
    
    m_Labels.clear();
    
    label  = "";
    quoted = false;
    add    = false;
    
    for (i = 0; i < value.length(); i++) {
      // quotes?
      if (value.charAt(i) == '"') {
	quoted = !quoted;
	if (!quoted)
	  add = true;
      }
      // comma
      else if ( (value.charAt(i) == ',') && (!quoted) ) {
	add = true;
      }
      // normal character
      else {
	label += value.charAt(i);
	// last character?
	if (i == value.length() - 1)
	  add = true;
      }
      
      if (add) {
	if (label.length() != 0)
	  m_Labels.add(label);
	label = "";
	add   = false;
      }
    }
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String sortTipText() {
    return "Whether to sort the labels alphabetically.";
  }

  /**
   * Gets whether the labels are sorted or not.
   *
   * @return 		true if the labels are sorted
   */
  public boolean getSort() {
    return m_Sort;
  }

  /**
   * Sets whether the labels are sorted.
   *
   * @param value	if true the labels are sorted
   */
  public void setSort(boolean value) {
    m_Sort = value;
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
   * Main method for testing and running this class.
   *
   * @param args 	should contain arguments to the filter: 
   * 			use -h for help
   */
  public static void main(String[] args) {
    runFilter(new AddValues(), args);
  }
}
