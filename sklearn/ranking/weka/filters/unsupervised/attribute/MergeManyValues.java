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
 *    MergeManyValues.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 * 
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
import weka.core.Range;
import weka.core.SingleIndex;
import weka.core.UnsupportedAttributeTypeException;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.StreamableFilter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Merges many values of a nominal attribute into one value.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C &lt;col&gt;
 *  Sets the attribute index
 *  (default: last)</pre>
 * 
 * <pre> -M &lt;label&gt;
 *  Sets the label of the newly merged classes
 *  (default: 'merged')</pre>
 * 
 * <pre> -R &lt;range&gt;
 *  Sets the merge range. 'first and 'last' are accepted as well.'
 *  E.g.: first-5,7,9,20-last
 *  (default: 1,2)</pre>
 * 
 <!-- options-end -->
 *
 * @author Kathryn Hempstalk (kah18 at cs.waikato.ac.nz) 
 * @version $Revision: 5987 $
 */
public class MergeManyValues 
  extends Filter
  implements UnsupervisedFilter, StreamableFilter, OptionHandler {

  /** for serialization */
  private static final long serialVersionUID = 4649332102154713625L;

  /** The attribute's index setting. */
  protected SingleIndex m_AttIndex = new SingleIndex("last"); 

  /** The first value's index setting. */
  protected String m_Label = "merged";

  /** The merge value's index setting. */
  protected Range m_MergeRange = new Range("1,2");

  /**
   * Returns a string describing this filter
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Merges many values of a nominal attribute into one value.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector newVector = new Vector();

    newVector.addElement(new Option(
	"\tSets the attribute index\n"
	+ "\t(default: last)",
	"C", 1, "-C <col>"));

    newVector.addElement(new Option(
	"\tSets the label of the newly merged classes\n"
	+ "\t(default: 'merged')",
	"M", 1, "-M <label>"));

    newVector.addElement(new Option(
	"\tSets the merge range. 'first and 'last' are accepted as well.'\n"
	+ "\tE.g.: first-5,7,9,20-last\n"
	+ "\t(default: 1,2)",
	"R", 1, "-R <range>"));

    return newVector.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -C &lt;col&gt;
   *  Sets the attribute index
   *  (default: last)</pre>
   * 
   * <pre> -M &lt;label&gt;
   *  Sets the label of the newly merged classes
   *  (default: 'merged')</pre>
   * 
   * <pre> -R &lt;range&gt;
   *  Sets the merge range. 'first and 'last' are accepted as well.'
   *  E.g.: first-5,7,9,20-last
   *  (default: 1,2)</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0) {
      setAttributeIndex(tmpStr);
    } else {
      setAttributeIndex("last");
    }

    tmpStr = Utils.getOption('L', options);
    if (tmpStr.length() != 0) {
      setLabel(tmpStr);
    } else {
      setLabel("merged");
    }

    tmpStr = Utils.getOption('R', options);
    if (tmpStr.length() != 0) {
      setMergeValueRange(tmpStr);
    } else {
      setMergeValueRange("1,2");
    }

    if (getInputFormat() != null) {
      setInputFormat(getInputFormat());
    }
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;
    
    result = new Vector<String>();
    
    result.add("-C");
    result.add(getAttributeIndex());
    
    result.add("-L");
    result.add(getLabel());

    result.add("-R"); 
    result.add(getMergeValueRange());

    return result.toArray(new String[result.size()]);
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();

    // attributes
    result.enableAllAttributes();
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
    result.enable(Capability.NO_CLASS);
    
    //RANKING BEGIN
    result.disable(Capability.RANKING);
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    //RANKING END

    return result;
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo 	an Instances object containing the input 
   * 				instance structure (any instances contained 
   * 				in the object are ignored - only the structure 
   * 				is required).
   * @return 			true if the outputFormat may be collected immediately
   * @throws Exception 		if the input format can't be set 
   * 				successfully
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {
    super.setInputFormat(instanceInfo);
    
    m_AttIndex.setUpper(instanceInfo.numAttributes() - 1);

    m_MergeRange.setUpper(instanceInfo.attribute(m_AttIndex.getIndex()).numValues() - 1);
    if ((instanceInfo.classIndex() > -1) && (instanceInfo.classIndex() == m_AttIndex.getIndex())) {
      throw new Exception("Cannot process class attribute.");
    }
    if (!instanceInfo.attribute(m_AttIndex.getIndex()).isNominal() && !instanceInfo.attribute(m_AttIndex.getIndex()).isRanking()) {
      throw new UnsupportedAttributeTypeException("Chosen attribute not nominal.");
    }
    if (instanceInfo.attribute(m_AttIndex.getIndex()).numValues() < 2) {
      throw new UnsupportedAttributeTypeException("Chosen attribute has less than " +
      "two values.");
    }

    setOutputFormat();
    return true;
  }

  /**
   * Set the output format. Takes the current average class values
   * and m_InputFormat and calls setOutputFormat(Instances) 
   * appropriately.
   */
  private void setOutputFormat() {
    Instances newData;
    FastVector newAtts, newVals;

    // Compute new attributes
    newAtts = new FastVector(getInputFormat().numAttributes());
    for (int j = 0; j < getInputFormat().numAttributes(); j++) {
      Attribute att = getInputFormat().attribute(j);
      if (j != m_AttIndex.getIndex()) {
	newAtts.addElement(att.copy());
      } else {

	// Compute list of attribute values	  
	newVals = new FastVector(att.numValues() - 1);
	for (int i = 0; i < att.numValues(); i++) {
	  boolean inMergeList = false;

	  if(att.value(i).equalsIgnoreCase(m_Label)){
	    //don't want to add this one.
	    inMergeList = true;		
	  }else{
	    inMergeList = m_MergeRange.isInRange(i);
	  }

	  if(!inMergeList){
	    //add it.
	    newVals.addElement(att.value(i));
	  }
	}
	newVals.addElement(m_Label);


	newAtts.addElement(new Attribute(att.name(), newVals));
      }
    }

    // Construct new header
    newData = new Instances(getInputFormat().relationName(), newAtts, 0);
    newData.setClassIndex(getInputFormat().classIndex());
    setOutputFormat(newData);
  }

  /**
   * Input an instance for filtering. The instance is processed
   * and made available for output immediately.
   *
   * @param instance 	the input instance
   * @return 		true if the filtered instance may now be
   * 			collected with output().
   * @throws IllegalStateException	if no input format has been set.
   */
  public boolean input(Instance instance) {
    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }

    Attribute att = getInputFormat().attribute(m_AttIndex.getIndex());
    FastVector newVals = new FastVector(att.numValues() - 1);
    for (int i = 0; i < att.numValues(); i++) {
      boolean inMergeList = false;

      if(att.value(i).equalsIgnoreCase(m_Label)){
	//don't want to add this one.
	inMergeList = true;		
      }else{
	inMergeList = m_MergeRange.isInRange(i);
      }

      if(!inMergeList){
	//add it.
	newVals.addElement(att.value(i));
      }
    }
    newVals.addElement(m_Label);

    Attribute temp = new Attribute(att.name(), newVals);

    Instance newInstance = (Instance)instance.copy();    
    if (!newInstance.isMissing(m_AttIndex.getIndex())) {
      String currValue = newInstance.stringValue(m_AttIndex.getIndex());
      if(temp.indexOfValue(currValue) == -1)
	newInstance.setValue(m_AttIndex.getIndex(), temp.indexOfValue(m_Label));
      else
	newInstance.setValue(m_AttIndex.getIndex(), temp.indexOfValue(currValue));
    }

    push(newInstance);
    return true;
  }

  /**
   * Returns the tip text for this property.
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
   * @return the index of the attribute
   */
  public String getAttributeIndex() {
    return m_AttIndex.getSingleIndex();
  }

  /**
   * Sets index of the attribute used.
   *
   * @param attIndex the index of the attribute
   */
  public void setAttributeIndex(String attIndex) {
    m_AttIndex.setSingleIndex(attIndex);
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String labelTipText() {
    return "The new label for the merged values.";
  }

  /**
   * Get the label for the new merged class.
   *
   * @return the label for the merged class.
   */
  public String getLabel() {
    return m_Label;
  }

  /**
   * Sets label of the merged class.
   *
   * @param alabel the new label.
   */
  public void setLabel(String alabel) {
    m_Label = alabel;
  }

  /**
   * Get the range of the merge values used.
   *
   * @return the range of the merge values
   */
  public String getMergeValueRange() {
    return m_MergeRange.getRanges();
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String mergeValueRangeTipText() {
    return "The range of values to merge.";
  }

  /**
   * Sets range of the merge values used.
   *
   * @param range the range of the merged values
   */
  public void setMergeValueRange(String range) {
    m_MergeRange.setRanges(range);
  }

  /**
   * Returns the revision string.
   *
   * @return The revision string.
   */
  public String getRevision(){
    return "$Revision: 5987 $";
  }

  /**
   * Main method for executing this filter.
   *
   * @param args 	use -h to display all options
   */
  public static void main(String[] args) {
    runFilter(new MergeManyValues(), args);
  }
}

