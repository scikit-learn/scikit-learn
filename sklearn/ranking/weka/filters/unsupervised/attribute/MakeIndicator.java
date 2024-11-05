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
 *    MakeIndicator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
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
import weka.core.RevisionUtils;
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
 * A filter that creates a new dataset with a boolean attribute replacing a nominal attribute.  In the new dataset, a value of 1 is assigned to an instance that exhibits a particular range of attribute values, a 0 to an instance that doesn't. The boolean attribute is coded as numeric by default.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C &lt;col&gt;
 *  Sets the attribute index.</pre>
 * 
 * <pre> -V &lt;index1,index2-index4,...&gt;
 *  Specify the list of values to indicate. First and last are
 *  valid indexes (default last)</pre>
 * 
 * <pre> -N &lt;index&gt;
 *  Set if new boolean attribute nominal.</pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz) 
 * @version $Revision: 5987 $
 */
public class MakeIndicator 
  extends Filter
  implements UnsupervisedFilter, StreamableFilter, OptionHandler {

  /** for serialization */
  static final long serialVersionUID = 766001176862773163L;
  
  /** The attribute's index setting. */
  private SingleIndex m_AttIndex = new SingleIndex("last"); 

  /** The value's index */
  private Range m_ValIndex;
  
  /** Make boolean attribute numeric. */
  private boolean m_Numeric = true;

  /**
   * Constructor
   */
  public MakeIndicator() {

      m_ValIndex = new Range("last");
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
   * @param instanceInfo an Instances object containing the input 
   * instance structure (any instances contained in the object are 
   * ignored - only the structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws UnsupportedAttributeTypeException the selecte attribute is not nominal
   * @throws UnsupportedAttributeTypeException the selecte attribute has fewer than two values.
   */
  public boolean setInputFormat(Instances instanceInfo) 
       throws Exception {

    super.setInputFormat(instanceInfo);
    m_AttIndex.setUpper(instanceInfo.numAttributes() - 1);
    m_ValIndex.setUpper(instanceInfo.attribute(m_AttIndex.
					       getIndex()).numValues() - 1);
    if (!instanceInfo.attribute(m_AttIndex.getIndex()).isNominal() && !instanceInfo.attribute(m_AttIndex.getIndex()).isRanking() ) {
      throw new UnsupportedAttributeTypeException("Chosen attribute not nominal.");
    }
    if (instanceInfo.attribute(m_AttIndex.getIndex()).numValues() < 2) {
      throw new UnsupportedAttributeTypeException("Chosen attribute has less " +
						  "than two values.");
    }
    setOutputFormat();
    return true;
  }

  /**
   * Input an instance for filtering. The instance is processed
   * and made available for output immediately.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been set.
   */
  public boolean input(Instance instance) {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    Instance newInstance = (Instance)instance.copy();
    if (!newInstance.isMissing(m_AttIndex.getIndex())) {
      if (m_ValIndex.isInRange((int)newInstance.value(m_AttIndex.getIndex()))) {
	newInstance.setValue(m_AttIndex.getIndex(), 1);
      } else {
	newInstance.setValue(m_AttIndex.getIndex(), 0);
      }
    }
    push(newInstance);
    return true;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(3);

    newVector.addElement(new Option(
              "\tSets the attribute index.",
              "C", 1, "-C <col>"));

    newVector.addElement(new Option(
              "\tSpecify the list of values to indicate. First and last are\n"+
              "\tvalid indexes (default last)",
              "V", 1, "-V <index1,index2-index4,...>"));
    newVector.addElement(new Option(
              "\tSet if new boolean attribute nominal.",
              "N", 0, "-N <index>"));

    return newVector.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -C &lt;col&gt;
   *  Sets the attribute index.</pre>
   * 
   * <pre> -V &lt;index1,index2-index4,...&gt;
   *  Specify the list of values to indicate. First and last are
   *  valid indexes (default last)</pre>
   * 
   * <pre> -N &lt;index&gt;
   *  Set if new boolean attribute nominal.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    String attIndex = Utils.getOption('C', options);
    if (attIndex.length() != 0) {
      setAttributeIndex(attIndex);
    } else {
      setAttributeIndex("last");
    }

    String valIndex = Utils.getOption('V', options);
    if (valIndex.length() != 0) {
      setValueIndices(valIndex);
    } else {
      setValueIndices("last");
    }

    setNumeric(!Utils.getFlag('N', options));

    if (getInputFormat() != null) {
      setInputFormat(getInputFormat());
    }
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [5];
    int current = 0;

    options[current++] = "-C";
    options[current++] = "" + (getAttributeIndex());
    options[current++] = "-V"; 
    options[current++] = getValueIndices();
    if (!getNumeric()) {
      options[current++] = "-N"; 
    }
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "A filter that creates a new dataset with a boolean attribute "
      + "replacing a nominal attribute.  In the new dataset, a value of 1 is "
      + "assigned to an instance that exhibits a particular range of attribute "
      + "values, a 0 to an instance that doesn't. The boolean attribute is "
      + "coded as numeric by default.";
  }

  /**
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String attributeIndexTipText() {

    return "Sets which attribute should be replaced by the indicator. This "
      + "attribute must be nominal.";
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
   * Get the range containing the indicator values.
   *
   * @return the range containing the indicator values
   */
  public Range getValueRange() {

    return m_ValIndex;
  }

  /**
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String valueIndicesTipText() {

    return "Specify range of nominal values to act on."
      + " This is a comma separated list of attribute indices (numbered from"
      + " 1), with \"first\" and \"last\" valid values. Specify an inclusive"
      + " range with \"-\". E.g: \"first-3,5,6-10,last\".";
  }

  /**
   * Get the indices of the indicator values.
   *
   * @return the indices of the indicator values
   */
  public String getValueIndices() {

    return m_ValIndex.getRanges();
  }

  /**
   * Sets indices of the indicator values.
   *
   * @param range the string representation of the indicator value indices
   * @see Range
   */
  public void setValueIndices(String range) {
    
    m_ValIndex.setRanges(range);
  }

  /**
   * Sets index of the indicator value.
   *
   * @param index the index of the indicator value
   */
  public void setValueIndex(int index) {

    setValueIndices("" +  (index + 1));
  }

  /**
   * Set which attributes are to be deleted (or kept if invert is true)
   *
   * @param indices an array containing indexes of attributes to select.
   * Since the array will typically come from a program, attributes are indexed
   * from 0.
   * @throws InvalidArgumentException if an invalid set of ranges is supplied
   */
  public void setValueIndicesArray(int [] indices) {
    
    setValueIndices(Range.indicesToRangeList(indices));
  }

  /**
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numericTipText() {

    return "Determines whether the output indicator attribute is numeric. If "
      + "this is set to false, the output attribute will be nominal.";
  }

  /**
   * Sets if the new Attribute is to be numeric.
   *
   * @param bool true if new Attribute is to be numeric
   */
  public void setNumeric(boolean bool) {

    m_Numeric = bool;
  }

  /**
   * Check if new attribute is to be numeric.
   *
   * @return true if new attribute is to be numeric
   */
  public boolean getNumeric() {

    return m_Numeric;
  }

  /**
   * Set the output format.
   */
  private void setOutputFormat() {
    
    Instances newData;
    FastVector newAtts, newVals;
      
    // Compute new attributes
    
    newAtts = new FastVector(getInputFormat().numAttributes());
    for (int j = 0; j < getInputFormat().numAttributes(); j++) {
      Attribute att = getInputFormat().attribute(j);
      if (j != m_AttIndex.getIndex()) {

	// We don't have to copy the attribute because the
	// attribute index remains unchanged.
	newAtts.addElement(att);
      } else {
	if (m_Numeric) {
	  newAtts.addElement(new Attribute(att.name()));
	} else {
          String vals;
          int [] sel = m_ValIndex.getSelection();
          if (sel.length == 1) {
            vals = att.value(sel[0]);
          } else {
            vals = m_ValIndex.getRanges().replace(',','_');
          }
	  newVals = new FastVector(2);
	  newVals.addElement("neg_" + vals);
	  newVals.addElement("pos_" + vals);
	  newAtts.addElement(new Attribute(att.name(), newVals));
	}
      }
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
    return RevisionUtils.extract("$Revision: 5987 $");
  }
 
  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: 
   * use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new MakeIndicator(), argv);
  }
}
