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
 * NumericToNominal.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.unsupervised.attribute;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.Range;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.SimpleBatchFilter;

import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * A filter for turning numeric attributes into nominal ones. Unlike discretization, it just takes all numeric values and adds them to the list of nominal values of that attribute. Useful after CSV imports, to enforce certain attributes to become nominal, e.g., the class attribute, containing values from 1 to 5.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -R &lt;col1,col2-col4,...&gt;
 *  Specifies list of columns to Discretize. First and last are valid indexes.
 *  (default: first-last)</pre>
 * 
 * <pre> -V
 *  Invert matching sense of column indexes.</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class NumericToNominal
  extends SimpleBatchFilter {

  /** for serialization */
  private static final long serialVersionUID = -6614630932899796239L;

  /** the maximum number of decimals to use */
  protected final static int MAX_DECIMALS = 6;
  
  /** Stores which columns to turn into nominals */
  protected Range m_Cols = new Range("first-last");

  /** The default columns to turn into nominals */
  protected String m_DefaultCols = "first-last";

  /**
   * Returns a string describing this filter
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "A filter for turning numeric attributes into nominal ones. Unlike "
      + "discretization, it just takes all numeric values and adds them to "
      + "the list of nominal values of that attribute. Useful after CSV "
      + "imports, to enforce certain attributes to become nominal, e.g., "
      + "the class attribute, containing values from 1 to 5.";
  }

  /**
   * Gets an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    result.addElement(new Option(
	"\tSpecifies list of columns to Discretize. First"
	+ " and last are valid indexes.\n"
	+ "\t(default: first-last)",
	"R", 1, "-R <col1,col2-col4,...>"));

    result.addElement(new Option(
	"\tInvert matching sense of column indexes.",
	"V", 0, "-V"));

    return result.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -R &lt;col1,col2-col4,...&gt;
   *  Specifies list of columns to Discretize. First and last are valid indexes.
   *  (default: first-last)</pre>
   * 
   * <pre> -V
   *  Invert matching sense of column indexes.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;

    super.setOptions(options);
    
    setInvertSelection(Utils.getFlag('V', options));

    tmpStr = Utils.getOption('R', options);
    if (tmpStr.length() != 0)
      setAttributeIndices(tmpStr);
    else
      setAttributeIndices(m_DefaultCols);

    if (getInputFormat() != null)
      setInputFormat(getInputFormat());
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int       i;
    Vector    result;
    String[]  options;

    result = new Vector();
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    if (!getAttributeIndices().equals("")) {
      result.add("-R");
      result.add(getAttributeIndices());
    }

    if (getInvertSelection())
      result.add("-V");

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String invertSelectionTipText() {
    return 
        "Set attribute selection mode. If false, only selected"
      + " (numeric) attributes in the range will be 'nominalized'; if"
      + " true, only non-selected attributes will be 'nominalized'.";
  }

  /**
   * Gets whether the supplied columns are to be worked on or the others.
   *
   * @return 		true if the supplied columns will be worked on
   */
  public boolean getInvertSelection() {
    return m_Cols.getInvert();
  }

  /**
   * Sets whether selected columns should be worked on or all the others apart
   * from these. If true all the other columns are considered for 
   * "nominalization".
   *
   * @param value 	the new invert setting
   */
  public void setInvertSelection(boolean value) {
    m_Cols.setInvert(value);
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String attributeIndicesTipText() {
    return "Specify range of attributes to act on."
      + " This is a comma separated list of attribute indices, with"
      + " \"first\" and \"last\" valid values. Specify an inclusive"
      + " range with \"-\". E.g: \"first-3,5,6-10,last\".";
  }

  /**
   * Gets the current range selection
   *
   * @return 		a string containing a comma separated list of ranges
   */
  public String getAttributeIndices() {
    return m_Cols.getRanges();
  }

  /**
   * Sets which attributes are to be "nominalized" (only numeric
   * attributes among the selection will be transformed).
   *
   * @param value 	a string representing the list of attributes. Since
   * 			the string will typically come from a user, attributes 
   * 			are indexed from 1. <br> eg: first-3,5,6-last
   * @throws IllegalArgumentException if an invalid range list is supplied 
   */
  public void setAttributeIndices(String value) {
    m_Cols.setRanges(value);
  }

  /**
   * Sets which attributes are to be transoformed to nominal. (only numeric
   * attributes among the selection will be transformed).
   *
   * @param value 	an array containing indexes of attributes to nominalize.
   * 			Since the array will typically come from a program, 
   * 			attributes are indexed from 0.
   * @throws IllegalArgumentException if an invalid set of ranges is supplied 
   */
  public void setAttributeIndicesArray(int[] value) {
    setAttributeIndices(Range.indicesToRangeList(value));
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
    
    Instances 	data;
    Instances	result;
    FastVector	atts;
    FastVector	values;
    HashSet	hash;
    int		i;
    int		n;
    boolean	isDate;
    Instance	inst;
    Vector	sorted;

    m_Cols.setUpper(inputFormat.numAttributes() - 1);
    data = new Instances(inputFormat);
    atts = new FastVector();
    for (i = 0; i < data.numAttributes(); i++) {
      if (!m_Cols.isInRange(i) || !data.attribute(i).isNumeric()) {
	atts.addElement(data.attribute(i));
	continue;
      }
      
      // date attribute?
      isDate = (data.attribute(i).type() == Attribute.DATE);
      
      // determine all available attribtues in dataset
      hash   = new HashSet();
      for (n = 0; n < data.numInstances(); n++) {
	inst = data.instance(n);
	if (inst.isMissing(i))
	  continue;
	
	if (isDate)
	  hash.add(inst.stringValue(i));
	else
	  hash.add(new Double(inst.value(i)));
      }
      
      // sort values
      sorted = new Vector();
      for (Object o: hash)
	sorted.add(o);
      Collections.sort(sorted);
      
      // create attribute from sorted values
      values = new FastVector();
      for (Object o: sorted) {
	if (isDate)
	  values.addElement(
	      o.toString());
	else
	  values.addElement(
	      Utils.doubleToString(((Double) o).doubleValue(), MAX_DECIMALS));
      }
      atts.addElement(new Attribute(data.attribute(i).name(), values));
    }
    
    result = new Instances(inputFormat.relationName(), atts, 0);
    result.setClassIndex(inputFormat.classIndex());
    
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
    Instances	result;
    int		i;
    int		n;
    double[]	values;
    String	value;
    Instance	inst;
    Instance	newInst;
    
    // we need the complete input data!
    if (!isFirstBatchDone())
      setOutputFormat(determineOutputFormat(getInputFormat()));
    
    result = new Instances(getOutputFormat());
    
    for (i = 0; i < instances.numInstances(); i++) {
      inst   = instances.instance(i);
      values = inst.toDoubleArray();
      
      for (n = 0; n < values.length; n++) {
	if (    !m_Cols.isInRange(n)
	     || !instances.attribute(n).isNumeric() 
	     || inst.isMissing(n) )
	  continue;

	// get index of value
	if (instances.attribute(n).type() == Attribute.DATE)
	  value = inst.stringValue(n);
	else
	  value = Utils.doubleToString(inst.value(n), MAX_DECIMALS);
	
	values[n] = result.attribute(n).indexOfValue(value);
      }
      
      // generate new instance
      if (inst instanceof SparseInstance)
	newInst = new SparseInstance(inst.weight(), values);
      //RANKING BEGIN
      if (inst instanceof PreferenceDenseInstance){
    	  PreferenceDenseInstance pdi = (PreferenceDenseInstance) inst;
    	  newInst = new PreferenceDenseInstance(inst.weight(), values, pdi.getHashMap());
      }
      //RANKING END
      else
	newInst = new DenseInstance(inst.weight(), values);
      
      // copy possible string, relational values
      newInst.setDataset(getOutputFormat());
      copyValues(newInst, false, inst.dataset(), getOutputFormat());
      
      result.add(newInst);
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
   * Runs the filter with the given parameters. Use -h to list options.
   * 
   * @param args	the commandline options
   */
  public static void main(String[] args) {
    runFilter(new NumericToNominal(), args);
  }
}
