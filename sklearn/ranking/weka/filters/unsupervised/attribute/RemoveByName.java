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
 * RemoveByName.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import java.util.Enumeration;
import java.util.Vector;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.filters.SimpleStreamFilter;

/** 
 <!-- globalinfo-start -->
 * Removes attributes based on a regular expression matched against their names.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns on output of debugging information.</pre>
 * 
 * <pre> -E &lt;regular expression&gt;
 *  The regular expression to match the attribute names against.
 *  (default: ^.*id$)</pre>
 * 
 * <pre> -V
 *  Flag for inverting the matching sense. If set, attributes are kept
 *  instead of deleted.
 *  (default: off)</pre>
 * 
 <!-- options-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6076 $
 */
public class RemoveByName
  extends SimpleStreamFilter {

  /** for serialization. */
  private static final long serialVersionUID = -3335106965521265631L;

  /** the default expression. */
  public final static String DEFAULT_EXPRESSION = "^.*id$";
  
  /** the regular expression for selecting the attributes by name. */
  protected String m_Expression = DEFAULT_EXPRESSION;
  
  /** whether to invert the matching sense. */
  protected boolean m_InvertSelection;

  /** the Remove filter used internally for removing the attributes. */
  protected Remove m_Remove;
  
  /**
   * Returns a string describing this classifier.
   *
   * @return      a description of the classifier suitable for
   *              displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Removes attributes based on a regular expression matched against "
      + "their names.";
  }

  /**
   * Gets an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector	result;
    Enumeration	enm;

    result = new Vector();

    enm = super.listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    result.addElement(new Option(
	"\tThe regular expression to match the attribute names against.\n"
	+ "\t(default: " + DEFAULT_EXPRESSION + ")",
	"E", 1, "-E <regular expression>"));

    result.addElement(new Option(
	"\tFlag for inverting the matching sense. If set, attributes are kept\n"
	+ "\tinstead of deleted.\n"
	+ "\t(default: off)",
	"V", 0, "-V"));

    return result.elements();
  }

  /**
   * returns the options of the current setup.
   *
   * @return      the current options
   */
  public String[] getOptions() {
    int       		i;
    Vector<String>	result;
    String[]		options;

    result = new Vector();
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-E");
    result.add("" + getExpression());

    if (getInvertSelection())
      result.add("-V");

    return (String[]) result.toArray(new String[result.size()]);	  
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
   * <pre> -E &lt;regular expression&gt;
   *  The regular expression to match the attribute names against.
   *  (default: ^.*id$)</pre>
   * 
   * <pre> -V
   *  Flag for inverting the matching sense. If set, attributes are kept
   *  instead of deleted.
   *  (default: off)</pre>
   * 
   <!-- options-end -->
   *
   * @param options	the options to use
   * @throws Exception	if the option setting fails
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;

    super.setOptions(options);

    tmpStr = Utils.getOption("E", options);
    if (tmpStr.length() != 0)
      setExpression(tmpStr);
    else
      setExpression(DEFAULT_EXPRESSION);

    setInvertSelection(Utils.getFlag("V", options));
  }

  /**
   * Sets the regular expression to match the attribute names against.
   *
   * @param value 	the regular expression
   */
  public void setExpression(String value) {
    m_Expression = value;
  }

  /**
   * Returns the regular expression in use.
   *
   * @return 		the regular expression
   */
  public String getExpression() {
    return m_Expression;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String expressionTipText() {
    return "The regular expression to match the attribute names against.";
  }

  /**
   * Set whether selected columns should be removed or kept. If true the 
   * selected columns are kept and unselected columns are deleted. If false
   * selected columns are deleted and unselected columns are kept.
   *
   * @param value 	the new invert setting
   */
  public void setInvertSelection(boolean value) {
    m_InvertSelection = value;
  }

  /**
   * Get whether the supplied columns are to be removed or kept.
   *
   * @return 		true if the supplied columns will be kept
   */
  public boolean getInvertSelection() {
    return m_InvertSelection;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String invertSelectionTipText() {
    return "Determines whether action is to select or delete."
      + " If set to true, only the specified attributes will be kept;"
      + " If set to false, specified attributes will be deleted.";
  }
  
  /**
   * Determines the output format based on the input format and returns 
   * this. In case the output format cannot be returned immediately, i.e.,
   * immediateOutputFormat() returns false, then this method will be called
   * from batchFinished().
   *
   * @param inputFormat	the input format to base the output format on
   * @return		the output format
   * @throws Exception	in case the determination goes wrong
   */
  protected Instances determineOutputFormat(Instances inputFormat) throws Exception {
    Vector<Integer>	indices;
    int[]		attributes;
    int			i;
    
    // determine indices
    indices = new Vector<Integer>();
    for (i = 0; i < inputFormat.numAttributes(); i++) {
      // skip class
      if (i == inputFormat.classIndex())
	continue;
      if (inputFormat.attribute(i).name().matches(m_Expression))
	indices.add(i);
    }
    attributes = new int[indices.size()];
    for (i = 0; i < indices.size(); i++)
      attributes[i] = indices.get(i);
    
    m_Remove = new Remove();
    m_Remove.setAttributeIndicesArray(attributes);
    m_Remove.setInvertSelection(getInvertSelection());
    m_Remove.setInputFormat(inputFormat);
    
    return m_Remove.getOutputFormat();
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result;
    
    result = new Remove().getCapabilities();
    result.setOwner(this);
    
    return result;
  }
  
  /**
   * processes the given instance (may change the provided instance) and
   * returns the modified version.
   *
   * @param instance    the instance to process
   * @return            the modified data
   * @throws Exception  in case the processing goes wrong
   */
  protected Instance process(Instance instance) throws Exception {
    m_Remove.input(instance);
    return m_Remove.output();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6076 $");
  }

  /**
   * runs the filter with the given arguments.
   *
   * @param args      the commandline arguments
   */
  public static void main(String[] args) {
    runFilter(new RemoveByName(), args);
  }
}
