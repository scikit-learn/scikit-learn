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
 * ClassAssigner.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.unsupervised.attribute;

import weka.core.Capabilities;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.SimpleStreamFilter;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Filter that can set and unset the class index.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns on output of debugging information.</pre>
 * 
 * <pre> -C &lt;num|first|last|0&gt;
 *  The index of the class attribute. Index starts with 1, 'first'
 *  and 'last' are accepted, '0' unsets the class index.
 *  (default: last)</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class ClassAssigner
  extends SimpleStreamFilter {

  /** for serialization. */
  private static final long serialVersionUID = 1775780193887394115L;

  /** use the first attribute as class. */
  public final static int FIRST = 0;
  
  /** use the last attribute as class. */
  public final static int LAST = -2;
  
  /** unset the class attribute. */
  public final static int UNSET = -1;
  
  /** the class index. */
  protected int m_ClassIndex = LAST;
  
  /**
   * Returns a string describing this classifier.
   * 
   * @return 		a description of the classifier suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Filter that can set and unset the class index.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
    Enumeration enm = super.listOptions();
    while (enm.hasMoreElements())
      result.add(enm.nextElement());
      
    result.addElement(new Option(
	"\tThe index of the class attribute. Index starts with 1, 'first'\n"
	+ "\tand 'last' are accepted, '0' unsets the class index.\n"
	+ "\t(default: last)",
	"C", 1, "-C <num|first|last|0>"));

    return result.elements();
  }

  /**
   * Parses a list of options for this object. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Turns on output of debugging information.</pre>
   * 
   * <pre> -C &lt;num|first|last|0&gt;
   *  The index of the class attribute. Index starts with 1, 'first'
   *  and 'last' are accepted, '0' unsets the class index.
   *  (default: last)</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String        tmpStr;

    super.setOptions(options);
    
    tmpStr = Utils.getOption("C", options);
    if (tmpStr.length() != 0)
      setClassIndex(tmpStr);
    else
      setClassIndex("last");
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector        result;
    String[]      options;
    int           i;

    result = new Vector();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-C");
    result.add(getClassIndex());

    return (String[]) result.toArray(new String[result.size()]);
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String classIndexTipText() {
    return 
        "The index of the class attribute, starts with 1, 'first' and 'last' "
      + "are accepted as well, '0' unsets the class index.";
  }

  /**
   * sets the class index.
   * 
   * @param value	the class index
   */
  public void setClassIndex(String value) {
    if (value.equalsIgnoreCase("first")) {
      m_ClassIndex = FIRST;
    }
    else if (value.equalsIgnoreCase("last")) {
      m_ClassIndex = LAST;
    }
    else if (value.equalsIgnoreCase("0")) {
      m_ClassIndex = UNSET;
    }
    else {
      try {
	m_ClassIndex = Integer.parseInt(value) - 1;
      }
      catch (Exception e) {
	System.err.println("Error parsing '" + value + "'!");
      }
    }
  }

  /**
   * returns the class index.
   * 
   * @return		the class index
   */
  public String getClassIndex() {
    if (m_ClassIndex == FIRST)
      return "first";
    else if (m_ClassIndex == LAST)
      return "last";
    else if (m_ClassIndex == UNSET)
      return "0";
    else
      return "" + (m_ClassIndex + 1);
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
    result.enable(Capability.NO_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }

  /**
   * Determines the output format based on the input format and returns 
   * this.
   *
   * @param inputFormat     the input format to base the output format on
   * @return                the output format
   * @throws Exception      in case the class index is invalid
   */
  protected Instances determineOutputFormat(Instances inputFormat)
      throws Exception {
    
    Instances result = new Instances(inputFormat, 0);
    
    if (m_ClassIndex == FIRST)
      result.setClassIndex(0);
    else if (m_ClassIndex == LAST)
      result.setClassIndex(result.numAttributes() - 1);
    else if (m_ClassIndex == UNSET)
      result.setClassIndex(-1);
    else
      result.setClassIndex(m_ClassIndex);
    
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
    return instance;
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
   * Main method for executing this class.
   *
   * @param args should contain arguments for the filter: use -h for help
   */
  public static void main(String[] args) {
    runFilter(new ClassAssigner(), args);
  }
}
