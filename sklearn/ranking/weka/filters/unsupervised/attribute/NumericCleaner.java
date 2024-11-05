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
 * NumericCleaner.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.unsupervised.attribute;

import weka.core.Capabilities;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.Range;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.SimpleStreamFilter;

import java.util.Enumeration;
import java.util.Vector;


/**
 <!-- globalinfo-start -->
 * A filter that 'cleanses' the numeric data from values that are too small, too big or very close to a certain value (e.g., 0) and sets these values to a pre-defined default.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns on output of debugging information.</pre>
 * 
 * <pre> -min &lt;double&gt;
 *  The minimum threshold. (default -Double.MAX_VALUE)</pre>
 * 
 * <pre> -min-default &lt;double&gt;
 *  The replacement for values smaller than the minimum threshold.
 *  (default -Double.MAX_VALUE)</pre>
 * 
 * <pre> -max &lt;double&gt;
 *  The maximum threshold. (default Double.MAX_VALUE)</pre>
 * 
 * <pre> -max-default &lt;double&gt;
 *  The replacement for values larger than the maximum threshold.
 *  (default Double.MAX_VALUE)</pre>
 * 
 * <pre> -closeto &lt;double&gt;
 *  The number values are checked for closeness. (default 0)</pre>
 * 
 * <pre> -closeto-default &lt;double&gt;
 *  The replacement for values that are close to '-closeto'.
 *  (default 0)</pre>
 * 
 * <pre> -closeto-tolerance &lt;double&gt;
 *  The tolerance below which numbers are considered being close to 
 *  to each other. (default 1E-6)</pre>
 * 
 * <pre> -decimals &lt;int&gt;
 *  The number of decimals to round to, -1 means no rounding at all.
 *  (default -1)</pre>
 * 
 * <pre> -R &lt;col1,col2,...&gt;
 *  The list of columns to cleanse, e.g., first-last or first-3,5-last.
 *  (default first-last)</pre>
 * 
 * <pre> -V
 *  Inverts the matching sense.</pre>
 * 
 * <pre> -include-class
 *  Whether to include the class in the cleansing.
 *  The class column will always be skipped, if this flag is not
 *  present. (default no)</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class NumericCleaner
  extends SimpleStreamFilter {

  /** for serialization */
  private static final long serialVersionUID = -352890679895066592L;

  /** the minimum threshold */
  protected double m_MinThreshold = -Double.MAX_VALUE;

  /** the minimum default replacement value */
  protected double m_MinDefault = -Double.MAX_VALUE;

  /** the maximum threshold */
  protected double m_MaxThreshold = Double.MAX_VALUE;

  /** the maximum default replacement value */
  protected double m_MaxDefault = Double.MAX_VALUE;

  /** the number the values are checked for closeness to */
  protected double m_CloseTo = 0;

  /** the default replacement value for numbers "close-to" */
  protected double m_CloseToDefault = 0;

  /** the tolerance distance, below which numbers are considered being "close-to" */
  protected double m_CloseToTolerance = 1E-6;

  /** Stores which columns to cleanse */
  protected Range m_Cols = new Range("first-last");

  /** whether to include the class attribute */
  protected boolean m_IncludeClass = false;
  
  /** the number of decimals to round to (-1 means no rounding) */
  protected int m_Decimals = -1;
  
  /**
   * Returns a string describing this filter.
   *
   * @return      a description of the filter suitable for
   *              displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "A filter that 'cleanses' the numeric data from values that are too "
      + "small, too big or very close to a certain value (e.g., 0) and sets "
      + "these values to a pre-defined default.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector        result;
    Enumeration   enm;

    result = new Vector();

    enm = super.listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    result.addElement(new Option(
	"\tThe minimum threshold. (default -Double.MAX_VALUE)",
	"min", 1, "-min <double>"));
    
    result.addElement(new Option(
	"\tThe replacement for values smaller than the minimum threshold.\n"
	+ "\t(default -Double.MAX_VALUE)",
	"min-default", 1, "-min-default <double>"));

    result.addElement(new Option(
	"\tThe maximum threshold. (default Double.MAX_VALUE)",
	"max", 1, "-max <double>"));
    
    result.addElement(new Option(
	"\tThe replacement for values larger than the maximum threshold.\n"
	+ "\t(default Double.MAX_VALUE)",
	"max-default", 1, "-max-default <double>"));

    result.addElement(new Option(
	"\tThe number values are checked for closeness. (default 0)",
	"closeto", 1, "-closeto <double>"));
    
    result.addElement(new Option(
	"\tThe replacement for values that are close to '-closeto'.\n"
	+ "\t(default 0)",
	"closeto-default", 1, "-closeto-default <double>"));
    
    result.addElement(new Option(
	"\tThe tolerance below which numbers are considered being close to \n"
	+ "\tto each other. (default 1E-6)",
	"closeto-tolerance", 1, "-closeto-tolerance <double>"));

    result.addElement(new Option(
	"\tThe number of decimals to round to, -1 means no rounding at all.\n"
	+ "\t(default -1)",
	"decimals", 1, "-decimals <int>"));
    
    result.addElement(new Option(
	"\tThe list of columns to cleanse, e.g., first-last or first-3,5-last.\n"
	+ "\t(default first-last)",
	"R", 1, "-R <col1,col2,...>"));

    result.addElement(new Option(
	"\tInverts the matching sense.",
	"V", 0, "-V"));

    result.addElement(new Option(
	"\tWhether to include the class in the cleansing.\n"
	+ "\tThe class column will always be skipped, if this flag is not\n"
	+ "\tpresent. (default no)",
	"include-class", 0, "-include-class"));

    return result.elements();
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

    result.add("-min"); 
    result.add("" + m_MinThreshold);

    result.add("-min-default"); 
    result.add("" + m_MinDefault);

    result.add("-max"); 
    result.add("" + m_MaxThreshold);

    result.add("-max-default"); 
    result.add("" + m_MaxDefault);

    result.add("-closeto"); 
    result.add("" + m_CloseTo);

    result.add("-closeto-default"); 
    result.add("" + m_CloseToDefault);
    
    result.add("-closeto-tolerance"); 
    result.add("" + m_CloseToTolerance);

    result.add("-R"); 
    result.add("" + m_Cols.getRanges());

    if (m_Cols.getInvert())
      result.add("-V");
    
    if (m_IncludeClass)
      result.add("-include-class"); 

    result.add("-decimals"); 
    result.add("" + getDecimals());

    return (String[]) result.toArray(new String[result.size()]);	  
  }	  

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Turns on output of debugging information.</pre>
   * 
   * <pre> -min &lt;double&gt;
   *  The minimum threshold. (default -Double.MAX_VALUE)</pre>
   * 
   * <pre> -min-default &lt;double&gt;
   *  The replacement for values smaller than the minimum threshold.
   *  (default -Double.MAX_VALUE)</pre>
   * 
   * <pre> -max &lt;double&gt;
   *  The maximum threshold. (default Double.MAX_VALUE)</pre>
   * 
   * <pre> -max-default &lt;double&gt;
   *  The replacement for values larger than the maximum threshold.
   *  (default Double.MAX_VALUE)</pre>
   * 
   * <pre> -closeto &lt;double&gt;
   *  The number values are checked for closeness. (default 0)</pre>
   * 
   * <pre> -closeto-default &lt;double&gt;
   *  The replacement for values that are close to '-closeto'.
   *  (default 0)</pre>
   * 
   * <pre> -closeto-tolerance &lt;double&gt;
   *  The tolerance below which numbers are considered being close to 
   *  to each other. (default 1E-6)</pre>
   * 
   * <pre> -decimals &lt;int&gt;
   *  The number of decimals to round to, -1 means no rounding at all.
   *  (default -1)</pre>
   * 
   * <pre> -R &lt;col1,col2,...&gt;
   *  The list of columns to cleanse, e.g., first-last or first-3,5-last.
   *  (default first-last)</pre>
   * 
   * <pre> -V
   *  Inverts the matching sense.</pre>
   * 
   * <pre> -include-class
   *  Whether to include the class in the cleansing.
   *  The class column will always be skipped, if this flag is not
   *  present. (default no)</pre>
   * 
   <!-- options-end -->
   * 
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported 
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;

    tmpStr = Utils.getOption("min", options);
    if (tmpStr.length() != 0)
      setMinThreshold(Double.parseDouble(tmpStr));
    else
      setMinThreshold(-Double.MAX_VALUE);
    
    tmpStr = Utils.getOption("min-default", options);
    if (tmpStr.length() != 0)
      setMinDefault(Double.parseDouble(tmpStr));
    else
      setMinDefault(-Double.MAX_VALUE);
    
    tmpStr = Utils.getOption("max", options);
    if (tmpStr.length() != 0)
      setMaxThreshold(Double.parseDouble(tmpStr));
    else
      setMaxThreshold(Double.MAX_VALUE);
    
    tmpStr = Utils.getOption("max-default", options);
    if (tmpStr.length() != 0)
      setMaxDefault(Double.parseDouble(tmpStr));
    else
      setMaxDefault(Double.MAX_VALUE);
    
    tmpStr = Utils.getOption("closeto", options);
    if (tmpStr.length() != 0)
      setCloseTo(Double.parseDouble(tmpStr));
    else
      setCloseTo(0);
    
    tmpStr = Utils.getOption("closeto-default", options);
    if (tmpStr.length() != 0)
      setCloseToDefault(Double.parseDouble(tmpStr));
    else
      setCloseToDefault(0);
    
    tmpStr = Utils.getOption("closeto-tolerance", options);
    if (tmpStr.length() != 0)
      setCloseToTolerance(Double.parseDouble(tmpStr));
    else
      setCloseToTolerance(1E-6);
    
    tmpStr = Utils.getOption("R", options);
    if (tmpStr.length() != 0)
      setAttributeIndices(tmpStr);
    else
      setAttributeIndices("first-last");
    
    setInvertSelection(Utils.getFlag("V", options));
    
    setIncludeClass(Utils.getFlag("include-class", options));

    tmpStr = Utils.getOption("decimals", options);
    if (tmpStr.length() != 0)
      setDecimals(Integer.parseInt(tmpStr));
    else
      setDecimals(-1);
    
    super.setOptions(options);
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

    m_Cols.setUpper(inputFormat.numAttributes() - 1);
    
    return new Instances(inputFormat);
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
    Instance		result;
    int			i;
    double		val;
    double		factor;
    
    result = (Instance) instance.copy();
    
    if (m_Decimals > -1)
      factor = StrictMath.pow(10, m_Decimals);
    else
      factor = 1;
    
    for (i = 0; i < result.numAttributes(); i++) {
      // only numeric attributes
      if (!result.attribute(i).isNumeric())
	continue;

      // out of range?
      if (!m_Cols.isInRange(i))
	continue;
      
      // skip class?
      if ( (result.classIndex() == i) && (!m_IncludeClass) )
	continue;
      
      // too small?
      if (result.value(i) < m_MinThreshold) {
	if (getDebug())
	  System.out.println("Too small: " + result.value(i) + " -> " + m_MinDefault);
	result.setValue(i, m_MinDefault);
      }
      // too big?
      else if (result.value(i) > m_MaxThreshold) {
	if (getDebug())
	  System.out.println("Too big: " + result.value(i) + " -> " + m_MaxDefault);
	result.setValue(i, m_MaxDefault);
      }
      // too close?
      else if (    (result.value(i) - m_CloseTo < m_CloseToTolerance) 
	        && (m_CloseTo - result.value(i) < m_CloseToTolerance) 
	        && (result.value(i) != m_CloseTo) ) {
	if (getDebug())
	  System.out.println("Too close: " + result.value(i) + " -> " + m_CloseToDefault);
	result.setValue(i, m_CloseToDefault);
      }
      
      // decimals?
      if (m_Decimals > -1) {
	val = result.value(i);
	val = StrictMath.round(val * factor) / factor;
	result.setValue(i, val);
      }
    }

    return result;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String minThresholdTipText() {
    return "The minimum threshold below values are replaced by a default.";
  }

  /**
   * Get the minimum threshold. 
   *
   * @return 		the minimum threshold.
   */
  public double getMinThreshold() {
    return m_MinThreshold;
  }

  /**
   * Set the minimum threshold. 
   *
   * @param value	the minimum threshold to use.
   */
  public void setMinThreshold(double value) {
    m_MinThreshold = value;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String minDefaultTipText() {
    return "The default value to replace values that are below the minimum threshold.";
  }

  /**
   * Get the minimum default. 
   *
   * @return 		the minimum default.
   */
  public double getMinDefault() {
    return m_MinDefault;
  }

  /**
   * Set the minimum default. 
   *
   * @param value	the minimum default to use.
   */
  public void setMinDefault(double value) {
    m_MinDefault = value;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String maxThresholdTipText() {
    return "The maximum threshold above values are replaced by a default.";
  }

  /**
   * Get the maximum threshold. 
   *
   * @return 		the maximum threshold.
   */
  public double getMaxThreshold() {
    return m_MaxThreshold;
  }

  /**
   * Set the maximum threshold. 
   *
   * @param value	the maximum threshold to use.
   */
  public void setMaxThreshold(double value) {
    m_MaxThreshold = value;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String maxDefaultTipText() {
    return "The default value to replace values that are above the maximum threshold.";
  }

  /**
   * Get the maximum default. 
   *
   * @return 		the maximum default.
   */
  public double getMaxDefault() {
    return m_MaxDefault;
  }

  /**
   * Set the naximum default. 
   *
   * @param value	the maximum default to use.
   */
  public void setMaxDefault(double value) {
    m_MaxDefault = value;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String closeToTipText() {
    return 
        "The number values are checked for whether they are too close to "
      + "and get replaced by a default.";
  }

  /**
   * Get the "close to" number.
   *
   * @return 		the "close to" number.
   */
  public double getCloseTo() {
    return m_CloseTo;
  }

  /**
   * Set the "close to" number.
   *
   * @param value	the number to use for checking closeness.
   */
  public void setCloseTo(double value) {
    m_CloseTo = value;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String closeToDefaultTipText() {
    return "The default value to replace values with that are too close.";
  }

  /**
   * Get the "close to" default.
   *
   * @return 		the "close to" default.
   */
  public double getCloseToDefault() {
    return m_CloseToDefault;
  }

  /**
   * Set the "close to" default. 
   *
   * @param value	the "close to" default to use.
   */
  public void setCloseToDefault(double value) {
    m_CloseToDefault = value;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String closeToToleranceTipText() {
    return "The value below which values are considered close to.";
  }

  /**
   * Get the "close to" Tolerance.
   *
   * @return 		the "close to" Tolerance.
   */
  public double getCloseToTolerance() {
    return m_CloseToTolerance;
  }

  /**
   * Set the "close to" Tolerance. 
   *
   * @param value	the "close to" Tolerance to use.
   */
  public void setCloseToTolerance(double value) {
    m_CloseToTolerance = value;
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String attributeIndicesTipText() {
    return "The selection of columns to use in the cleansing processs, first and last are valid indices.";
  }

  /**
   * Gets the selection of the columns, e.g., first-last or first-3,5-last
   *
   * @return 		the selected indices
   */
  public String getAttributeIndices() {
    return m_Cols.getRanges();
  }

  /**
   * Sets the columns to use, e.g., first-last or first-3,5-last
   *
   * @param value 	the columns to use
   */
  public void setAttributeIndices(String value) {
    m_Cols.setRanges(value);
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String invertSelectionTipText() {
    return "If enabled the selection of the columns is inverted.";
  }

  /**
   * Gets whether the selection of the columns is inverted
   *
   * @return 		true if the selection is inverted
   */
  public boolean getInvertSelection() {
    return m_Cols.getInvert();
  }

  /**
   * Sets whether the selection of the indices is inverted or not
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
  public String includeClassTipText() {
    return "If disabled, the class attribute will be always left out of the cleaning process.";
  }

  /**
   * Gets whether the class is included in the cleaning process or always 
   * skipped.
   *
   * @return 		true if the class can be considered for cleaning.
   */
  public boolean getIncludeClass() {
    return m_IncludeClass;
  }

  /**
   * Sets whether the class can be cleaned, too.
   *
   * @param value	true if the class can be cleansed, too
   */
  public void setIncludeClass(boolean value) {
    m_IncludeClass = value;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String decimalsTipText() {
    return "The number of decimals to round to, -1 means no rounding at all.";
  }

  /**
   * Get the number of decimals to round to. 
   *
   * @return 		the number of decimals.
   */
  public int getDecimals() {
    return m_Decimals;
  }

  /**
   * Set the number of decimals to round to.
   *
   * @param value	the number of decimals.
   */
  public void setDecimals(int value) {
    m_Decimals = value;
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
   * Runs the filter from commandline, use "-h" to see all options.
   * 
   * @param args the commandline options for the filter
   */
  public static void main(String[] args) {
    runFilter(new NumericCleaner(), args);
  }
}
