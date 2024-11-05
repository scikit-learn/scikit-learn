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
 * InterquartileRange.java
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
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.SimpleBatchFilter;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * A filter for detecting outliers and extreme values based on interquartile ranges. The filter skips the class attribute.<br/>
 * <br/>
 * Outliers:<br/>
 *   Q3 + OF*IQR &lt; x &lt;= Q3 + EVF*IQR<br/>
 *   or<br/>
 *   Q1 - EVF*IQR &lt;= x &lt; Q1 - OF*IQR<br/>
 * <br/>
 * Extreme values:<br/>
 *   x &gt; Q3 + EVF*IQR<br/>
 *   or<br/>
 *   x &lt; Q1 - EVF*IQR<br/>
 * <br/>
 * Key:<br/>
 *   Q1  = 25% quartile<br/>
 *   Q3  = 75% quartile<br/>
 *   IQR = Interquartile Range, difference between Q1 and Q3<br/>
 *   OF  = Outlier Factor<br/>
 *   EVF = Extreme Value Factor
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns on output of debugging information.</pre>
 * 
 * <pre> -R &lt;col1,col2-col4,...&gt;
 *  Specifies list of columns to base outlier/extreme value detection
 *  on. If an instance is considered in at least one of those
 *  attributes an outlier/extreme value, it is tagged accordingly.
 *  'first' and 'last' are valid indexes.
 *  (default none)</pre>
 * 
 * <pre> -O &lt;num&gt;
 *  The factor for outlier detection.
 *  (default: 3)</pre>
 * 
 * <pre> -E &lt;num&gt;
 *  The factor for extreme values detection.
 *  (default: 2*Outlier Factor)</pre>
 * 
 * <pre> -E-as-O
 *  Tags extreme values also as outliers.
 *  (default: off)</pre>
 * 
 * <pre> -P
 *  Generates Outlier/ExtremeValue pair for each numeric attribute in
 *  the range, not just a single indicator pair for all the attributes.
 *  (default: off)</pre>
 * 
 * <pre> -M
 *  Generates an additional attribute 'Offset' per Outlier/ExtremeValue
 *  pair that contains the multiplier that the value is off the median.
 *     value = median + 'multiplier' * IQR
 * Note: implicitely sets '-P'. (default: off)</pre>
 * 
 <!-- options-end -->
 * 
 * Thanks to Dale for a few brainstorming sessions.
 *
 * @author  Dale Fletcher (dale at cs dot waikato dot ac dot nz)
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class InterquartileRange
  extends SimpleBatchFilter {

  /** for serialization */
  private static final long serialVersionUID = -227879653639723030L;

  /** indicator for non-numeric attributes */
  public final static int NON_NUMERIC = -1;
  
  /** the attribute range to work on */
  protected Range m_Attributes = new Range("first-last");
  
  /** the generated indices (only for performance reasons) */
  protected int[] m_AttributeIndices = null;

  /** the factor for detecting outliers */
  protected double m_OutlierFactor = 3;
  
  /** the factor for detecting extreme values, by default 2*m_OutlierFactor */
  protected double m_ExtremeValuesFactor = 2*m_OutlierFactor;
  
  /** whether extreme values are also tagged as outliers */
  protected boolean m_ExtremeValuesAsOutliers = false;

  /** the upper extreme value threshold (= Q3 + EVF*IQR) */
  protected double[] m_UpperExtremeValue = null;

  /** the upper outlier threshold (= Q3 + OF*IQR) */
  protected double[] m_UpperOutlier = null;

  /** the lower outlier threshold (= Q1 - OF*IQR) */
  protected double[] m_LowerOutlier = null;

  /** the interquartile range  */
  protected double[] m_IQR = null;

  /** the median  */
  protected double[] m_Median = null;

  /** the lower extreme value threshold (= Q1 - EVF*IQR) */
  protected double[] m_LowerExtremeValue = null;
  
  /** whether to generate Outlier/ExtremeValue attributes for each attribute
   * instead of a general one */
  protected boolean m_DetectionPerAttribute = false;

  /** the position of the outlier attribute */
  protected int[] m_OutlierAttributePosition = null;

  /** whether to add another attribute called "Offset", that lists the 
   * 'multiplier' by which the outlier/extreme value is away from the median,
   * i.e., value = median + 'multiplier' * IQR <br/>
   * automatically enables m_DetectionPerAttribute!
   */
  protected boolean m_OutputOffsetMultiplier = false;
  
  /**
   * Returns a string describing this filter
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "A filter for detecting outliers and extreme values based on "
      + "interquartile ranges. The filter skips the class attribute.\n\n"
      + "Outliers:\n"
      + "  Q3 + OF*IQR < x <= Q3 + EVF*IQR\n"
      + "  or\n"
      + "  Q1 - EVF*IQR <= x < Q1 - OF*IQR\n"
      + "\n"
      + "Extreme values:\n"
      + "  x > Q3 + EVF*IQR\n"
      + "  or\n"
      + "  x < Q1 - EVF*IQR\n"
      + "\n"
      + "Key:\n"
      + "  Q1  = 25% quartile\n"
      + "  Q3  = 75% quartile\n"
      + "  IQR = Interquartile Range, difference between Q1 and Q3\n"
      + "  OF  = Outlier Factor\n"
      + "  EVF = Extreme Value Factor";
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
	"\tSpecifies list of columns to base outlier/extreme value detection\n"
	+ "\ton. If an instance is considered in at least one of those\n"
	+ "\tattributes an outlier/extreme value, it is tagged accordingly.\n"
	+ " 'first' and 'last' are valid indexes.\n"
	+ "\t(default none)",
	"R", 1, "-R <col1,col2-col4,...>"));

    result.addElement(new Option(
        "\tThe factor for outlier detection.\n"
	+ "\t(default: 3)",
        "O", 1, "-O <num>"));

    result.addElement(new Option(
        "\tThe factor for extreme values detection.\n"
	+ "\t(default: 2*Outlier Factor)",
        "E", 1, "-E <num>"));

    result.addElement(new Option(
        "\tTags extreme values also as outliers.\n"
	+ "\t(default: off)",
        "E-as-O", 0, "-E-as-O"));

    result.addElement(new Option(
        "\tGenerates Outlier/ExtremeValue pair for each numeric attribute in\n"
	+ "\tthe range, not just a single indicator pair for all the attributes.\n"
	+ "\t(default: off)",
        "P", 0, "-P"));

    result.addElement(new Option(
        "\tGenerates an additional attribute 'Offset' per Outlier/ExtremeValue\n"
	+ "\tpair that contains the multiplier that the value is off the median.\n"
	+ "\t   value = median + 'multiplier' * IQR\n"
	+ "Note: implicitely sets '-P'."
	+ "\t(default: off)",
        "M", 0, "-M"));

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
   * <pre> -R &lt;col1,col2-col4,...&gt;
   *  Specifies list of columns to base outlier/extreme value detection
   *  on. If an instance is considered in at least one of those
   *  attributes an outlier/extreme value, it is tagged accordingly.
   *  'first' and 'last' are valid indexes.
   *  (default none)</pre>
   * 
   * <pre> -O &lt;num&gt;
   *  The factor for outlier detection.
   *  (default: 3)</pre>
   * 
   * <pre> -E &lt;num&gt;
   *  The factor for extreme values detection.
   *  (default: 2*Outlier Factor)</pre>
   * 
   * <pre> -E-as-O
   *  Tags extreme values also as outliers.
   *  (default: off)</pre>
   * 
   * <pre> -P
   *  Generates Outlier/ExtremeValue pair for each numeric attribute in
   *  the range, not just a single indicator pair for all the attributes.
   *  (default: off)</pre>
   * 
   * <pre> -M
   *  Generates an additional attribute 'Offset' per Outlier/ExtremeValue
   *  pair that contains the multiplier that the value is off the median.
   *     value = median + 'multiplier' * IQR
   * Note: implicitely sets '-P'. (default: off)</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String        tmpStr;

    super.setOptions(options);

    tmpStr = Utils.getOption("R", options);
    if (tmpStr.length() != 0)
      setAttributeIndices(tmpStr);
    else
      setAttributeIndices("first-last");

    tmpStr = Utils.getOption("O", options);
    if (tmpStr.length() != 0)
      setOutlierFactor(Double.parseDouble(tmpStr));
    else
      setOutlierFactor(3);

    tmpStr = Utils.getOption("E", options);
    if (tmpStr.length() != 0)
      setExtremeValuesFactor(Double.parseDouble(tmpStr));
    else
      setExtremeValuesFactor(2*getOutlierFactor());
    
    setExtremeValuesAsOutliers(Utils.getFlag("E-as-O", options));
    
    setDetectionPerAttribute(Utils.getFlag("P", options));

    setOutputOffsetMultiplier(Utils.getFlag("M", options));
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

    result.add("-R");
    if (!getAttributeIndices().equals(""))
      result.add(getAttributeIndices());
    else
      result.add("first-last");
    
    result.add("-O");
    result.add("" + getOutlierFactor());

    result.add("-E");
    result.add("" + getExtremeValuesFactor());

    if (getExtremeValuesAsOutliers())
      result.add("-E-as-O");
    
    if (getDetectionPerAttribute())
      result.add("-P");
    
    if (getOutputOffsetMultiplier())
      result.add("-M");
    
    return (String[]) result.toArray(new String[result.size()]);
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String attributeIndicesTipText() {
    return 
        "Specify range of attributes to act on; "
      + " this is a comma separated list of attribute indices, with"
      + " \"first\" and \"last\" valid values; specify an inclusive"
      + " range with \"-\", eg: \"first-3,5,6-10,last\".";
  }

  /**
   * Gets the current range selection
   *
   * @return 		a string containing a comma separated list of ranges
   */
  public String getAttributeIndices() {
    return m_Attributes.getRanges();
  }

  /**
   * Sets which attributes are to be used for interquartile calculations and
   * outlier/extreme value detection (only numeric attributes among the 
   * selection will be used).
   *
   * @param value 	a string representing the list of attributes. Since
   * 			the string will typically come from a user, attributes 
   * 			are indexed from 1. <br> eg: first-3,5,6-last
   * @throws IllegalArgumentException if an invalid range list is supplied 
   */
  public void setAttributeIndices(String value) {
    m_Attributes.setRanges(value);
  }

  /**
   * Sets which attributes are to be used for interquartile calculations and
   * outlier/extreme value detection (only numeric attributes among the 
   * selection will be used).
   *
   * @param value 	an array containing indexes of attributes to work on.
   * 			Since the array will typically come from a program, 
   * 			attributes are indexed from 0.
   * @throws IllegalArgumentException if an invalid set of ranges is supplied 
   */
  public void setAttributeIndicesArray(int[] value) {
    setAttributeIndices(Range.indicesToRangeList(value));
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String outlierFactorTipText() {
    return "The factor for determining the thresholds for outliers.";
  }

  /**
   * Sets the factor for determining the thresholds for outliers.
   *
   * @param value 	the factor.
   */
  public void setOutlierFactor(double value) {
    if (value >= getExtremeValuesFactor())
      System.err.println("OutlierFactor must be smaller than ExtremeValueFactor");
    else
      m_OutlierFactor = value;
  }

  /**
   * Gets the factor for determining the thresholds for outliers.
   *
   * @return 		the factor.
   */
  public double getOutlierFactor() {
    return m_OutlierFactor;
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String extremeValuesFactorTipText() {
    return "The factor for determining the thresholds for extreme values.";
  }

  /**
   * Sets the factor for determining the thresholds for extreme values.
   *
   * @param value 	the factor.
   */
  public void setExtremeValuesFactor(double value) {
    if (value <= getOutlierFactor())
      System.err.println("ExtremeValuesFactor must be greater than OutlierFactor!");
    else
      m_ExtremeValuesFactor = value;
  }

  /**
   * Gets the factor for determining the thresholds for extreme values.
   *
   * @return 		the factor.
   */
  public double getExtremeValuesFactor() {
    return m_ExtremeValuesFactor;
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String extremeValuesAsOutliersTipText() {
    return "Whether to tag extreme values also as outliers.";
  }

  /**
   * Set whether extreme values are also tagged as outliers.
   *
   * @param value 	whether or not to tag extreme values also as outliers.
   */
  public void setExtremeValuesAsOutliers(boolean value) {
    m_ExtremeValuesAsOutliers = value;
  }

  /**
   * Get whether extreme values are also tagged as outliers.
   *
   * @return 		true if extreme values are also tagged as outliers.
   */
  public boolean getExtremeValuesAsOutliers() {
    return m_ExtremeValuesAsOutliers;
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String detectionPerAttributeTipText() {
    return 
        "Generates Outlier/ExtremeValue attribute pair for each numeric "
      + "attribute, not just a single pair for all numeric attributes together.";
  }

  /**
   * Set whether an Outlier/ExtremeValue attribute pair is generated for 
   * each numeric attribute ("true") or just one pair for all numeric 
   * attributes together ("false").
   *
   * @param value 	whether or not to generate indicator attribute pairs 
   * 			for each numeric attribute.
   */
  public void setDetectionPerAttribute(boolean value) {
    m_DetectionPerAttribute = value;
    if (!m_DetectionPerAttribute)
      m_OutputOffsetMultiplier = false;
  }

  /**
   * Gets whether an Outlier/ExtremeValue attribute pair is generated for 
   * each numeric attribute ("true") or just one pair for all numeric 
   * attributes together ("false").
   *
   * @return 		true if indicator attribute pairs are generated for
   * 			each numeric attribute.
   */
  public boolean getDetectionPerAttribute() {
    return m_DetectionPerAttribute;
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String outputOffsetMultiplierTipText() {
    return 
        "Generates an additional attribute 'Offset' that contains the "
      + "multiplier the value is off the median: "
      + "value = median + 'multiplier' * IQR";
  }

  /**
   * Set whether an additional attribute "Offset" is generated per 
   * Outlier/ExtremeValue attribute pair that lists the multiplier the value
   * is off the median: value = median + 'multiplier' * IQR.
   *
   * @param value 	whether or not to generate the additional attribute.
   */
  public void setOutputOffsetMultiplier(boolean value) {
    m_OutputOffsetMultiplier = value;
    if (m_OutputOffsetMultiplier)
      m_DetectionPerAttribute = true;
  }

  /**
   * Gets whether an additional attribute "Offset" is generated per 
   * Outlier/ExtremeValue attribute pair that lists the multiplier the value
   * is off the median: value = median + 'multiplier' * IQR.
   *
   * @return 		true if the additional attribute is generated.
   */
  public boolean getOutputOffsetMultiplier() {
    return m_OutputOffsetMultiplier;
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
   * hasImmediateOutputFormat() returns false, then this method will called
   * from batchFinished() after the call of preprocess(Instances), in which,
   * e.g., statistics for the actual processing step can be gathered.
   *
   * @param inputFormat     the input format to base the output format on
   * @return                the output format
   * @throws Exception      in case the determination goes wrong
   * @see                   #hasImmediateOutputFormat()
   * @see                   #batchFinished()
   */
  protected Instances determineOutputFormat(Instances inputFormat)
      throws Exception {
    
    FastVector		atts;
    FastVector		values;
    Instances		result;
    int			i;

    // attributes must be numeric
    m_Attributes.setUpper(inputFormat.numAttributes() - 1);
    m_AttributeIndices = m_Attributes.getSelection();
    for (i = 0; i < m_AttributeIndices.length; i++) {
      // ignore class
      if (m_AttributeIndices[i] == inputFormat.classIndex()) {
	m_AttributeIndices[i] = NON_NUMERIC;
	continue;
      }
      // not numeric -> ignore it
      if (!inputFormat.attribute(m_AttributeIndices[i]).isNumeric())
	m_AttributeIndices[i] = NON_NUMERIC;
    }
    
    // get old attributes
    atts = new FastVector();
    for (i = 0; i < inputFormat.numAttributes(); i++)
      atts.addElement(inputFormat.attribute(i));
    
    if (!getDetectionPerAttribute()) {
      m_OutlierAttributePosition    = new int[1];
      m_OutlierAttributePosition[0] = atts.size();
      
      // add 2 new attributes
      values = new FastVector();
      values.addElement("no");
      values.addElement("yes");
      atts.addElement(new Attribute("Outlier", values));
      
      values = new FastVector();
      values.addElement("no");
      values.addElement("yes");
      atts.addElement(new Attribute("ExtremeValue", values));
    }
    else {
      m_OutlierAttributePosition = new int[m_AttributeIndices.length];
      
      for (i = 0; i < m_AttributeIndices.length; i++) {
	if (m_AttributeIndices[i] == NON_NUMERIC)
	  continue;
	
	m_OutlierAttributePosition[i] = atts.size();

	// add new attributes
	values = new FastVector();
	values.addElement("no");
	values.addElement("yes");
	atts.addElement(
	    new Attribute(
		inputFormat.attribute(
		    m_AttributeIndices[i]).name() + "_Outlier", values));
	
	values = new FastVector();
	values.addElement("no");
	values.addElement("yes");
	atts.addElement(
	    new Attribute(
		inputFormat.attribute(
		    m_AttributeIndices[i]).name() + "_ExtremeValue", values));

	if (getOutputOffsetMultiplier())
	  atts.addElement(
	      new Attribute(
		  inputFormat.attribute(
		      m_AttributeIndices[i]).name() + "_Offset"));
      }
    }

    // generate header
    result = new Instances(inputFormat.relationName(), atts, 0);
    result.setClassIndex(inputFormat.classIndex());
    
    return result;
  }

  /**
   * computes the thresholds for outliers and extreme values
   * 
   * @param instances	the data to work on
   */
  protected void computeThresholds(Instances instances) {
    int		i;
    double[]	values;
    int[]	sortedIndices;
    int		half;
    int		quarter;
    double	q1;
    double	q2;
    double	q3;
    
    m_UpperExtremeValue = new double[m_AttributeIndices.length];
    m_UpperOutlier      = new double[m_AttributeIndices.length];
    m_LowerOutlier      = new double[m_AttributeIndices.length];
    m_LowerExtremeValue = new double[m_AttributeIndices.length];
    m_Median            = new double[m_AttributeIndices.length];
    m_IQR               = new double[m_AttributeIndices.length];
    
    for (i = 0; i < m_AttributeIndices.length; i++) {
      // non-numeric attribute?
      if (m_AttributeIndices[i] == NON_NUMERIC)
	continue;
      
      // sort attribute data
      values        = instances.attributeToDoubleArray(m_AttributeIndices[i]);
      sortedIndices = Utils.sort(values);
      
      // determine indices
      half    = sortedIndices.length / 2;
      quarter = half / 2;
      
      if (sortedIndices.length % 2 == 1) {
	q2 = values[sortedIndices[half]];
      }
      else {
	q2 = (values[sortedIndices[half]] + values[sortedIndices[half + 1]]) / 2;
      }
      
      if (half % 2 == 1) {
	q1 = values[sortedIndices[quarter]];
	q3 = values[sortedIndices[sortedIndices.length - quarter - 1]];
      }
      else {
	q1 = (values[sortedIndices[quarter]] + values[sortedIndices[quarter + 1]]) / 2;
	q3 = (values[sortedIndices[sortedIndices.length - quarter - 1]] + values[sortedIndices[sortedIndices.length - quarter]]) / 2;
      }
      
      // determine thresholds and other values
      m_Median[i]            = q2;
      m_IQR[i]               = q3 - q1;
      m_UpperExtremeValue[i] = q3 + getExtremeValuesFactor() * m_IQR[i];
      m_UpperOutlier[i]      = q3 + getOutlierFactor()       * m_IQR[i];
      m_LowerOutlier[i]      = q1 - getOutlierFactor()       * m_IQR[i];
      m_LowerExtremeValue[i] = q1 - getExtremeValuesFactor() * m_IQR[i];
    }
  }
  
  /**
   * returns whether the instance has an outlier in the specified attribute 
   * or not
   * 
   * @param inst	the instance to test
   * @param index	the attribute index
   * @return		true if the instance is an outlier
   */
  protected boolean isOutlier(Instance inst, int index) {
    boolean	result;
    double	value;

    value  = inst.value(m_AttributeIndices[index]);
    result =    ((m_UpperOutlier[index]      <  value) && (value <= m_UpperExtremeValue[index]))
             || ((m_LowerExtremeValue[index] <= value) && (value <  m_LowerOutlier[index]));
    
    return result;
  }
  
  /**
   * returns whether the instance is an outlier or not
   * 
   * @param inst	the instance to test
   * @return		true if the instance is an outlier
   */
  protected boolean isOutlier(Instance inst) {
    boolean	result;
    int		i;

    result = false;
    
    for (i = 0; i < m_AttributeIndices.length; i++) {
      // non-numeric attribute?
      if (m_AttributeIndices[i] == NON_NUMERIC)
	continue;

      result = isOutlier(inst, i);
      
      if (result)
	break;
    }
    
    return result;
  }
  
  /**
   * returns whether the instance has an extreme value in the specified 
   * attribute or not
   * 
   * @param inst	the instance to test
   * @param index	the attribute index
   * @return		true if the instance is an extreme value
   */
  protected boolean isExtremeValue(Instance inst, int index) {
    boolean	result;
    double	value;

    value  = inst.value(m_AttributeIndices[index]);
    result =    (value > m_UpperExtremeValue[index]) 
             || (value < m_LowerExtremeValue[index]);
      
    return result;
  }
  
  /**
   * returns whether the instance is an extreme value or not
   * 
   * @param inst	the instance to test
   * @return		true if the instance is an extreme value
   */
  protected boolean isExtremeValue(Instance inst) {
    boolean	result;
    int		i;

    result = false;
    
    for (i = 0; i < m_AttributeIndices.length; i++) {
      // non-numeric attribute?
      if (m_AttributeIndices[i] == NON_NUMERIC)
	continue;
      
      result = isExtremeValue(inst, i);
      
      if (result)
	break;
    }
    
    return result;
  }
  
  /**
   * returns the mulitplier of the IQR the instance is off the median for this
   * particular attribute.
   * 
   * @param inst	the instance to test
   * @param index	the attribute index
   * @return		the multiplier
   */
  protected double calculateMultiplier(Instance inst, int index) {
    double	result;
    double	value;

    value  = inst.value(m_AttributeIndices[index]);
    result = (value - m_Median[index]) / m_IQR[index];
      
    return result;
  }
  
  /**
   * Processes the given data (may change the provided dataset) and returns
   * the modified version. This method is called in batchFinished().
   * This implementation only calls process(Instance) for each instance
   * in the given dataset.
   *
   * @param instances   the data to process
   * @return            the modified data
   * @throws Exception  in case the processing goes wrong
   * @see               #batchFinished()
   */
  protected Instances process(Instances instances) throws Exception {
    Instances	result;
    Instance	instOld;
    Instance	instNew;
    int		i;
    int		n;
    double[]	values;
    int		numAttNew;
    int		numAttOld;
    
    if (!isFirstBatchDone())
      computeThresholds(instances);
    
    result    = getOutputFormat();
    numAttOld = instances.numAttributes();
    numAttNew = result.numAttributes();
    
    for (n = 0; n < instances.numInstances(); n++) {
      instOld = instances.instance(n);
      values  = new double[numAttNew];
      System.arraycopy(instOld.toDoubleArray(), 0, values, 0, numAttOld);
      
      // generate new instance
      //RANKING BEGIN
      if(instOld instanceof PreferenceDenseInstance){
    	  PreferenceDenseInstance pdi = (PreferenceDenseInstance) instOld;
    	  instNew = new PreferenceDenseInstance(1.0, values, pdi.getHashMap());
      }
      else
    	  instNew = new DenseInstance(1.0, values);
      //RANKING END
      instNew.setDataset(result);

      // per attribute?
      if (!getDetectionPerAttribute()) {
	// outlier?
	if (isOutlier(instOld))
	  instNew.setValue(m_OutlierAttributePosition[0], 1);
	// extreme value?
	if (isExtremeValue(instOld)) {
	  instNew.setValue(m_OutlierAttributePosition[0] + 1, 1);
	  // tag extreme values also as outliers?
	  if (getExtremeValuesAsOutliers())
	    instNew.setValue(m_OutlierAttributePosition[0], 1);
	}
      }
      else {
	for (i = 0; i < m_AttributeIndices.length; i++) {
	  // non-numeric attribute?
	  if (m_AttributeIndices[i] == NON_NUMERIC)
	    continue;
	  
	  // outlier?
	  if (isOutlier(instOld, m_AttributeIndices[i]))
	    instNew.setValue(m_OutlierAttributePosition[i], 1);
	  // extreme value?
	  if (isExtremeValue(instOld, m_AttributeIndices[i])) {
	    instNew.setValue(m_OutlierAttributePosition[i] + 1, 1);
	    // tag extreme values also as outliers?
	    if (getExtremeValuesAsOutliers())
	      instNew.setValue(m_OutlierAttributePosition[i], 1);
	  }
	  // add multiplier?
	  if (getOutputOffsetMultiplier())
	    instNew.setValue(
		m_OutlierAttributePosition[i] + 2, 
		calculateMultiplier(instOld, m_AttributeIndices[i]));
	}
      }
      
      // copy possible strings, relational values...
      copyValues(instNew, false, instOld.dataset(), getOutputFormat());
      
      // add to output
      result.add(instNew);
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
   * Main method for testing this class.
   *
   * @param args should contain arguments to the filter: use -h for help
   */
  public static void main(String[] args) {
    runFilter(new InterquartileRange(), args);
  }
}
