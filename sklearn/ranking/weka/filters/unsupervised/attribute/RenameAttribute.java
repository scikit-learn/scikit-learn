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
 *    RenameAttribute.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Vector;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.Range;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.SimpleStreamFilter;

/**
 <!-- globalinfo-start -->
 * This filter is used for renaming attribute names.<br/>
 * Regular expressions can be used in the matching and replacing.<br/>
 * See Javadoc of java.util.regex.Pattern class for more information:<br/>
 * http://java.sun.com/javase/6/docs/api/java/util/regex/Pattern.html
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -find &lt;regexp&gt;
 *  The regular expression that the attribute names must match.
 *  (default: ([\s\S]+))</pre>
 * 
 * <pre> -replace &lt;regexp&gt;
 *  The regular expression to replace matching attributes with.
 *  (default: $0)</pre>
 * 
 * <pre> -all
 *  Replaces all occurrences instead of just the first.
 *  (default: only first occurrence)</pre>
 * 
 * <pre> -R &lt;range&gt;
 *  The attribute range to work on.
 * This is a comma separated list of attribute indices, with "first" and "last" valid values.
 *  Specify an inclusive range with "-".
 *  E.g: "first-3,5,6-10,last".
 *  (default: first-last)</pre>
 * 
 * <pre> -V
 *  Inverts the attribute selection range.
 *  (default: off)</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6108 $
 */
public class RenameAttribute
  extends SimpleStreamFilter {

  /** for serialization. */
  private static final long serialVersionUID = 4216491776378279596L;

  /** the regular expression that the attribute names have to match. */
  protected String m_Find = "([\\s\\S]+)";
  
  /** the regular expression to replace the attribute name with. */
  protected String m_Replace = "$0";
  
  /** the attribute range to work on. */
  protected Range m_AttributeIndices = new Range("first-last");

  /** whether to replace all occurrences or just the first. */
  protected boolean m_ReplaceAll = false;
  
  /**
   * Returns a string describing this filter.
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "This filter is used for renaming attribute names.\n"
      + "Regular expressions can be used in the matching and replacing.\n"
      + "See Javadoc of java.util.regex.Pattern class for more information:\n"
      + "http://java.sun.com/javase/6/docs/api/java/util/regex/Pattern.html";
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
	"\tThe regular expression that the attribute names must match.\n"
	+ "\t(default: ([\\s\\S]+))",
	"find", 1, "-find <regexp>"));

    result.addElement(new Option(
	"\tThe regular expression to replace matching attributes with.\n"
	+ "\t(default: $0)",
	"replace", 1, "-replace <regexp>"));

    result.addElement(new Option(
	"\tReplaces all occurrences instead of just the first.\n"
	+ "\t(default: only first occurrence)",
	"all", 0, "-all"));

    result.addElement(new Option(
	"\tThe attribute range to work on.\n"
	+ "This is a comma separated list of attribute indices, with "
        + "\"first\" and \"last\" valid values.\n"
        + "\tSpecify an inclusive range with \"-\".\n"
        + "\tE.g: \"first-3,5,6-10,last\".\n"
	+ "\t(default: first-last)",
	"R", 1, "-R <range>"));

    result.addElement(new Option(
	"\tInverts the attribute selection range.\n"
	+ "\t(default: off)",
	"V", 0, "-V"));

    return result.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -find &lt;regexp&gt;
   *  The regular expression that the attribute names must match.
   *  (default: ([\s\S]+))</pre>
   * 
   * <pre> -replace &lt;regexp&gt;
   *  The regular expression to replace matching attributes with.
   *  (default: $0)</pre>
   * 
   * <pre> -all
   *  Replaces all occurrences instead of just the first.
   *  (default: only first occurrence)</pre>
   * 
   * <pre> -R &lt;range&gt;
   *  The attribute range to work on.
   * This is a comma separated list of attribute indices, with "first" and "last" valid values.
   *  Specify an inclusive range with "-".
   *  E.g: "first-3,5,6-10,last".
   *  (default: first-last)</pre>
   * 
   * <pre> -V
   *  Inverts the attribute selection range.
   *  (default: off)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption("find", options);
    if (tmpStr.length() != 0)
      setFind(tmpStr);
    else
      setFind("([\\s\\S]+)");
    
    tmpStr = Utils.getOption("replace", options);
    if (tmpStr.length() != 0)
      setReplace(tmpStr);
    else
      setReplace("$0");

    setReplaceAll(Utils.getFlag("all", options));
    
    tmpStr = Utils.getOption("R", options);
    if (tmpStr.length() != 0)
      setAttributeIndices(tmpStr);
    else
      setAttributeIndices("first-last");

    setInvertSelection(Utils.getFlag("V", options));
    
    if (getInputFormat() != null)
      setInputFormat(getInputFormat());
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;

    result = new Vector<String>(Arrays.asList(super.getOptions()));
    
    result.add("-find");
    result.add(getFind());
    
    result.add("-replace");
    result.add(getReplace());
    
    if (getReplaceAll())
      result.add("-all");
    
    result.add("-R");
    result.add(getAttributeIndices());
    
    if (getInvertSelection())
      result.add("-V");

    return result.toArray(new String[result.size()]);
  }

  /**
   * Sets the regular expression that the attribute names must match.
   *
   * @param value 	the regular expression
   */
  public void setFind(String value) {
    m_Find = value;
  }

  /**
   * Returns the current regular expression for .
   *
   * @return 		a string containing a comma separated list of ranges
   */
  public String getFind() {
    return m_Find;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String findTipText() {
    return "The regular expression that the attribute names must match.";
  }

  /**
   * Sets the regular expression to replace matching attribute names with.
   *
   * @param value 	the regular expression
   */
  public void setReplace(String value) {
    m_Replace = value;
  }

  /**
   * Returns the regular expression to replace matching attribute names with.
   *
   * @return 		the regular expression
   */
  public String getReplace() {
    return m_Replace;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String replaceTipText() {
    return 
        "The regular expression to use for replacing the matching attribute "
      + "names with.";
  }

  /**
   * Sets whether to replace all occurrences or just the first one.
   *
   * @param value 	if true then all occurrences are replace
   */
  public void setReplaceAll(boolean value) {
    m_ReplaceAll = value;
  }

  /**
   * Returns whether all occurrences are replaced or just the first one.
   *
   * @return 		true if all occurrences are replaced
   */
  public boolean getReplaceAll() {
    return m_ReplaceAll;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String replaceAllTipText() {
    return 
        "If set to true, then all occurrences of the match will be replaced; "
      + "otherwise only the first.";
  }

  /**
   * Sets which attributes are to be acted on.
   *
   * @param value 	a string representing the list of attributes. Since
   * 			the string will typically come from a user, attributes 
   * 			are indexed from1. <br/>
   * 			eg: first-3,5,6-last
   */
  public void setAttributeIndices(String value) {
    m_AttributeIndices.setRanges(value);
  }

  /**
   * Gets the current range selection.
   *
   * @return a string containing a comma separated list of ranges
   */
  public String getAttributeIndices() {
    return m_AttributeIndices.getRanges();
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String attributeIndicesTipText() {
    return 
        "Specify range of attributes to act on; "
      + "this is a comma separated list of attribute indices, with "
      + "\"first\" and \"last\" valid values; specify an inclusive "
      + "range with \"-\"; eg: \"first-3,5,6-10,last\".";
  }

  /**
   * Sets whether to invert the selection of the attributes.
   *
   * @param value 	if true then the selection is inverted
   */
  public void setInvertSelection(boolean value) {
    m_AttributeIndices.setInvert(value);
  }

  /**
   * Gets whether to invert the selection of the attributes.
   *
   * @return 		true if the selection is inverted
   */
  public boolean getInvertSelection() {
    return m_AttributeIndices.getInvert();
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String invertSelectionTipText() {
    return 
        "If set to true, the selection will be inverted; eg: the attribute "
      + "indices '2-4' then mean everything apart from '2-4'.";
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
   */
  protected Instances determineOutputFormat(Instances inputFormat) throws Exception {
    Instances			result;
    Attribute			att;
    ArrayList<Attribute>	atts;
    int				i;
    
    m_AttributeIndices.setUpper(inputFormat.numAttributes() - 1);
    
    // generate new header
    atts = new ArrayList<Attribute>();
    for (i = 0; i < inputFormat.numAttributes(); i++) {
      att = inputFormat.attribute(i);
      if (m_AttributeIndices.isInRange(i)) {
	if (m_ReplaceAll)
	  atts.add(att.copy(att.name().replaceAll(m_Find, m_Replace)));
	else
	  atts.add(att.copy(att.name().replaceFirst(m_Find, m_Replace)));
      }
      else {
	atts.add((Attribute) att.copy());
      }
    }
    result = new Instances(inputFormat.relationName(), atts, 0);
    result.setClassIndex(inputFormat.classIndex());
    
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
    return (Instance) instance.copy();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6108 $");
  }

  /**
   * Main method for executing this filter.
   *
   * @param args 	the arguments to the filter: use -h for help
   */
  public static void main(String[] args) {
    runFilter(new RenameAttribute(), args);
  }
}
