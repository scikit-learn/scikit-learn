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
 * SortLabels.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
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
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.SimpleStreamFilter;

import java.io.Serializable;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * A simple filter for sorting the labels of nominal attributes.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns on output of debugging information.</pre>
 * 
 * <pre> -R &lt;index1,index2-index4,...&gt;
 *  Specify list of string attributes to convert to words.
 *  (default: select all relational attributes)</pre>
 * 
 * <pre> -V
 *  Inverts the matching sense of the selection.</pre>
 * 
 * <pre> -S &lt;CASE|NON-CASE&gt;
 *  Determines the type of sorting:
 *  CASE = Case-sensitive
 *  NON-CASE = Case-insensitive
 *  (default: CASE)</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class SortLabels
  extends SimpleStreamFilter {

  /** for serialization. */
  private static final long serialVersionUID = 7815204879694105691L;
  
  /**
   * Represents a case-sensitive comparator for two strings.
   *
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5987 $
   */
  public static class CaseSensitiveComparator
    implements Comparator, Serializable {

    /** for serialization. */
    private static final long serialVersionUID = 7071450356783873277L;

    /**
     * compares the two strings, returns -1 if o1 is smaller than o2, 0
     * if equal and +1 if greater.
     * 
     * @param o1	the first string to compare
     * @param o2	the second string to compare
     * @return		returns -1 if o1 is smaller than o2, 0 if equal and +1 
     *			if greater
     */
    public int compare(Object o1, Object o2) {
      String	s1;
      String	s2;
      
      if ((o1 == null) && (o2 == null))
	return 0;
      else if (o1 == null)
	return -1;
      else if (o2 == null)
	return +1;
      
      s1 = (String) o1;
      s2 = (String) o2;
      
      return s1.compareTo(s2);
    }
  }
  
  /**
   * Represents a case-insensitive comparator for two strings.
   *
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5987 $
   */
  public static class CaseInsensitiveComparator
    implements Comparator, Serializable {

    /** for serialization. */
    private static final long serialVersionUID = -4515292733342486066L;

    /**
     * compares the two strings, returns -1 if o1 is smaller than o2, 0
     * if equal and +1 if greater.
     * 
     * @param o1	the first string to compare
     * @param o2	the second string to compare
     * @return		returns -1 if o1 is smaller than o2, 0 if equal and +1 
     *			if greater
     */
    public int compare(Object o1, Object o2) {
      String	s1;
      String	s2;
      
      if ((o1 == null) && (o2 == null))
	return 0;
      else if (o1 == null)
	return -1;
      else if (o2 == null)
	return +1;
      
      s1 = (String) o1;
      s2 = (String) o2;
      
      return s1.toLowerCase().compareTo(s2.toLowerCase());
    }
  }
  
  /** sorts the strings case-sensitive. */
  public final static int SORT_CASESENSITIVE = 0;
  
  /** sorts the strings case-insensitive. */
  public final static int SORT_CASEINSENSITIVE = 1;
  
  /** Tag allowing selection of sort type. */
  public final static Tag[] TAGS_SORTTYPE = {
    new Tag(SORT_CASESENSITIVE, "case", "Case-sensitive"),
    new Tag(SORT_CASEINSENSITIVE, "non-case", "Case-insensitive")
  };
  
  /** the range of attributes to process (only relational ones will be processed). */
  protected Range m_AttributeIndices = new Range("first-last");

  /** the new order for the labels. */
  protected int[][] m_NewOrder = null;

  /** the sort type. */
  protected int m_SortType = SORT_CASEINSENSITIVE;
  
  /** the comparator to use for sorting. */
  protected Comparator m_Comparator = new CaseSensitiveComparator();
  
  /**
   * Returns a string describing this filter.
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "A simple filter for sorting the labels of nominal attributes.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector        	result;
    Enumeration   	en;
    String		desc;
    int			i;
    SelectedTag		tag;

    result = new Vector();

    en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());

    result.addElement(new Option(
	"\tSpecify list of string attributes to convert to words.\n"
	+ "\t(default: select all relational attributes)",
	"R", 1, "-R <index1,index2-index4,...>"));

    result.addElement(new Option(
	"\tInverts the matching sense of the selection.",
	"V", 0, "-V"));

    desc  = "";
    for (i = 0; i < TAGS_SORTTYPE.length; i++) {
      tag = new SelectedTag(TAGS_SORTTYPE[i].getID(), TAGS_SORTTYPE);
      desc  +=   "\t" + tag.getSelectedTag().getIDStr() 
      	       + " = " + tag.getSelectedTag().getReadable()
      	       + "\n";
    }
    result.addElement(new Option(
	"\tDetermines the type of sorting:\n"
	+ desc
	+ "\t(default: " + new SelectedTag(SORT_CASESENSITIVE, TAGS_SORTTYPE) + ")",
	"S", 1, "-S " + Tag.toOptionList(TAGS_SORTTYPE)));

    return result.elements();
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
   * <pre> -R &lt;index1,index2-index4,...&gt;
   *  Specify list of string attributes to convert to words.
   *  (default: select all relational attributes)</pre>
   * 
   * <pre> -V
   *  Inverts the matching sense of the selection.</pre>
   * 
   * <pre> -S &lt;CASE|NON-CASE&gt;
   *  Determines the type of sorting:
   *  CASE = Case-sensitive
   *  NON-CASE = Case-insensitive
   *  (default: CASE)</pre>
   * 
   <!-- options-end -->
   *
   * @param options	the options to use
   * @throws Exception	if setting of options fails
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;

    tmpStr = Utils.getOption('R', options);
    if (tmpStr.length() != 0)
      setAttributeIndices(tmpStr);
    else
      setAttributeIndices("first-last");

    setInvertSelection(Utils.getFlag('V', options));

    tmpStr = Utils.getOption('S', options);
    if (tmpStr.length() != 0)
      setSortType(new SelectedTag(tmpStr, TAGS_SORTTYPE));
    else
      setSortType(new SelectedTag(SORT_CASESENSITIVE, TAGS_SORTTYPE));

    super.setOptions(options);
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int			i;
    Vector<String>	result;
    String[]		options;

    result = new Vector<String>();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-R"); 
    result.add(getAttributeIndices().getRanges());

    if (getInvertSelection())
      result.add("-V");

    result.add("-S");
    result.add("" + getSortType());
    
    return result.toArray(new String[result.size()]);	  
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
      + "\"first\" and \"last\" valid values; Specify an inclusive "
      + "range with \"-\"; eg: \"first-3,5,6-10,last\".";
  }
  
  /**
   * Set the range of attributes to process.
   *
   * @param value 	the new range.
   */
  public void setAttributeIndices(String value) {
    m_AttributeIndices = new Range(value);
  }

  /**
   * Gets the current selected attributes.
   *
   * @return 		current selection.
   */
  public Range getAttributeIndices() {
    return m_AttributeIndices;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String invertSelectionTipText() {
    return 
        "Set attribute selection mode. If false, only selected "
      + "attributes in the range will be worked on; if "
      + "true, only non-selected attributes will be processed.";
  }

  /**
   * Sets whether selected columns should be processed or skipped.
   *
   * @param value 	the new invert setting
   */
  public void setInvertSelection(boolean value) {
    m_AttributeIndices.setInvert(value);
  }

  /**
   * Gets whether the supplied columns are to be processed or skipped.
   *
   * @return 		true if the supplied columns will be kept
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
  public String sortTypeTipText() {
    return "The type of sorting to use.";
  }

  /**
   * Sets the sort type to be used.
   *
   * @param type 	the type of sorting
   */
  public void setSortType(SelectedTag type) {
    if (type.getTags() == TAGS_SORTTYPE) {
      m_SortType = type.getSelectedTag().getID();
      
      if (m_SortType == SORT_CASESENSITIVE)
	m_Comparator = new CaseSensitiveComparator();
      else if (m_SortType == SORT_CASEINSENSITIVE)
	m_Comparator = new CaseInsensitiveComparator();
      else
	throw new IllegalStateException("Unhandled sort type '" + type + "'!");
    }
  }

  /**
   * Gets the sort type to be used.
   *
   * @return 		the sort type
   */
  public SelectedTag getSortType() {
    return new SelectedTag(m_SortType, TAGS_SORTTYPE);
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
    
    //RANKING BEGIN
    result.disable(Capability.RANKING);
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    //RANKING END
    
    return result;
  }
  
  /**
   * Determines the output format based on the input format and returns 
   * this. 
   *
   * @param inputFormat     the input format to base the output format on
   * @return                the output format
   * @throws Exception      in case the determination goes wrong
   */
  protected Instances determineOutputFormat(Instances inputFormat) throws Exception {
    Instances		result;
    Attribute		att;
    Attribute		attSorted;
    FastVector		atts;
    FastVector		values;
    Vector<String>	sorted;
    int			i;
    int			n;
    
    m_AttributeIndices.setUpper(inputFormat.numAttributes() - 1);
    
    // determine sorted indices
    atts       = new FastVector();
    m_NewOrder = new int[inputFormat.numAttributes()][];
    for (i = 0; i < inputFormat.numAttributes(); i++) {
      att = inputFormat.attribute(i);
      if (!att.isNominal() && !att.isRanking() || !m_AttributeIndices.isInRange(i)) {
	m_NewOrder[i] = new int[0];
	atts.addElement(inputFormat.attribute(i).copy());
	continue;
      }

      // sort labels
      sorted = new Vector<String>();
      for (n = 0; n < att.numValues(); n++)
	sorted.add(att.value(n));
      Collections.sort(sorted, m_Comparator);
      
      // determine new indices
      m_NewOrder[i] = new int[att.numValues()];
      values        = new FastVector();
      for (n = 0; n < att.numValues(); n++) {
	m_NewOrder[i][n] = sorted.indexOf(att.value(n));
	values.addElement(sorted.get(n));
      }
      attSorted = new Attribute(att.name(), values);
      attSorted.setWeight(att.weight());
      atts.addElement(attSorted);
    }
    
    // generate new header
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
    Instance	result;
    Attribute	att;
    double[]	values;
    int		i;

    // adjust indices
    values = new double[instance.numAttributes()];
    for (i = 0; i < instance.numAttributes(); i++) {
      att = instance.attribute(i);
      if (!att.isNominal() && !att.isRanking() || !m_AttributeIndices.isInRange(i) || instance.isMissing(i))
	values[i] = instance.value(i);
      else
	values[i] = m_NewOrder[i][(int) instance.value(i)];
    }

    // create new instance
    result = new DenseInstance(instance.weight(), values);
    
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
   * runs the filter with the given arguments.
   *
   * @param args      the commandline arguments
   */
  public static void main(String[] args) {
    runFilter(new SortLabels(), args);
  }
}
