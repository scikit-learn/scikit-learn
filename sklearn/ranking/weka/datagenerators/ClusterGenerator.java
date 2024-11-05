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
 * ClusterGenerator.java
 * Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.datagenerators;

import weka.core.Option;
import weka.core.Range;
import weka.core.Utils;
import weka.datagenerators.DataGenerator;

import java.util.Enumeration;
import java.util.Vector;

/** 
 * Abstract class for cluster data generators. <p/>
 *
 * Example usage as the main of a datagenerator called RandomGenerator:
 * <pre>
 * public static void main(String[] args) {
 *   try {
 *     DataGenerator.makeData(new RandomGenerator(), args);
 *   } 
 *   catch (Exception e) {
 *     e.printStackTrace();
 *     System.err.println(e.getMessage());
 *   }
 * }
 * </pre>
 * <p/>
 *
 * @author Gabi Schmidberger (gabi@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.6 $
 */
public abstract class ClusterGenerator 
  extends DataGenerator {

  /** for serialization */
  private static final long serialVersionUID = 6131722618472046365L;

  /** Number of attribute the dataset should have */
  protected int m_NumAttributes;

  /** class flag  */
  protected boolean m_ClassFlag = false;

  /** Stores which columns are boolean (default numeric) */
  protected Range m_booleanCols;

  /** Stores which columns are nominal (default numeric)  */
  protected Range m_nominalCols;

  /**
   * initializes the generator 
   */
  public ClusterGenerator() {
    super();

    setNumAttributes(defaultNumAttributes());
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = enumToVector(super.listOptions());

    result.addElement(new Option(
          "\tThe number of attributes (default " 
          + defaultNumAttributes() + ").",
          "a", 1, "-a <num>"));

    result.addElement(new Option(
        "\tClass Flag, if set, the cluster is listed in extra attribute.",
        "c", 0, "-c"));
    
    result.addElement(new Option(
        "\tThe indices for boolean attributes.",
        "b", 1, "-b <range>"));
    
    result.addElement(new Option(
        "\tThe indices for nominal attributes.",
        "m", 1, "-m <range>"));
    
    return result.elements();
  }

  /**
   * Sets the options.
   *
   * @param options the options 
   * @throws Exception if invalid option
   */
  public void setOptions(String[] options) throws Exception { 
    String        tmpStr;
   
    super.setOptions(options);

    tmpStr = Utils.getOption('a', options);
    if (tmpStr.length() != 0)
      setNumAttributes(Integer.parseInt(tmpStr));
    else
      setNumAttributes(defaultNumAttributes());

    setClassFlag(Utils.getFlag('c', options));

    tmpStr = Utils.getOption('b', options);
    setBooleanIndices(tmpStr);
    m_booleanCols.setUpper(getNumAttributes());

    tmpStr = Utils.getOption('m', options);
    setNominalIndices(tmpStr);
    m_nominalCols.setUpper(getNumAttributes());

    // check indices
    tmpStr = checkIndices();
    if (tmpStr.length() > 0)
      throw new IllegalArgumentException(tmpStr);
  }
  
  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector        result;
    String[]      options;
    int           i;
    
    result  = new Vector();
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    result.add("-a");
    result.add("" + getNumAttributes());

    if (getClassFlag())
      result.add("-c");
    
    if (!getBooleanCols().toString().equalsIgnoreCase("empty")) {
      result.add("-b");
      result.add("" + getBooleanCols());
    }
    
    if (!getNominalCols().toString().equalsIgnoreCase("empty")) {
      result.add("-m");
      result.add("" + getNominalCols());
    }
    
    return (String[]) result.toArray(new String[result.size()]);
  }

  /**
   * returns the default number of attributes
   * 
   * @return the default number of attributes
   */
  protected int defaultNumAttributes() {
    return 10;
  }

  /**
   * Sets the number of attributes the dataset should have.
   * @param numAttributes the new number of attributes
   */
  public void setNumAttributes(int numAttributes) {
    m_NumAttributes = numAttributes;
    getBooleanCols().setUpper(getNumAttributes());
    getNominalCols().setUpper(getNumAttributes());
  }

  /**
   * Gets the number of attributes that should be produced.
   * @return the number of attributes that should be produced
   */
  public int getNumAttributes() { 
    return m_NumAttributes; 
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String numAttributesTipText() {
    return "The number of attributes the generated data will contain.";
  }

  /**
   * Sets the class flag, if class flag is set, 
   * the cluster is listed as class atrribute in an extra attribute.
   * @param classFlag the new class flag
   */
  public void setClassFlag(boolean classFlag) { 
    m_ClassFlag = classFlag; 
  }

  /**
   * Gets the class flag.
   * @return the class flag 
   */
  public boolean getClassFlag() {
    return m_ClassFlag; 
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String classFlagTipText() {
    return "If set to TRUE, lists the cluster as an extra attribute.";
  }

  /**
   * Sets which attributes are boolean 
   * @param rangeList a string representing the list of attributes. Since
   * the string will typically come from a user, attributes are indexed from
   * 1. <br/>
   * eg: first-3,5,6-last
   * @throws IllegalArgumentException if an invalid range list is supplied 
   */
  public void setBooleanIndices(String rangeList) {
    m_booleanCols.setRanges(rangeList);
  }

  /**
   * Sets which attributes are boolean.
   * @param value the range to use
   */
  public void setBooleanCols(Range value) {
    m_booleanCols.setRanges(value.getRanges());
  }

  /**
   * returns the range of boolean attributes.
   * 
   * @return the range of boolean attributes
   */
  public Range getBooleanCols() {
    if (m_booleanCols == null)
      m_booleanCols = new Range();

    return m_booleanCols;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String booleanColsTipText() {
    return "The range of attributes that are generated as boolean ones.";
  }

  /**
   * Sets which attributes are nominal
   * @param rangeList a string representing the list of attributes. Since
   * the string will typically come from a user, attributes are indexed from
   * 1. <br/>
   * eg: first-3,5,6-last
   * @throws IllegalArgumentException if an invalid range list is supplied 
   */
  public void setNominalIndices(String rangeList) {
    m_nominalCols.setRanges(rangeList);
  }

  /**
   * Sets which attributes are nominal.
   * @param value the range to use
   */
  public void setNominalCols(Range value) {
    m_nominalCols.setRanges(value.getRanges());
  }

  /**
   * returns the range of nominal attributes
   * 
   * @return the range of nominal attributes
   */
  public Range getNominalCols() {
    if (m_nominalCols == null)
      m_nominalCols = new Range();

    return m_nominalCols;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String nominalColsTipText() {
    return "The range of attributes to generate as nominal ones.";
  }

  /**
   * check if attribute types are not contradicting
   * 
   * @return empty string if no problem, otherwise error message
   */
  protected String checkIndices() {
    for (int i = 1; i < getNumAttributes() + 1; i++) {
      m_booleanCols.isInRange(i);
      if (m_booleanCols.isInRange(i) && m_nominalCols.isInRange(i)) {
	return   "Error in attribute type: Attribute " 
               + i + " is set boolean and nominal.";
      }
    } 
    return "";
  }
}


