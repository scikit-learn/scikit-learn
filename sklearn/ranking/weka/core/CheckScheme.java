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
 * CheckScheme.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.util.Enumeration;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.Vector;

import weka.core.labelranking.PreferenceAttribute;

/**
 * Abstract general class for testing schemes in Weka. Derived classes are
 * also used for JUnit tests.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 * @see TestInstances
 */
public abstract class CheckScheme
  extends Check {
  
  /** a class for postprocessing the test-data */
  public static class PostProcessor
    implements RevisionHandler {
    
    /**
     * Provides a hook for derived classes to further modify the data. Currently,
     * the data is just passed through.
     * 
     * @param data	the data to process
     * @return		the processed data
     */
    public Instances process(Instances data) {
      return data;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 5953 $");
    }
  }
  
  /** The number of instances in the datasets */
  protected int m_NumInstances = 20;
  
  /** the number of nominal attributes */
  protected int m_NumNominal = 2;
  
  /** the number of numeric attributes */
  protected int m_NumNumeric = 1;
  
  /** the number of string attributes */
  protected int m_NumString = 1;
  
  /** the number of date attributes */
  protected int m_NumDate = 1;
  
  /** the number of relational attributes */
  protected int m_NumRelational = 1;
  
  /** the number of instances in relational attributes (applies also for bags
   * in multi-instance) */
  protected int m_NumInstancesRelational = 10;
  
  /** for generating String attributes/classes */
  protected String[] m_Words = TestInstances.DEFAULT_WORDS;
  
  /** for generating String attributes/classes */
  protected String m_WordSeparators = TestInstances.DEFAULT_SEPARATORS;
  
  /** for post-processing the data even further */
  protected PostProcessor m_PostProcessor = null;
  
  /** whether classpath problems occurred */
  protected boolean m_ClasspathProblems = false;
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> result = new Vector<Option>();
    
    Enumeration en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement((Option)en.nextElement());
    
    result.addElement(new Option(
        "\tThe number of instances in the datasets (default 20).",
        "N", 1, "-N <num>"));

    result.addElement(new Option(
        "\tThe number of nominal attributes (default 2).",
        "nominal", 1, "-nominal <num>"));
    
    result.addElement(new Option(
        "\tThe number of values for nominal attributes (default 1).",
        "nominal-values", 1, "-nominal-values <num>"));
    
    result.addElement(new Option(
        "\tThe number of numeric attributes (default 1).",
        "numeric", 1, "-numeric <num>"));
    
    result.addElement(new Option(
        "\tThe number of string attributes (default 1).",
        "string", 1, "-string <num>"));
    
    result.addElement(new Option(
        "\tThe number of date attributes (default 1).",
        "date", 1, "-date <num>"));
    
    result.addElement(new Option(
        "\tThe number of relational attributes (default 1).",
        "relational", 1, "-relational <num>"));
    
    result.addElement(new Option(
        "\tThe number of instances in relational/bag attributes (default 10).",
        "num-instances-relational", 1, "-num-instances-relational <num>"));
    
    result.addElement(new Option(
        "\tThe words to use in string attributes.",
        "words", 1, "-words <comma-separated-list>"));
    
    result.addElement(new Option(
        "\tThe word separators to use in string attributes.",
        "word-separators", 1, "-word-separators <chars>"));
    
    return result.elements();
  }
  
  /**
   * Parses a given list of options. 
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      tmpStr;
    
    super.setOptions(options);
    
    tmpStr = Utils.getOption('N', options);
    if (tmpStr.length() != 0)
      setNumInstances(Integer.parseInt(tmpStr));
    else
      setNumInstances(20);
    
    tmpStr = Utils.getOption("nominal", options);
    if (tmpStr.length() != 0)
      setNumNominal(Integer.parseInt(tmpStr));
    else
      setNumNominal(2);
    
    tmpStr = Utils.getOption("numeric", options);
    if (tmpStr.length() != 0)
      setNumNumeric(Integer.parseInt(tmpStr));
    else
      setNumNumeric(1);
    
    tmpStr = Utils.getOption("string", options);
    if (tmpStr.length() != 0)
      setNumString(Integer.parseInt(tmpStr));
    else
      setNumString(1);
    
    tmpStr = Utils.getOption("date", options);
    if (tmpStr.length() != 0)
      setNumDate(Integer.parseInt(tmpStr));
    else
      setNumDate(1);
    
    tmpStr = Utils.getOption("relational", options);
    if (tmpStr.length() != 0)
      setNumRelational(Integer.parseInt(tmpStr));
    else
      setNumRelational(1);
    
    tmpStr = Utils.getOption("num-instances-relational", options);
    if (tmpStr.length() != 0)
      setNumInstancesRelational(Integer.parseInt(tmpStr));
    else
      setNumInstancesRelational(10);
    
    tmpStr = Utils.getOption("words", options);
    if (tmpStr.length() != 0)
      setWords(tmpStr);
    else
      setWords(new TestInstances().getWords());
    
    if (Utils.getOptionPos("word-separators", options) > -1) {
      tmpStr = Utils.getOption("word-separators", options);
      setWordSeparators(tmpStr);
    }
    else {
      setWordSeparators(TestInstances.DEFAULT_SEPARATORS);
    }
  }
  
  /**
   * Gets the current settings of the CheckClassifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>        result;
    String[]      options;
    int           i;
    
    result = new Vector<String>();
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    result.add("-N");
    result.add("" + getNumInstances());
    
    result.add("-nominal");
    result.add("" + getNumNominal());
    
    result.add("-numeric");
    result.add("" + getNumNumeric());
    
    result.add("-string");
    result.add("" + getNumString());
    
    result.add("-date");
    result.add("" + getNumDate());
    
    result.add("-relational");
    result.add("" + getNumRelational());
    
    result.add("-words");
    result.add("" + getWords());
    
    result.add("-word-separators");
    result.add("" + getWordSeparators());
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * sets the PostProcessor to use
   * 
   * @param value	the new PostProcessor
   * @see #m_PostProcessor
   */
  public void setPostProcessor(PostProcessor value) {
    m_PostProcessor = value;
  }
  
  /**
   * returns the current PostProcessor, can be null
   * 
   * @return		the current PostProcessor
   */
  public PostProcessor getPostProcessor() {
    return m_PostProcessor;
  }
  
  /**
   * returns TRUE if the classifier returned a "not in classpath" Exception
   * 
   * @return	true if CLASSPATH problems occurred
   */
  public boolean hasClasspathProblems() {
    return m_ClasspathProblems;
  }
  
  /**
   * Begin the tests, reporting results to System.out
   */
  public abstract void doTests();
  
  /**
   * Sets the number of instances to use in the datasets (some classifiers
   * might require more instances).
   *
   * @param value the number of instances to use
   */
  public void setNumInstances(int value) {
    m_NumInstances = value;
  }
  
  /**
   * Gets the current number of instances to use for the datasets.
   *
   * @return the number of instances
   */
  public int getNumInstances() {
    return m_NumInstances;
  }
  
  /**
   * sets the number of nominal attributes
   * 
   * @param value	the number of nominal attributes
   */
  public void setNumNominal(int value) {
    m_NumNominal = value;
  }
  
  /**
   * returns the current number of nominal attributes
   * 
   * @return 		the number of nominal attributes
   */
  public int getNumNominal() {
    return m_NumNominal;
  }
  
  /**
   * sets the number of numeric attributes
   * 
   * @param value 	the number of numeric attributes
   */
  public void setNumNumeric(int value) {
    m_NumNumeric = value;
  }
  
  /**
   * returns the current number of numeric attributes
   * 
   * @return 		the number of numeric attributes
   */
  public int getNumNumeric() {
    return m_NumNumeric;
  }
  
  /**
   * sets the number of string attributes
   * 
   * @param value 	the number of string attributes
   */
  public void setNumString(int value) {
    m_NumString = value;
  }
  
  /**
   * returns the current number of string attributes
   * 
   * @return 		the number of string attributes
   */
  public int getNumString() {
    return m_NumString;
  }
  
  /**
   * sets the number of data attributes
   * 
   * @param value	the number of date attributes
   */
  public void setNumDate(int value) {
    m_NumDate = value;
  }
  
  /**
   * returns the current number of date attributes
   * 
   * @return		the number of date attributes
   */
  public int getNumDate() {
    return m_NumDate;
  }
  
  /**
   * sets the number of relational attributes
   * 
   * @param value	the number of relational attributes
   */
  public void setNumRelational(int value) {
    m_NumRelational = value;
  }
  
  /**
   * returns the current number of relational attributes
   * 
   * @return		the number of relational attributes
   */
  public int getNumRelational() {
    return m_NumRelational;
  }
  
  /**
   * sets the number of instances in relational/bag attributes to produce
   * 
   * @param value	the number of instances
   */
  public void setNumInstancesRelational(int value) {
    m_NumInstancesRelational = value;
  }
  
  /**
   * returns the current number of instances in relational/bag attributes to produce
   * 
   * @return		the number of instances
   */
  public int getNumInstancesRelational() {
    return m_NumInstancesRelational;
  }

  /**
   * turns the comma-separated list into an array
   * 
   * @param value	the list to process
   * @return		the list as array
   */
  protected static String[] listToArray(String value) {
    StringTokenizer	tok;
    Vector<String>		list;
    
    list = new Vector<String>();
    tok = new StringTokenizer(value, ",");
    while (tok.hasMoreTokens())
      list.add(tok.nextToken());
    
    return (String[]) list.toArray(new String[list.size()]);
  }
  
  /**
   * turns the array into a comma-separated list
   * 
   * @param value	the array to process
   * @return		the array as list
   */
  protected static String arrayToList(String[] value) {
    String	result;
    int		i;
    
    result = "";
    
    for (i = 0; i < value.length; i++) {
      if (i > 0)
	result += ",";
      result += value[i];
    }
    
    return result;
  }
  
  /**
   * returns a string representation of the attribute type
   * 
   * @param type	the attribute type to get a string rerpresentation for
   * @return		the string representation
   */
  public static String attributeTypeToString(int type) {
    String	result;
    
    switch (type) {
      case Attribute.NUMERIC:
	result = "numeric";
	break;
	
      case Attribute.NOMINAL:
	result = "nominal";
	break;
	
      case Attribute.STRING:
	result = "string";
	break;
	
      case Attribute.DATE:
	result = "date";
	break;
	
      case Attribute.RELATIONAL:
	result = "relational";
	break;
	
      case PreferenceAttribute.RANKING:
    result = "ranking";
    break;

      default:
	result = "???";
    }
    
    return result;
  }
  
  /**
   * Sets the comma-separated list of words to use for generating strings. The
   * list must contain at least 2 words, otherwise an exception will be thrown.
   * 
   * @param value			the list of words
   * @throws IllegalArgumentException	if not at least 2 words are provided
   */
  public void setWords(String value) {
    if (listToArray(value).length < 2)
      throw new IllegalArgumentException("At least 2 words must be provided!");
    
    m_Words = listToArray(value);
  }
  
  /**
   * returns the words used for assembling strings in a comma-separated list.
   * 
   * @return		the words as comma-separated list
   */
  public String getWords() {
    return arrayToList(m_Words);
  }

  /**
   * sets the word separators (chars) to use for assembling strings.
   * 
   * @param value	the characters to use as separators
   */
  public void setWordSeparators(String value) {
    m_WordSeparators = value;
  }
  
  /**
   * returns the word separators (chars) to use for assembling strings.
   * 
   * @return		the current separators
   */
  public String getWordSeparators() {
    return m_WordSeparators;
  }
  
  /**
   * Compare two datasets to see if they differ.
   *
   * @param data1 one set of instances
   * @param data2 the other set of instances
   * @throws Exception if the datasets differ
   */
  protected void compareDatasets(Instances data1, Instances data2)
    throws Exception {
    
    if (!data2.equalHeaders(data1)) {
      throw new Exception("header has been modified\n" + data2.equalHeadersMsg(data1));
    }
    if (!(data2.numInstances() == data1.numInstances())) {
      throw new Exception("number of instances has changed");
    }
    for (int i = 0; i < data2.numInstances(); i++) {
      Instance orig = data1.instance(i);
      Instance copy = data2.instance(i);
      for (int j = 0; j < orig.numAttributes(); j++) {
        if (orig.isMissing(j)) {
          if (!copy.isMissing(j)) {
            throw new Exception("instances have changed");
          }
        } else if (orig.value(j) != copy.value(j)) {
          throw new Exception("instances have changed");
        }
        if (orig.weight() != copy.weight()) {
          throw new Exception("instance weights have changed");
        }	  
      }
    }
  }
  
  /**
   * Add missing values to a dataset.
   *
   * @param data the instances to add missing values to
   * @param level the level of missing values to add (if positive, this
   * is the probability that a value will be set to missing, if negative
   * all but one value will be set to missing (not yet implemented))
   * @param predictorMissing if true, predictor attributes will be modified
   * @param classMissing if true, the class attribute will be modified
   */
  protected void addMissing(Instances data, int level,
      boolean predictorMissing, boolean classMissing) {
    
    int classIndex = data.classIndex();
    Random random = new Random(1);
    for (int i = 0; i < data.numInstances(); i++) {
      Instance current = data.instance(i);
      for (int j = 0; j < data.numAttributes(); j++) {
        if (((j == classIndex) && classMissing) ||
            ((j != classIndex) && predictorMissing)) {
          if (Math.abs(random.nextInt()) % 100 < level)
            current.setMissing(j);
        }
      }
    }
  }
  
  /**
   * Provides a hook for derived classes to further modify the data. 
   * 
   * @param data	the data to process
   * @return		the processed data
   * @see #m_PostProcessor
   */
  protected Instances process(Instances data) {
    if (getPostProcessor() == null)
      return data;
    else
      return getPostProcessor().process(data);
  }
}
