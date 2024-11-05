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
 * TestInstances.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceAttribute;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.ArrayList;

/**
 * Generates artificial datasets for testing. In case of Multi-Instance data
 * the settings for the number of attributes applies to the data inside
 * the bag. Originally based on code from the CheckClassifier.<p/>
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -relation &lt;name&gt;
 *  The name of the data set.</pre>
 * 
 * <pre> -seed &lt;num&gt;
 *  The seed value.</pre>
 * 
 * <pre> -num-instances &lt;num&gt;
 *  The number of instances in the datasets (default 20).</pre>
 * 
 * <pre> -class-type &lt;num&gt;
 *  The class type, see constants in weka.core.Attribute
 *  (default 1=nominal).</pre>
 * 
 * <pre> -class-values &lt;num&gt;
 *  The number of classes to generate (for nominal classes only)
 *  (default 2).</pre>
 * 
 * <pre> -class-index &lt;num&gt;
 *  The class index, with -1=last, (default -1).</pre>
 * 
 * <pre> -no-class
 *  Doesn't include a class attribute in the output.</pre>
 * 
 * <pre> -nominal &lt;num&gt;
 *  The number of nominal attributes (default 1).</pre>
 * 
 * <pre> -nominal-values &lt;num&gt;
 *  The number of values for nominal attributes (default 2).</pre>
 * 
 * <pre> -numeric &lt;num&gt;
 *  The number of numeric attributes (default 0).</pre>
 * 
 * <pre> -string &lt;num&gt;
 *  The number of string attributes (default 0).</pre>
 * 
 * <pre> -words &lt;comma-separated-list&gt;
 *  The words to use in string attributes.</pre>
 * 
 * <pre> -word-separators &lt;chars&gt;
 *  The word separators to use in string attributes.</pre>
 * 
 * <pre> -date &lt;num&gt;
 *  The number of date attributes (default 0).</pre>
 * 
 * <pre> -relational &lt;num&gt;
 *  The number of relational attributes (default 0).</pre>
 * 
 * <pre> -relational-nominal &lt;num&gt;
 *  The number of nominal attributes in a rel. attribute (default 1).</pre>
 * 
 * <pre> -relational-nominal-values &lt;num&gt;
 *  The number of values for nominal attributes in a rel. attribute (default 2).</pre>
 * 
 * <pre> -relational-numeric &lt;num&gt;
 *  The number of numeric attributes in a rel. attribute (default 0).</pre>
 * 
 * <pre> -relational-string &lt;num&gt;
 *  The number of string attributes in a rel. attribute (default 0).</pre>
 * 
 * <pre> -relational-date &lt;num&gt;
 *  The number of date attributes in a rel. attribute (default 0).</pre>
 * 
 * <pre> -num-instances-relational &lt;num&gt;
 *  The number of instances in relational/bag attributes (default 10).</pre>
 * 
 * <pre> -multi-instance
 *  Generates multi-instance data.</pre>
 * 
 * <pre> -W &lt;classname&gt;
 *  The Capabilities handler to base the dataset on.
 *  The other parameters can be used to override the ones
 *  determined from the handler. Additional parameters for
 *  handler can be passed on after the '--'.</pre>
 * 
 <!-- options-end -->
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6326 $
 * @see weka.classifiers.CheckClassifier
 */
public class TestInstances 
  implements Cloneable, Serializable, OptionHandler, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -6263968936330390469L;

  /** can be used for settting the class attribute index to last 
   * @see #setClassIndex(int) */
  public final static int CLASS_IS_LAST = -1;
  
  /** can be used to avoid generating a class attribute
   * @see #setClassIndex(int) */
  public final static int NO_CLASS = -2;

  /** the default list of words used in strings */
  public final static String[] DEFAULT_WORDS = {"The", "quick", "brown", "fox", "jumps", "over", "the", "lazy", "dog"};

  /** the default word separators used in strings */
  public final static String DEFAULT_SEPARATORS = " ";
  
  /** for generating String attributes/classes */
  protected String[] m_Words = DEFAULT_WORDS;
  
  /** for generating String attributes/classes */
  protected String m_WordSeparators = DEFAULT_SEPARATORS;
  
  /** the name of the relation */
  protected String m_Relation = "Testdata";
  
  /** the seed value */
  protected int m_Seed = 1;
  
  /** the random number generator */
  protected Random m_Random = new Random(m_Seed);
  
  /** the number of instances */
  protected int m_NumInstances = 20;
  
  /** the class type */
  protected int m_ClassType = Attribute.NOMINAL;
  
  /** the number of classes (in case of NOMINAL class) */
  protected int m_NumClasses = 2;

  /** the class index (-1 is LAST, -2 means no class) 
   * @see #CLASS_IS_LAST
   * @see #NO_CLASS */
  protected int m_ClassIndex = CLASS_IS_LAST;
  
  /** the number of nominal attributes */
  protected int m_NumNominal = 1;
  
  /** the number of values for nominal attributes */
  protected int m_NumNominalValues = 2;
  
  /** the number of numeric attributes */
  protected int m_NumNumeric = 0;
  
  /** the number of string attributes */
  protected int m_NumString = 0;
  
  /** the number of date attributes */
  protected int m_NumDate = 0;
  
  /** the number of relational attributes */
  protected int m_NumRelational = 0;
  
  /** the number of nominal attributes in a relational attribute */
  protected int m_NumRelationalNominal = 1;
  
  /** the number of values for nominal attributes in relational attributes */
  protected int m_NumRelationalNominalValues = 2;
  
  /** the number of numeric attributes in a relational attribute */
  protected int m_NumRelationalNumeric = 0;
  
  /** the number of string attributes in a relational attribute */
  protected int m_NumRelationalString = 0;
  
  /** the number of date attributes in a relational attribute */
  protected int m_NumRelationalDate = 0;
  
  /** whether to generate Multi-Instance data or not */
  protected boolean m_MultiInstance = false;
  
  /** the number of instances in relational attributes (applies also for bags
   * in multi-instance) */
  protected int m_NumInstancesRelational = 10;
  
  /** the format of the multi-instance data */
  protected Instances[] m_RelationalFormat = null;
  
  /** the format of the multi-instance data of the class */
  protected Instances m_RelationalClassFormat = null;
  
  /** the generated data */
  protected Instances m_Data = null;
  
  /** the CapabilitiesHandler to get the Capabilities from */
  protected CapabilitiesHandler m_Handler = null;

  /**
   * the default constructor
   */
  public TestInstances() {
    super();
    
    setRelation("Testdata");
    setSeed(1);
    setNumInstances(20);
    setClassType(Attribute.NOMINAL);
    setNumClasses(2);
    setClassIndex(CLASS_IS_LAST);
    setNumNominal(1);
    setNumNominalValues(2);
    setNumNumeric(0);
    setNumString(0);
    setNumDate(0);
    setNumRelational(0);
    setNumRelationalNominal(1);
    setNumRelationalNominalValues(2);
    setNumRelationalNumeric(0);
    setNumRelationalString(0);
    setNumRelationalDate(0);
    setNumInstancesRelational(10);
    setMultiInstance(false);
    setWords(arrayToList(DEFAULT_WORDS));
    setWordSeparators(DEFAULT_SEPARATORS);
  }
  
  /**
   * creates a clone of the current object
   * 
   * @return		a clone of the current object
   */
  public Object clone() {
    TestInstances     result;
    
    result = new TestInstances();
    result.assign(this);
    
    return result;
  }
  
  /**
   * updates itself with all the settings from the given TestInstances
   * object
   * 
   * @param t		the object to get the settings from
   */
  public void assign(TestInstances t) {
    setRelation(t.getRelation());
    setSeed(t.getSeed());
    setNumInstances(t.getNumInstances());
    setClassType(t.getClassType());
    setNumClasses(t.getNumClasses());
    setClassIndex(t.getClassIndex());
    setNumNominal(t.getNumNominal());
    setNumNominalValues(t.getNumNominalValues());
    setNumNumeric(t.getNumNumeric());
    setNumString(t.getNumString());
    setNumDate(t.getNumDate());
    setNumRelational(t.getNumRelational());
    setNumRelationalNominal(t.getNumRelationalNominal());
    setNumRelationalNominalValues(t.getNumRelationalNominalValues());
    setNumRelationalNumeric(t.getNumRelationalNumeric());
    setNumRelationalString(t.getNumRelationalString());
    setNumRelationalDate(t.getNumRelationalDate());
    setMultiInstance(t.getMultiInstance());
    for (int i = 0; i < t.getNumRelational(); i++)
      setRelationalFormat(i, t.getRelationalFormat(i));
    setRelationalClassFormat(t.getRelationalClassFormat());
    setNumInstancesRelational(t.getNumInstancesRelational());
    setWords(t.getWords());
    setWordSeparators(t.getWordSeparators());
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> result = new Vector<Option>();
    
    result.add(new Option(
        "\tThe name of the data set.",
        "relation", 1, "-relation <name>"));
    
    result.add(new Option(
        "\tThe seed value.",
        "seed", 1, "-seed <num>"));
    
    result.add(new Option(
        "\tThe number of instances in the datasets (default 20).",
        "num-instances", 1, "-num-instances <num>"));
    
    result.add(new Option(
        "\tThe class type, see constants in weka.core.Attribute\n"
	+ "\t(default 1=nominal).",
        "class-type", 1, "-class-type <num>"));
    
    result.add(new Option(
        "\tThe number of classes to generate (for nominal classes only)\n"
	+ "\t(default 2).",
        "class-values", 1, "-class-values <num>"));
    
    result.add(new Option(
        "\tThe class index, with -1=last, (default -1).",
        "class-index", 1, "-class-index <num>"));
    
    result.add(new Option(
        "\tDoesn't include a class attribute in the output.",
        "no-class", 0, "-no-class"));

    result.add(new Option(
        "\tThe number of nominal attributes (default 1).",
        "nominal", 1, "-nominal <num>"));
    
    result.add(new Option(
        "\tThe number of values for nominal attributes (default 2).",
        "nominal-values", 1, "-nominal-values <num>"));
    
    result.add(new Option(
        "\tThe number of numeric attributes (default 0).",
        "numeric", 1, "-numeric <num>"));
    
    result.add(new Option(
        "\tThe number of string attributes (default 0).",
        "string", 1, "-string <num>"));
    
    result.add(new Option(
        "\tThe words to use in string attributes.",
        "words", 1, "-words <comma-separated-list>"));
    
    result.add(new Option(
        "\tThe word separators to use in string attributes.",
        "word-separators", 1, "-word-separators <chars>"));
    
    result.add(new Option(
        "\tThe number of date attributes (default 0).",
        "date", 1, "-date <num>"));
    
    result.add(new Option(
        "\tThe number of relational attributes (default 0).",
        "relational", 1, "-relational <num>"));
    
    result.add(new Option(
        "\tThe number of nominal attributes in a rel. attribute (default 1).",
        "relational-nominal", 1, "-relational-nominal <num>"));
    
    result.add(new Option(
        "\tThe number of values for nominal attributes in a rel. attribute (default 2).",
        "relational-nominal-values", 1, "-relational-nominal-values <num>"));
    
    result.add(new Option(
        "\tThe number of numeric attributes in a rel. attribute (default 0).",
        "relational-numeric", 1, "-relational-numeric <num>"));
    
    result.add(new Option(
        "\tThe number of string attributes in a rel. attribute (default 0).",
        "relational-string", 1, "-relational-string <num>"));
    
    result.add(new Option(
        "\tThe number of date attributes in a rel. attribute (default 0).",
        "relational-date", 1, "-relational-date <num>"));
    
    result.add(new Option(
        "\tThe number of instances in relational/bag attributes (default 10).",
        "num-instances-relational", 1, "-num-instances-relational <num>"));
    
    result.add(new Option(
        "\tGenerates multi-instance data.",
        "multi-instance", 0, "-multi-instance"));

    result.add(new Option(
        "\tThe Capabilities handler to base the dataset on.\n"
	+ "\tThe other parameters can be used to override the ones\n"
	+ "\tdetermined from the handler. Additional parameters for\n"
	+ "\thandler can be passed on after the '--'.",
        "W", 1, "-W <classname>"));
    
    return result.elements();
  }
  
  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -relation &lt;name&gt;
   *  The name of the data set.</pre>
   * 
   * <pre> -seed &lt;num&gt;
   *  The seed value.</pre>
   * 
   * <pre> -num-instances &lt;num&gt;
   *  The number of instances in the datasets (default 20).</pre>
   * 
   * <pre> -class-type &lt;num&gt;
   *  The class type, see constants in weka.core.Attribute
   *  (default 1=nominal).</pre>
   * 
   * <pre> -class-values &lt;num&gt;
   *  The number of classes to generate (for nominal classes only)
   *  (default 2).</pre>
   * 
   * <pre> -class-index &lt;num&gt;
   *  The class index, with -1=last, (default -1).</pre>
   * 
   * <pre> -no-class
   *  Doesn't include a class attribute in the output.</pre>
   * 
   * <pre> -nominal &lt;num&gt;
   *  The number of nominal attributes (default 1).</pre>
   * 
   * <pre> -nominal-values &lt;num&gt;
   *  The number of values for nominal attributes (default 2).</pre>
   * 
   * <pre> -numeric &lt;num&gt;
   *  The number of numeric attributes (default 0).</pre>
   * 
   * <pre> -string &lt;num&gt;
   *  The number of string attributes (default 0).</pre>
   * 
   * <pre> -words &lt;comma-separated-list&gt;
   *  The words to use in string attributes.</pre>
   * 
   * <pre> -word-separators &lt;chars&gt;
   *  The word separators to use in string attributes.</pre>
   * 
   * <pre> -date &lt;num&gt;
   *  The number of date attributes (default 0).</pre>
   * 
   * <pre> -relational &lt;num&gt;
   *  The number of relational attributes (default 0).</pre>
   * 
   * <pre> -relational-nominal &lt;num&gt;
   *  The number of nominal attributes in a rel. attribute (default 1).</pre>
   * 
   * <pre> -relational-nominal-values &lt;num&gt;
   *  The number of values for nominal attributes in a rel. attribute (default 2).</pre>
   * 
   * <pre> -relational-numeric &lt;num&gt;
   *  The number of numeric attributes in a rel. attribute (default 0).</pre>
   * 
   * <pre> -relational-string &lt;num&gt;
   *  The number of string attributes in a rel. attribute (default 0).</pre>
   * 
   * <pre> -relational-date &lt;num&gt;
   *  The number of date attributes in a rel. attribute (default 0).</pre>
   * 
   * <pre> -num-instances-relational &lt;num&gt;
   *  The number of instances in relational/bag attributes (default 10).</pre>
   * 
   * <pre> -multi-instance
   *  Generates multi-instance data.</pre>
   * 
   * <pre> -W &lt;classname&gt;
   *  The Capabilities handler to base the dataset on.
   *  The other parameters can be used to override the ones
   *  determined from the handler. Additional parameters for
   *  handler can be passed on after the '--'.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      		tmpStr;
    Class			cls;
    CapabilitiesHandler		handler;
    boolean			initialized;
    
    initialized = false;

    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() > 0) {
      cls = Class.forName(tmpStr);
      if (ClassDiscovery.hasInterface(CapabilitiesHandler.class, cls)) {
	initialized = true;
	handler = (CapabilitiesHandler) cls.newInstance();
	if (handler instanceof OptionHandler)
	  ((OptionHandler) handler).setOptions(Utils.partitionOptions(options));
	setHandler(handler);
	// initialize
	this.assign(forCapabilities(handler.getCapabilities()));
      }
      else {
	throw new IllegalArgumentException("Class '" + tmpStr + "' is not a CapabilitiesHandler!");
      }
    }

    tmpStr = Utils.getOption("relation", options);
    if (tmpStr.length() != 0)
      setRelation(tmpStr);
    else if (!initialized)
      setRelation("Testdata");
    
    tmpStr = Utils.getOption("seed", options);
    if (tmpStr.length() != 0)
      setSeed(Integer.parseInt(tmpStr));
    else if (!initialized)
      setSeed(1);
    
    tmpStr = Utils.getOption("num-instances", options);
    if (tmpStr.length() != 0)
      setNumInstances(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumInstances(20);
    
    setNoClass(Utils.getFlag("no-class", options));
    
    if (!getNoClass()) {
      tmpStr = Utils.getOption("class-type", options);
      if (tmpStr.length() != 0)
        setClassType(Integer.parseInt(tmpStr));
      else if (!initialized)
        setClassType(Attribute.NOMINAL);

      tmpStr = Utils.getOption("class-values", options);
      if (tmpStr.length() != 0)
        setNumClasses(Integer.parseInt(tmpStr));
      else if (!initialized)
        setNumClasses(2);

      tmpStr = Utils.getOption("class-index", options);
      if (tmpStr.length() != 0)
        setClassIndex(Integer.parseInt(tmpStr));
      else if (!initialized)
        setClassIndex(-1);
    }
    
    tmpStr = Utils.getOption("nominal", options);
    if (tmpStr.length() != 0)
      setNumNominal(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumNominal(1);
    
    tmpStr = Utils.getOption("nominal-values", options);
    if (tmpStr.length() != 0)
      setNumNominalValues(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumNominalValues(2);
    
    tmpStr = Utils.getOption("numeric", options);
    if (tmpStr.length() != 0)
      setNumNumeric(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumNumeric(0);
    
    tmpStr = Utils.getOption("string", options);
    if (tmpStr.length() != 0)
      setNumString(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumString(0);
    
    tmpStr = Utils.getOption("words", options);
    if (tmpStr.length() != 0)
      setWords(tmpStr);
    else if (!initialized)
      setWords(arrayToList(DEFAULT_WORDS));
    
    if (Utils.getOptionPos("word-separators", options) > -1) {
      tmpStr = Utils.getOption("word-separators", options);
      setWordSeparators(tmpStr);
    }
    else if (!initialized) {
      setWordSeparators(DEFAULT_SEPARATORS);
    }
    
    tmpStr = Utils.getOption("date", options);
    if (tmpStr.length() != 0)
      setNumDate(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumDate(0);
    
    tmpStr = Utils.getOption("relational", options);
    if (tmpStr.length() != 0)
      setNumRelational(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumRelational(0);
    
    tmpStr = Utils.getOption("relational-nominal", options);
    if (tmpStr.length() != 0)
      setNumRelationalNominal(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumRelationalNominal(1);
    
    tmpStr = Utils.getOption("relational-nominal-values", options);
    if (tmpStr.length() != 0)
      setNumRelationalNominalValues(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumRelationalNominalValues(2);
    
    tmpStr = Utils.getOption("relational-numeric", options);
    if (tmpStr.length() != 0)
      setNumRelationalNumeric(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumRelationalNumeric(0);
    
    tmpStr = Utils.getOption("relational-string", options);
    if (tmpStr.length() != 0)
      setNumRelationalString(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumRelationalString(0);
    
    tmpStr = Utils.getOption("num-instances-relational", options);
    if (tmpStr.length() != 0)
      setNumInstancesRelational(Integer.parseInt(tmpStr));
    else if (!initialized)
      setNumInstancesRelational(10);
    
    if (!initialized)
      setMultiInstance(Utils.getFlag("multi-instance", options));
  }
  
  /**
   * Gets the current settings of this object.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String> 	result;
    String[]	options;
    int		i;

    result = new Vector<String>();
    
    result.add("-relation");
    result.add(getRelation());
    
    result.add("-seed");
    result.add("" + getSeed());
    
    result.add("-num-instances");
    result.add("" + getNumInstances());
    
    if (getNoClass()) {
      result.add("-no-class");
    }
    else {
      result.add("-class-type");
      result.add("" + getClassType());
      
      result.add("-class-values");
      result.add("" + getNumClasses());
      
      result.add("-class-index");
      result.add("" + getClassIndex());
    }
    
    result.add("-nominal");
    result.add("" + getNumNominal());
    
    result.add("-nominal-values");
    result.add("" + getNumNominalValues());
    
    result.add("-numeric");
    result.add("" + getNumNumeric());
    
    result.add("-string");
    result.add("" + getNumString());
    
    result.add("-words");
    result.add("" + getWords());
    
    result.add("-word-separators");
    result.add("" + getWordSeparators());
    
    result.add("-date");
    result.add("" + getNumDate());
    
    result.add("-relational");
    result.add("" + getNumRelational());
    
    result.add("-relational-nominal");
    result.add("" + getNumRelationalNominal());
    
    result.add("-relational-nominal-values");
    result.add("" + getNumRelationalNominalValues());
    
    result.add("-relational-numeric");
    result.add("" + getNumRelationalNumeric());
    
    result.add("-relational-string");
    result.add("" + getNumRelationalString());
    
    result.add("-relational-date");
    result.add("" + getNumRelationalDate());
    
    result.add("-num-instances-relational");
    result.add("" + getNumInstancesRelational());

    if (getMultiInstance())
      result.add("-multi-instance");

    if (getHandler() != null) {
      result.add("-W");
      result.add(getHandler().getClass().getName());
      if (getHandler() instanceof OptionHandler) {
	result.add("--");
	options = ((OptionHandler) getHandler()).getOptions();
	for (i = 0; i < options.length; i++)
	  result.add(options[i]);
      }
    }
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * sets the name of the relation
   * 
   * @param value 	the name of the relation
   */
  public void setRelation(String value) {
    m_Relation = value;
  }
  
  /**
   * returns the current name of the relation
   * 
   * @return 		the name of the relation
   */
  public String getRelation() {
    return m_Relation;
  }
  
  /**
   * sets the seed value for the random number generator
   * 
   * @param value	the seed
   */
  public void setSeed(int value) {
    m_Seed   = value;
    m_Random = new Random(m_Seed);
  }
  
  /**
   * returns the current seed value
   * 
   * @return 		the seed value
   */
  public int getSeed() {
    return m_Seed;
  }
  
  /**
   * sets the number of instances to produce
   * 
   * @param value	the number of instances
   */
  public void setNumInstances(int value) {
    m_NumInstances = value;
  }
  
  /**
   * returns the current number of instances to produce
   * 
   * @return the number of instances
   */
  public int getNumInstances() {
    return m_NumInstances;
  }
  
  /**
   * sets the class attribute type
   * 
   * @param value	the class attribute type
   */
  public void setClassType(int value) {
    m_ClassType = value;
    m_RelationalClassFormat = null;
  }
  
  /**
   * returns the current class type
   * 
   * @return the class attribute type
   */
  public int getClassType() {
    return m_ClassType;
  }
  
  /**
   * sets the number of classes
   * 
   * @param value	the number of classes
   */
  public void setNumClasses(int value) {
    m_NumClasses = value;
  }
  
  /**
   * returns the current number of classes
   * 
   * @return the number of classes
   */
  public int getNumClasses() {
    return m_NumClasses;
  }
  
  /**
   * sets the class index (0-based)
   * 
   * @param value	the class index 
   * @see #CLASS_IS_LAST
   * @see #NO_CLASS
   */
  public void setClassIndex(int value) {
    m_ClassIndex = value;
  }
  
  /**
   * returns the current class index (0-based), -1 is last attribute
   * 
   * @return 		the class index
   * @see #CLASS_IS_LAST
   * @see #NO_CLASS
   */
  public int getClassIndex() {
    return m_ClassIndex;
  }
  
  /**
   * whether to have no class, e.g., for clusterers; otherwise the class
   * attribute index is set to last
   * 
   * @param value whether to have no class
   * @see #CLASS_IS_LAST
   * @see #NO_CLASS
   */
  public void setNoClass(boolean value) {
    if (value)
      setClassIndex(NO_CLASS);
    else
      setClassIndex(CLASS_IS_LAST);
  }
  
  /**
   * whether no class attribute is generated
   * 
   * @return 		true if no class attribute is generated
   */
  public boolean getNoClass() {
    return (getClassIndex() == NO_CLASS);
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
   * sets the number of values for nominal attributes
   * 
   * @param value	the number of values
   */
  public void setNumNominalValues(int value) {
    m_NumNominalValues = value;
  }
  
  /**
   * returns the current number of values for nominal attributes
   * 
   * @return 		the number of values
   */
  public int getNumNominalValues() {
    return m_NumNominalValues;
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
   * sets the number of date attributes
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
    m_RelationalFormat = new Instances[value];
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
   * sets the number of nominal attributes in a relational attribute
   * 
   * @param value	the number of nominal attributes
   */
  public void setNumRelationalNominal(int value) {
    m_NumRelationalNominal = value;
  }
  
  /**
   * returns the current number of nominal attributes in a relational attribute
   * 
   * @return 		the number of nominal attributes
   */
  public int getNumRelationalNominal() {
    return m_NumRelationalNominal;
  }
  
  /**
   * sets the number of values for nominal attributes in a relational attribute
   * 
   * @param value	the number of values
   */
  public void setNumRelationalNominalValues(int value) {
    m_NumRelationalNominalValues = value;
  }
  
  /**
   * returns the current number of values for nominal attributes in a relational attribute
   * 
   * @return 		the number of values
   */
  public int getNumRelationalNominalValues() {
    return m_NumRelationalNominalValues;
  }
  
  /**
   * sets the number of numeric attributes in a relational attribute
   * 
   * @param value 	the number of numeric attributes
   */
  public void setNumRelationalNumeric(int value) {
    m_NumRelationalNumeric = value;
  }
  
  /**
   * returns the current number of numeric attributes in a relational attribute
   * 
   * @return 		the number of numeric attributes
   */
  public int getNumRelationalNumeric() {
    return m_NumRelationalNumeric;
  }
  
  /**
   * sets the number of string attributes in a relational attribute
   * 
   * @param value 	the number of string attributes
   */
  public void setNumRelationalString(int value) {
    m_NumRelationalString = value;
  }
  
  /**
   * returns the current number of string attributes in a relational attribute
   * 
   * @return 		the number of string attributes
   */
  public int getNumRelationalString() {
    return m_NumRelationalString;
  }
  
  /**
   * sets the number of date attributes in a relational attribute
   * 
   * @param value	the number of date attributes
   */
  public void setNumRelationalDate(int value) {
    m_NumRelationalDate = value;
  }
  
  /**
   * returns the current number of date attributes in a relational attribute
   * 
   * @return		the number of date attributes
   */
  public int getNumRelationalDate() {
    return m_NumRelationalDate;
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
   * sets whether multi-instance data should be generated (with a fixed
   * data structure)
   * 
   * @param value	whether multi-instance data is generated
   */
  public void setMultiInstance(boolean value) {
    m_MultiInstance = value;
  }
  
  /**
   * Gets whether multi-instance data (with a fixed structure) is generated
   * 
   * @return		true if multi-instance data is generated
   */
  public boolean getMultiInstance() {
    return m_MultiInstance;
  }
  
  /**
   * sets the structure for the bags for the relational attribute
   * 
   * @param index       the index of the relational attribute
   * @param value       the new structure
   */
  public void setRelationalFormat(int index, Instances value) {
    if (value != null)
      m_RelationalFormat[index] = new Instances(value, 0);
    else
      m_RelationalFormat[index] = null;
  }
  
  /**
   * returns the format for the specified relational attribute, can be null
   * 
   * @param index       the index of the relational attribute
   * @return            the current structure
   */
  public Instances getRelationalFormat(int index) {
    return m_RelationalFormat[index];
  }

  /**
   * sets the structure for the relational class attribute
   * 
   * @param value	the structure for the relational attribute
   */
  public void setRelationalClassFormat(Instances value) {
    if (value != null)
      m_RelationalClassFormat = new Instances(value, 0);
    else
      m_RelationalClassFormat = null;
  }
  
  /**
   * returns the current strcuture of the relational class attribute, can
   * be null
   * 
   * @return 		the relational structure of the class attribute
   */
  public Instances getRelationalClassFormat() {
    return m_RelationalClassFormat;
  }
  
  /**
   * returns the overall number of attributes (incl. class, if that is also
   * generated)
   * 
   * @return 		the overall number of attributes
   */
  public int getNumAttributes() {
    int		result;
    
    result = m_NumNominal + m_NumNumeric + m_NumString + m_NumDate + m_NumRelational;
    
    if (!getNoClass())
      result++;
      
    return result;
  }
  
  /**
   * returns the current dataset, can be null
   * 
   * @return		the current dataset
   */
  public Instances getData() {
    return m_Data;
  }
  
  /**
   * sets the Capabilities handler to generate the data for
   * 
   * @param value	the handler to generate the data for
   */
  public void setHandler(CapabilitiesHandler value) {
    m_Handler = value;
  }
  
  /**
   * returns the current set CapabilitiesHandler to generate the dataset
   * for, can be null
   * 
   * @return		the handler to generate the data for
   */
  public CapabilitiesHandler getHandler() {
    return m_Handler;
  }
  
  /**
   * creates a new Attribute of the given type
   * 
   * @param index the index of the current attribute (0-based)
   * @param attType the attribute type (NUMERIC, NOMINAL, etc.)
   * @return the configured attribute
   * @throws Exception if something goes wrong, e.g., an unknown attribute type
   * 
   * @see Attribute#type()
   * @see #CLASS_IS_LAST
   * @see #NO_CLASS
   */
  protected Attribute generateAttribute(int index, int attType, 
      String namePrefix) throws Exception {
    Attribute     result;
    String        name;
    int           valIndex;
    int           nomCount;
    String        prefix;
    
    result = null;

    // determine name and start-index
    if (index == CLASS_IS_LAST) {
      valIndex = 0;
      name     = "Class";
      prefix   = "class";
      nomCount = getNumClasses();
    }
    else {
      valIndex = index;
      nomCount = getNumNominalValues();
      prefix   = "att" + (valIndex + 1) + "val";
      
      switch (attType) {
        case Attribute.NOMINAL:
          name = "Nominal" + (valIndex + 1);
          break;
          
        case PreferenceAttribute.RANKING:
          name = "Ranking" + (valIndex + 1);
          break;
          
        case Attribute.NUMERIC:
          name = "Numeric" + (valIndex + 1);
          break;
          
        case Attribute.STRING:
          name = "String" + (valIndex + 1);
          break;
          
        case Attribute.DATE:
          name = "Date" + (valIndex + 1);
          break;
          
        case Attribute.RELATIONAL:
          name = "Relational" + (valIndex + 1);
          break;
          
        default:
          throw new IllegalArgumentException("Attribute type '" + attType + "' unknown!");
      }
    }
    
    switch (attType) {
      case Attribute.NOMINAL:
        ArrayList<String> nomStrings = new ArrayList<String>(valIndex + 1);
        for (int j = 0; j < nomCount; j++)
          nomStrings.add(prefix + (j + 1));
        result = new Attribute(namePrefix + name, nomStrings);
        break;
        
      case PreferenceAttribute.RANKING:
          nomStrings = new ArrayList<String>(valIndex + 1);
          for (int j = 0; j < nomCount; j++)
            nomStrings.add(prefix + (j + 1));
          result = new PreferenceAttribute(namePrefix + name, nomStrings);
          break;
        
      case Attribute.NUMERIC:
        result = new Attribute(namePrefix + name);
        break;
        
      case Attribute.STRING:
        result = new Attribute(namePrefix + name, (ArrayList<String>) null);
        break;
        
      case Attribute.DATE:
        result = new Attribute(namePrefix + name, "yyyy-mm-dd");
        break;
        
      case Attribute.RELATIONAL:
        Instances rel;
        if (index == CLASS_IS_LAST)
          rel = getRelationalClassFormat();
        else
          rel = getRelationalFormat(index);
        
        if (rel == null) {
          TestInstances dataset = new TestInstances();
          dataset.setNumNominal(getNumRelationalNominal());
          dataset.setNumNominalValues(getNumRelationalNominalValues());
          dataset.setNumNumeric(getNumRelationalNumeric());
          dataset.setNumString(getNumRelationalString());
          dataset.setNumDate(getNumRelationalDate());
          dataset.setNumInstances(0);
          dataset.setClassType(Attribute.NOMINAL);  // dummy to avoid endless recursion, will be deleted anyway
          rel = new Instances(dataset.generate());
          if (!getNoClass()) {
            int clsIndex = rel.classIndex();
            rel.setClassIndex(-1);
            rel.deleteAttributeAt(clsIndex);
          }
        }
        result = new Attribute(namePrefix + name, rel);
        break;
        
      default:
        throw new IllegalArgumentException("Attribute type '" + attType + "' unknown!");
    }
    
    return result;
  }
  
  /**
   * Generates the class value
   * 
   * @param data  	the dataset to work on
   * @return      	the new class value
   * @throws Exception 	if something goes wrong
   */
  protected double generateClassValue(Instances data) throws Exception {
    double result = Double.NaN;
    
    switch (m_ClassType) {
      case Attribute.NUMERIC:
        result = m_Random.nextFloat() * 0.25
            + Math.abs(m_Random.nextInt())
            % Math.max(2, m_NumNominal);
        break;
        
      case Attribute.NOMINAL:
      case PreferenceAttribute.RANKING:
        result = Math.abs(m_Random.nextInt()) % data.numClasses();
        break;
        
      case Attribute.STRING:
        String str = "";
        for (int n = 0; n < m_Words.length; n++) {
          if ( (n > 0) && (m_WordSeparators.length() != 0) )
            str += m_WordSeparators.charAt(m_Random.nextInt(m_WordSeparators.length()));
          str += m_Words[m_Random.nextInt(m_Words.length)];
        }
        result = data.classAttribute().addStringValue(str);
        break;
        
      case Attribute.DATE:
        result = data.classAttribute().parseDate(
                (2000 + m_Random.nextInt(100)) + "-01-01");
        break;
        
      case Attribute.RELATIONAL:
        if (getRelationalClassFormat() != null) {
          result = data.classAttribute().addRelation(getRelationalClassFormat());
        }
        else {
          TestInstances dataset = new TestInstances();
          dataset.setNumNominal(getNumRelationalNominal());
          dataset.setNumNominalValues(getNumRelationalNominalValues());
          dataset.setNumNumeric(getNumRelationalNumeric());
          dataset.setNumString(getNumRelationalString());
          dataset.setNumDate(getNumRelationalDate());
          dataset.setNumInstances(getNumInstancesRelational());
          dataset.setClassType(Attribute.NOMINAL);  // dummy to avoid endless recursion, will be deleted anyway
          Instances rel = new Instances(dataset.generate());
          int clsIndex = rel.classIndex();
          rel.setClassIndex(-1);
          rel.deleteAttributeAt(clsIndex);
          result = data.classAttribute().addRelation(rel);
        }
        break;
    }
    
    return result;
  }
  
  /**
   * Generates a new value for the specified attribute. The classValue
   * might be used in the process.
   * 
   * @param data          the dataset to work on
   * @param index         the index of the attribute
   * @param classVal      the class value for the current instance, might be 
   *                      used in the calculation
   * @return              the new attribute value
   * @throws Exception    if something goes wrong
   */
  protected double generateAttributeValue(Instances data, int index, double classVal) throws Exception {
    double result = Double.NaN;
    
    switch (data.attribute(index).type()) {
      case Attribute.NUMERIC:
        result = classVal * 4 + m_Random.nextFloat() * 1 - 0.5;
        break;
        
      case Attribute.NOMINAL:
      case PreferenceAttribute.RANKING:
        if (m_Random.nextFloat() < 0.2) {
          result = Math.abs(m_Random.nextInt())
          % data.attribute(index).numValues();
        } else {
          result = ((int)classVal) % data.attribute(index).numValues();
        }
	//result = m_Random.nextInt(data.attribute(index).numValues());
        break;
        
      case Attribute.STRING:
        String str = "";
        for (int n = 0; n < m_Words.length; n++) {
          if ( (n > 0) && (m_WordSeparators.length() != 0) )
            str += m_WordSeparators.charAt(m_Random.nextInt(m_WordSeparators.length()));
          str += m_Words[m_Random.nextInt(m_Words.length)];
        }
        result = data.attribute(index).addStringValue(str);
        break;
        
      case Attribute.DATE:
        result = data.attribute(index).parseDate(
                (2000 + m_Random.nextInt(100)) + "-01-01");
        break;
        
      case Attribute.RELATIONAL:
        Instances rel = new Instances(data.attribute(index).relation(), 0);
        for (int n = 0; n < getNumInstancesRelational(); n++) {
          Instance inst = new DenseInstance(rel.numAttributes());
          inst.setDataset(data);
          for (int i = 0; i < rel.numAttributes(); i++) {
            inst.setValue(i, generateAttributeValue(rel, i, 0));
          }
          rel.add(inst);
        }
        result = data.attribute(index).addRelation(rel);
        break;
    }
    
    return result;
  }
  
  /**
   * Generates a new dataset
   * 
   * @return the generated data
   * @throws Exception if something goes wrong
   */
  public Instances generate() throws Exception {
    return generate("");
  }
  
  /**
   * generates a new dataset.
   *
   * @param namePrefix the prefix to add to the name of an attribute
   * @return 		the generated data
   * @throws Exception	if something goes wrong
   */
  public Instances generate(String namePrefix) throws Exception {
    if (getMultiInstance()) {
      TestInstances bag = (TestInstances) this.clone();
      bag.setMultiInstance(false);
      bag.setNumInstances(0);
      bag.setSeed(m_Random.nextInt());
      Instances bagFormat = bag.generate("bagAtt_");
      bagFormat.setClassIndex(-1);
      bagFormat.deleteAttributeAt(bagFormat.numAttributes() - 1);

      // generate multi-instance structure
      TestInstances structure = new TestInstances();
      structure.setSeed(m_Random.nextInt());
      structure.setNumNominal(1);
      structure.setNumRelational(1);
      structure.setRelationalFormat(0, bagFormat);
      structure.setClassType(getClassType());
      structure.setNumClasses(getNumClasses());
      structure.setRelationalClassFormat(getRelationalClassFormat());
      structure.setNumInstances(getNumInstances());
      m_Data = structure.generate();
      
      // generate bags
      bag.setNumInstances(getNumInstancesRelational());
      for (int i = 0; i < getNumInstances(); i++) {
        bag.setSeed(m_Random.nextInt());
        Instances bagData = new Instances(bag.generate("bagAtt_"));
        bagData.setClassIndex(-1);
        bagData.deleteAttributeAt(bagData.numAttributes() - 1);
        double val = m_Data.attribute(1).addRelation(bagData);
        m_Data.instance(i).setValue(1, val);
      }
    }
    else {
      // initialize
      int clsIndex = m_ClassIndex;
      if (clsIndex == CLASS_IS_LAST)
        clsIndex = getNumAttributes() - 1;

      // generate attributes
      ArrayList<Attribute> attributes = new ArrayList<Attribute>(getNumAttributes());
      // Add Nominal attributes
      for (int i = 0; i < getNumNominal(); i++)
        attributes.add(generateAttribute(i, Attribute.NOMINAL, namePrefix));
      
      // Add m_Numeric attributes
      for (int i = 0; i < getNumNumeric(); i++)
        attributes.add(generateAttribute(i, Attribute.NUMERIC, namePrefix));
      
      // Add some String attributes...
      for (int i = 0; i < getNumString(); i++)
        attributes.add(generateAttribute(i, Attribute.STRING, namePrefix));
      
      // Add some Date attributes...
      for (int i = 0; i < getNumDate(); i++)
        attributes.add(generateAttribute(i, Attribute.DATE, namePrefix));
      
      // Add some Relational attributes...
      for (int i = 0; i < getNumRelational(); i++)
        attributes.add(generateAttribute(i, Attribute.RELATIONAL, namePrefix));
      
      // Add class attribute
      if (clsIndex != NO_CLASS)
	attributes.add(clsIndex, generateAttribute(CLASS_IS_LAST, getClassType(), namePrefix));
      
      m_Data = new Instances(getRelation(), attributes, getNumInstances());
      m_Data.setClassIndex(clsIndex);

      // generate instances
      for (int i = 0; i < getNumInstances(); i++) {
        Instance current = new DenseInstance(getNumAttributes());
        current.setDataset(m_Data);

        // class
        double classVal;
        if (clsIndex != NO_CLASS) {
          classVal = generateClassValue(m_Data);
          current.setClassValue(classVal);
        }
        else {
          classVal = m_Random.nextFloat();
        }
        
        // other attributes
        for (int n = 0; n < getNumAttributes(); n++) {
          if (clsIndex == n)
            continue;
          
          current.setValue(n, generateAttributeValue(m_Data, n, classVal));
        }
        
        m_Data.add(current);
      }
    }

    if (m_Data.classIndex() == NO_CLASS)
      m_Data.setClassIndex(-1);
    
    return getData();
  }
  
  /**
   * returns a TestInstances instance setup already for the the given
   * capabilities.
   * 
   * @param c		the capabilities to base the TestInstances on
   * @return		the configured TestInstances object
   */
  public static TestInstances forCapabilities(Capabilities c) {
    TestInstances	result;
    
    result = new TestInstances();
    
    // multi-instance?
    if (c.getOwner() instanceof MultiInstanceCapabilitiesHandler) {
      Capabilities multi = (Capabilities) ((MultiInstanceCapabilitiesHandler) c.getOwner()).getMultiInstanceCapabilities().clone();
      multi.setOwner(null);  // otherwise recursive!
      result = forCapabilities(multi);
      result.setMultiInstance(true);
    }
    else  {
      // class
      if (c.handles(Capability.NO_CLASS))
	result.setClassIndex(NO_CLASS);
      //RANKING BEGIN
      else if (c.handles(Capability.RANKING))
    result.setClassType(PreferenceAttribute.RANKING);
      //RANKING END
      else if (c.handles(Capability.NOMINAL_CLASS))
	result.setClassType(Attribute.NOMINAL);
      else if (c.handles(Capability.BINARY_CLASS))
	result.setClassType(Attribute.NOMINAL);
      else if (c.handles(Capability.NUMERIC_CLASS))
	result.setClassType(Attribute.NUMERIC);
      else if (c.handles(Capability.DATE_CLASS))
	result.setClassType(Attribute.DATE);
      else if (c.handles(Capability.STRING_CLASS))
	result.setClassType(Attribute.STRING);
      else if (c.handles(Capability.RELATIONAL_CLASS))
	result.setClassType(Attribute.RELATIONAL);

      // # of classes
      if (c.handles(Capability.UNARY_CLASS))
	result.setNumClasses(1);
      if (c.handles(Capability.BINARY_CLASS))
	result.setNumClasses(2);
      if (c.handles(Capability.NOMINAL_CLASS))
	result.setNumClasses(4);
      
      // attributes
      if (c.handles(Capability.NOMINAL_ATTRIBUTES)) {
	result.setNumNominal(1);
	result.setNumRelationalNominal(1);
      }
      else {
	result.setNumNominal(0);
	result.setNumRelationalNominal(0);
      }

      if (c.handles(Capability.NUMERIC_ATTRIBUTES)) {
	result.setNumNumeric(1);
	result.setNumRelationalNumeric(1);
      }
      else {
	result.setNumNumeric(0);
	result.setNumRelationalNumeric(0);
      }

      if (c.handles(Capability.DATE_ATTRIBUTES)) {
	result.setNumDate(1);
	result.setNumRelationalDate(1);
      }
      else {
	result.setNumDate(0);
	result.setNumRelationalDate(0);
      }
      
      if (c.handles(Capability.STRING_ATTRIBUTES)) {
	result.setNumString(1);
	result.setNumRelationalString(1);
      }
      else {
	result.setNumString(0);
	result.setNumRelationalString(0);
      }
      
      if (c.handles(Capability.RELATIONAL_ATTRIBUTES))
	result.setNumRelational(1);
      else
	result.setNumRelational(0);
    }
    
    return result;
  }
  
  /**
   * returns a string representation of the object
   * 
   * @return		a string representation of the object
   */
  public String toString() {
    String	result;
    
    result = "";
    result += "Relation: " + getRelation() + "\n";
    result += "Seed: " + getSeed() + "\n";
    result += "# Instances: " + getNumInstances() + "\n";
    result += "ClassType: " + getClassType() + "\n";
    result += "# Classes: " + getNumClasses() + "\n";
    result += "Class index: " + getClassIndex() + "\n";
    result += "# Nominal: " +     getNumNominal() + "\n";
    result += "# Nominal values: " + getNumNominalValues() + "\n";
    result += "# Numeric: " + getNumNumeric() + "\n";
    result += "# String: " + getNumString() + "\n";
    result += "# Date: " + getNumDate() + "\n";
    result += "# Relational: " + getNumRelational() + "\n";
    result += "  - # Nominal: " +     getNumRelationalNominal() + "\n";
    result += "  - # Nominal values: " + getNumRelationalNominalValues() + "\n";
    result += "  - # Numeric: " + getNumRelationalNumeric() + "\n";
    result += "  - # String: " + getNumRelationalString() + "\n";
    result += "  - # Date: " + getNumRelationalDate() + "\n";
    result += "  - # Instances: " + getNumInstancesRelational() + "\n";
    result += "Multi-Instance: " + getMultiInstance() + "\n";
    result += "Words: " + getWords() + "\n";
    result += "Word separators: " + getWordSeparators() + "\n";
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6326 $");
  }
  
  /**
   * for running the class from commandline, prints the generated data
   * to stdout
   * 
   * @param args	the commandline parameters
   * @throws Exception	if something goes wrong
   */
  public static void main(String[] args) throws Exception {
    TestInstances inst;
    
    inst = new TestInstances();

    // help requested?
    if (Utils.getFlag("h", args) || Utils.getFlag("help", args)) {
      StringBuffer result = new StringBuffer();
      result.append("\nTest data generator options:\n\n");

      result.append("-h|-help\n\tprints this help\n");
      
      Enumeration enm = inst.listOptions();
      while (enm.hasMoreElements()) {
        Option option = (Option) enm.nextElement();
        result.append(option.synopsis() + "\n" + option.description() + "\n");
      }

      System.out.println(result);
      System.exit(0);
    }
    
    // generate data
    inst.setOptions(args);
    System.out.println(inst.generate());
  }
}
