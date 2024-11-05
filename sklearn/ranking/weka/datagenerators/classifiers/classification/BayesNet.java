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
 * BayesNet.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.datagenerators.classifiers.classification;

import weka.classifiers.bayes.net.BayesNetGenerator;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.datagenerators.ClassificationGenerator;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Generates random instances based on a Bayes network.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -h
 *  Prints this help.</pre>
 * 
 * <pre> -o &lt;file&gt;
 *  The name of the output file, otherwise the generated data is
 *  printed to stdout.</pre>
 * 
 * <pre> -r &lt;name&gt;
 *  The name of the relation.</pre>
 * 
 * <pre> -d
 *  Whether to print debug informations.</pre>
 * 
 * <pre> -S
 *  The seed for random function (default 1)</pre>
 * 
 * <pre> -n &lt;num&gt;
 *  The number of examples to generate (default 100)</pre>
 * 
 * <pre> -A &lt;num&gt;
 *  The number of arcs to use. (default 20)</pre>
 * 
 * <pre> -C &lt;num&gt;
 *  The cardinality of the attributes and the class. (default 2)</pre>
 * 
 <!-- options-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 * @see BayesNetGenerator
 */

public class BayesNet
  extends ClassificationGenerator {
  
  /** for serialization */
  static final long serialVersionUID = -796118162379901512L;
  
  /** the bayesian net generator, that produces the actual data */
  protected BayesNetGenerator m_Generator;

  /**
   * initializes the generator
   */
  public BayesNet() {
    super();

    setNumAttributes(defaultNumAttributes());
    setNumArcs(defaultNumArcs());
    setCardinality(defaultCardinality());
  }
  
  /**
   * Returns a string describing this data generator.
   *
   * @return a description of the data generator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Generates random instances based on a Bayes network.";
  }

 /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {
    Vector result = enumToVector(super.listOptions());

    result.add(new Option(
              "\tThe number of arcs to use. (default " 
              + defaultNumArcs() + ")",
              "A", 1, "-A <num>"));

    result.add(new Option(
              "\tThe cardinality of the attributes and the class. (default " 
              + defaultCardinality() + ")",
              "C", 1, "-C <num>"));

    return result.elements();
  }

  /**
   * Parses a list of options for this object. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -h
   *  Prints this help.</pre>
   * 
   * <pre> -o &lt;file&gt;
   *  The name of the output file, otherwise the generated data is
   *  printed to stdout.</pre>
   * 
   * <pre> -r &lt;name&gt;
   *  The name of the relation.</pre>
   * 
   * <pre> -d
   *  Whether to print debug informations.</pre>
   * 
   * <pre> -S
   *  The seed for random function (default 1)</pre>
   * 
   * <pre> -n &lt;num&gt;
   *  The number of examples to generate (default 100)</pre>
   * 
   * <pre> -A &lt;num&gt;
   *  The number of arcs to use. (default 20)</pre>
   * 
   * <pre> -C &lt;num&gt;
   *  The cardinality of the attributes and the class. (default 2)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String        tmpStr;
    Vector        list;

    super.setOptions(options);

    list = new Vector();

    list.add("-N");
    list.add("" + getNumAttributes());

    list.add("-M");
    list.add("" + getNumExamples());
    
    list.add("-S");
    list.add("" + getSeed());
    
    list.add("-A");
    tmpStr = Utils.getOption('A', options);
    if (tmpStr.length() != 0)
      list.add(tmpStr);
    else
      list.add("" + defaultNumArcs());

    list.add("-C");
    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0)
      list.add(tmpStr);
    else
      list.add("" + defaultCardinality());

    setGeneratorOptions(list);
  }

  /**
   * Gets the current settings of the datagenerator.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector        result;
    String[]      options;
    int           i;
    
    result  = new Vector();
    options = removeBlacklist(super.getOptions());
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    // determine options from generator
    options = getGenerator().getOptions();

    try {
      result.add("-A");
      result.add(Utils.getOption('A', options));
    }
    catch (Exception e) {
      e.printStackTrace();
    }

    try {
      result.add("-C");
      result.add(Utils.getOption('C', options));
    }
    catch (Exception e) {
      e.printStackTrace();
    }
    
    return (String[]) result.toArray(new String[result.size()]);
  }

  /**
   * sets the given options of the BayesNetGenerator
   * 
   * @param generator the generator to set the options for
   * @param options the options to set
   */
  protected void setGeneratorOptions(
      BayesNetGenerator generator, Vector options) {

    try {
      generator.setOptions(
          (String[]) options.toArray(new String[options.size()]));
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * returns the actual datagenerator
   * 
   * @return the actual datagenerator
   */
  protected BayesNetGenerator getGenerator() {
    if (m_Generator == null)
      m_Generator = new BayesNetGenerator();

    return m_Generator;
  }

  /**
   * sets the given options of the BayesNetGenerator
   * 
   * @param options the options to set
   */
  protected void setGeneratorOptions(Vector options) {
    setGeneratorOptions(getGenerator(), options);
  }

  /**
   * sets a specific option/value of the generator (option must be w/o
   * then '-')
   * @param generator       the generator to set the option for
   * @param option          the option to set
   * @param value           the new value for the option
   */
  protected void setGeneratorOption( BayesNetGenerator generator, 
                                     String option, String value ) {

    String[]      options;
    Vector        list;
    int           i;

    try {
      // get options and remove specific option
      options = generator.getOptions();
      Utils.getOption(option, options);

      // add option and set the new options
      list = new Vector();
      for (i = 0; i < options.length; i++) {
        if (options[i].length() != 0)
          list.add(options[i]);
      }
      list.add("-" + option);
      list.add(value);
      setGeneratorOptions(generator, list);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * sets a specific option/value of the generator (option must be w/o
   * then '-')
   * @param option          the option to set
   * @param value           the new value for the option
   */
  protected void setGeneratorOption(String option, String value) {
    setGeneratorOption(getGenerator(), option, value);
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
    setGeneratorOption("N", "" + numAttributes);
  }

  /**
   * Gets the number of attributes that should be produced.
   * @return the number of attributes that should be produced
   */
  public int getNumAttributes() { 
    int       result;
    
    result = -1;
    try {
      result = Integer.parseInt(
          Utils.getOption('N', getGenerator().getOptions()));
    }
    catch (Exception e) {
      e.printStackTrace();
      result = -1;
    }

    return result;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String numAttributesTipText() {
    return "The number of attributes the generated data will contain (including class attribute), ie the number of nodes in the bayesian net.";
  }

  /**
   * returns the default cardinality
   * 
   * @return the default cardinality
   */
  protected int defaultCardinality() {
    return 2;
  }

  /**
   * Sets the cardinality of the attributes (incl class attribute)
   * @param value the cardinality
   */
  public void setCardinality(int value) { 
    setGeneratorOption("C", "" + value);
  }

  /**
   * Gets the cardinality of the attributes (incl class attribute)
   * @return the cardinality of the attributes
   */
  public int getCardinality() { 
    int       result;
    
    result = -1;
    try {
      result = Integer.parseInt(
          Utils.getOption('C', getGenerator().getOptions()));
    }
    catch (Exception e) {
      e.printStackTrace();
      result = -1;
    }

    return result;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String cardinalityTipText() {
    return "The cardinality of the attributes, incl the class attribute.";
  }

  /**
   * returns the default number of arcs
   * 
   * @return the default number of arcs
   */
  protected int defaultNumArcs() {
    return 20;
  }

  /**
   * Sets the number of arcs for the bayesian net
   * @param value the number of arcs
   */
  public void setNumArcs(int value) {
    int       nodes;
    int       minArcs;
    int       maxArcs;

    nodes   = getNumAttributes();
    minArcs = nodes - 1;
    maxArcs = nodes * (nodes - 1) / 2;
    
    if (value > maxArcs)
      throw new IllegalArgumentException(
          "Number of arcs should be at most nodes * (nodes - 1) / 2 = " 
          + maxArcs + " instead of " + value + " (nodes = numAttributes)!");
    else if (value < minArcs)
      throw new IllegalArgumentException(
          "Number of arcs should be at least (nodes - 1) = " + minArcs 
          + " instead of " + value + " (nodes = numAttributes)!");
    else
      setGeneratorOption("A", "" + value);
  }

  /**
   * Gets the number of arcs for the bayesian net
   * @return the number of arcs
   */
  public int getNumArcs() { 
    int       result;
    
    result = -1;
    try {
      result = Integer.parseInt(
          Utils.getOption('A', getGenerator().getOptions()));
    }
    catch (Exception e) {
      e.printStackTrace();
      result = -1;
    }

    return result;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String numArcsTipText() {
    return "The number of arcs in the bayesian net, at most: n * (n - 1) / 2 and at least: (n - 1); with n = numAttributes";
  }

  /**
   * Sets the number of examples, given by option.
   * @param numExamples the new number of examples
   */
  public void setNumExamples(int numExamples) { 
    super.setNumExamples(numExamples);
    setGeneratorOption("M", "" + numExamples);
  }

  /**
   * Gets the number of examples, given by option.
   * @return the number of examples, given by option 
   */
  public int getNumExamples() { 
    int       result;
    
    result = -1;
    try {
      result = Integer.parseInt(
          Utils.getOption('M', getGenerator().getOptions()));
    }
    catch (Exception e) {
      e.printStackTrace();
      result = -1;
    }

    return result;
  }

  /**
   * Return if single mode is set for the given data generator
   * mode depends on option setting and or generator type.
   * 
   * @return single mode flag
   * @throws Exception if mode is not set yet
   */
  public boolean getSingleModeFlag() throws Exception {
    return false;
  }

  /**
   * Initializes the format for the dataset produced. 
   * Must be called before the generateExample or generateExamples
   * methods are used.
   * Re-initializes the random number generator with the given seed.
   *
   * @return the format for the dataset 
   * @throws Exception if the generating of the format failed
   * @see  #getSeed()
   */
  public Instances defineDataFormat() throws Exception {
    BayesNetGenerator   bng;

    bng = new BayesNetGenerator();
    bng.setOptions(getGenerator().getOptions());
    setGeneratorOption(bng, "M", "1");
    bng.generateRandomNetwork();
    bng.generateInstances();
    bng.m_Instances.renameAttribute(0, "class");
    bng.m_Instances.setRelationName(getRelationNameToUse());
    
    return bng.m_Instances;
  }

  /**
   * Generates one example of the dataset. 
   *
   * @return the generated example
   * @throws Exception if the format of the dataset is not yet defined
   * @throws Exception if the generator only works with generateExamples
   * which means in non single mode
   */
  public Instance generateExample() throws Exception {
    throw new Exception("Cannot generate examples one-by-one!");
  }

  /**
   * Generates all examples of the dataset. Re-initializes the random number
   * generator with the given seed, before generating instances.
   *
   * @return the generated dataset
   * @throws Exception if the format of the dataset is not yet defined
   * @throws Exception if the generator only works with generateExample,
   * which means in single mode
   * @see   #getSeed()
   */
  public Instances generateExamples() throws Exception {
    getGenerator().setOptions(getGenerator().getOptions());
    getGenerator().generateRandomNetwork();
    getGenerator().generateInstances();
    getGenerator().m_Instances.renameAttribute(0, "class");
    getGenerator().m_Instances.setRelationName(getRelationNameToUse());
    
    return getGenerator().m_Instances;
  }

  /**
   * Generates a comment string that documentates the data generator.
   * By default this string is added at the beginning of the produced output
   * as ARFF file type, next after the options.
   * 
   * @return string contains info about the generated rules
   */
  public String generateStart () {
    return "";
  }

  /**
   * Generates a comment string that documentats the data generator.
   * By default this string is added at the end of theproduces output
   * as ARFF file type.
   * 
   * @return string contains info about the generated rules
   * @throws Exception if the generating of the documentaion fails
   */
  public String generateFinished() throws Exception {
    return "";
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
   * @param args should contain arguments for the data producer: 
   */
  public static void main(String[] args) {
    runDataGenerator(new BayesNet(), args);
  }
}
