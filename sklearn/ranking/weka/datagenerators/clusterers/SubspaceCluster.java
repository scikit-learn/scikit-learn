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
 *    SubspaceCluster.java
 *    Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.datagenerators.clusterers;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.Range;
import weka.core.RevisionUtils;
import weka.core.Tag;
import weka.core.Utils;
import weka.datagenerators.ClusterDefinition;
import weka.datagenerators.ClusterGenerator;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * A data generator that produces data points in hyperrectangular subspace clusters.
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
 * <pre> -a &lt;num&gt;
 *  The number of attributes (default 1).</pre>
 * 
 * <pre> -c
 *  Class Flag, if set, the cluster is listed in extra attribute.</pre>
 * 
 * <pre> -b &lt;range&gt;
 *  The indices for boolean attributes.</pre>
 * 
 * <pre> -m &lt;range&gt;
 *  The indices for nominal attributes.</pre>
 * 
 * <pre> -P &lt;num&gt;
 *  The noise rate in percent (default 0.0).
 *  Can be between 0% and 30%. (Remark: The original 
 *  algorithm only allows noise up to 10%.)</pre>
 * 
 * <pre> -C &lt;cluster-definition&gt;
 *  A cluster definition of class 'SubspaceClusterDefinition'
 *  (definition needs to be quoted to be recognized as 
 *  a single argument).</pre>
 * 
 * <pre> 
 * Options specific to weka.datagenerators.clusterers.SubspaceClusterDefinition:
 * </pre>
 * 
 * <pre> -A &lt;range&gt;
 *  Generates randomly distributed instances in the cluster.</pre>
 * 
 * <pre> -U &lt;range&gt;
 *  Generates uniformly distributed instances in the cluster.</pre>
 * 
 * <pre> -G &lt;range&gt;
 *  Generates gaussian distributed instances in the cluster.</pre>
 * 
 * <pre> -D &lt;num&gt;,&lt;num&gt;
 *  The attribute min/max (-A and -U) or mean/stddev (-G) for
 *  the cluster.</pre>
 * 
 * <pre> -N &lt;num&gt;..&lt;num&gt;
 *  The range of number of instances per cluster (default 1..50).</pre>
 * 
 * <pre> -I
 *  Uses integer instead of continuous values (default continuous).</pre>
 * 
 <!-- options-end -->
 *
 * @author Gabi Schmidberger (gabi@cs.waikato.ac.nz)
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $ 
 */
public class SubspaceCluster 
  extends ClusterGenerator {

  /** for serialization */
  static final long serialVersionUID = -3454999858505621128L;
  
  /** noise rate in percent (option P,  between 0 and 30)*/ 
  protected double m_NoiseRate;

  /** cluster list */
  protected ClusterDefinition[] m_Clusters;

  /** if nominal, store number of values */
  protected int[] m_numValues;

  /** store global min values */
  protected double[] m_globalMinValue;

  /** store global max values */
  protected double[] m_globalMaxValue;

  /** cluster type: uniform/random */
  public static final int UNIFORM_RANDOM = 0;  
  /** cluster type: total uniform */
  public static final int TOTAL_UNIFORM = 1;
  /** cluster type: gaussian */
  public static final int GAUSSIAN = 2;
  /** the tags for the cluster types */
  public static final Tag[] TAGS_CLUSTERTYPE = {
    new Tag(UNIFORM_RANDOM, "uniform/random"),
    new Tag(TOTAL_UNIFORM,  "total uniform"),
    new Tag(GAUSSIAN,       "gaussian")
  };

  /** cluster subtype: continuous */
  public static final int CONTINUOUS = 0;
  /** cluster subtype: integer */
  public static final int INTEGER = 1;
  /** the tags for the cluster types */
  public static final Tag[] TAGS_CLUSTERSUBTYPE = {
    new Tag(CONTINUOUS, "continuous"),
    new Tag(INTEGER,    "integer")
  };

  /**
   * initializes the generator, sets the number of clusters to 0, since user
   * has to specify them explicitly
   */
  public SubspaceCluster() {
    super();

    setNoiseRate(defaultNoiseRate());
  }

  /**
   * Returns a string describing this data generator.
   *
   * @return a description of the data generator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "A data generator that produces data points in "
      + "hyperrectangular subspace clusters.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {
    Vector result = enumToVector(super.listOptions());

    result.addElement(new Option(
          "\tThe noise rate in percent (default " 
          + defaultNoiseRate() + ").\n"
          + "\tCan be between 0% and 30%. (Remark: The original \n"
          + "\talgorithm only allows noise up to 10%.)",
          "P", 1, "-P <num>"));

    result.addElement(new Option(
          "\tA cluster definition of class '" 
	  + SubspaceClusterDefinition.class.getName().replaceAll(".*\\.", "") + "'\n"
	  + "\t(definition needs to be quoted to be recognized as \n"
	  + "\ta single argument).",
          "C", 1, "-C <cluster-definition>"));

    result.addElement(new Option(
	      "", "", 0, 
	      "\nOptions specific to " 
	      + SubspaceClusterDefinition.class.getName() + ":"));

    result.addAll(
        enumToVector(new SubspaceClusterDefinition(this).listOptions()));

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
   * <pre> -a &lt;num&gt;
   *  The number of attributes (default 1).</pre>
   * 
   * <pre> -c
   *  Class Flag, if set, the cluster is listed in extra attribute.</pre>
   * 
   * <pre> -b &lt;range&gt;
   *  The indices for boolean attributes.</pre>
   * 
   * <pre> -m &lt;range&gt;
   *  The indices for nominal attributes.</pre>
   * 
   * <pre> -P &lt;num&gt;
   *  The noise rate in percent (default 0.0).
   *  Can be between 0% and 30%. (Remark: The original 
   *  algorithm only allows noise up to 10%.)</pre>
   * 
   * <pre> -C &lt;cluster-definition&gt;
   *  A cluster definition of class 'SubspaceClusterDefinition'
   *  (definition needs to be quoted to be recognized as 
   *  a single argument).</pre>
   * 
   * <pre> 
   * Options specific to weka.datagenerators.clusterers.SubspaceClusterDefinition:
   * </pre>
   * 
   * <pre> -A &lt;range&gt;
   *  Generates randomly distributed instances in the cluster.</pre>
   * 
   * <pre> -U &lt;range&gt;
   *  Generates uniformly distributed instances in the cluster.</pre>
   * 
   * <pre> -G &lt;range&gt;
   *  Generates gaussian distributed instances in the cluster.</pre>
   * 
   * <pre> -D &lt;num&gt;,&lt;num&gt;
   *  The attribute min/max (-A and -U) or mean/stddev (-G) for
   *  the cluster.</pre>
   * 
   * <pre> -N &lt;num&gt;..&lt;num&gt;
   *  The range of number of instances per cluster (default 1..50).</pre>
   * 
   * <pre> -I
   *  Uses integer instead of continuous values (default continuous).</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String                      tmpStr;
    SubspaceClusterDefinition   cl;
    Vector                      list;
    int                         clCount;

    super.setOptions(options);

    m_numValues = new int[getNumAttributes()];
    // numValues might be changed by a cluster definition
    // (only relevant for nominal data)
    for (int i = 0; i < getNumAttributes(); i++)
      m_numValues[i] = 1;

    tmpStr = Utils.getOption('P', options);
    if (tmpStr.length() != 0)
      setNoiseRate(Double.parseDouble(tmpStr));
    else
      setNoiseRate(defaultNoiseRate());

    // cluster definitions
    list = new Vector();
    
    clCount = 0;
    do {
      tmpStr = Utils.getOption('C', options);
      if (tmpStr.length() != 0) {
        clCount++;
        cl = new SubspaceClusterDefinition(this);
        cl.setOptions(Utils.splitOptions(tmpStr));
        list.add(cl);
      }
    }
    while (tmpStr.length() != 0);

    m_Clusters = (ClusterDefinition[]) 
                    list.toArray(new ClusterDefinition[list.size()]);
    // in case no cluster definition was provided, make sure that there's at
    // least one definition present -> see getClusters()
    getClusters();
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
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-P"); 
    result.add("" + getNoiseRate());

    for (i = 0; i < getClusters().length; i++)  {
      result.add("-C");
      result.add(Utils.joinOptions(getClusters()[i].getOptions()));
    }

    return (String[]) result.toArray(new String[result.size()]);
  }

  /**
   * returns the current cluster definitions, if necessary initializes them
   * 
   * @return the current cluster definitions
   */
  protected ClusterDefinition[] getClusters() {
    if ( (m_Clusters == null) || (m_Clusters.length == 0) ) {
      if (m_Clusters != null)
        System.out.println("NOTE: at least 1 cluster definition is necessary, " 
            + "created default one.");
      m_Clusters = new ClusterDefinition[]{new SubspaceClusterDefinition(this)};
    }

    return m_Clusters;
  }

  /**
   * returns the default number of attributes
   * 
   * @return the default number of attributes
   */
  protected int defaultNumAttributes() {
    return 1;
  }

  /**
   * Sets the number of attributes the dataset should have.
   * @param numAttributes the new number of attributes
   */
  public void setNumAttributes(int numAttributes) {
    super.setNumAttributes(numAttributes);
    m_numValues = new int[getNumAttributes()];
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String numAttributesTipText() {
    return "The number of attributes the generated data will contain (Note: they must be covered by the cluster definitions!)";
  }

  /**
   * returns the default noise rate
   * 
   * @return the default noise rate
   */
  protected double defaultNoiseRate() {
    return 0.0;
  }

  /**
   * Gets the percentage of noise set.
   *
   * @return the percentage of noise set
   */
  public double getNoiseRate() { 
    return m_NoiseRate; 
  }

  /**
   * Sets the percentage of noise set.
   *
   * @param newNoiseRate new percentage of noise 
   */
  public void setNoiseRate(double newNoiseRate) {
    m_NoiseRate = newNoiseRate;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String noiseRateTipText() {
    return "The noise rate to use.";
  }

  /**
   * returns the currently set clusters
   * 
   * @return the currently set clusters
   */
  public ClusterDefinition[] getClusterDefinitions() {
    return getClusters();
  }

  /**
   * sets the clusters to use
   * 
   * @param value the clusters do use
   * @throws Exception if clusters are not the correct class
   */
  public void setClusterDefinitions(ClusterDefinition[] value) 
    throws Exception {

    String      indexStr;
    
    indexStr   = "";
    m_Clusters = value;
    for (int i = 0; i < getClusters().length; i++) {
      if (!(getClusters()[i] instanceof SubspaceClusterDefinition)) {
        if (indexStr.length() != 0)
          indexStr += ",";
        indexStr += "" + (i+1);
      }
      getClusters()[i].setParent(this);
      getClusters()[i].setOptions(getClusters()[i].getOptions()); // for initializing!
    }

    // any wrong classes encountered?
    if (indexStr.length() != 0)
      throw new Exception("These cluster definitions are not '" 
          + SubspaceClusterDefinition.class.getName() + "': " + indexStr);
  }

  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String clusterDefinitionsTipText() {
    return "The clusters to use.";
  }

  /**
   * Checks, whether all attributes are covered by cluster definitions and 
   * returns TRUE in that case.
   * 
   * @return whether all attributes are covered
   */
  protected boolean checkCoverage() {
    int         i;
    int         n;
    int[]       count;
    Range       r;
    String      attrIndex;
    SubspaceClusterDefinition  cl;
    
    // check whether all the attributes are covered
    count = new int[getNumAttributes()];
    for (i = 0; i < getNumAttributes(); i++) {
      for (n = 0; n < getClusters().length; n++) {
        cl = (SubspaceClusterDefinition) getClusters()[n];
        r  = new Range(cl.getAttrIndexRange());
        r.setUpper(getNumAttributes());
        if (r.isInRange(i))
          count[i]++;
      }
    }

    // list all indices that are not covered
    attrIndex = "";
    for (i = 0; i < count.length; i++) {
      if (count[i] == 0) {
        if (attrIndex.length() != 0)
          attrIndex += ",";
        attrIndex += (i+1);
      }
    }

    if (attrIndex.length() != 0)
      throw new IllegalArgumentException(
          "The following attributes are not covered by a cluster "
          + "definition: " + attrIndex + "\n");

    return true;
  }

  /**
   * Gets the single mode flag.
   *
   * @return true if methode generateExample can be used.
   */
  public boolean getSingleModeFlag() { 
    return false; 
  }

  /**
   * Initializes the format for the dataset produced. 
   *
   * @return the output data format
   * @throws Exception data format could not be defined 
   */

  public Instances defineDataFormat() throws Exception {

    // initialize
    setOptions(getOptions());

    checkCoverage();

    Random random = new Random (getSeed());
    setRandom(random);
    Instances dataset;
    FastVector attributes = new FastVector(3);
    Attribute attribute;
    boolean classFlag = getClassFlag();

    FastVector classValues = null;
    if (classFlag) 
      classValues = new FastVector(getClusters().length);     
    FastVector boolValues = new FastVector(2);
    boolValues.addElement("false");
    boolValues.addElement("true");
    FastVector nomValues = null;

    // define dataset
    for (int i = 0; i < getNumAttributes(); i++) {
      // define boolean attribute
      if (m_booleanCols.isInRange(i)) {
        attribute = new Attribute("B" + i, boolValues);
      } 
      else if (m_nominalCols.isInRange(i)) {
        // define nominal attribute
        nomValues = new FastVector(m_numValues[i]);
        for (int j = 0; j < m_numValues[i]; j++)
          nomValues.addElement("value-" + j);
        attribute = new Attribute("N" + i, nomValues);
      } 
      else {
        // numerical attribute
        attribute = new Attribute("X" + i); 
      }
      attributes.addElement(attribute);
    }

    if (classFlag) {
      for (int i = 0; i < getClusters().length; i++)
        classValues.addElement("c" + i);
      attribute = new Attribute ("class", classValues); 
      attributes.addElement(attribute);
    }

    dataset = new Instances(getRelationNameToUse(), attributes, 0);
    if (classFlag) 
      dataset.setClassIndex(m_NumAttributes);

    // set dataset format of this class
    Instances format = new Instances(dataset, 0);
    setDatasetFormat(format);

    for (int i = 0; i < getClusters().length; i++) {
      SubspaceClusterDefinition cl = (SubspaceClusterDefinition) getClusters()[i];
      cl.setNumInstances(random);
      cl.setParent(this);
    }

    return dataset; 
  }

  /**
   * Returns true if attribute is boolean
   *@param index of the attribute
   *@return true if the attribute is boolean
   */
  public boolean isBoolean(int index) {
    return m_booleanCols.isInRange(index); 
  }

  /**
   * Returns true if attribute is nominal
   *@param index of the attribute
   *@return true if the attribute is nominal
   */
  public boolean isNominal(int index) {
    return m_nominalCols.isInRange(index);
  }

  /**
   * returns array that stores the number of values for a nominal attribute.
   * 
   * @return the array that stores the number of values for a nominal attribute
   */
  public int[] getNumValues() {
    return m_numValues;
  }

  /**
   * Generate an example of the dataset. 
   * @return the instance generated
   * @throws Exception if format not defined or generating <br/>
   * examples one by one is not possible, because voting is chosen
   */

  public Instance generateExample() throws Exception {
    throw new Exception("Examples cannot be generated one by one.");
  }

  /**
   * Generate all examples of the dataset. 
   * @return the instance generated
   * @throws Exception if format not defined 
   */

  public Instances generateExamples() throws Exception {
    Instances format = getDatasetFormat();
    Instance example = null;

    if (format == null) 
      throw new Exception("Dataset format not defined.");

    // generate examples for one cluster after another
    for (int cNum = 0; cNum < getClusters().length; cNum++) {
      SubspaceClusterDefinition cl  = (SubspaceClusterDefinition) getClusters()[cNum];

      //get the number of instances to create
      int instNum = cl.getNumInstances();

      //class value is c + cluster number
      String cName = "c" + cNum;

      switch (cl.getClusterType().getSelectedTag().getID()) {
        case (UNIFORM_RANDOM):
          for (int i = 0; i < instNum; i++) {
            // generate example
            example = generateExample(format, getRandom(), cl, cName);
            if (example != null)
              format.add(example);
          }
          break;
        case (TOTAL_UNIFORM):
          // generate examples
          if (!cl.isInteger())
            generateUniformExamples(format, instNum, cl, cName);
          else
            generateUniformIntegerExamples(format, instNum, cl, cName);
          break;
        case (GAUSSIAN):
          // generate examples
          generateGaussianExamples(format, instNum, getRandom(), cl, cName);
          break;
      }
    }

    return format;
  }

  /**
   * Generate an example of the dataset. 
   * 
   * @param format the dataset format
   * @param randomG the random number generator to use
   * @param cl the cluster definition
   * @param cName the class value
   * @return the generated instance
   */
  private Instance generateExample(
      Instances format, Random randomG, SubspaceClusterDefinition cl, 
      String cName) {

    boolean makeInteger = cl.isInteger();
    int num = -1;
    Instance example = null;
    int numAtts = m_NumAttributes;
    if (getClassFlag()) numAtts++;

    example = new DenseInstance(numAtts);
    example.setDataset(format);
    boolean[] attributes = cl.getAttributes();
    double[] minValue = cl.getMinValue();
    double[] maxValue = cl.getMaxValue();
    double value;

    int clusterI = -1;
    for (int i = 0; i < m_NumAttributes; i++) {
      if (attributes[i]) {
        clusterI++;
        num++;
        // boolean  or nominal attribute
        if (isBoolean(i) || isNominal(i)) {

          if (minValue[clusterI] == maxValue[clusterI]) {
            value = minValue[clusterI];
          } 
          else {
            int numValues = (int)(maxValue[clusterI] - minValue[clusterI] + 1.0);
            value = randomG.nextInt(numValues);
            value += minValue[clusterI];
          }
        } 
        else {
          // numeric attribute
          value = randomG.nextDouble() * 
            (maxValue[num] - minValue[num]) + minValue[num];
          if (makeInteger)
            value = Math.round(value);
        }
        example.setValue(i, value);
      } 
      else {
        example.setMissing(i);
      }
    }

    if (getClassFlag())
      example.setClassValue(cName);

    return example; 
  }

  /**
   * Generate examples for a uniform cluster dataset. 
   * 
   * @param format the dataset format
   * @param numInstances the number of instances to generator
   * @param cl the cluster definition
   * @param cName the class value
   */
  private void generateUniformExamples(
      Instances format, int numInstances, SubspaceClusterDefinition cl, 
      String cName) {

    Instance example = null;
    int numAtts = m_NumAttributes;
    if (getClassFlag()) numAtts++;

    example = new DenseInstance(numAtts);
    example.setDataset(format);
    boolean[] attributes = cl.getAttributes();
    double[] minValue = cl.getMinValue();
    double[] maxValue = cl.getMaxValue();
    double[] diff = new double[minValue.length];

    for (int i = 0; i < minValue.length; i++)
      diff[i] = (maxValue[i] - minValue[i]);

    for (int j = 0; j < numInstances; j++) {
      int num = -1;
      for (int i = 0; i < m_NumAttributes; i++) {
        if (attributes[i]) {
          num++;
          double value = minValue[num] + (diff[num] * (double)((double)j / (double)(numInstances - 1)));
          example.setValue(i, value);
        } 
        else {
          example.setMissing(i);
        }
      }
      if (getClassFlag())
        example.setClassValue(cName);
      format.add(example);
    }
  }

  /**
   * Generate examples for a uniform cluster dataset. 
   * 
   * @param format the dataset format
   * @param numInstances the number of instances to generator
   * @param cl the cluster definition
   * @param cName the class value
   */
  private void generateUniformIntegerExamples(
      Instances format, int numInstances, SubspaceClusterDefinition cl, 
      String cName) {

    Instance example = null;
    int numAtts = m_NumAttributes;
    if (getClassFlag()) numAtts++;

    example = new DenseInstance(numAtts);
    example.setDataset(format);
    boolean[] attributes = cl.getAttributes();
    double[] minValue = cl.getMinValue();
    double[] maxValue = cl.getMaxValue();
    int[] minInt = new int[minValue.length];
    int[] maxInt = new int[maxValue.length];
    int[] intValue = new int[maxValue.length];
    int[] numInt = new int[minValue.length];

    int num = 1;
    for (int i = 0; i < minValue.length; i++) {
      minInt[i] = (int)Math.ceil(minValue[i]);
      maxInt[i] = (int)Math.floor(maxValue[i]);
      numInt[i] = (maxInt[i] - minInt[i] + 1);
      num = num * numInt[i];
    }
    int numEach = numInstances / num;
    int rest = numInstances - numEach * num;

    // initialize with smallest values combination
    for (int i = 0; i < m_NumAttributes; i++) {
      if (attributes[i]) {
        example.setValue(i, (double)minInt[i]);
        intValue[i] = minInt[i];
      } 
      else {
        example.setMissing(i);
      }
    }
    if (getClassFlag())
      example.setClassValue(cName);
    int added = 0;
    int attr = 0;
    // do while not added all
    do {
      // add all for one value combination
      for (int k = 0; k < numEach; k++) {
        format.add(example);
        example = (Instance) example.copy();
        added++;
      }
      if (rest > 0) {
        format.add(example);
        example = (Instance) example.copy();
        added++;
        rest--;
      }

      if (added >= numInstances) break;
      // switch to the next value combination
      boolean done = false;
      do {
        if (attributes[attr] && (intValue[attr] + 1 <= maxInt[attr])) {
          intValue[attr]++;
          done = true;
        } 
        else {
          attr++;
        }
      } while (!done);

      example.setValue(attr, (double)intValue[attr]);
    } while (added < numInstances);
  }

  /**
   * Generate examples for a uniform cluster dataset. 
   * 
   * @param format the dataset format
   * @param numInstances the number of instances to generate
   * @param random the random number generator
   * @param cl the cluster definition
   * @param cName the class value
   */
  private void generateGaussianExamples(
      Instances format, int numInstances, Random random, 
      SubspaceClusterDefinition cl, String cName) {

    boolean makeInteger = cl.isInteger();
    Instance example = null;
    int numAtts = m_NumAttributes;
    if (getClassFlag()) numAtts++;

    example = new DenseInstance(numAtts);
    example.setDataset(format);
    boolean[] attributes = cl.getAttributes();
    double[] meanValue = cl.getMeanValue();
    double[] stddevValue = cl.getStddevValue();

    for (int j = 0; j < numInstances; j++) {
      int num = -1;
      for (int i = 0; i < m_NumAttributes; i++) {
        if (attributes[i]) {
          num++;
          double value = meanValue[num] + (random.nextGaussian() * stddevValue[num]);
          if (makeInteger)
            value = Math.round(value);
          example.setValue(i, value);
        } 
        else {
          example.setMissing(i);
        }
      }
      if (getClassFlag())
        example.setClassValue(cName);
      format.add(example);
    }
  }

  /**
   * Compiles documentation about the data generation after
   * the generation process
   *
   * @return string with additional information about generated dataset
   * @throws Exception no input structure has been defined
   */
  public String generateFinished() throws Exception {
    return "";
  }

  /**
   * Compiles documentation about the data generation before
   * the generation process
   *
   * @return string with additional information 
   */
  public String generateStart() {
    StringBuffer docu = new StringBuffer();

    int sumInst = 0;
    for (int cNum = 0; cNum < getClusters().length; cNum++) {
      SubspaceClusterDefinition cl  = (SubspaceClusterDefinition) getClusters()[cNum];
      docu.append("%\n");
      docu.append("% Cluster: c"+ cNum + "   ");
      switch (cl.getClusterType().getSelectedTag().getID()) {
        case UNIFORM_RANDOM: 
          docu.append("Uniform Random");
          break;
        case TOTAL_UNIFORM: 
          docu.append("Total Random");
          break;
        case GAUSSIAN: 
          docu.append("Gaussian");
          break;
      }
      if (cl.isInteger()) {
        docu.append(" / INTEGER");
      }

      docu.append("\n% ----------------------------------------------\n");
      docu.append("%"+cl.attributesToString());

      docu.append("\n% Number of Instances:            "  + cl.getInstNums() + "\n");
      docu.append(  "% Generated Number of Instances:  "  + cl.getNumInstances() + "\n");
      sumInst += cl.getNumInstances();
        }
    docu.append("%\n% ----------------------------------------------\n"); 
    docu.append("% Total Number of Instances: " + sumInst + "\n");
    docu.append("%                            in " + getClusters().length + " Cluster(s)\n%");

    return docu.toString();
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
   * @param args should contain arguments for the data producer: 
   */
  public static void main(String[] args) {
    runDataGenerator(new SubspaceCluster(), args);
  }
}
