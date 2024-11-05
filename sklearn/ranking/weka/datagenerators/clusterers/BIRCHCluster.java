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
 *    BIRCHCluster.java
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
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.datagenerators.ClusterGenerator;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Cluster data generator designed for the BIRCH System<br/>
 * <br/>
 * Dataset is generated with instances in K clusters.<br/>
 * Instances are 2-d data points.<br/>
 * Each cluster is characterized by the number of data points in itits radius and its center. The location of the cluster centers isdetermined by the pattern parameter. Three patterns are currentlysupported grid, sine and random.<br/>
 * <br/>
 * For more information refer to:<br/>
 * <br/>
 * Tian Zhang, Raghu Ramakrishnan, Miron Livny: BIRCH: An Efficient Data Clustering Method for Very Large Databases. In: ACM SIGMOD International Conference on Management of Data, 103-114, 1996.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Zhang1996,
 *    author = {Tian Zhang and Raghu Ramakrishnan and Miron Livny},
 *    booktitle = {ACM SIGMOD International Conference on Management of Data},
 *    pages = {103-114},
 *    publisher = {ACM Press},
 *    title = {BIRCH: An Efficient Data Clustering Method for Very Large Databases},
 *    year = {1996}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
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
 *  The number of attributes (default 10).</pre>
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
 * <pre> -k &lt;num&gt;
 *  The number of clusters (default 4)</pre>
 * 
 * <pre> -G
 *  Set pattern to grid (default is random).
 *  This flag cannot be used at the same time as flag I.
 *  The pattern is random, if neither flag G nor flag I is set.</pre>
 * 
 * <pre> -I
 *  Set pattern to sine (default is random).
 *  This flag cannot be used at the same time as flag I.
 *  The pattern is random, if neither flag G nor flag I is set.</pre>
 * 
 * <pre> -N &lt;num&gt;..&lt;num&gt;
 *  The range of number of instances per cluster (default 1..50).
 *  Lower number must be between 0 and 2500,
 *  upper number must be between 50 and 2500.</pre>
 * 
 * <pre> -R &lt;num&gt;..&lt;num&gt;
 *  The range of radius per cluster (default 0.1..1.4142135623730951).
 *  Lower number must be between 0 and SQRT(2), 
 *  upper number must be between SQRT(2) and SQRT(32).</pre>
 * 
 * <pre> -M &lt;num&gt;
 *  The distance multiplier (default 4.0).</pre>
 * 
 * <pre> -C &lt;num&gt;
 *  The number of cycles (default 4).</pre>
 * 
 * <pre> -O
 *  Flag for input order is ORDERED. If flag is not set then 
 *  input order is RANDOMIZED. RANDOMIZED is currently not 
 *  implemented, therefore is the input order always ORDERED.</pre>
 * 
 * <pre> -P &lt;num&gt;
 *  The noise rate in percent (default 0.0).
 *  Can be between 0% and 30%. (Remark: The original 
 *  algorithm only allows noise up to 10%.)</pre>
 * 
 <!-- options-end -->
 *
 * @author Gabi Schmidberger (gabi@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $ 
 */
public class BIRCHCluster 
  extends ClusterGenerator
  implements TechnicalInformationHandler {

  /** for serialization */
  static final long serialVersionUID = -334820527230755027L;
  
  /** Number of Clusters the dataset should have */
  protected int m_NumClusters;

  /** minimal number of instances per cluster (option N)*/ 
  private int m_MinInstNum;
  
  /** maximal number of instances per cluster (option N)*/ 
  private int m_MaxInstNum;
  
  /** minimum radius (option R)*/ 
  private double m_MinRadius;
  
  /** maximum radius (option R)*/ 
  private double m_MaxRadius;
  
  /**  Constant set for choice of pattern. (option G)*/
  public static final int GRID = 0;
  /**  Constant set for choice of pattern. (option I)*/
  public static final int SINE = 1;
  /**  Constant set for choice of pattern. (default)*/
  public static final int RANDOM = 2;
  /** the pattern tags */
  public static final Tag[] TAGS_PATTERN = {
    new Tag(GRID,   "Grid"),
    new Tag(SINE,   "Sine"),
    new Tag(RANDOM, "Random")
  };
  
  /** pattern (changed with options G or S)*/ 
  private int m_Pattern;
  
  /** distance multiplier (option M)*/
  private double m_DistMult;

  /** number of cycles (option C)*/
  private int m_NumCycles;

  /**  Constant set for input order (option O)*/
  public static final int ORDERED = 0;
  /**  Constant set for input order (default)*/
  public static final int RANDOMIZED = 1;
  /** the input order tags */
  public static final Tag[] TAGS_INPUTORDER = {
    new Tag(ORDERED,    "ordered"),
    new Tag(RANDOMIZED, "randomized")
  };

  /** input order (changed with option O)*/ 
  private int m_InputOrder;

  /** noise rate in percent (option P,  between 0 and 30)*/ 
  private double m_NoiseRate;

  /** cluster list */
  private FastVector m_ClusterList;

  // following are used for pattern is GRID
  /** grid size*/
  private int m_GridSize;

  /** grid width*/
  private double m_GridWidth;
  
  /**
   * class to represent cluster
   */
  private class Cluster 
    implements Serializable, RevisionHandler {

    /** for serialization */
    static final long serialVersionUID = -8336901069823498140L;    
    
    /** number of instances for this cluster */
    private int m_InstNum;

    /** radius of cluster
     *   variance is radius ** 2 / 2 */
    private double m_Radius;

    /** center of cluster = array of Double values */
    private double[] m_Center;

    /**
     * Constructor, used for pattern = RANDOM
     *
     * @param instNum the number of instances
     * @param radius radius of the cluster
     * @param random the random number generator to use 
     */
    private Cluster(int instNum, double radius, Random random) {
      m_InstNum = instNum;
      m_Radius = radius;
      m_Center = new double[getNumAttributes()];
      for (int i = 0; i < getNumAttributes(); i++) {
	m_Center[i] = random.nextDouble() * (double) m_NumClusters;
      }
    }

    /**
     * Constructor, used for pattern = GRID
     *
     * @param instNum the number of instances
     * @param radius radius of the cluster
     * @param gridVector vector for grid positions
     * @param gridWidth factor for grid position
     */
      // center is defined in the constructor of cluster
    private Cluster(int instNum,
		    double radius,
		    int[] gridVector,
		    double gridWidth) {
      m_InstNum = instNum;
      m_Radius = radius;
      m_Center = new double[getNumAttributes()];
      for (int i = 0; i < getNumAttributes(); i++) {
	m_Center[i] = ((double) gridVector[i] + 1.0) * gridWidth;
      }
      
    }
   
    /**
     * returns the number of instances
     * 
     * @return the number of instances
     */
    private int getInstNum() { 
      return m_InstNum; 
    }
    
    /**
     * returns the radius
     * 
     * @return the radius
     */
    private double getRadius() { 
      return m_Radius; 
    }
    
    /**
     * returns the variance
     * 
     * @return the variance
     */
    private double getVariance() { 
      return Math.pow(m_Radius, 2.0) / 2.0; 
    }
    
    /**
     * returns the standard deviation
     * 
     * @return the standard deviation
     */
    private double getStdDev() { 
      return (m_Radius / Math.pow(2.0, 0.5)); 
    }
    
    /**
     * returns the centers
     * 
     * @return the centers
     */
    private double[] getCenter() { 
      return m_Center; 
    }
    
    /**
     * returns the center value for a given dimension
     * 
     * @param dimension the dimension to return the center for
     * @return the center value for the given dimension
     * @throws Exception if dimension invalid
     */
    private double getCenterValue(int dimension) throws Exception {
      if (dimension >= m_Center.length)
	throw new Exception("Current system has only " +
			    m_Center.length + " dimensions.");
      return m_Center[dimension];
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 5987 $");
    }
  } // end class Cluster

  /**
   * class to represent Vector for placement of the center in space
   */
  private class GridVector 
    implements Serializable, RevisionHandler {

    /** for serialization */
    static final long serialVersionUID = -1900309948991039522L;
    
    /** array of integer */
    private int[] m_GridVector;

    /**  one higher then the highest possible integer value
     *  in any of the integers in the gridvector */
    private int m_Base;

    /** size of vector */
    private int m_Size;

    /**
     * Constructor
     *
     * @param numDim number of dimensions = number of attributes
     * @param base is one higher then the highest possible integer value
     * in any of the integers in the gridvector
     */
    private GridVector(int numDim, int base) {
      m_Size = numDim;
      m_Base = base;
      m_GridVector = new int [numDim];
      for (int i = 0; i < numDim; i++)
	m_GridVector[i] = 0;
    }

    /**
     * returns the integer array
     *
     * @return the integer array
     */
    private int[] getGridVector() {
      return m_GridVector;
    }

    /**
     * Overflow has occurred when integer is zero.
     *
     *@param digit the input integer
     *@return true if digit is 0
     */
    private boolean overflow(int digit) {
      return (digit == 0);
    }

    /**
     * Adds one to integer and sets to zero, if new value was
     * equal m_Base.
     *
     *@param digit the input integer
     *@return new integer object
     */
    private int addOne(int digit) {
      int value = digit + 1;
      if (value >= m_Base) value = 0;
      return value;
    }

    /**
     * add 1 to vector
     */
    private void addOne() {
      m_GridVector[0] = addOne(m_GridVector[0]);
      int i = 1;
      while (overflow(m_GridVector[i - 1]) && i < m_Size) {
        m_GridVector[i] = addOne(m_GridVector[i]);
	i++;
      }
	
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 5987 $");
    }
  } // end class GridVector
  
  /**
   * initializes the generator with default values
   */
  public BIRCHCluster() {
    super();

    setNumClusters(defaultNumClusters());
    setMinInstNum(defaultMinInstNum());
    setMaxInstNum(defaultMaxInstNum());
    setMinRadius(defaultMinRadius());
    setMaxRadius(defaultMaxRadius());
    setPattern(defaultPattern());
    setDistMult(defaultDistMult());
    setNumCycles(defaultNumCycles());
    setInputOrder(defaultInputOrder());
    setNoiseRate(defaultNoiseRate());
  }
  
  /**
   * Returns a string describing this data generator.
   *
   * @return a description of the data generator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Cluster data generator designed for the BIRCH System\n\n"
      + "Dataset is generated with instances in K clusters.\n"
      + "Instances are 2-d data points.\n"
      + "Each cluster is characterized by the number of data points in it"
      + "its radius and its center. The location of the cluster centers is"
      + "determined by the pattern parameter. Three patterns are currently"
      + "supported grid, sine and random.\n\n"
      + "For more information refer to:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Tian Zhang and Raghu Ramakrishnan and Miron Livny");
    result.setValue(Field.TITLE, "BIRCH: An Efficient Data Clustering Method for Very Large Databases");
    result.setValue(Field.BOOKTITLE, "ACM SIGMOD International Conference on Management of Data");
    result.setValue(Field.YEAR, "1996");
    result.setValue(Field.PAGES, "103-114");
    result.setValue(Field.PUBLISHER, "ACM Press");
    
    return result;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {
    Vector result = enumToVector(super.listOptions());

    result.addElement(new Option(
          "\tThe number of clusters (default "
          + defaultNumClusters() + ")",
          "k", 1, "-k <num>"));

    result.addElement(new Option(
          "\tSet pattern to grid (default is random).\n"
	  + "\tThis flag cannot be used at the same time as flag I.\n"
	  + "\tThe pattern is random, if neither flag G nor flag I is set.",
          "G", 0, "-G"));

    result.addElement(new Option(
          "\tSet pattern to sine (default is random).\n"
	  + "\tThis flag cannot be used at the same time as flag I.\n"
	  + "\tThe pattern is random, if neither flag G nor flag I is set.",
          "I", 0, "-I"));

    result.addElement(new Option(
          "\tThe range of number of instances per cluster (default "
          + defaultMinInstNum() + ".." + defaultMaxInstNum() + ").\n"
          + "\tLower number must be between 0 and 2500,\n"
          + "\tupper number must be between 50 and 2500.",
          "N", 1, "-N <num>..<num>"));

    result.addElement(new Option(
          "\tThe range of radius per cluster (default "
          + defaultMinRadius() + ".." + defaultMaxRadius() + ").\n"
          + "\tLower number must be between 0 and SQRT(2), \n"
          + "\tupper number must be between SQRT(2) and SQRT(32).",
          "R", 1, "-R <num>..<num>"));

    result.addElement(new Option(
          "\tThe distance multiplier (default " 
          + defaultDistMult() + ").",
          "M", 1, "-M <num>"));

    result.addElement(new Option(
          "\tThe number of cycles (default "
          + defaultNumCycles() + ").",
          "C", 1, "-C <num>"));

    result.addElement(new Option(
  	  "\tFlag for input order is ORDERED. If flag is not set then \n"
	  + "\tinput order is RANDOMIZED. RANDOMIZED is currently not \n"
	  + "\timplemented, therefore is the input order always ORDERED.",
          "O", 0, "-O"));

    result.addElement(new Option(
          "\tThe noise rate in percent (default " 
          + defaultNoiseRate() + ").\n"
          + "\tCan be between 0% and 30%. (Remark: The original \n"
          + "\talgorithm only allows noise up to 10%.)",
          "P", 1, "-P <num>"));

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
   *  The number of attributes (default 10).</pre>
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
   * <pre> -k &lt;num&gt;
   *  The number of clusters (default 4)</pre>
   * 
   * <pre> -G
   *  Set pattern to grid (default is random).
   *  This flag cannot be used at the same time as flag I.
   *  The pattern is random, if neither flag G nor flag I is set.</pre>
   * 
   * <pre> -I
   *  Set pattern to sine (default is random).
   *  This flag cannot be used at the same time as flag I.
   *  The pattern is random, if neither flag G nor flag I is set.</pre>
   * 
   * <pre> -N &lt;num&gt;..&lt;num&gt;
   *  The range of number of instances per cluster (default 1..50).
   *  Lower number must be between 0 and 2500,
   *  upper number must be between 50 and 2500.</pre>
   * 
   * <pre> -R &lt;num&gt;..&lt;num&gt;
   *  The range of radius per cluster (default 0.1..1.4142135623730951).
   *  Lower number must be between 0 and SQRT(2), 
   *  upper number must be between SQRT(2) and SQRT(32).</pre>
   * 
   * <pre> -M &lt;num&gt;
   *  The distance multiplier (default 4.0).</pre>
   * 
   * <pre> -C &lt;num&gt;
   *  The number of cycles (default 4).</pre>
   * 
   * <pre> -O
   *  Flag for input order is ORDERED. If flag is not set then 
   *  input order is RANDOMIZED. RANDOMIZED is currently not 
   *  implemented, therefore is the input order always ORDERED.</pre>
   * 
   * <pre> -P &lt;num&gt;
   *  The noise rate in percent (default 0.0).
   *  Can be between 0% and 30%. (Remark: The original 
   *  algorithm only allows noise up to 10%.)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String        tmpStr;
   
    super.setOptions(options);

    tmpStr = Utils.getOption('k', options);
    if (tmpStr.length() != 0)
      setNumClusters(Integer.parseInt(tmpStr));
    else
      setNumClusters(defaultNumClusters());

    tmpStr = Utils.getOption('N', options);
    if (tmpStr.length() != 0)
      setInstNums(tmpStr);
    else
      setInstNums(defaultMinInstNum() + ".." + defaultMaxInstNum());
    
    tmpStr = Utils.getOption('R', options);
    if (tmpStr.length() != 0)
      setRadiuses(tmpStr);
    else
      setRadiuses(defaultMinRadius() + ".." + defaultMaxRadius());

    boolean grid = Utils.getFlag('G', options);
    boolean sine = Utils.getFlag('I', options);

    if (grid && sine)
      throw new Exception("Flags -G and -I can only be set mutually exclusiv.");

    setPattern(new SelectedTag(RANDOM, TAGS_PATTERN));
    if (grid)
      setPattern(new SelectedTag(GRID, TAGS_PATTERN));
    if (sine)
      setPattern(new SelectedTag(SINE, TAGS_PATTERN));

    tmpStr= Utils.getOption('M', options);
    if (tmpStr.length() != 0) {
      if (!grid)
	throw new Exception("Option M can only be used with GRID pattern.");
      setDistMult(Double.parseDouble(tmpStr));
    }
    else {
      setDistMult(defaultDistMult());
    }

    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0) {
      if (!sine)
	throw new Exception("Option C can only be used with SINE pattern.");
      setNumCycles(Integer.parseInt(tmpStr));
    } 
    else {
      setNumCycles(defaultNumCycles());
    }

    if (Utils.getFlag('O', options))
      setInputOrder(new SelectedTag(ORDERED, TAGS_INPUTORDER));
    else
      setInputOrder(defaultInputOrder());

    tmpStr = Utils.getOption('P', options);
    if (tmpStr.length() != 0)
      setNoiseRate(Double.parseDouble(tmpStr));
    else
      setNoiseRate(defaultNoiseRate());
  }

  /**
   * Gets the current settings of the datagenerator BIRCHCluster.
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
    
    result.add("-k");
    result.add("" + getNumClusters());
    
    result.add("-N"); 
    result.add("" + getInstNums());

    result.add("-R"); 
    result.add("" + getRadiuses());

    if (m_Pattern == GRID) {
      result.add("-G");
      
      result.add("-M"); 
      result.add("" + getDistMult());
    }

    if (m_Pattern == SINE) {
      result.add("-I");
      
      result.add("-C"); 
      result.add("" + getNumCycles());
    }

    if (getOrderedFlag())
      result.add("-O");

    result.add("-P"); 
    result.add("" + getNoiseRate());
    
    return (String[]) result.toArray(new String[result.size()]);
  }

  /**
   * returns the default number of clusters
   * 
   * @return the default number of clusters
   */
  protected int defaultNumClusters() {
    return 4;
  }

  /**
   * Sets the number of clusters the dataset should have.
   * @param numClusters the new number of clusters
   */
  public void setNumClusters(int numClusters) { 
    m_NumClusters = numClusters; 
  }

  /**
   * Gets the number of clusters the dataset should have.
   * @return the number of clusters the dataset should have
   */
  public int getNumClusters() { 
    return m_NumClusters; 
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String numClustersTipText() {
    return "The number of clusters to generate.";
  }

  /**
   * Sets the upper and lower boundary for instances per cluster.
   *
   * @param fromTo  the string containing the upper and lower boundary for
   *                instances per cluster separated by ..
   */
  protected void setInstNums(String fromTo) {
    int i = fromTo.indexOf("..");
    String from = fromTo.substring(0, i);
    setMinInstNum(Integer.parseInt(from));
    String to = fromTo.substring(i + 2, fromTo.length());
    setMaxInstNum(Integer.parseInt(to));
  }
  
  /**
   * Gets the upper and lower boundary for instances per cluster.
   *
   * @return the string containing the upper and lower boundary for
   * instances per cluster separated by ..
   */
  protected String getInstNums() {
    String fromTo = "" 
                    + getMinInstNum() + ".."
                    + getMaxInstNum();
    return fromTo;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  protected String instNumsTipText() {
    return "The upper and lowet boundary for instances per cluster.";
  }

  /**
   * returns the default min number of instances
   * 
   * @return the default min number of instances
   */
  protected int defaultMinInstNum() {
    return 1;
  }

  /**
   * Gets the lower boundary for instances per cluster.
   *
   * @return the the lower boundary for instances per cluster
   */
  public int getMinInstNum() { 
    return m_MinInstNum; 
  }
  
  /**
   * Sets the lower boundary for instances per cluster.
   *
   * @param newMinInstNum new lower boundary for instances per cluster
   */
  public void setMinInstNum(int newMinInstNum) {
    m_MinInstNum = newMinInstNum;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String minInstNumTipText() {
    return "The lower boundary for instances per cluster.";
  }

  /**
   * returns the default max number of instances
   * 
   * @return the default max number of instances
   */
  protected int defaultMaxInstNum() {
    return 50;
  }

  /**
   * Gets the upper boundary for instances per cluster.
   *
   * @return the upper boundary for instances per cluster
   */
  public int getMaxInstNum() { 
    return m_MaxInstNum; 
  }
  
  /**
   * Sets the upper boundary for instances per cluster.
   *
   * @param newMaxInstNum new upper boundary for instances per cluster
   */
  public void setMaxInstNum(int newMaxInstNum) {
    m_MaxInstNum = newMaxInstNum;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String maxInstNumTipText() {
    return "The upper boundary for instances per cluster.";
  }

  /**
   * Sets the upper and lower boundary for the radius of the clusters.
   *
   * @param fromTo the string containing the upper and lower boundary for
   * the radius  of the clusters, separated by ..
   */
  protected void setRadiuses(String fromTo) {
    int i = fromTo.indexOf("..");
    String from = fromTo.substring(0, i);
    setMinRadius(Double.valueOf(from).doubleValue());
    String to = fromTo.substring(i + 2, fromTo.length());
    setMaxRadius(Double.valueOf(to).doubleValue());
  }

  /**
   * Gets the upper and lower boundary for the radius of the clusters.
   *
   * @return the string containing the upper and lower boundary for
   * the radius  of the clusters, separated by ..
   */
  protected String getRadiuses() {
    String fromTo = "" 
                    + Utils.doubleToString(getMinRadius(), 2) + ".."
                    + Utils.doubleToString(getMaxRadius(), 2);
    return fromTo;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  protected String radiusesTipText() {
    return "The upper and lower boundary for the radius of the clusters.";
  }

  /**
   * returns the default min radius
   * 
   * @return the default min radius
   */
  protected double defaultMinRadius() {
    return 0.1;
  }

  /**
   * Gets the lower boundary for the radiuses of the clusters.
   *
   * @return the lower boundary for the radiuses of the clusters
   */
  public double getMinRadius() { 
    return m_MinRadius; 
  }
  
  /**
   * Sets the lower boundary for the radiuses of the clusters.
   *
   * @param newMinRadius new lower boundary for the radiuses of the clusters
   */
  public void setMinRadius(double newMinRadius) {
    m_MinRadius = newMinRadius;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String minRadiusTipText() {
    return "The lower boundary for the radius of the clusters.";
  }

  /**
   * returns the default max radius
   * 
   * @return the default max radius
   */
  protected double defaultMaxRadius() {
    return Math.sqrt(2.0);
  }

  /**
   * Gets the upper boundary for the radiuses of the clusters.
   *
   * @return the upper boundary for the radiuses of the clusters
   */
  public double getMaxRadius() { 
    return m_MaxRadius; 
  }
  
  /**
   * Sets the upper boundary for the radiuses of the clusters.
   *
   * @param newMaxRadius new upper boundary for the radiuses of the clusters
   */
  public void setMaxRadius(double newMaxRadius) {
    m_MaxRadius = newMaxRadius;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String maxRadiusTipText() {
    return "The upper boundary for the radius of the clusters.";
  }

  /**
   * returns the default pattern
   * 
   * @return the default pattern
   */
  protected SelectedTag defaultPattern() {
    return new SelectedTag(RANDOM, TAGS_PATTERN);
  }
  
  /**
   * Gets the pattern type.
   *
   * @return the current pattern type
   */
  public SelectedTag getPattern() { 
    return new SelectedTag(m_Pattern, TAGS_PATTERN);
  }

  /**
   * Sets the pattern type.
   *
   * @param value new pattern type 
   */
  public void setPattern(SelectedTag value) {
    if (value.getTags() == TAGS_PATTERN)
      m_Pattern = value.getSelectedTag().getID();
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String patternTipText() {
    return "The pattern for generating the data.";
  }

  /**
   * returns the default distance multiplier
   * 
   * @return the default distance multiplier
   */
  protected double defaultDistMult() {
    return 4.0;
  }

  /**
   * Gets the distance multiplier.
   *
   * @return the distance multiplier
   */
  public double getDistMult() { 
    return m_DistMult; 
  }
  
  /**
   * Sets the distance multiplier.
   *
   * @param newDistMult new distance multiplier
   */
  public void setDistMult(double newDistMult) {
    m_DistMult = newDistMult;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String distMultTipText() {
    return "The distance multiplier (in combination with the 'Grid' pattern).";
  }

  /**
   * returns the default number of cycles
   * 
   * @return the default number of cycles
   */
  protected int defaultNumCycles() {
    return 4;
  }

  /**
   * Gets the number of cycles.
   *
   * @return the number of cycles
   */
  public int getNumCycles() { 
    return m_NumCycles; 
  }
  
  /**
   * Sets the the number of cycles.
   *
   * @param newNumCycles new number of cycles 
   */
  public void setNumCycles(int newNumCycles) {
    m_NumCycles = newNumCycles;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String numCyclesTipText() {
    return "The number of cycles to use (in combination with the 'Sine' pattern).";
  }

  /**
   * returns the default input order
   * 
   * @return the default input order
   */
  protected SelectedTag defaultInputOrder() {
    return new SelectedTag(ORDERED, TAGS_INPUTORDER);  // TODO: the only one that is currently implemented, normally RANDOMIZED
  }

  /**
   * Gets the input order.
   *
   * @return the current input order
   */
  public SelectedTag getInputOrder() { 
    return new SelectedTag(m_InputOrder, TAGS_INPUTORDER);
  }

  /**
   * Sets the input order.
   *
   * @param value new input order 
   */
  public void setInputOrder(SelectedTag value) {
    if (value.getTags() == TAGS_INPUTORDER)
      m_InputOrder = value.getSelectedTag().getID();
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String inputOrderTipText() {
    return "The input order to use.";
  }

  /**
   * Gets the ordered flag (option O).
   *
   * @return true if ordered flag is set
   */
  public boolean getOrderedFlag() { 
    return m_InputOrder == ORDERED; 
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
    Random random = new Random (getSeed());
    setRandom(random);
    Instances dataset;
    FastVector attributes = new FastVector(3);
    Attribute attribute;
    boolean classFlag = getClassFlag();
    
    FastVector classValues = null;
    if (classFlag) classValues = new FastVector (m_NumClusters);     

    // define dataset
    for (int i = 0; i < getNumAttributes(); i++) {
      attribute = new Attribute("X" + i); 
      attributes.addElement(attribute);
    }
    
    if (classFlag) {
      for (int i = 0; i < m_NumClusters; i++)
	classValues.addElement("c" + i);
      attribute = new Attribute ("class", classValues); 
      attributes.addElement(attribute);
    }

    dataset = new Instances(getRelationNameToUse(), attributes, 0);
    if (classFlag) 
      dataset.setClassIndex(getNumAttributes());

    // set dataset format of this class
    Instances format = new Instances(dataset, 0);
    setDatasetFormat(format);

    m_ClusterList = defineClusters(random);

    //System.out.println("dataset" + dataset.numAttributes());
    return dataset; 
  }

  /**
   * Generate an example of the dataset. 
   * @return the instance generated
   * @throws Exception if format not defined or generating <br/>
   * examples one by one is not possible, because voting is chosen
   */

  public Instance generateExample() throws Exception {
    throw new Exception("Examples cannot be generated" +
                                           " one by one.");
  }

  /**
   * Generate all examples of the dataset. 
   * @return the instance generated
   * @throws Exception if format not defined 
   */

  public Instances generateExamples() throws Exception {
    Random random = getRandom();
    Instances data = getDatasetFormat();
    if (data == null) throw new Exception("Dataset format not defined.");

    // generate examples
    if (getOrderedFlag())
      data = generateExamples(random, data);
    else
      throw new Exception("RANDOMIZED is not yet implemented.");
  
    return (data);
  }

  /**
   * Generate all examples of the dataset. 
   * 
   * @param random the random number generator to use
   * @param format the dataset format
   * @return the instance generated
   * @throws Exception if format not defined
   */
  public Instances generateExamples(Random random,
				    Instances format) throws Exception {
    Instance example = null;
    
    if (format == null) 
      throw new Exception("Dataset format not defined.");

    // generate examples for one cluster after another
    int cNum = 0;
    for (Enumeration enm = m_ClusterList.elements();
	 enm.hasMoreElements(); cNum++) {
      Cluster cl  = (Cluster) enm.nextElement();
      double stdDev = cl.getStdDev();
      int instNum = cl.getInstNum();
      double[] center = cl.getCenter();
      String cName = "c" + cNum;

      for (int i = 0; i < instNum; i++) {
	// generate example
	example = generateInstance(
                    format, random, stdDev, center, cName);
       
	if (example != null)
	  example.setDataset(format);
	format.add(example);
      }
    }

    return (format);
  }

  /**
   * Generate an example of the dataset. 
   * 
   * @param format the dataset format
   * @param randomG the random number generator
   * @param stdDev the standard deviation to use
   * @param center the centers
   * @param cName the class value
   * @return the instance generated
   * examples one by one is not possible, because voting is chosen
   */
  private Instance generateInstance (Instances format,
				     Random randomG,
				     double stdDev,
				     double[] center,
				     String cName) {
    Instance example;
    int numAtts = getNumAttributes();
    if (getClassFlag()) 
      numAtts++;

    example = new DenseInstance(numAtts);
    example.setDataset(format);
        
    for (int i = 0; i < getNumAttributes(); i++)
      example.setValue(i, randomG.nextGaussian() * stdDev + center[i]); 
    
    if (getClassFlag())
      example.setClassValue(cName);

    return example; 
  }

 /**
   * Defines the clusters 
   *
   * @param random random number generator
   * @return the cluster definitions
   * @throws Exception if defining fails
   */
  private FastVector defineClusters(Random random)
   throws Exception {

    if (m_Pattern == GRID)
      return defineClustersGRID(random);
    else
      return defineClustersRANDOM(random);
  }

  /**
   * Defines the clusters if pattern is GRID
   *
   * @param random random number generator
   * @return the defined clusters for GRID
   * @throws Exception if something goes wrong
   */
  private FastVector defineClustersGRID(Random random)
    throws Exception {

    FastVector clusters = new FastVector(m_NumClusters);
    double diffInstNum = (double) (m_MaxInstNum - m_MinInstNum);
    double minInstNum = (double) m_MinInstNum;
    double diffRadius = m_MaxRadius - m_MinRadius;
    Cluster cluster;

    // compute gridsize
    double gs = Math.pow(m_NumClusters, 1.0 / getNumAttributes());
    
    if (gs - ((double) ((int) gs))  > 0.0) {
      m_GridSize = (int) (gs + 1.0);
    } else { m_GridSize = (int) gs; }

    // compute gridwidth
    m_GridWidth = ((m_MaxRadius + m_MinRadius) / 2) * m_DistMult;

    //System.out.println("GridSize= " + m_GridSize);
    //System.out.println("GridWidth= " + m_GridWidth);
    
    // initialize gridvector with zeros
    GridVector gv = new GridVector(getNumAttributes(), m_GridSize);

    for (int i = 0; i < m_NumClusters; i++) {
      int instNum = (int) (random.nextDouble() * diffInstNum
                                   + minInstNum);
      double radius = (random.nextDouble() * diffRadius) + m_MinRadius;

      // center is defined in the constructor of cluster
      cluster = new Cluster(instNum, radius,
			    gv.getGridVector(), m_GridWidth);
      clusters.addElement((Object) cluster);
      gv.addOne();
    }
    return clusters;
  }

 /**
   * Defines the clusters if pattern is RANDOM
   *
   * @param random random number generator
   * @return the cluster definitions
   * @throws Exception if something goes wrong
   */
  private FastVector defineClustersRANDOM(Random random)
    throws Exception {

    FastVector clusters = new FastVector(m_NumClusters);
    double diffInstNum = (double) (m_MaxInstNum - m_MinInstNum);
    double minInstNum = (double) m_MinInstNum;
    double diffRadius = m_MaxRadius - m_MinRadius;
    Cluster cluster;

    for (int i = 0; i < m_NumClusters; i++) {
      int instNum = (int) (random.nextDouble() * diffInstNum
                                   + minInstNum);
      double radius = (random.nextDouble() * diffRadius) + m_MinRadius;

      // center is defined in the constructor of cluster
      cluster = new Cluster(instNum, radius, random);
      clusters.addElement((Object) cluster);
    }
    return clusters;
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
    int cNum = 0;
    for (Enumeration enm = m_ClusterList.elements();
	 enm.hasMoreElements(); cNum++) {
      Cluster cl  = (Cluster) enm.nextElement();
      docu.append("%\n");
      docu.append("% Cluster: c"+ cNum + "\n");
      docu.append("% ----------------------------------------------\n");
      docu.append("% StandardDeviation: "
		  + Utils.doubleToString(cl.getStdDev(), 2) + "\n");
      docu.append("% Number of instances: "
		  + cl.getInstNum() + "\n");
      sumInst += cl.getInstNum();
      double[] center = cl.getCenter();
      docu.append("% "); 
      for (int i = 0; i < center.length - 1; i++) {
        docu.append(Utils.doubleToString(center[i], 2) + ", ");
      }
      docu.append(Utils.doubleToString(center[center.length - 1], 2) + "\n");
    }
    docu.append("%\n% ----------------------------------------------\n"); 
    docu.append("% Total number of instances: " + sumInst + "\n");
    docu.append("%                            in " + cNum + " clusters\n");
    docu.append("% Pattern chosen           : ");
    if (m_Pattern == GRID) 
      docu.append(
          "GRID, " + "distance multiplier = " 
          + Utils.doubleToString(m_DistMult, 2) + "\n");
    else if (m_Pattern == SINE) 
      docu.append("SINE\n");
    else
      docu.append("RANDOM\n");

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
    runDataGenerator(new BIRCHCluster(), args);
  }
}

