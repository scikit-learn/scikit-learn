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
 *    SimpleKMeans.java
 *    Copyright (C) 2000-2010 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.clusterers;

import weka.classifiers.rules.DecisionTableHashKey;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.DistanceFunction;
import weka.core.EuclideanDistance;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.ManhattanDistance;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

import java.util.Enumeration;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Cluster data using the k means algorithm. Can use either the Euclidean distance (default) or the Manhattan distance. If the Manhattan distance is used, then centroids are computed as the component-wise median rather than mean.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -N &lt;num&gt;
 *  number of clusters.
 *  (default 2).</pre>
 * 
 * <pre> -V
 *  Display std. deviations for centroids.
 * </pre>
 * 
 * <pre> -M
 *  Replace missing values with mean/mode.
 * </pre>
 * 
 * <pre> -A &lt;classname and options&gt;
 *  Distance function to use.
 *  (default: weka.core.EuclideanDistance)</pre>
 * 
 * <pre> -I &lt;num&gt;
 *  Maximum number of iterations.
 * </pre>
 * 
 * <pre> -O
 *  Preserve order of instances.
 * </pre>
 * 
 * <pre> -fast
 *  Enables faster distance calculations, using cut-off values.
 *  Disables the calculation/output of squared errors/distances.
 * </pre>
 * 
 * <pre> -S &lt;num&gt;
 *  Random number seed.
 *  (default 10)</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6450 $
 * @see RandomizableClusterer
 */
public class SimpleKMeans
  extends RandomizableClusterer 
  implements NumberOfClustersRequestable, WeightedInstancesHandler {

  /** for serialization. */
  static final long serialVersionUID = -3235809600124455376L;
  
  /**
   * replace missing values in training instances.
   */
  private ReplaceMissingValues m_ReplaceMissingFilter;

  /**
   * number of clusters to generate.
   */
  private int m_NumClusters = 2;

  /**
   * holds the cluster centroids.
   */
  private Instances m_ClusterCentroids;

  /**
   * Holds the standard deviations of the numeric attributes in each cluster.
   */
  private Instances m_ClusterStdDevs;

  
  /**
   * For each cluster, holds the frequency counts for the values of each 
   * nominal attribute.
   */
  private int[][][] m_ClusterNominalCounts;
  private int[][] m_ClusterMissingCounts;
  
  /**
   * Stats on the full data set for comparison purposes.
   * In case the attribute is numeric the value is the mean if is 
   * being used the Euclidian distance or the median if Manhattan distance
   * and if the attribute is nominal then it's mode is saved.
   */
  private double[] m_FullMeansOrMediansOrModes;
  private double[] m_FullStdDevs;
  private int[][] m_FullNominalCounts;
  private int[] m_FullMissingCounts;

  /**
   * Display standard deviations for numeric atts.
   */
  private boolean m_displayStdDevs;

  /**
   * Replace missing values globally?
   */
  private boolean m_dontReplaceMissing = false;

  /**
   * The number of instances in each cluster.
   */
  private int[] m_ClusterSizes;

  /**
   * Maximum number of iterations to be executed.
   */
  private int m_MaxIterations = 500;

  /**
   * Keep track of the number of iterations completed before convergence.
   */
  private int m_Iterations = 0;

  /**
   * Holds the squared errors for all clusters.
   */
  private double[] m_squaredErrors;

  /** the distance function used. */
  protected DistanceFunction m_DistanceFunction = new EuclideanDistance();

  /**
   * Preserve order of instances.
   */
  private boolean m_PreserveOrder = false;
	
  /**
   * Assignments obtained.
   */
  protected int[] m_Assignments = null;
	
  /** whether to use fast calculation of distances (using a cut-off). */
  protected boolean m_FastDistanceCalc = false;
  
  /**
   * the default constructor.
   */
  public SimpleKMeans() {
    super();
    
    m_SeedDefault = 10;
    setSeed(m_SeedDefault);
  }
  
  /**
   * Returns a string describing this clusterer.
   * @return a description of the evaluator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Cluster data using the k means algorithm. Can use either "
      + "the Euclidean distance (default) or the Manhattan distance."
      + " If the Manhattan distance is used, then centroids are computed "
      + "as the component-wise median rather than mean.";
  }

  /**
   * Returns default capabilities of the clusterer.
   *
   * @return      the capabilities of this clusterer
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();
    result.enable(Capability.NO_CLASS);

    // attributes
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    return result;
  }

  /**
   * Generates a clusterer. Has to initialize all fields of the clusterer
   * that are not being set via options.
   *
   * @param data set of instances serving as training data 
   * @throws Exception if the clusterer has not been 
   * generated successfully
   */
  public void buildClusterer(Instances data) throws Exception {

    // can clusterer handle the data?
    getCapabilities().testWithFail(data);

    m_Iterations = 0;

    m_ReplaceMissingFilter = new ReplaceMissingValues();
    Instances instances = new Instances(data);
				
    instances.setClassIndex(-1);
    if (!m_dontReplaceMissing) {
      m_ReplaceMissingFilter.setInputFormat(instances);
      instances = Filter.useFilter(instances, m_ReplaceMissingFilter);
    }

    m_FullMissingCounts = new int[instances.numAttributes()];
    if (m_displayStdDevs) {
      m_FullStdDevs = new double[instances.numAttributes()];
    }
    m_FullNominalCounts = new int[instances.numAttributes()][0];
		
    m_FullMeansOrMediansOrModes = moveCentroid(0, instances, false);
    for (int i = 0; i < instances.numAttributes(); i++) {
      m_FullMissingCounts[i] = instances.attributeStats(i).missingCount;
      if (instances.attribute(i).isNumeric()) {
        if (m_displayStdDevs) {
          m_FullStdDevs[i] = Math.sqrt(instances.variance(i));
        }
        if (m_FullMissingCounts[i] == instances.numInstances()) {
          m_FullMeansOrMediansOrModes[i] = Double.NaN; // mark missing as mean
        }
      } else {
        m_FullNominalCounts[i] = instances.attributeStats(i).nominalCounts;
        if (m_FullMissingCounts[i] 
            > m_FullNominalCounts[i][Utils.maxIndex(m_FullNominalCounts[i])]) {
          m_FullMeansOrMediansOrModes[i] = -1; // mark missing as most common value
        }
      }
    }

    m_ClusterCentroids = new Instances(instances, m_NumClusters);
    int[] clusterAssignments = new int [instances.numInstances()];

    if (m_PreserveOrder)
      m_Assignments = clusterAssignments;
		
    m_DistanceFunction.setInstances(instances);
    
    Random RandomO = new Random(getSeed());
    int instIndex;
    HashMap initC = new HashMap();
    DecisionTableHashKey hk = null;

    Instances initInstances = null;
    if (m_PreserveOrder)
      initInstances = new Instances(instances);
    else
      initInstances = instances;
		
    for (int j = initInstances.numInstances() - 1; j >= 0; j--) {
      instIndex = RandomO.nextInt(j+1);
      hk = new DecisionTableHashKey(initInstances.instance(instIndex),
                                    initInstances.numAttributes(), true);
      if (!initC.containsKey(hk)) {
        m_ClusterCentroids.add(initInstances.instance(instIndex));
	initC.put(hk, null);
      }
      initInstances.swap(j, instIndex);
      
      if (m_ClusterCentroids.numInstances() == m_NumClusters) {
	break;
      }
    }

    m_NumClusters = m_ClusterCentroids.numInstances();
    
    //removing reference
    initInstances = null;
		
    int i;
    boolean converged = false;
    int emptyClusterCount;
    Instances[] tempI = new Instances[m_NumClusters];
    m_squaredErrors = new double [m_NumClusters];
    m_ClusterNominalCounts = new int [m_NumClusters][instances.numAttributes()][0];
    m_ClusterMissingCounts = new int[m_NumClusters][instances.numAttributes()];
    while (!converged) {
      emptyClusterCount = 0;
      m_Iterations++;
      converged = true;
      for (i = 0; i < instances.numInstances(); i++) {
	Instance toCluster = instances.instance(i);
	int newC = clusterProcessedInstance(toCluster, false, true);
	if (newC != clusterAssignments[i]) {
	  converged = false;
	}
	clusterAssignments[i] = newC;
      }
      
      // update centroids
      m_ClusterCentroids = new Instances(instances, m_NumClusters);
      for (i = 0; i < m_NumClusters; i++) {
	tempI[i] = new Instances(instances, 0);
      }
      for (i = 0; i < instances.numInstances(); i++) {
	tempI[clusterAssignments[i]].add(instances.instance(i));
      }
      for (i = 0; i < m_NumClusters; i++) {
	if (tempI[i].numInstances() == 0) {
	  // empty cluster
	  emptyClusterCount++;
	} else {
          moveCentroid( i, tempI[i], true  );					
	}
      }

      if (emptyClusterCount > 0) {
	m_NumClusters -= emptyClusterCount;
        if (converged) {
          Instances[] t = new Instances[m_NumClusters];
          int index = 0;
          for (int k = 0; k < tempI.length; k++) {
            if (tempI[k].numInstances() > 0) {
              t[index++] = tempI[k];
            }
          }
          tempI = t;
        } else {
          tempI = new Instances[m_NumClusters];
        }
      }
			
      if (m_Iterations == m_MaxIterations)
        converged = true;
			
      if (!converged) {
	m_ClusterNominalCounts = new int [m_NumClusters][instances.numAttributes()][0];
      }
    }

    // calculate errors
    if (!m_FastDistanceCalc) {
      for (i = 0; i < instances.numInstances(); i++) {
	clusterProcessedInstance(instances.instance(i), true, false);
      }
    }
    
    if (m_displayStdDevs) {
      m_ClusterStdDevs = new Instances(instances, m_NumClusters);
    }
    m_ClusterSizes = new int [m_NumClusters];
    for (i = 0; i < m_NumClusters; i++) {
      if (m_displayStdDevs) {
        double[] vals2 = new double[instances.numAttributes()];
        for (int j = 0; j < instances.numAttributes(); j++) {
          if (instances.attribute(j).isNumeric()) {
            vals2[j] = Math.sqrt(tempI[i].variance(j));
          } else {
            vals2[j] = Utils.missingValue();
          }	
        }    
        m_ClusterStdDevs.add(new DenseInstance(1.0, vals2));
      }
      m_ClusterSizes[i] = tempI[i].numInstances();
    }
  }

  /**
   * Move the centroid to it's new coordinates. Generate the centroid coordinates based 
   * on it's  members (objects assigned to the cluster of the centroid) and the distance 
   * function being used.
   * @param centroidIndex index of the centroid which the coordinates will be computed
   * @param members the objects that are assigned to the cluster of this centroid
   * @param updateClusterInfo if the method is supposed to update the m_Cluster arrays
   * @return the centroid coordinates
   */
  protected double[] moveCentroid(int centroidIndex, Instances members, boolean updateClusterInfo) {
    double[] vals = new double[members.numAttributes()];
		
    //used only for Manhattan Distance
    Instances sortedMembers = null;
    int middle = 0;
    boolean dataIsEven = false;
		
    if (m_DistanceFunction instanceof ManhattanDistance) {
      middle = (members.numInstances()-1)/2;
      dataIsEven = ((members.numInstances()%2)==0);
      if (m_PreserveOrder) {
        sortedMembers = members;
      }else{
        sortedMembers = new Instances(members);
      }
    }
		
    for (int j = 0; j < members.numAttributes(); j++) {						
			
      //in case of Euclidian distance the centroid is the mean point
      //in case of Manhattan distance the centroid is the median point
      //in both cases, if the attribute is nominal, the centroid is the mode
      if (m_DistanceFunction instanceof EuclideanDistance ||
         members.attribute(j).isNominal())
        {													
          vals[j] = members.meanOrMode(j);
        }else if (m_DistanceFunction instanceof ManhattanDistance) {
        //singleton special case
        if (members.numInstances() == 1) {
          vals[j] = members.instance(0).value(j);
        }else{
          sortedMembers.kthSmallestValue(j, middle+1);
          vals[j] = sortedMembers.instance(middle).value(j);
          if ( dataIsEven ) {						
            sortedMembers.kthSmallestValue(j, middle+2);						
            vals[j] = (vals[j]+sortedMembers.instance(middle+1).value(j))/2;
          }
        }
      }	
			
      if (updateClusterInfo) {
        m_ClusterMissingCounts[centroidIndex][j] = members.attributeStats(j).missingCount;
        m_ClusterNominalCounts[centroidIndex][j] = members.attributeStats(j).nominalCounts;
        if (members.attribute(j).isNominal()) {
          if (m_ClusterMissingCounts[centroidIndex][j] >  
              m_ClusterNominalCounts[centroidIndex][j][Utils.maxIndex(m_ClusterNominalCounts[centroidIndex][j])]) 
            {
              vals[j] = Utils.missingValue(); // mark mode as missing
            }
        } else {
          if (m_ClusterMissingCounts[centroidIndex][j] == members.numInstances()) {
            vals[j] = Utils.missingValue(); // mark mean as missing
          }
        }
      }
    }
    if (updateClusterInfo)
      m_ClusterCentroids.add(new DenseInstance(1.0, vals));
    return vals;
  }
	
  /**
   * clusters an instance that has been through the filters.
   *
   * @param instance the instance to assign a cluster to
   * @param updateErrors if true, update the within clusters sum of errors
   * @param useFastDistCalc whether to use the fast distance calculation or not
   * @return a cluster number
   */
  private int clusterProcessedInstance(Instance instance, boolean updateErrors, boolean useFastDistCalc) {
    double minDist = Integer.MAX_VALUE;
    int bestCluster = 0;
    for (int i = 0; i < m_NumClusters; i++) {
      double dist;
      if (useFastDistCalc)
	dist = m_DistanceFunction.distance(instance, m_ClusterCentroids.instance(i), minDist);
      else
	dist = m_DistanceFunction.distance(instance, m_ClusterCentroids.instance(i));
      if (dist < minDist) {
	minDist = dist;
	bestCluster = i;
      }
    }
    if (updateErrors) {
      if (m_DistanceFunction instanceof EuclideanDistance) {
        //Euclidean distance to Squared Euclidean distance
        minDist *= minDist;
      }
      m_squaredErrors[bestCluster] += minDist;
    }
    return bestCluster;
  }

  /**
   * Classifies a given instance.
   *
   * @param instance the instance to be assigned to a cluster
   * @return the number of the assigned cluster as an interger
   * if the class is enumerated, otherwise the predicted value
   * @throws Exception if instance could not be classified
   * successfully
   */
  public int clusterInstance(Instance instance) throws Exception {
    Instance inst = null;
    if (!m_dontReplaceMissing) {
      m_ReplaceMissingFilter.input(instance);
      m_ReplaceMissingFilter.batchFinished();
      inst = m_ReplaceMissingFilter.output();
    } else {
      inst = instance;
    }

    return clusterProcessedInstance(inst, false, true);
  }

  /**
   * Returns the number of clusters.
   *
   * @return the number of clusters generated for a training dataset.
   * @throws Exception if number of clusters could not be returned
   * successfully
   */
  public int numberOfClusters() throws Exception {
    return m_NumClusters;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions () {
    Vector result = new Vector();

    result.addElement(new Option(
                                 "\tnumber of clusters.\n"
                                 + "\t(default 2).", 
                                 "N", 1, "-N <num>"));
    result.addElement(new Option(
                                 "\tDisplay std. deviations for centroids.\n", 
                                 "V", 0, "-V"));
    result.addElement(new Option(
                                 "\tReplace missing values with mean/mode.\n", 
                                 "M", 0, "-M"));

    result.add(new Option(
                          "\tDistance function to use.\n"
                          + "\t(default: weka.core.EuclideanDistance)",
                          "A", 1,"-A <classname and options>"));
		
    result.add(new Option(
                          "\tMaximum number of iterations.\n",
                          "I",1,"-I <num>"));

    result.addElement(new Option(
                                 "\tPreserve order of instances.\n", 
                                 "O", 0, "-O"));

    result.addElement(new Option(
                                 "\tEnables faster distance calculations, using cut-off values.\n"
                                 + "\tDisables the calculation/output of squared errors/distances.\n", 
                                 "fast", 0, "-fast"));
		
    Enumeration en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());

    return  result.elements();
  }

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numClustersTipText() {
    return "set number of clusters";
  }

  /**
   * set the number of clusters to generate.
   *
   * @param n the number of clusters to generate
   * @throws Exception if number of clusters is negative
   */
  public void setNumClusters(int n) throws Exception {
    if (n <= 0) {
      throw new Exception("Number of clusters must be > 0");
    }
    m_NumClusters = n;
  }

  /**
   * gets the number of clusters to generate.
   *
   * @return the number of clusters to generate
   */
  public int getNumClusters() {
    return m_NumClusters;
  }

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String maxIterationsTipText() {
    return "set maximum number of iterations";
  }

  /**
   * set the maximum number of iterations to be executed.
   *
   * @param n the maximum number of iterations
   * @throws Exception if maximum number of iteration is smaller than 1
   */
  public void setMaxIterations(int n) throws Exception {
    if (n <= 0) {
      throw new Exception("Maximum number of iterations must be > 0");
    }
    m_MaxIterations = n;
  }

  /**
   * gets the number of maximum iterations to be executed.
   *
   * @return the number of clusters to generate
   */
  public int getMaxIterations() {
    return m_MaxIterations;
  }
	

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String displayStdDevsTipText() {
    return "Display std deviations of numeric attributes "
      + "and counts of nominal attributes.";
  }

  /**
   * Sets whether standard deviations and nominal count.
   * Should be displayed in the clustering output.
   *
   * @param stdD true if std. devs and counts should be 
   * displayed
   */
  public void setDisplayStdDevs(boolean stdD) {
    m_displayStdDevs = stdD;
  }

  /**
   * Gets whether standard deviations and nominal count.
   * Should be displayed in the clustering output.
   *
   * @return true if std. devs and counts should be 
   * displayed
   */
  public boolean getDisplayStdDevs() {
    return m_displayStdDevs;
  }

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String dontReplaceMissingValuesTipText() {
    return "Replace missing values globally with mean/mode.";
  }

  /**
   * Sets whether missing values are to be replaced.
   *
   * @param r true if missing values are to be
   * replaced
   */
  public void setDontReplaceMissingValues(boolean r) {
    m_dontReplaceMissing = r;
  }

  /**
   * Gets whether missing values are to be replaced.
   *
   * @return true if missing values are to be
   * replaced
   */
  public boolean getDontReplaceMissingValues() {
    return m_dontReplaceMissing;
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   *         		displaying in the explorer/experimenter gui
   */
  public String distanceFunctionTipText() {
    return "The distance function to use for instances comparison " +
      "(default: weka.core.EuclideanDistance). ";
  }

  /**
   * returns the distance function currently in use.
   * 
   * @return the distance function
   */
  public DistanceFunction getDistanceFunction() {
    return m_DistanceFunction;
  }

  /**
   * sets the distance function to use for instance comparison.
   * 
   * @param df the new distance function to use
   * @throws Exception if instances cannot be processed
   */
  public void setDistanceFunction(DistanceFunction df) throws Exception {
    if (!(df instanceof EuclideanDistance) && 
	!(df instanceof ManhattanDistance))      {
      throw new Exception("SimpleKMeans currently only supports the Euclidean and Manhattan distances.");
    }
    m_DistanceFunction = df;
  }	

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String preserveInstancesOrderTipText() {
    return "Preserve order of instances.";
  }

  /**
   * Sets whether order of instances must be preserved.
   *
   * @param r true if missing values are to be
   * replaced
   */
  public void setPreserveInstancesOrder(boolean r) {
    m_PreserveOrder = r;
  }

  /**
   * Gets whether order of instances must be preserved.
   *
   * @return true if missing values are to be
   * replaced
   */
  public boolean getPreserveInstancesOrder() {
    return m_PreserveOrder;
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String fastDistanceCalcTipText() {
    return 
        "Uses cut-off values for speeding up distance calculation, but "
      + "suppresses also the calculation and output of the within cluster sum "
      + "of squared errors/sum of distances.";
  }

  /**
   * Sets whether to use faster distance calculation.
   *
   * @param value 	true if faster calculation to be used
   */
  public void setFastDistanceCalc(boolean value) {
    m_FastDistanceCalc = value;
  }

  /**
   * Gets whether to use faster distance calculation.
   *
   * @return 		true if faster calculation is used
   */
  public boolean getFastDistanceCalc() {
    return m_FastDistanceCalc;
  }
	
  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -N &lt;num&gt;
   *  number of clusters.
   *  (default 2).</pre>
   * 
   * <pre> -V
   *  Display std. deviations for centroids.
   * </pre>
   * 
   * <pre> -M
   *  Replace missing values with mean/mode.
   * </pre>
   * 
   * <pre> -A &lt;classname and options&gt;
   *  Distance function to use.
   *  (default: weka.core.EuclideanDistance)</pre>
   * 
   * <pre> -I &lt;num&gt;
   *  Maximum number of iterations.
   * </pre>
   * 
   * <pre> -O
   *  Preserve order of instances.
   * </pre>
   * 
   * <pre> -fast
   *  Enables faster distance calculations, using cut-off values.
   *  Disables the calculation/output of squared errors/distances.
   * </pre>
   * 
   * <pre> -S &lt;num&gt;
   *  Random number seed.
   *  (default 10)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions (String[] options)
    throws Exception {

    m_displayStdDevs = Utils.getFlag("V", options);
    m_dontReplaceMissing = Utils.getFlag("M", options);

    String optionString = Utils.getOption('N', options);

    if (optionString.length() != 0) {
      setNumClusters(Integer.parseInt(optionString));
    }
    
    optionString = Utils.getOption("I", options);
    if (optionString.length() != 0) {
      setMaxIterations(Integer.parseInt(optionString));
    }
		
    String distFunctionClass = Utils.getOption('A', options);
    if (distFunctionClass.length() != 0) {
      String distFunctionClassSpec[] = Utils.splitOptions(distFunctionClass);
      if (distFunctionClassSpec.length == 0) { 
        throw new Exception("Invalid DistanceFunction specification string."); 
      }
      String className = distFunctionClassSpec[0];
      distFunctionClassSpec[0] = "";

      setDistanceFunction( (DistanceFunction)
                           Utils.forName( DistanceFunction.class, 
                                          className, distFunctionClassSpec) );
    }
    else {
      setDistanceFunction(new EuclideanDistance());
    }
		
    m_PreserveOrder = Utils.getFlag("O", options);

    m_FastDistanceCalc = Utils.getFlag("fast", options);
    
    super.setOptions(options);
  }

  /**
   * Gets the current settings of SimpleKMeans.
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String[] getOptions () {
    int       	i;
    Vector    	result;
    String[]  	options;

    result = new Vector();

    if (m_displayStdDevs) {
      result.add("-V");
    }

    if (m_dontReplaceMissing) {
      result.add("-M");
    }

    result.add("-N");
    result.add("" + getNumClusters());

    result.add("-A");
    result.add((m_DistanceFunction.getClass().getName() + " " +
                Utils.joinOptions(m_DistanceFunction.getOptions())).trim());
		
    result.add("-I");
    result.add(""+ getMaxIterations());

    if (m_PreserveOrder) {
      result.add("-O");
    }
		
    if (m_FastDistanceCalc) {
      result.add("-fast");
    }
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * return a string describing this clusterer.
   *
   * @return a description of the clusterer as a string
   */
  public String toString() {
    if (m_ClusterCentroids == null) {
      return "No clusterer built yet!";
    }

    int maxWidth = 0;
    int maxAttWidth = 0;
    boolean containsNumeric = false;
    for (int i = 0; i < m_NumClusters; i++) {
      for (int j = 0 ;j < m_ClusterCentroids.numAttributes(); j++) {
        if (m_ClusterCentroids.attribute(j).name().length() > maxAttWidth) {
          maxAttWidth = m_ClusterCentroids.attribute(j).name().length();
        }
	if (m_ClusterCentroids.attribute(j).isNumeric()) {
          containsNumeric = true;
	  double width = Math.log(Math.abs(m_ClusterCentroids.instance(i).value(j))) /
	    Math.log(10.0);
          //          System.err.println(m_ClusterCentroids.instance(i).value(j)+" "+width);
          if (width < 0) {
            width = 1;
          }
          // decimal + # decimal places + 1
	  width += 6.0;
	  if ((int)width > maxWidth) {
	    maxWidth = (int)width;
	  }
	}
      }
    }

    for (int i = 0; i < m_ClusterCentroids.numAttributes(); i++) {
      if (m_ClusterCentroids.attribute(i).isNominal()) {
        Attribute a = m_ClusterCentroids.attribute(i);
        for (int j = 0; j < m_ClusterCentroids.numInstances(); j++) {
          String val = a.value((int)m_ClusterCentroids.instance(j).value(i));
          if (val.length() > maxWidth) {
            maxWidth = val.length();
          }
        }
        for (int j = 0; j < a.numValues(); j++) {
          String val = a.value(j) + " ";
          if (val.length() > maxAttWidth) {
            maxAttWidth = val.length();
          }
        }
      }
    }

    if (m_displayStdDevs) {
      // check for maximum width of maximum frequency count
      for (int i = 0; i < m_ClusterCentroids.numAttributes(); i++) {
        if (m_ClusterCentroids.attribute(i).isNominal()) {
          int maxV = Utils.maxIndex(m_FullNominalCounts[i]);
          /*          int percent = (int)((double)m_FullNominalCounts[i][maxV] /
                      Utils.sum(m_ClusterSizes) * 100.0); */
          int percent = 6; // max percent width (100%)
          String nomV = "" + m_FullNominalCounts[i][maxV];
          //            + " (" + percent + "%)";
          if (nomV.length() + percent > maxWidth) {
            maxWidth = nomV.length() + 1;
          }
        }
      }
    }

    // check for size of cluster sizes
    for (int i = 0; i < m_ClusterSizes.length; i++) {
      String size = "(" + m_ClusterSizes[i] + ")";
      if (size.length() > maxWidth) {
        maxWidth = size.length();
      }
    }
    
    if (m_displayStdDevs && maxAttWidth < "missing".length()) {
      maxAttWidth = "missing".length();
    }
    
    String plusMinus = "+/-";
    maxAttWidth += 2;
    if (m_displayStdDevs && containsNumeric) {
      maxWidth += plusMinus.length();
    }
    if (maxAttWidth < "Attribute".length() + 2) {
      maxAttWidth = "Attribute".length() + 2;
    }

    if (maxWidth < "Full Data".length()) {
      maxWidth = "Full Data".length() + 1;
    }

    if (maxWidth < "missing".length()) {
      maxWidth = "missing".length() + 1;
    }
    
    StringBuffer temp = new StringBuffer();
    temp.append("\nkMeans\n======\n");
    temp.append("\nNumber of iterations: " + m_Iterations);
		
    if (!m_FastDistanceCalc) {
      temp.append("\n");
      if (m_DistanceFunction instanceof EuclideanDistance) {
	temp.append("Within cluster sum of squared errors: " + Utils.sum(m_squaredErrors));
      }else{
	temp.append("Sum of within cluster distances: " + Utils.sum(m_squaredErrors));
      }
    }
		
    if (!m_dontReplaceMissing) {
      temp.append("\nMissing values globally replaced with mean/mode");
    }

    temp.append("\n\nCluster centroids:\n");
    temp.append(pad("Cluster#", " ", (maxAttWidth + (maxWidth * 2 + 2)) - "Cluster#".length(), true));

    temp.append("\n");
    temp.append(pad("Attribute", " ", maxAttWidth - "Attribute".length(), false));

    
    temp.append(pad("Full Data", " ", maxWidth + 1 - "Full Data".length(), true));

    // cluster numbers
    for (int i = 0; i < m_NumClusters; i++) {
      String clustNum = "" + i;
      temp.append(pad(clustNum, " ", maxWidth + 1 - clustNum.length(), true));
    }
    temp.append("\n");

    // cluster sizes
    String cSize = "(" + Utils.sum(m_ClusterSizes) + ")";
    temp.append(pad(cSize, " ", maxAttWidth + maxWidth + 1 - cSize.length(), true));
    for (int i = 0; i < m_NumClusters; i++) {
      cSize = "(" + m_ClusterSizes[i] + ")";
      temp.append(pad(cSize, " ",maxWidth + 1 - cSize.length(), true));
    }
    temp.append("\n");

    temp.append(pad("", "=", maxAttWidth + 
                    (maxWidth * (m_ClusterCentroids.numInstances()+1) 
                     + m_ClusterCentroids.numInstances() + 1), true));
    temp.append("\n");

    for (int i = 0; i < m_ClusterCentroids.numAttributes(); i++) {
      String attName = m_ClusterCentroids.attribute(i).name();
      temp.append(attName);
      for (int j = 0; j < maxAttWidth - attName.length(); j++) {
        temp.append(" ");
      }

      String strVal;
      String valMeanMode;
      // full data
      if (m_ClusterCentroids.attribute(i).isNominal()) {
        if (m_FullMeansOrMediansOrModes[i] == -1) { // missing
          valMeanMode = pad("missing", " ", maxWidth + 1 - "missing".length(), true);
        } else {
          valMeanMode = 
            pad((strVal = m_ClusterCentroids.attribute(i).value((int)m_FullMeansOrMediansOrModes[i])),
                " ", maxWidth + 1 - strVal.length(), true);
        }
      } else {
        if (Double.isNaN(m_FullMeansOrMediansOrModes[i])) {
          valMeanMode = pad("missing", " ", maxWidth + 1 - "missing".length(), true);
        } else {
          valMeanMode =  pad((strVal = Utils.doubleToString(m_FullMeansOrMediansOrModes[i],
                                                            maxWidth,4).trim()), 
                             " ", maxWidth + 1 - strVal.length(), true);
        }
      }
      temp.append(valMeanMode);

      for (int j = 0; j < m_NumClusters; j++) {
        if (m_ClusterCentroids.attribute(i).isNominal()) {
          if (m_ClusterCentroids.instance(j).isMissing(i)) {
            valMeanMode = pad("missing", " ", maxWidth + 1 - "missing".length(), true);
          } else {
            valMeanMode = 
              pad((strVal = m_ClusterCentroids.attribute(i).value((int)m_ClusterCentroids.instance(j).value(i))),
                  " ", maxWidth + 1 - strVal.length(), true);
          }
        } else {
          if (m_ClusterCentroids.instance(j).isMissing(i)) {
            valMeanMode = pad("missing", " ", maxWidth + 1 - "missing".length(), true);
          } else {
            valMeanMode = pad((strVal = Utils.doubleToString(m_ClusterCentroids.instance(j).value(i),
                                                             maxWidth,4).trim()), 
                              " ", maxWidth + 1 - strVal.length(), true);
          }
        }
        temp.append(valMeanMode);
      }
      temp.append("\n");

      if (m_displayStdDevs) {
        // Std devs/max nominal
        String stdDevVal = "";

        if (m_ClusterCentroids.attribute(i).isNominal()) {
          // Do the values of the nominal attribute
          Attribute a = m_ClusterCentroids.attribute(i);
          for (int j = 0; j < a.numValues(); j++) {
            // full data
            String val = "  " + a.value(j);
            temp.append(pad(val, " ", maxAttWidth + 1 - val.length(), false));
            int count = m_FullNominalCounts[i][j];
            int percent = (int)((double)m_FullNominalCounts[i][j] /
                                Utils.sum(m_ClusterSizes) * 100.0);
            String percentS = "" + percent + "%)";
            percentS = pad(percentS, " ", 5 - percentS.length(), true);
            stdDevVal = "" + count + " (" + percentS;
            stdDevVal = 
              pad(stdDevVal, " ", maxWidth + 1 - stdDevVal.length(), true);
            temp.append(stdDevVal);

            // Clusters
            for (int k = 0; k < m_NumClusters; k++) {
              count = m_ClusterNominalCounts[k][i][j];
              percent = (int)((double)m_ClusterNominalCounts[k][i][j] /
                              m_ClusterSizes[k] * 100.0);
              percentS = "" + percent + "%)";
              percentS = pad(percentS, " ", 5 - percentS.length(), true);
              stdDevVal = "" + count + " (" + percentS;
              stdDevVal = 
                pad(stdDevVal, " ", maxWidth + 1 - stdDevVal.length(), true);
              temp.append(stdDevVal);
            }
            temp.append("\n");
          }
          // missing (if any)
          if (m_FullMissingCounts[i] > 0) {
            // Full data
            temp.append(pad("  missing", " ", maxAttWidth + 1 - "  missing".length(), false));
            int count = m_FullMissingCounts[i];
            int percent = (int)((double)m_FullMissingCounts[i] /
                                Utils.sum(m_ClusterSizes) * 100.0);
            String percentS = "" + percent + "%)";
            percentS = pad(percentS, " ", 5 - percentS.length(), true);
            stdDevVal = "" + count + " (" + percentS;
            stdDevVal = 
              pad(stdDevVal, " ", maxWidth + 1 - stdDevVal.length(), true);
            temp.append(stdDevVal);
           
            // Clusters
            for (int k = 0; k < m_NumClusters; k++) {
              count = m_ClusterMissingCounts[k][i];
              percent = (int)((double)m_ClusterMissingCounts[k][i] /
                              m_ClusterSizes[k] * 100.0);
              percentS = "" + percent + "%)";
              percentS = pad(percentS, " ", 5 - percentS.length(), true);
              stdDevVal = "" + count + " (" + percentS;
              stdDevVal = 
                pad(stdDevVal, " ", maxWidth + 1 - stdDevVal.length(), true);
              temp.append(stdDevVal);
            }

            temp.append("\n");
          }

          temp.append("\n");
        } else {
          // Full data
          if (Double.isNaN(m_FullMeansOrMediansOrModes[i])) {
            stdDevVal = pad("--", " ", maxAttWidth + maxWidth + 1 - 2, true);
          } else {
            stdDevVal = pad((strVal = plusMinus 
                             + Utils.doubleToString(m_FullStdDevs[i],
                                                    maxWidth,4).trim()), 
                            " ", maxWidth + maxAttWidth + 1 - strVal.length(), true);
          }
          temp.append(stdDevVal);

          // Clusters
          for (int j = 0; j < m_NumClusters; j++) {
            if (m_ClusterCentroids.instance(j).isMissing(i)) {
              stdDevVal = pad("--", " ", maxWidth + 1 - 2, true);
            } else {
              stdDevVal = 
                pad((strVal = plusMinus 
                     + Utils.doubleToString(m_ClusterStdDevs.instance(j).value(i),
                                            maxWidth,4).trim()), 
                    " ", maxWidth + 1 - strVal.length(), true);
            }
            temp.append(stdDevVal);
          }
          temp.append("\n\n");
        }
      }
    }

    temp.append("\n\n");
    return temp.toString();
  }

  private String pad(String source, String padChar, 
                     int length, boolean leftPad) {
    StringBuffer temp = new StringBuffer();

    if (leftPad) {
      for (int i = 0; i< length; i++) {
        temp.append(padChar);
      }
      temp.append(source);
    } else {
      temp.append(source);
      for (int i = 0; i< length; i++) {
        temp.append(padChar);
      }
    }
    return temp.toString();
  }

  /**
   * Gets the the cluster centroids.
   * 
   * @return		the cluster centroids
   */
  public Instances getClusterCentroids() {
    return m_ClusterCentroids;
  }

  /**
   * Gets the standard deviations of the numeric attributes in each cluster.
   * 
   * @return		the standard deviations of the numeric attributes 
   * 			in each cluster
   */
  public Instances getClusterStandardDevs() {
    return m_ClusterStdDevs;
  }

  /**
   * Returns for each cluster the frequency counts for the values of each 
   * nominal attribute.
   * 
   * @return		the counts
   */
  public int[][][] getClusterNominalCounts() {
    return m_ClusterNominalCounts;
  }

  /**
   * Gets the squared error for all clusters.
   * 
   * @return		the squared error, NaN if fast distance calculation is
   * 			used
   * @see		#m_FastDistanceCalc
   */
  public double getSquaredError() {
    if (m_FastDistanceCalc)
      return Double.NaN;
    else
      return Utils.sum(m_squaredErrors);
  }

  /**
   * Gets the number of instances in each cluster.
   * 
   * @return		The number of instances in each cluster
   */
  public int[] getClusterSizes() {
    return m_ClusterSizes;
  }
  
  /**
   * Gets the assignments for each instance.
   * @return Array of indexes of the centroid assigned to each instance
   * @throws Exception if order of instances wasn't preserved or no assignments were made
   */
  public int[] getAssignments() throws Exception{
    if (!m_PreserveOrder) {
      throw new Exception("The assignments are only available when order of instances is preserved (-O)");
    }
    if (m_Assignments == null) {
      throw new Exception("No assignments made.");
    }
    return m_Assignments;
  }
	
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6450 $");
  }

  /**
   * Main method for executing this class.
   *
   * @param args 	use -h to list all parameters
   */
  public static void main (String[] args) {
    runClusterer(new SimpleKMeans(), args);
  }
}
