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
 *    EM.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.clusterers;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Attribute;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.estimators.DiscreteEstimator;
import weka.estimators.Estimator;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Simple EM (expectation maximisation) class.<br/>
 * <br/>
 * EM assigns a probability distribution to each instance which indicates the probability of it belonging to each of the clusters. EM can decide how many clusters to create by cross validation, or you may specify apriori how many clusters to generate.<br/>
 * <br/>
 * The cross validation performed to determine the number of clusters is done in the following steps:<br/>
 * 1. the number of clusters is set to 1<br/>
 * 2. the training set is split randomly into 10 folds.<br/>
 * 3. EM is performed 10 times using the 10 folds the usual CV way.<br/>
 * 4. the loglikelihood is averaged over all 10 results.<br/>
 * 5. if loglikelihood has increased the number of clusters is increased by 1 and the program continues at step 2. <br/>
 * <br/>
 * The number of folds is fixed to 10, as long as the number of instances in the training set is not smaller 10. If this is the case the number of folds is set equal to the number of instances.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -N &lt;num&gt;
 *  number of clusters. If omitted or -1 specified, then 
 *  cross validation is used to select the number of clusters.</pre>
 * 
 * <pre> -I &lt;num&gt;
 *  max iterations.
 * (default 100)</pre>
 * 
 * <pre> -V
 *  verbose.</pre>
 * 
 * <pre> -M &lt;num&gt;
 *  minimum allowable standard deviation for normal density
 *  computation
 *  (default 1e-6)</pre>
 * 
 * <pre> -O
 *  Display model in old format (good when there are many clusters)
 * </pre>
 * 
 * <pre> -S &lt;num&gt;
 *  Random number seed.
 *  (default 100)</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6298 $
 */
public class EM
  extends RandomizableDensityBasedClusterer
  implements NumberOfClustersRequestable, WeightedInstancesHandler {

  /** for serialization */
  static final long serialVersionUID = 8348181483812829475L;
  
  /** hold the discrete estimators for each cluster */
  private Estimator m_model[][];

  /** hold the normal estimators for each cluster */
  private double m_modelNormal[][][];

  /** default minimum standard deviation */
  private double m_minStdDev = 1e-6;

  private double [] m_minStdDevPerAtt;

  /** hold the weights of each instance for each cluster */
  private double m_weights[][];

  /** the prior probabilities for clusters */
  private double m_priors[];

  /** the loglikelihood of the data */
  private double m_loglikely;

  /** training instances */
  private Instances m_theInstances = null;

  /** number of clusters selected by the user or cross validation */
  private int m_num_clusters;

  /** the initial number of clusters requested by the user--- -1 if
      xval is to be used to find the number of clusters */
  private int m_initialNumClusters;

  /** number of attributes */
  private int m_num_attribs;

  /** number of training instances */
  private int m_num_instances;

  /** maximum iterations to perform */
  private int m_max_iterations;

  /** attribute min values */
  private double [] m_minValues;

  /** attribute max values */
  private double [] m_maxValues;

  /** random number generator */
  private Random m_rr;

  /** Verbose? */
  private boolean m_verbose;

 /** globally replace missing values */
  private ReplaceMissingValues m_replaceMissing;

  /** display model output in old-style format */
  private boolean m_displayModelInOldFormat;

  /**
   * Returns a string describing this clusterer
   * @return a description of the evaluator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return
        "Simple EM (expectation maximisation) class.\n\n"
      + "EM assigns a probability distribution to each instance which "
      + "indicates the probability of it belonging to each of the clusters. "
      + "EM can decide how many clusters to create by cross validation, or you "
      + "may specify apriori how many clusters to generate.\n\n"
      + "The cross validation performed to determine the number of clusters "
      + "is done in the following steps:\n"
      + "1. the number of clusters is set to 1\n"
      + "2. the training set is split randomly into 10 folds.\n"
      + "3. EM is performed 10 times using the 10 folds the usual CV way.\n"
      + "4. the loglikelihood is averaged over all 10 results.\n"
      + "5. if loglikelihood has increased the number of clusters is increased "
      + "by 1 and the program continues at step 2. \n\n"
      + "The number of folds is fixed to 10, as long as the number of "
      + "instances in the training set is not smaller 10. If this is the case "
      + "the number of folds is set equal to the number of instances.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions () {
    Vector result = new Vector();
    
    result.addElement(new Option(
	"\tnumber of clusters. If omitted or -1 specified, then \n"
	+ "\tcross validation is used to select the number of clusters.", 
	"N", 1, "-N <num>"));

    result.addElement(new Option(
	"\tmax iterations."
	+ "\n(default 100)", 
	"I", 1, "-I <num>"));
    
    result.addElement(new Option(
	"\tverbose.",
	"V", 0, "-V"));
    
    result.addElement(new Option(
	"\tminimum allowable standard deviation for normal density\n"
	+ "\tcomputation\n"
	+ "\t(default 1e-6)",
	"M",1,"-M <num>"));

    result.addElement(
              new Option("\tDisplay model in old format (good when there are "
                         + "many clusters)\n",
                         "O", 0, "-O"));

    Enumeration en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());
    
    return  result.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -N &lt;num&gt;
   *  number of clusters. If omitted or -1 specified, then 
   *  cross validation is used to select the number of clusters.</pre>
   * 
   * <pre> -I &lt;num&gt;
   *  max iterations.
   * (default 100)</pre>
   * 
   * <pre> -V
   *  verbose.</pre>
   * 
   * <pre> -M &lt;num&gt;
   *  minimum allowable standard deviation for normal density
   *  computation
   *  (default 1e-6)</pre>
   * 
   * <pre> -O
   *  Display model in old format (good when there are many clusters)
   * </pre>
   * 
   * <pre> -S &lt;num&gt;
   *  Random number seed.
   *  (default 100)</pre>
   * 
   <!-- options-end -->
   * 
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions (String[] options)
    throws Exception {
    resetOptions();
    setDebug(Utils.getFlag('V', options));
    String optionString = Utils.getOption('I', options);

    if (optionString.length() != 0) {
      setMaxIterations(Integer.parseInt(optionString));
    }

    optionString = Utils.getOption('N', options);
    if (optionString.length() != 0) {
      setNumClusters(Integer.parseInt(optionString));
    }

    optionString = Utils.getOption('M', options);
    if (optionString.length() != 0) {
      setMinStdDev((new Double(optionString)).doubleValue());
    }

    setDisplayModelInOldFormat(Utils.getFlag('O', options));
    
    super.setOptions(options);
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String displayModelInOldFormatTipText() {
    return "Use old format for model output. The old format is "
      + "better when there are many clusters. The new format "
      + "is better when there are fewer clusters and many attributes.";
  }

  /**
   * Set whether to display model output in the old, original
   * format.
   *
   * @param d true if model ouput is to be shown in the old format
   */
  public void setDisplayModelInOldFormat(boolean d) {
    m_displayModelInOldFormat = d;
  }

  /**
   * Get whether to display model output in the old, original
   * format.
   *
   * @return true if model ouput is to be shown in the old format
   */
  public boolean getDisplayModelInOldFormat() {
    return m_displayModelInOldFormat;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String minStdDevTipText() {
    return "set minimum allowable standard deviation";
  }

  /**
   * Set the minimum value for standard deviation when calculating
   * normal density. Reducing this value can help prevent arithmetic
   * overflow resulting from multiplying large densities (arising from small
   * standard deviations) when there are many singleton or near singleton
   * values.
   * @param m minimum value for standard deviation
   */
  public void setMinStdDev(double m) {
    m_minStdDev = m;
  }

  public void setMinStdDevPerAtt(double [] m) {
    m_minStdDevPerAtt = m;
  }

  /**
   * Get the minimum allowable standard deviation.
   * @return the minumum allowable standard deviation
   */
  public double getMinStdDev() {
    return m_minStdDev;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numClustersTipText() {
    return "set number of clusters. -1 to select number of clusters "
      +"automatically by cross validation.";
  }

  /**
   * Set the number of clusters (-1 to select by CV).
   *
   * @param n the number of clusters
   * @throws Exception if n is 0
   */
  public void setNumClusters (int n)
    throws Exception {
    
    if (n == 0) {
      throw  new Exception("Number of clusters must be > 0. (or -1 to " 
			   + "select by cross validation).");
    }

    if (n < 0) {
      m_num_clusters = -1;
      m_initialNumClusters = -1;
    }
    else {
      m_num_clusters = n;
      m_initialNumClusters = n;
    }
  }


  /**
   * Get the number of clusters
   *
   * @return the number of clusters.
   */
  public int getNumClusters () {
    return  m_initialNumClusters;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String maxIterationsTipText() {
    return "maximum number of iterations";
  }

  /**
   * Set the maximum number of iterations to perform
   *
   * @param i the number of iterations
   * @throws Exception if i is less than 1
   */
  public void setMaxIterations (int i)
    throws Exception {
    if (i < 1) {
      throw  new Exception("Maximum number of iterations must be > 0!");
    }

    m_max_iterations = i;
  }


  /**
   * Get the maximum number of iterations
   *
   * @return the number of iterations
   */
  public int getMaxIterations () {
    return  m_max_iterations;
  }

  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String debugTipText() {
    return "If set to true, clusterer may output additional info to " +
      "the console.";
  }


  /**
   * Set debug mode - verbose output
   *
   * @param v true for verbose output
   */
  public void setDebug (boolean v) {
    m_verbose = v;
  }


  /**
   * Get debug mode
   *
   * @return true if debug mode is set
   */
  public boolean getDebug () {
    return  m_verbose;
  }


  /**
   * Gets the current settings of EM.
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String[] getOptions () {
    int       	i;
    Vector    	result;
    String[]  	options;

    result = new Vector();

    result.add("-I");
    result.add("" + m_max_iterations);
    result.add("-N");
    result.add("" + getNumClusters());
    result.add("-M");
    result.add("" + getMinStdDev());
    if (m_displayModelInOldFormat) {
      result.add("-O");
    }

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * Initialise estimators and storage.
   *
   * @param inst the instances
   * @throws Exception if initialization fails
   **/
  private void EM_Init (Instances inst)
    throws Exception {
    int i, j, k;

    // run k means 10 times and choose best solution
    SimpleKMeans bestK = null;
    double bestSqE = Double.MAX_VALUE;
    for (i = 0; i < 10; i++) {
      SimpleKMeans sk = new SimpleKMeans();
      sk.setSeed(m_rr.nextInt());
      sk.setNumClusters(m_num_clusters);
      sk.setDisplayStdDevs(true);
      sk.buildClusterer(inst);
      if (sk.getSquaredError() < bestSqE) {
	bestSqE = sk.getSquaredError();
	bestK = sk;
      }
    }
    
    // initialize with best k-means solution
    m_num_clusters = bestK.numberOfClusters();
    m_weights = new double[inst.numInstances()][m_num_clusters];
    m_model = new DiscreteEstimator[m_num_clusters][m_num_attribs];
    m_modelNormal = new double[m_num_clusters][m_num_attribs][3];
    m_priors = new double[m_num_clusters];
    Instances centers = bestK.getClusterCentroids();
    Instances stdD = bestK.getClusterStandardDevs();
    int [][][] nominalCounts = bestK.getClusterNominalCounts();
    int [] clusterSizes = bestK.getClusterSizes();

    for (i = 0; i < m_num_clusters; i++) {
      Instance center = centers.instance(i);
      for (j = 0; j < m_num_attribs; j++) {
	if (inst.attribute(j).isNominal()) {
	  m_model[i][j] = new DiscreteEstimator(m_theInstances.
						attribute(j).numValues()
						, true);
	  for (k = 0; k < inst.attribute(j).numValues(); k++) {
	    m_model[i][j].addValue(k, nominalCounts[i][j][k]);
	  }
	} else {
	  double minStdD = (m_minStdDevPerAtt != null)
	    ? m_minStdDevPerAtt[j]
	    : m_minStdDev;
	  double mean = (center.isMissing(j))
	    ? inst.meanOrMode(j)
	    : center.value(j);
	  m_modelNormal[i][j][0] = mean;
	  double stdv = (stdD.instance(i).isMissing(j))
	    ? ((m_maxValues[j] - m_minValues[j]) / (2 * m_num_clusters))
	    : stdD.instance(i).value(j);
	  if (stdv < minStdD) {
	    stdv = inst.attributeStats(j).numericStats.stdDev;
            if (Double.isInfinite(stdv)) {
              stdv = minStdD;
            }
	    if (stdv < minStdD) {
	      stdv = minStdD;
	    }
	  }
	  if (stdv <= 0) {
	    stdv = m_minStdDev;
	  }

	  m_modelNormal[i][j][1] = stdv;
	  m_modelNormal[i][j][2] = 1.0;
	}
      } 
    }    
    
    
    for (j = 0; j < m_num_clusters; j++) {
      //      m_priors[j] += 1.0;
      m_priors[j] = clusterSizes[j];
    }
    Utils.normalize(m_priors);
  }


  /**
   * calculate prior probabilites for the clusters
   *
   * @param inst the instances
   * @throws Exception if priors can't be calculated
   **/
  private void estimate_priors (Instances inst)
    throws Exception {

    for (int i = 0; i < m_num_clusters; i++) {
      m_priors[i] = 0.0;
    }

    for (int i = 0; i < inst.numInstances(); i++) {
      for (int j = 0; j < m_num_clusters; j++) {
        m_priors[j] += inst.instance(i).weight() * m_weights[i][j];
      }
    }

    Utils.normalize(m_priors);
  }


  /** Constant for normal distribution. */
  private static double m_normConst = Math.log(Math.sqrt(2*Math.PI));

  /**
   * Density function of normal distribution.
   * @param x input value
   * @param mean mean of distribution
   * @param stdDev standard deviation of distribution
   * @return the density
   */
  private double logNormalDens (double x, double mean, double stdDev) {

    double diff = x - mean;
    //    System.err.println("x: "+x+" mean: "+mean+" diff: "+diff+" stdv: "+stdDev);
    //    System.err.println("diff*diff/(2*stdv*stdv): "+ (diff * diff / (2 * stdDev * stdDev)));
    
    return - (diff * diff / (2 * stdDev * stdDev))  - m_normConst - Math.log(stdDev);
  }

  /**
   * New probability estimators for an iteration
   */
  private void new_estimators () {
    for (int i = 0; i < m_num_clusters; i++) {
      for (int j = 0; j < m_num_attribs; j++) {
        if (m_theInstances.attribute(j).isNominal()) {
          m_model[i][j] = new DiscreteEstimator(m_theInstances.
						attribute(j).numValues()
						, true);
        }
        else {
          m_modelNormal[i][j][0] = m_modelNormal[i][j][1] = 
	    m_modelNormal[i][j][2] = 0.0;
        }
      }
    }
  }


  /**
   * The M step of the EM algorithm.
   * @param inst the training instances
   * @throws Exception if something goes wrong
   */
  private void M (Instances inst)
    throws Exception {

    int i, j, l;

    new_estimators();
    estimate_priors(inst);

    for (i = 0; i < m_num_clusters; i++) {
      for (j = 0; j < m_num_attribs; j++) {
        for (l = 0; l < inst.numInstances(); l++) {
	  Instance in = inst.instance(l);
          if (!in.isMissing(j)) {
            if (inst.attribute(j).isNominal()) {
              m_model[i][j].addValue(in.value(j), 
				     in.weight() * m_weights[l][i]);
            }
            else {
              m_modelNormal[i][j][0] += (in.value(j) * in.weight() *
					 m_weights[l][i]);
              m_modelNormal[i][j][2] += in.weight() * m_weights[l][i];
              m_modelNormal[i][j][1] += (in.value(j) * 
					 in.value(j) * in.weight() * m_weights[l][i]);
            }
          }
        }
      }
    }
    
    // calcualte mean and std deviation for numeric attributes
    for (j = 0; j < m_num_attribs; j++) {
      if (!inst.attribute(j).isNominal()) {
        for (i = 0; i < m_num_clusters; i++) {
          if (m_modelNormal[i][j][2] <= 0) {
            m_modelNormal[i][j][1] = Double.MAX_VALUE;
	    //	    m_modelNormal[i][j][0] = 0;
	    m_modelNormal[i][j][0] = m_minStdDev;
          } else {
	      
	    // variance
	    m_modelNormal[i][j][1] = (m_modelNormal[i][j][1] - 
				      (m_modelNormal[i][j][0] * 
				       m_modelNormal[i][j][0] / 
				       m_modelNormal[i][j][2])) / 
	      (m_modelNormal[i][j][2]);
	    
	    if (m_modelNormal[i][j][1] < 0) {
	      m_modelNormal[i][j][1] = 0;
	    }
	    
	    // std dev      
	    double minStdD = (m_minStdDevPerAtt != null)
	    ? m_minStdDevPerAtt[j]
	    : m_minStdDev;

	    m_modelNormal[i][j][1] = Math.sqrt(m_modelNormal[i][j][1]);              

	    if ((m_modelNormal[i][j][1] <= minStdD)) {
	      m_modelNormal[i][j][1] = inst.attributeStats(j).numericStats.stdDev;
	      if ((m_modelNormal[i][j][1] <= minStdD)) {
		m_modelNormal[i][j][1] = minStdD;
	      }
	    }
	    if ((m_modelNormal[i][j][1] <= 0)) {
	      m_modelNormal[i][j][1] = m_minStdDev;
	    }
            if (Double.isInfinite(m_modelNormal[i][j][1])) {
              m_modelNormal[i][j][1] = m_minStdDev;
            }
	    
	    // mean
	    m_modelNormal[i][j][0] /= m_modelNormal[i][j][2];
	  }
        }
      }
    }
  }

  /**
   * The E step of the EM algorithm. Estimate cluster membership 
   * probabilities.
   *
   * @param inst the training instances
   * @param change_weights whether to change the weights
   * @return the average log likelihood
   * @throws Exception if computation fails
   */
  private double E (Instances inst, boolean change_weights)
    throws Exception {

    double loglk = 0.0, sOW = 0.0;

    for (int l = 0; l < inst.numInstances(); l++) {

      Instance in = inst.instance(l);

      loglk += in.weight() * logDensityForInstance(in);
      sOW += in.weight();

      if (change_weights) {
	m_weights[l] = distributionForInstance(in);
      }
    }
    
    // reestimate priors
    /*if (change_weights) {
      estimate_priors(inst);
    }*/
    return  loglk / sOW;
  }
  
  
  /**
   * Constructor.
   *
   **/
  public EM () {
    super();
    
    m_SeedDefault = 100;
    resetOptions();
  }


  /**
   * Reset to default options
   */
  protected void resetOptions () {
    m_minStdDev = 1e-6;
    m_max_iterations = 100;
    m_Seed = m_SeedDefault;
    m_num_clusters = -1;
    m_initialNumClusters = -1;
    m_verbose = false;
  }

  /**
   * Return the normal distributions for the cluster models
   *
   * @return a <code>double[][][]</code> value
   */
  public double [][][] getClusterModelsNumericAtts() {
    return m_modelNormal;
  }

  /**
   * Return the priors for the clusters
   *
   * @return a <code>double[]</code> value
   */
  public double [] getClusterPriors() {
    return m_priors;
  }

  /**
   * Outputs the generated clusters into a string.
   * 
   * @return the clusterer in string representation
   */
  public String toString() {
    if (m_displayModelInOldFormat) {
      return toStringOriginal();
    }

    if (m_priors == null) {
      return "No clusterer built yet!";
    }
    StringBuffer temp = new StringBuffer();
    temp.append("\nEM\n==\n");
    if (m_initialNumClusters == -1) {
      temp.append("\nNumber of clusters selected by cross validation: "
		  +m_num_clusters+"\n");
    } else {
      temp.append("\nNumber of clusters: " + m_num_clusters + "\n");
    }

    int maxWidth = 0;
    int maxAttWidth = 0;
    boolean containsKernel = false;
    
    // set up max widths
    // attributes
    for (int i = 0; i < m_num_attribs; i++) {
      Attribute a = m_theInstances.attribute(i);
      if (a.name().length() > maxAttWidth) {
        maxAttWidth = m_theInstances.attribute(i).name().length();
      }
      if (a.isNominal()) {
        // check values
        for (int j = 0; j < a.numValues(); j++) {
          String val = a.value(j) + "  ";
          if (val.length() > maxAttWidth) {
            maxAttWidth = val.length();
          }
        }
      }
    }

    for (int i = 0; i < m_num_clusters; i++) {
      for (int j = 0; j < m_num_attribs; j++) {
        if (m_theInstances.attribute(j).isNumeric()) {
          // check mean and std. dev. against maxWidth
          double mean = Math.log(Math.abs(m_modelNormal[i][j][0])) / Math.log(10.0);
          double stdD = Math.log(Math.abs(m_modelNormal[i][j][1])) / Math.log(10.0);
          double width = (mean > stdD)
            ? mean
            : stdD;
          if (width < 0) {
            width = 1;
          }
          // decimal + # decimal places + 1
          width += 6.0;
          if ((int)width > maxWidth) {
            maxWidth = (int)width;
          }
        } else {
          // nominal distributions
          DiscreteEstimator d = (DiscreteEstimator)m_model[i][j];
          for (int k = 0; k < d.getNumSymbols(); k++) {
            String size = Utils.doubleToString(d.getCount(k), maxWidth, 4).trim();
            if (size.length() > maxWidth) {
              maxWidth = size.length();
            }
          }
          int sum = 
            Utils.doubleToString(d.getSumOfCounts(), maxWidth, 4).trim().length();
          if (sum > maxWidth) {
            maxWidth = sum;
          }
        }
      }
    }

    if (maxAttWidth < "Attribute".length()) {
      maxAttWidth = "Attribute".length();
    }    
    
    maxAttWidth += 2;

    temp.append("\n\n");
    temp.append(pad("Cluster", " ", 
                    (maxAttWidth + maxWidth + 1) - "Cluster".length(), 
                    true));
    
    temp.append("\n");
    temp.append(pad("Attribute", " ", maxAttWidth - "Attribute".length(), false));

    // cluster #'s
    for (int i = 0; i < m_num_clusters; i++) {
      String classL = "" + i;
      temp.append(pad(classL, " ", maxWidth + 1 - classL.length(), true));
    }
    temp.append("\n");

    // cluster priors
    temp.append(pad("", " ", maxAttWidth, true));
    for (int i = 0; i < m_num_clusters; i++) {
      String priorP = Utils.doubleToString(m_priors[i], maxWidth, 2).trim();
      priorP = "(" + priorP + ")";
      temp.append(pad(priorP, " ", maxWidth + 1 - priorP.length(), true));
    }

    temp.append("\n");
    temp.append(pad("", "=", maxAttWidth + 
                    (maxWidth * m_num_clusters)
                    + m_num_clusters + 1, true));
    temp.append("\n");

    for (int i = 0; i < m_num_attribs; i++) {
      String attName = m_theInstances.attribute(i).name();
      temp.append(attName + "\n");

      if (m_theInstances.attribute(i).isNumeric()) {
        String meanL = "  mean";
        temp.append(pad(meanL, " ", maxAttWidth + 1 - meanL.length(), false));
        for (int j = 0; j < m_num_clusters; j++) {
          // means
          String mean = 
            Utils.doubleToString(m_modelNormal[j][i][0], maxWidth, 4).trim();
          temp.append(pad(mean, " ", maxWidth + 1 - mean.length(), true));
        }
        temp.append("\n");            
        // now do std deviations
        String stdDevL = "  std. dev.";
        temp.append(pad(stdDevL, " ", maxAttWidth + 1 - stdDevL.length(), false));
        for (int j = 0; j < m_num_clusters; j++) {
          String stdDev = 
            Utils.doubleToString(m_modelNormal[j][i][1], maxWidth, 4).trim();
          temp.append(pad(stdDev, " ", maxWidth + 1 - stdDev.length(), true));
        }
        temp.append("\n\n");
      } else {
        Attribute a = m_theInstances.attribute(i);
        for (int j = 0; j < a.numValues(); j++) {
          String val = "  " + a.value(j);
          temp.append(pad(val, " ", maxAttWidth + 1 - val.length(), false));
          for (int k = 0; k < m_num_clusters; k++) {
            DiscreteEstimator d = (DiscreteEstimator)m_model[k][i];
            String count = Utils.doubleToString(d.getCount(j), maxWidth, 4).trim();
            temp.append(pad(count, " ", maxWidth + 1 - count.length(), true));
          }
          temp.append("\n");
        }
        // do the totals
        String total = "  [total]";
        temp.append(pad(total, " ", maxAttWidth + 1 - total.length(), false));
        for (int k = 0; k < m_num_clusters; k++) {
          DiscreteEstimator d = (DiscreteEstimator)m_model[k][i];
          String count = 
            Utils.doubleToString(d.getSumOfCounts(), maxWidth, 4).trim();
            temp.append(pad(count, " ", maxWidth + 1 - count.length(), true));
        }
        temp.append("\n");        
      }
    }

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
   * Outputs the generated clusters into a string.
   * 
   * @return the clusterer in string representation
   */
  protected String toStringOriginal () {
    if (m_priors == null) {
      return "No clusterer built yet!";
    }
    StringBuffer temp = new StringBuffer();
    temp.append("\nEM\n==\n");
    if (m_initialNumClusters == -1) {
      temp.append("\nNumber of clusters selected by cross validation: "
		  +m_num_clusters+"\n");
    } else {
      temp.append("\nNumber of clusters: " + m_num_clusters + "\n");
    }

    for (int j = 0; j < m_num_clusters; j++) {
      temp.append("\nCluster: " + j + " Prior probability: " 
		  + Utils.doubleToString(m_priors[j], 4) + "\n\n");

      for (int i = 0; i < m_num_attribs; i++) {
        temp.append("Attribute: " + m_theInstances.attribute(i).name() + "\n");

        if (m_theInstances.attribute(i).isNominal()) {
          if (m_model[j][i] != null) {
            temp.append(m_model[j][i].toString());
          }
        }
        else {
          temp.append("Normal Distribution. Mean = " 
		      + Utils.doubleToString(m_modelNormal[j][i][0], 4) 
		      + " StdDev = " 
		      + Utils.doubleToString(m_modelNormal[j][i][1], 4) 
		      + "\n");
        }
      }
    }

    return  temp.toString();
  }


  /**
   * verbose output for debugging
   * @param inst the training instances
   */
  private void EM_Report (Instances inst) {
    int i, j, l, m;
    System.out.println("======================================");

    for (j = 0; j < m_num_clusters; j++) {
      for (i = 0; i < m_num_attribs; i++) {
	System.out.println("Clust: " + j + " att: " + i + "\n");

	if (m_theInstances.attribute(i).isNominal()) {
	  if (m_model[j][i] != null) {
	    System.out.println(m_model[j][i].toString());
	  }
	}
	else {
	  System.out.println("Normal Distribution. Mean = " 
			     + Utils.doubleToString(m_modelNormal[j][i][0]
						    , 8, 4) 
			     + " StandardDev = " 
			     + Utils.doubleToString(m_modelNormal[j][i][1]
						    , 8, 4) 
			     + " WeightSum = " 
			     + Utils.doubleToString(m_modelNormal[j][i][2]
						    , 8, 4));
	}
      }
    }
    
    for (l = 0; l < inst.numInstances(); l++) {
      m = Utils.maxIndex(m_weights[l]);
      System.out.print("Inst " + Utils.doubleToString((double)l, 5, 0) 
		       + " Class " + m + "\t");
      for (j = 0; j < m_num_clusters; j++) {
	System.out.print(Utils.doubleToString(m_weights[l][j], 7, 5) + "  ");
      }
      System.out.println();
    }
  }


  /**
   * estimate the number of clusters by cross validation on the training
   * data.
   *
   * @throws Exception if something goes wrong
   */
  private void CVClusters ()
    throws Exception {
    double CVLogLikely = -Double.MAX_VALUE;
    double templl, tll;
    boolean CVincreased = true;
    m_num_clusters = 1;
    int num_clusters = m_num_clusters;
    int i;
    Random cvr;
    Instances trainCopy;
    int numFolds = (m_theInstances.numInstances() < 10) 
      ? m_theInstances.numInstances() 
      : 10;

    boolean ok = true;
    int seed = getSeed();
    int restartCount = 0;
    CLUSTER_SEARCH: while (CVincreased) {
      // theInstances.stratify(10);
        
      CVincreased = false;
      cvr = new Random(getSeed());
      trainCopy = new Instances(m_theInstances);
      trainCopy.randomize(cvr);
      templl = 0.0;
      for (i = 0; i < numFolds; i++) {
	Instances cvTrain = trainCopy.trainCV(numFolds, i, cvr);
	if (num_clusters > cvTrain.numInstances()) {
	  break CLUSTER_SEARCH;
	}
	Instances cvTest = trainCopy.testCV(numFolds, i);
	m_rr = new Random(seed);
        for (int z=0; z<10; z++) m_rr.nextDouble();
	m_num_clusters = num_clusters;
	EM_Init(cvTrain);
	try {
	  iterate(cvTrain, false);
	} catch (Exception ex) {
	  // catch any problems - i.e. empty clusters occuring
	  ex.printStackTrace();
          //          System.err.println("Restarting after CV training failure ("+num_clusters+" clusters");
          seed++;
          restartCount++;
          ok = false;
          if (restartCount > 5) {
            break CLUSTER_SEARCH;
          }
	  break;
	}
        try {
          tll = E(cvTest, false);
        } catch (Exception ex) {
          // catch any problems - i.e. empty clusters occuring
          //          ex.printStackTrace();
          ex.printStackTrace();
          //          System.err.println("Restarting after CV testing failure ("+num_clusters+" clusters");
          //          throw new Exception(ex); 
          seed++;
          restartCount++;
          ok = false;
          if (restartCount > 5) {
            break CLUSTER_SEARCH;
          }
          break;
        }

	if (m_verbose) {
	  System.out.println("# clust: " + num_clusters + " Fold: " + i 
			     + " Loglikely: " + tll);
	}
	templl += tll;
      }

      if (ok) {
        restartCount = 0;
        seed = getSeed();
        templl /= (double)numFolds;
        
        if (m_verbose) {
          System.out.println("===================================" 
                             + "==============\n# clust: " 
                             + num_clusters 
                             + " Mean Loglikely: " 
                             + templl 
                             + "\n================================" 
                             + "=================");
        }
        
        if (templl > CVLogLikely) {
          CVLogLikely = templl;
          CVincreased = true;
          num_clusters++;
        }
      }
    }

    if (m_verbose) {
      System.out.println("Number of clusters: " + (num_clusters - 1));
    }

    m_num_clusters = num_clusters - 1;
  }


  /**
   * Returns the number of clusters.
   *
   * @return the number of clusters generated for a training dataset.
   * @throws Exception if number of clusters could not be returned
   * successfully
   */
  public int numberOfClusters ()
    throws Exception {
    if (m_num_clusters == -1) {
      throw  new Exception("Haven't generated any clusters!");
    }

    return  m_num_clusters;
  }

 /**
  * Updates the minimum and maximum values for all the attributes
  * based on a new instance.
  *
  * @param instance the new instance
  */
  private void updateMinMax(Instance instance) {
    
    for (int j = 0; j < m_theInstances.numAttributes(); j++) {
      if (!instance.isMissing(j)) {
	if (Double.isNaN(m_minValues[j])) {
	  m_minValues[j] = instance.value(j);
	  m_maxValues[j] = instance.value(j);
	} else {
	  if (instance.value(j) < m_minValues[j]) {
	    m_minValues[j] = instance.value(j);
	  } else {
	    if (instance.value(j) > m_maxValues[j]) {
	      m_maxValues[j] = instance.value(j);
	    }
	  }
	}
      }
    }
  }

  /**
   * Returns default capabilities of the clusterer (i.e., the ones of 
   * SimpleKMeans).
   *
   * @return      the capabilities of this clusterer
   */
  public Capabilities getCapabilities() {
    Capabilities result = new SimpleKMeans().getCapabilities();
    result.setOwner(this);
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
  public void buildClusterer (Instances data)
    throws Exception {
    
    // can clusterer handle the data?
    getCapabilities().testWithFail(data);

    m_replaceMissing = new ReplaceMissingValues();
    Instances instances = new Instances(data);
    instances.setClassIndex(-1);
    m_replaceMissing.setInputFormat(instances);
    data = weka.filters.Filter.useFilter(instances, m_replaceMissing);
    instances = null;
    
    m_theInstances = data;

    // calculate min and max values for attributes
    m_minValues = new double [m_theInstances.numAttributes()];
    m_maxValues = new double [m_theInstances.numAttributes()];
    for (int i = 0; i < m_theInstances.numAttributes(); i++) {
      m_minValues[i] = m_maxValues[i] = Double.NaN;
    }
    for (int i = 0; i < m_theInstances.numInstances(); i++) {
      updateMinMax(m_theInstances.instance(i));
    }

    doEM();
    
    // save memory
    m_theInstances = new Instances(m_theInstances,0);
  }

  /**
   * Returns the cluster priors.
   * 
   * @return the cluster priors
   */
  public double[] clusterPriors() {

    double[] n = new double[m_priors.length];
  
    System.arraycopy(m_priors, 0, n, 0, n.length);
    return n;
  }

  /**
   * Computes the log of the conditional density (per cluster) for a given instance.
   * 
   * @param inst the instance to compute the density for
   * @return an array containing the estimated densities
   * @throws Exception if the density could not be computed
   * successfully
   */
  public double[] logDensityPerClusterForInstance(Instance inst) throws Exception {

    int i, j;
    double logprob;
    double[] wghts = new double[m_num_clusters];
    
    m_replaceMissing.input(inst);
    inst = m_replaceMissing.output();

    for (i = 0; i < m_num_clusters; i++) {
      //      System.err.println("Cluster : "+i);
      logprob = 0.0;

      for (j = 0; j < m_num_attribs; j++) {
	if (!inst.isMissing(j)) {
	  if (inst.attribute(j).isNominal()) {
	    logprob += Math.log(m_model[i][j].getProbability(inst.value(j)));
	  }
	  else { // numeric attribute
	    logprob += logNormalDens(inst.value(j), 
				     m_modelNormal[i][j][0], 
				     m_modelNormal[i][j][1]);
	    /*	    System.err.println(logNormalDens(inst.value(j), 
				     m_modelNormal[i][j][0], 
				     m_modelNormal[i][j][1]) + " "); */
	  }
	}
      }
      //      System.err.println("");

      wghts[i] = logprob;
    }
    return  wghts;
  }


  /**
   * Perform the EM algorithm
   * 
   * @throws Exception if something goes wrong
   */
  private void doEM ()
    throws Exception {
    
    if (m_verbose) {
      System.out.println("Seed: " + getSeed());
    }

    m_rr = new Random(getSeed());

    // throw away numbers to avoid problem of similar initial numbers
    // from a similar seed
    for (int i=0; i<10; i++) m_rr.nextDouble();

    m_num_instances = m_theInstances.numInstances();
    m_num_attribs = m_theInstances.numAttributes();

    if (m_verbose) {
      System.out.println("Number of instances: " 
			 + m_num_instances 
			 + "\nNumber of atts: " 
			 + m_num_attribs 
			 + "\n");
    }

    // setDefaultStdDevs(theInstances);
    // cross validate to determine number of clusters?
    if (m_initialNumClusters == -1) {
      if (m_theInstances.numInstances() > 9) {
	CVClusters();
	m_rr = new Random(getSeed());
	for (int i=0; i<10; i++) m_rr.nextDouble();
      } else {
	m_num_clusters = 1;
      }
    }

    // fit full training set
    EM_Init(m_theInstances);
    m_loglikely = iterate(m_theInstances, m_verbose);
  }


  /**
   * iterates the E and M steps until the log likelihood of the data
   * converges.
   *
   * @param inst the training instances.
   * @param report be verbose.
   * @return the log likelihood of the data
   * @throws Exception if something goes wrong
   */
  private double iterate (Instances inst, boolean report)
    throws Exception {
    
    int i;
    double llkold = 0.0;
    double llk = 0.0;

    if (report) {
      EM_Report(inst);
    }

    boolean ok = false;
    int seed = getSeed();
    int restartCount = 0;
    while (!ok) {
      try {
        for (i = 0; i < m_max_iterations; i++) {
          llkold = llk;
          llk = E(inst, true);
          
          if (report) {
            System.out.println("Loglikely: " + llk);
          }
          
          if (i > 0) {
            if ((llk - llkold) < 1e-6) {
              break;
            }
          }
          M(inst);
        }
        ok = true;
      } catch (Exception ex) {
        //        System.err.println("Restarting after training failure");
        ex.printStackTrace();
        seed++;
        restartCount++;
        m_rr = new Random(seed);
        for (int z = 0; z < 10; z++) {
          m_rr.nextDouble(); m_rr.nextInt();
        }
        if (restartCount > 5) {
          //          System.err.println("Reducing the number of clusters");
          m_num_clusters--;
          restartCount = 0;
        }
        EM_Init(m_theInstances);
      }
    }
      
    if (report) {
      EM_Report(inst);
    }

    return  llk;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6298 $");
  }

  // ============
  // Test method.
  // ============
  /**
   * Main method for testing this class.
   *
   * @param argv should contain the following arguments: <p>
   * -t training file [-T test file] [-N number of clusters] [-S random seed]
   */
  public static void main (String[] argv) {
    runClusterer(new EM(), argv);
  }
}

