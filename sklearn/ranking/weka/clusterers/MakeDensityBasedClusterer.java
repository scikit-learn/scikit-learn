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
 *    MakeDensityBasedClusterer.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.clusterers;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.estimators.DiscreteEstimator;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class for wrapping a Clusterer to make it return a distribution and density. Fits normal distributions and discrete distributions within each cluster produced by the wrapped clusterer. Supports the NumberOfClustersRequestable interface only if the wrapped Clusterer does.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -M &lt;num&gt;
 *  minimum allowable standard deviation for normal density computation 
 *  (default 1e-6)</pre>
 * 
 * <pre> -W &lt;clusterer name&gt;
 *  Clusterer to wrap.
 *  (default weka.clusterers.SimpleKMeans)</pre>
 * 
 * <pre> 
 * Options specific to clusterer weka.clusterers.SimpleKMeans:
 * </pre>
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
 * <pre> -S &lt;num&gt;
 *  Random number seed.
 *  (default 10)</pre>
 * 
 <!-- options-end -->
 * 
 * Options after "--" are passed on to the base clusterer.
 *
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5488 $
 */
public class MakeDensityBasedClusterer 
  extends AbstractDensityBasedClusterer
  implements NumberOfClustersRequestable, 
	     OptionHandler, 
	     WeightedInstancesHandler {

  /** for serialization */
  static final long serialVersionUID = -5643302427972186631L;
  
  /** holds training instances header information */
  private Instances m_theInstances;
  /** prior probabilities for the fitted clusters */
  private double [] m_priors;
  /** normal distributions fitted to each numeric attribute in each cluster */
  private double [][][] m_modelNormal;
  /** discrete distributions fitted to each discrete attribute in each cluster */
  private DiscreteEstimator [][] m_model;
  /** default minimum standard deviation */
  private double m_minStdDev = 1e-6;
  /** The clusterer being wrapped */
  private Clusterer m_wrappedClusterer = new weka.clusterers.SimpleKMeans();
  /** globally replace missing values */
  private ReplaceMissingValues m_replaceMissing;

  /**
   * Default constructor.
   * 
   */  
  public MakeDensityBasedClusterer() {
    super();
  }
   
  /**
   * Contructs a MakeDensityBasedClusterer wrapping a given Clusterer.
   * 
   * @param toWrap the clusterer to wrap around
   */    
  public MakeDensityBasedClusterer(Clusterer toWrap) {

    setClusterer(toWrap);
  }
  
  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Class for wrapping a Clusterer to make it return a distribution "
      + "and density. Fits normal distributions and discrete distributions "
      + "within each cluster produced by the wrapped clusterer. Supports the "
      + "NumberOfClustersRequestable interface only if the wrapped Clusterer "
      + "does.";
  }

  /**
   * String describing default clusterer.
   * 
   * @return 		the default clusterer classname
   */
  protected String defaultClustererString() {
    return SimpleKMeans.class.getName();
  }

  /**
   * Set the number of clusters to generate.
   *
   * @param n the number of clusters to generate
   * @throws Exception if the wrapped clusterer has not been set, or if
   * the wrapped clusterer does not implement this facility.
   */
  public void setNumClusters(int n) throws Exception {
    if (m_wrappedClusterer == null) {
      throw new Exception("Can't set the number of clusters to generate - "
			  +"no clusterer has been set yet.");
    }
    if (!(m_wrappedClusterer instanceof NumberOfClustersRequestable)) {
      throw new Exception("Can't set the number of clusters to generate - "
			  +"wrapped clusterer does not support this facility.");
    }

    ((NumberOfClustersRequestable)m_wrappedClusterer).setNumClusters(n);
  }

  /**
   * Returns default capabilities of the clusterer (i.e., of the wrapper
   * clusterer).
   *
   * @return      the capabilities of this clusterer
   */
  public Capabilities getCapabilities() {
    if (m_wrappedClusterer != null) {
      return m_wrappedClusterer.getCapabilities();
    }
    Capabilities result = super.getCapabilities();
    result.disableAll();
    result.enable(Capability.NO_CLASS);
    
    return result;
  }
  
  /**
   * Builds a clusterer for a set of instances.
   *
   * @param data the instances to train the clusterer with
   * @throws Exception if the clusterer hasn't been set or something goes wrong
   */  
  public void buildClusterer(Instances data) throws Exception {
    // can clusterer handle the data?
    getCapabilities().testWithFail(data);

    m_replaceMissing = new ReplaceMissingValues();
    m_replaceMissing.setInputFormat(data);
    data = weka.filters.Filter.useFilter(data, m_replaceMissing);

    m_theInstances = new Instances(data, 0);
    if (m_wrappedClusterer == null) {
      throw new Exception("No clusterer has been set");
    }
    m_wrappedClusterer.buildClusterer(data);
    m_model = 
       new DiscreteEstimator[m_wrappedClusterer.numberOfClusters()][data.numAttributes()];
    m_modelNormal = 
      new double[m_wrappedClusterer.numberOfClusters()][data.numAttributes()][2];
    double[][] weights =  new double[m_wrappedClusterer.numberOfClusters()][data.numAttributes()];
    m_priors = new double[m_wrappedClusterer.numberOfClusters()]; 
     for (int i = 0; i < m_wrappedClusterer.numberOfClusters(); i++) {
       m_priors[i] = 1.0; // laplace correction
       for (int j = 0; j < data.numAttributes(); j++) {
	 if (data.attribute(j).isNominal()) {
	   m_model[i][j] = new DiscreteEstimator(data.attribute(j).numValues(),
						 true);
	 }
       }
     }
     
     Instance inst = null;

     // Compute mean, etc.
     int[] clusterIndex = new int[data.numInstances()];
     for (int i = 0; i < data.numInstances(); i++) {
       inst = data.instance(i);
       int cluster = m_wrappedClusterer.clusterInstance(inst);
       m_priors[cluster] += inst.weight();
       for (int j = 0; j < data.numAttributes(); j++) {
	 if (!inst.isMissing(j)) {
	   if (data.attribute(j).isNominal()) {
	     m_model[cluster][j].addValue(inst.value(j),inst.weight());
	   } else {
	     m_modelNormal[cluster][j][0] += inst.weight() * inst.value(j);
	     weights[cluster][j] += inst.weight();
	   }
	 }
       }
       clusterIndex[i] = cluster;
     }

     for (int j = 0; j < data.numAttributes(); j++) {
       if (data.attribute(j).isNumeric()) {
	 for (int i = 0; i < m_wrappedClusterer.numberOfClusters(); i++) {	   
	   if (weights[i][j] > 0) {
	     m_modelNormal[i][j][0] /= weights[i][j];
	   }
	 }
       }
     }

     // Compute standard deviations
     for (int i = 0; i < data.numInstances(); i++) {
       inst = data.instance(i);
       for (int j = 0; j < data.numAttributes(); j++) {
	 if (!inst.isMissing(j)) {
	   if (data.attribute(j).isNumeric()) {
	     double diff = m_modelNormal[clusterIndex[i]][j][0] - inst.value(j);
	     m_modelNormal[clusterIndex[i]][j][1] += inst.weight() * diff * diff;
	   }
	 }
       }
     }

     for (int j = 0; j < data.numAttributes(); j++) {
       if (data.attribute(j).isNumeric()) {
	 for (int i = 0; i < m_wrappedClusterer.numberOfClusters(); i++) {	   
	   if (weights[i][j] > 0) {
	     m_modelNormal[i][j][1] = 
	       Math.sqrt(m_modelNormal[i][j][1] / weights[i][j]);
	   } else if (weights[i][j] <= 0) {
	     m_modelNormal[i][j][1] = Double.MAX_VALUE;
	   }
	   if (m_modelNormal[i][j][1] <= m_minStdDev) {
	     m_modelNormal[i][j][1] = data.attributeStats(j).numericStats.stdDev;
	     if (m_modelNormal[i][j][1] <= m_minStdDev) {
	       m_modelNormal[i][j][1] = m_minStdDev;
	     }
	   }
	 }
       }
     }
     
     Utils.normalize(m_priors);
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
    double[] wghts = new double[m_wrappedClusterer.numberOfClusters()];
    
    m_replaceMissing.input(inst);
    inst = m_replaceMissing.output();

    for (i = 0; i < m_wrappedClusterer.numberOfClusters(); i++) {
      logprob = 0;
      for (j = 0; j < inst.numAttributes(); j++) {
	if (!inst.isMissing(j)) {
	  if (inst.attribute(j).isNominal()) {
	    logprob += Math.log(m_model[i][j].getProbability(inst.value(j)));
	  } else { // numeric attribute
	    logprob += logNormalDens(inst.value(j), 
				     m_modelNormal[i][j][0], 
				     m_modelNormal[i][j][1]);
	  }
	}
      }
      wghts[i] = logprob;
    }
    return  wghts;
  }

  /** Constant for normal distribution. */
  private static double m_normConst = 0.5 * Math.log(2 * Math.PI);

  /**
   * Density function of normal distribution.
   * @param x input value
   * @param mean mean of distribution
   * @param stdDev standard deviation of distribution
   * @return the density
   */
  private double logNormalDens (double x, double mean, double stdDev) {

    double diff = x - mean;
    
    return - (diff * diff / (2 * stdDev * stdDev))  - m_normConst - Math.log(stdDev);
  }
  
  /**
   * Returns the number of clusters.
   *
   * @return the number of clusters generated for a training dataset.
   * @throws Exception if number of clusters could not be returned successfully
   */
  public int numberOfClusters() throws Exception {

    return m_wrappedClusterer.numberOfClusters();
  }

  /**
   * Returns a description of the clusterer.
   *
   * @return a string containing a description of the clusterer
   */
  public String toString() {
    if (m_priors == null) {
      return "No clusterer built yet!";
    }

    StringBuffer text = new StringBuffer();
    text.append("MakeDensityBasedClusterer: \n\nWrapped clusterer: " 
		+ m_wrappedClusterer.toString());

    text.append("\nFitted estimators (with ML estimates of variance):\n");
    
    for (int j = 0; j < m_priors.length; j++) {
      text.append("\nCluster: " + j + " Prior probability: " 
		  + Utils.doubleToString(m_priors[j], 4) + "\n\n");
      
      for (int i = 0; i < m_model[0].length; i++) {
        text.append("Attribute: " + m_theInstances.attribute(i).name() + "\n");
	
        if (m_theInstances.attribute(i).isNominal()) {
          if (m_model[j][i] != null) {
            text.append(m_model[j][i].toString());
          }
        }
        else {
          text.append("Normal Distribution. Mean = " 
		      + Utils.doubleToString(m_modelNormal[j][i][0], 4) 
		      + " StdDev = " 
		      + Utils.doubleToString(m_modelNormal[j][i][1], 4) 
		      + "\n");
        }
      }
    }

    return  text.toString();
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String clustererTipText() {
    return "the clusterer to wrap";
  }

  /**
   * Sets the clusterer to wrap.
   *
   * @param toWrap the clusterer
   */
  public void setClusterer(Clusterer toWrap) {

    m_wrappedClusterer = toWrap;
  }

  /**
   * Gets the clusterer being wrapped.
   *
   * @return the clusterer
   */
  public Clusterer getClusterer() {

    return m_wrappedClusterer;
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

  /**
   * Get the minimum allowable standard deviation.
   * @return the minumum allowable standard deviation
   */
  public double getMinStdDev() {
    return m_minStdDev;
  }

  /**
   * Returns an enumeration describing the available options..
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    result.addElement(new Option(
	"\tminimum allowable standard deviation for normal density computation "
	+"\n\t(default 1e-6)"
	,"M",1,"-M <num>"));
	
    result.addElement(new Option(
	"\tClusterer to wrap.\n"
	+ "\t(default " + defaultClustererString() + ")",
	"W", 1,"-W <clusterer name>"));

    if ((m_wrappedClusterer != null) &&
	(m_wrappedClusterer instanceof OptionHandler)) {
      result.addElement(new Option(
	  "",
	  "", 0, "\nOptions specific to clusterer "
	  + m_wrappedClusterer.getClass().getName() + ":"));
      Enumeration enu = ((OptionHandler)m_wrappedClusterer).listOptions();
      while (enu.hasMoreElements()) {
	result.addElement(enu.nextElement());
      }
    }
    
    return result.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -M &lt;num&gt;
   *  minimum allowable standard deviation for normal density computation 
   *  (default 1e-6)</pre>
   * 
   * <pre> -W &lt;clusterer name&gt;
   *  Clusterer to wrap.
   *  (default weka.clusterers.SimpleKMeans)</pre>
   * 
   * <pre> 
   * Options specific to clusterer weka.clusterers.SimpleKMeans:
   * </pre>
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
   * <pre> -S &lt;num&gt;
   *  Random number seed.
   *  (default 10)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    String optionString = Utils.getOption('M', options);
    if (optionString.length() != 0)
      setMinStdDev((new Double(optionString)).doubleValue());
    else
      setMinStdDev(1e-6);
     
    String wString = Utils.getOption('W', options);
    if (wString.length() == 0)
      wString = defaultClustererString();
    setClusterer(AbstractClusterer.forName(wString, Utils.partitionOptions(options)));
  }

  /**
   * Gets the current settings of the clusterer.
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String[] getOptions() {

    String [] clustererOptions = new String [0];
    if ((m_wrappedClusterer != null) &&
	(m_wrappedClusterer instanceof OptionHandler)) {
      clustererOptions = ((OptionHandler)m_wrappedClusterer).getOptions();
    }
    String [] options = new String [clustererOptions.length + 5];
    int current = 0;

    options[current++] = "-M";
    options[current++] = ""+getMinStdDev();

    if (getClusterer() != null) {
      options[current++] = "-W";
      options[current++] = getClusterer().getClass().getName();
    }
    options[current++] = "--";

    System.arraycopy(clustererOptions, 0, options, current, 
		     clustererOptions.length);
    current += clustererOptions.length;
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5488 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv the options
   */
  public static void main(String [] argv) {
    runClusterer(new MakeDensityBasedClusterer(), argv);
  }
}

