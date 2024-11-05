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
 *   KDDataGenerator.java
 *   Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.boundaryvisualizer;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;

import java.io.Serializable;
import java.util.Random;

/**
 * KDDataGenerator. Class that uses kernels to generate new random
 * instances based on a supplied set of instances.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.7 $
 * @since 1.0
 * @see DataGenerator
 * @see Serializable
 */
public class KDDataGenerator
  implements DataGenerator, Serializable {

  /** for serialization */
  private static final long serialVersionUID = -958573275606402792L;

  /** the instances to use */
  private Instances m_instances;

  /** standard deviations of the normal distributions for numeric attributes in
   * each KD estimator */
  private double [] m_standardDeviations;

  /** global means or modes to use for missing values */
  private double [] m_globalMeansOrModes;

  /** minimum standard deviation for numeric attributes */
  private double m_minStdDev = 1e-5;

  /** Laplace correction for discrete distributions */
  private double m_laplaceConst = 1.0;

  /** random number seed */
  private int m_seed = 1;

  /** random number generator */
  private Random m_random;

  /** which dimensions to use for computing a weight for each generated
   * instance */
  private boolean [] m_weightingDimensions;
  
  /** the values for the weighting dimensions to use for computing the weight 
   * for the next instance to be generated */
  private double [] m_weightingValues;

  private static double m_normConst = Math.sqrt(2*Math.PI);

  /** Number of neighbours to use for kernel bandwidth */
  private int m_kernelBandwidth = 3;

  /** standard deviations for numeric attributes computed from the 
   * m_kernelBandwidth nearest neighbours for each kernel. */
  private double [][] m_kernelParams;

  /** The minimum values for numeric attributes. */
  protected double [] m_Min;
  
  /** The maximum values for numeric attributes. */
  protected double [] m_Max;

  /**
   * Initialize the generator using the supplied instances
   *
   * @param inputInstances the instances to use as the basis of the kernels
   * @throws Exception if an error occurs
   */
  public void buildGenerator(Instances inputInstances) throws Exception {
    m_random = new Random(m_seed);
    
    m_instances = inputInstances;
    m_standardDeviations = new double [m_instances.numAttributes()];
    m_globalMeansOrModes = new double [m_instances.numAttributes()];
    if (m_weightingDimensions == null) {
      m_weightingDimensions = new boolean[m_instances.numAttributes()];
    }
    /*    for (int i = 0; i < m_instances.numAttributes(); i++) {
      if (i != m_instances.classIndex()) {
	if (m_instances.attribute(i).isNumeric()) {
	  // global standard deviations
	  double var = m_instances.variance(i);
	  if (var == 0) {
	    var = m_minStdDev;
	  } else {
	    var = Math.sqrt(var);
	    //  heuristic to take into account # instances and dimensions
	    double adjust = Math.pow((double) m_instances.numInstances(), 
				     1.0 / m_instances.numAttributes());
	    //	  double adjust = m_instances.numInstances();
	    var /= adjust;
	  }
	  m_standardDeviations[i] = var;
	} else {
	  m_globalMeansOrModes[i] = m_instances.meanOrMode(i);
	}
      }
      } */
    for (int i = 0; i < m_instances.numAttributes(); i++) {
      if (i != m_instances.classIndex()) {
	m_globalMeansOrModes[i] = m_instances.meanOrMode(i);
      }
    }

    m_kernelParams = 
      new double [m_instances.numInstances()][m_instances.numAttributes()];
    computeParams();
  }

  public double [] getWeights() {

    double [] weights = new double[m_instances.numInstances()];

    for (int k = 0; k < m_instances.numInstances(); k++) {
      double weight = 1;
      for (int i = 0; i < m_instances.numAttributes(); i++) {
	if (m_weightingDimensions[i]) {
	  double mean = 0;
	  if (!m_instances.instance(k).isMissing(i)) {
	    mean = m_instances.instance(k).value(i);
	  } else {
	    mean = m_globalMeansOrModes[i];
	  }
	  double wm = 1.0;
	  
	  //	    wm = normalDens(m_weightingValues[i], mean, m_standardDeviations[i]);
	  wm = normalDens(m_weightingValues[i], mean, 
			  m_kernelParams[k][i]);
	  
	  weight *= wm;
	}
      }
      weights[k] = weight;
    }
    return weights;
  }

  /**
   * Return a cumulative distribution from a discrete distribution
   *
   * @param dist the distribution to use
   * @return the cumulative distribution
   */
  private double [] computeCumulativeDistribution(double [] dist) {

    double [] cumDist = new double[dist.length];
    double sum = 0;
    for (int i = 0; i < dist.length; i++) {
      sum += dist[i];
      cumDist[i] = sum;
    }
    
    return cumDist;
  }

  /**
   * Generates a new instance using one kernel estimator. Each successive
   * call to this method incremets the index of the kernel to use.
   *
   * @return the new random instance
   * @throws Exception if an error occurs
   */
  public double [][] generateInstances(int [] indices) throws Exception {
    
    double [][] values = new double[m_instances.numInstances()][];

    for (int k = 0; k < indices.length; k++) {
      values[indices[k]] = new double[m_instances.numAttributes()];
      for (int i = 0; i < m_instances.numAttributes(); i++) {
	if ((!m_weightingDimensions[i]) && (i != m_instances.classIndex())) {
	  if (m_instances.attribute(i).isNumeric()) {
	    double mean = 0;
	    double val = m_random.nextGaussian();
	    if (!m_instances.instance(indices[k]).isMissing(i)) {
	      mean = m_instances.instance(indices[k]).value(i);
	    } else {
	      mean = m_globalMeansOrModes[i];
	    }
	    
	    val *= m_kernelParams[indices[k]][i];
	    val += mean;

	    values[indices[k]][i] = val;
	  } else {
	    // nominal attribute
	    double [] dist = new double[m_instances.attribute(i).numValues()];
	    for (int j = 0; j < dist.length; j++) {
	      dist[j] = m_laplaceConst;
	    }
	    if (!m_instances.instance(indices[k]).isMissing(i)) {
	      dist[(int)m_instances.instance(indices[k]).value(i)]++;
	    } else {
	      dist[(int)m_globalMeansOrModes[i]]++;
	    }
	    Utils.normalize(dist);
	    double [] cumDist = computeCumulativeDistribution(dist);
	    double randomVal = m_random.nextDouble();
	    int instVal = 0;
	    for (int j = 0; j < cumDist.length; j++) {
	      if (randomVal <= cumDist[j]) {
		instVal = j;
		break;
	      }
	    }
	    values[indices[k]][i] = (double)instVal;
	  }
	}
      }
    }
    return values;
  }

  /**
   * Density function of normal distribution.
   * @param x input value
   * @param mean mean of distribution
   * @param stdDev standard deviation of distribution
   */
  private double normalDens (double x, double mean, double stdDev) {
    double diff = x - mean;
   
    return  (1/(m_normConst*stdDev))*Math.exp(-(diff*diff/(2*stdDev*stdDev)));
  }

  /**
   * Set which dimensions to use when computing a weight for the next
   * instance to generate
   *
   * @param dims an array of booleans indicating which dimensions to use
   */
  public void setWeightingDimensions(boolean [] dims) {
    m_weightingDimensions = dims;
  }

  /**
   * Set the values for the weighting dimensions to be used when computing
   * the weight for the next instance to be generated
   *
   * @param vals an array of doubles containing the values of the
   * weighting dimensions (corresponding to the entries that are set to
   * true throw setWeightingDimensions)
   */
  public void setWeightingValues(double [] vals) {
    m_weightingValues = vals;
  }

  /**
   * Return the number of kernels (there is one per training instance)
   *
   * @return the number of kernels
   */
  public int getNumGeneratingModels() {
    if (m_instances != null) {
      return m_instances.numInstances();
    }
    return 0;
  }

  /**
   * Set the kernel bandwidth (number of nearest neighbours to cover)
   *
   * @param kb an <code>int</code> value
   */
  public void setKernelBandwidth(int kb) {
    m_kernelBandwidth = kb;
  }

  /**
   * Get the kernel bandwidth
   *
   * @return an <code>int</code> value
   */
  public int getKernelBandwidth() {
    return m_kernelBandwidth;
  } 

  /**
   * Initializes a new random number generator using the
   * supplied seed.
   *
   * @param seed an <code>int</code> value
   */
  public void setSeed(int seed) {
    m_seed = seed;
    m_random = new Random(m_seed);
  }

  /**
   * Calculates the distance between two instances
   *
   * @param test the first instance
   * @param train the second instance
   * @return the distance between the two given instances, between 0 and 1
   */          
  private double distance(Instance first, Instance second) {  

    double diff, distance = 0;

    for(int i = 0; i < m_instances.numAttributes(); i++) { 
      if (i == m_instances.classIndex()) {
	continue;
      }
      double firstVal = m_globalMeansOrModes[i];
      double secondVal = m_globalMeansOrModes[i];

      switch (m_instances.attribute(i).type()) {
      case Attribute.NUMERIC:
	// If attribute is numeric
	if (!first.isMissing(i)) {
	  firstVal = first.value(i);
	}
	
	if (!second.isMissing(i)) {
	  secondVal = second.value(i);
	}

	diff = norm(firstVal,i) - norm(secondVal,i);

	break;
      default:
	diff = 0;
	break;
      }
      distance += diff * diff;
    }
    return Math.sqrt(distance);
  }

  /**
   * Normalizes a given value of a numeric attribute.
   *
   * @param x the value to be normalized
   * @param i the attribute's index
   */
  private double norm(double x,int i) {
    
    if (Double.isNaN(m_Min[i]) || Utils.eq(m_Max[i], m_Min[i])) {
      return 0;
    } else {
      return (x - m_Min[i]) / (m_Max[i] - m_Min[i]);
    }
  }

  /**
   * Updates the minimum and maximum values for all the attributes
   * based on a new instance.
   *
   * @param instance the new instance
   */
  private void updateMinMax(Instance instance) {  

    for (int j = 0; j < m_instances.numAttributes(); j++) {
      if (!instance.isMissing(j)) {
	if (Double.isNaN(m_Min[j])) {
	  m_Min[j] = instance.value(j);
	  m_Max[j] = instance.value(j);
	} else if (instance.value(j) < m_Min[j]) {
	  m_Min[j] = instance.value(j);
	} else if (instance.value(j) > m_Max[j]) {
	  m_Max[j] = instance.value(j);
	}
      }
    }
  }

  private void computeParams() throws Exception {
    // Calculate the minimum and maximum values
    m_Min = new double [m_instances.numAttributes()];
    m_Max = new double [m_instances.numAttributes()];
    for (int i = 0; i < m_instances.numAttributes(); i++) {
      m_Min[i] = m_Max[i] = Double.NaN;
    }
    for (int i = 0; i < m_instances.numInstances(); i++) {
      updateMinMax(m_instances.instance(i));
    }

    double [] distances = new double[m_instances.numInstances()];
    for (int i = 0; i < m_instances.numInstances(); i++) {
      Instance current = m_instances.instance(i);
      for (int j = 0; j < m_instances.numInstances(); j++) {
	distances[j] = distance(current, m_instances.instance(j));
      }
      int [] sorted = Utils.sort(distances);
      int k = m_kernelBandwidth;
      double bandwidth = distances[sorted[k]];

      // Check for bandwidth zero
      if (bandwidth <= 0) {
	for (int j = k + 1; j < sorted.length; j++) {
	  if (distances[sorted[j]] > bandwidth) {
	    bandwidth = distances[sorted[j]];
	    break;
	  }
	}
	if (bandwidth <= 0) {
	  throw new Exception("All training instances coincide with "
			      +"test instance!");
	}
      }
      for (int j = 0; j < m_instances.numAttributes(); j++) {
	if ((m_Max[j] - m_Min[j]) > 0) {
	  m_kernelParams[i][j] = bandwidth * (m_Max[j] - m_Min[j]);
	}
      }
    }
  }
}


