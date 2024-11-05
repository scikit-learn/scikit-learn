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
 *    UnivariateEqualFrequencyEstimator.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import java.util.Random;
import java.util.Collection;
import java.util.Set;
import java.util.Map;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.ArrayList;
import java.util.Arrays;

import weka.core.Statistics;
import weka.core.Utils;

/**
 * Simple histogram density estimator. Uses equal-frequency histograms
 * based on the specified number of bins (default: 10).
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5680 $
 */
public class UnivariateEqualFrequencyHistogramEstimator implements UnivariateDensityEstimator,
                                                     UnivariateIntervalEstimator {

  /** The collection used to store the weighted values. */
  protected TreeMap<Double, Double> m_TM = new TreeMap<Double, Double>();

  /** The interval boundaries. */
  protected double[] m_Boundaries = null;

  /** The weight of each interval. */
  protected double[] m_Weights = null;

  /** The weighted sum of values */
  protected double m_WeightedSum = 0;

  /** The weighted sum of squared values */
  protected double m_WeightedSumSquared = 0;

  /** The total sum of weights. */
  protected double m_SumOfWeights = 0;

  /** The number of bins to use. */
  protected int m_NumBins = 10;

  /** The current bandwidth (only computed when needed) */
  protected double m_Width = Double.MAX_VALUE;

  /** The exponent to use in computation of bandwidth (default: -0.25) */
  protected double m_Exponent = -0.25;

  /** The minimum allowed value of the kernel width (default: 1.0E-6) */
  protected double m_MinWidth = 1.0E-6;

  /** Constant for Gaussian density. */
  public static final double CONST = - 0.5 * Math.log(2 * Math.PI);

  /** The number of intervals used to approximate prediction interval. */
  protected int m_NumIntervals = 1000;

  /** Whether boundaries are updated or only weights. */
  protected boolean m_UpdateWeightsOnly = false;

  /**
   * Gets the number of bins 
   *
   * @return the number of bins.
   */
  public int getNumBins() {

    return m_NumBins;
  }

  /**
   * Sets the number of bins 
   *
   * @param numBins the number of bins
   */
  public void setNumBins(int numBins) {

    m_NumBins = numBins;
  }

  /**
   * Triggers construction of estimator based on current data
   * and then initializes the statistics.
   */
  public void initializeStatistics() {

    updateBoundariesAndOrWeights();

    m_TM = new TreeMap<Double, Double>();
    m_WeightedSum = 0;
    m_WeightedSumSquared = 0;
    m_SumOfWeights = 0;
    m_Weights = null;
  }    

  /**
   * Sets whether only weights should be udpated.
   */
  public void setUpdateWeightsOnly(boolean flag) {

    m_UpdateWeightsOnly = flag;
  }

  /**
   * Gets whether only weights should be udpated.*
   */
  public boolean getUpdateWeightsOnly() {

    return m_UpdateWeightsOnly;
  }

  /**
   * Adds a value to the density estimator.
   *
   * @param value the value to add
   * @param weight the weight of the value
   */
  public void addValue(double value, double weight) {

    // Add data point to collection
    m_WeightedSum += value * weight;
    m_WeightedSumSquared += value * value * weight;
    m_SumOfWeights += weight;
    if (m_TM.get(value) == null) {
      m_TM.put(value, weight);
    } else {
      m_TM.put(value, m_TM.get(value) + weight);
    }

    // Make sure estimator is updated
    if (!getUpdateWeightsOnly()) {
      m_Boundaries = null;
    }
    m_Weights = null;
  }

  /**
   * Updates the boundaries if necessary.
   */
  protected void updateBoundariesAndOrWeights() {

    // Do we need to update?
    if (m_Weights != null) {
      return;
    }

    // Update widths for cases that are out of bounds,
    // using same code as in kernel estimator

    // First, compute variance for scaling
    double mean = m_WeightedSum / m_SumOfWeights;
    double variance = m_WeightedSumSquared / m_SumOfWeights - mean * mean;
    if (variance < 0) {
      variance = 0;
    }
    
    // Compute kernel bandwidth
    m_Width = Math.sqrt(variance) * Math.pow(m_SumOfWeights, m_Exponent);
    
    if (m_Width <= m_MinWidth) {
      m_Width = m_MinWidth;
    }
    
    // Do we need to update weights only
    if (getUpdateWeightsOnly()) {
      updateWeightsOnly();
    } else {
      updateBoundariesAndWeights();
    }
  }
   
  /**
   * Updates the weights only.
   */
  protected void updateWeightsOnly() throws IllegalArgumentException {

    // Get values and keys from tree map
    Iterator<Map.Entry<Double,Double>> itr = m_TM.entrySet().iterator();
    int j = 1;
    m_Weights = new double[m_Boundaries.length - 1];
    while(itr.hasNext()) {
      Map.Entry<Double,Double> entry = itr.next();
      double value = entry.getKey();
      double weight = entry.getValue();
      if ((value < m_Boundaries[0]) || (value > m_Boundaries[m_Boundaries.length - 1])) {
        throw new IllegalArgumentException("Out-of-range value during weight update");
      }
      while (value > m_Boundaries[j]) {
        j++;
      }
      m_Weights[j - 1] += weight;
    }
  }

  /**
   * Updates the boundaries and weights.
   */
  protected void updateBoundariesAndWeights() {

    // Get values and keys from tree map
    double[] values = new double[m_TM.size()];
    double[] weights = new double[m_TM.size()];
    Iterator<Map.Entry<Double,Double>> itr = m_TM.entrySet().iterator();
    int j = 0;
    while(itr.hasNext()) {
      Map.Entry<Double,Double> entry = itr.next();
      values[j] = entry.getKey();
      weights[j] = entry.getValue();
      j++;
    }

    double freq = m_SumOfWeights / m_NumBins;
    double[] cutPoints = new double[m_NumBins - 1];
    double[] binWeights = new double[m_NumBins];
    double sumOfWeights = m_SumOfWeights;

    // Compute break points
    double weightSumSoFar = 0, lastWeightSum = 0;
    int cpindex = 0, lastIndex = -1;
    for (int i = 0; i < values.length - 1; i++) {

      // Update weight statistics
      weightSumSoFar += weights[i];
      sumOfWeights -= weights[i];

      // Have we passed the ideal size?
      if (weightSumSoFar >= freq) {

        // Is this break point worse than the last one?
        if (((freq - lastWeightSum) < (weightSumSoFar - freq)) && (lastIndex != -1)) {
          cutPoints[cpindex] = (values[lastIndex] + values[lastIndex + 1]) / 2;
          weightSumSoFar -= lastWeightSum;
          binWeights[cpindex] = lastWeightSum;
          lastWeightSum = weightSumSoFar;
          lastIndex = i;
        } else {
          cutPoints[cpindex] = (values[i] + values[i + 1]) / 2;
          binWeights[cpindex] = weightSumSoFar;
          weightSumSoFar = 0;
          lastWeightSum = 0;
          lastIndex = -1;
        }
        cpindex++;
        freq = (sumOfWeights + weightSumSoFar) / ((cutPoints.length + 1) - cpindex);
      } else {
        lastIndex = i;
        lastWeightSum = weightSumSoFar;
      }
    }

    // Check whether there was another possibility for a cut point
    if ((cpindex < cutPoints.length) && (lastIndex != -1)) {
      cutPoints[cpindex] = (values[lastIndex] + values[lastIndex + 1]) / 2;      
      binWeights[cpindex] = lastWeightSum;
      cpindex++;
      binWeights[cpindex] = weightSumSoFar - lastWeightSum;
    } else {
      binWeights[cpindex] = weightSumSoFar;
    }

    // Did we find any cutpoints?
    if (cpindex == 0) {
      m_Boundaries = null;
      m_Weights = null;
    } else {

      // Need to add weight of last data point to right-most bin
      binWeights[cpindex] += weights[values.length - 1];

      // Copy over boundaries and weights
      m_Boundaries = new double[cpindex + 2];
      m_Boundaries[0] = m_TM.firstKey();
      m_Boundaries[cpindex + 1] = m_TM.lastKey();
      System.arraycopy(cutPoints, 0, m_Boundaries, 1, cpindex);
      m_Weights = new double[cpindex + 1];
      System.arraycopy(binWeights, 0, m_Weights, 0, cpindex + 1);
    }
  }
   

  /**
   * Returns the interval for the given confidence value. 
   * 
   * @param conf the confidence value in the interval [0, 1]
   * @return the interval
   */
  public double[][] predictIntervals(double conf) {

    // Update the bandwidth
    updateBoundariesAndOrWeights();

    // Compute minimum and maximum value, and delta
    double val = Statistics.normalInverse(1.0 - (1.0 - conf) / 2);
    double min = m_TM.firstKey() - val * m_Width;
    double max = m_TM.lastKey() + val * m_Width;
    double delta = (max - min) / m_NumIntervals;

    // Create array with estimated probabilities
    double[] probabilities = new double[m_NumIntervals];
    double leftVal = Math.exp(logDensity(min));
    for (int i = 0; i < m_NumIntervals; i++) {
      double rightVal = Math.exp(logDensity(min + (i + 1) * delta));
      probabilities[i] = 0.5 * (leftVal + rightVal) * delta;
      leftVal = rightVal;
    }

    // Sort array based on area of bin estimates
    int[] sortedIndices = Utils.sort(probabilities);

    // Mark the intervals to use
    double sum = 0;
    boolean[] toUse = new boolean[probabilities.length];
    int k = 0;
    while ((sum < conf) && (k < toUse.length)){
      toUse[sortedIndices[toUse.length - (k + 1)]] = true;
      sum += probabilities[sortedIndices[toUse.length - (k + 1)]];
      k++;
    }

    // Don't need probabilities anymore
    probabilities = null;

    // Create final list of intervals
    ArrayList<double[]> intervals = new ArrayList<double[]>();

    // The current interval
    double[] interval = null;
    
    // Iterate through kernels
    boolean haveStartedInterval = false;
    for (int i = 0; i < m_NumIntervals; i++) {

      // Should the current bin be used?
      if (toUse[i]) {

        // Do we need to create a new interval?
        if (haveStartedInterval == false) {
          haveStartedInterval = true;
          interval = new double[2];
          interval[0] = min + i * delta;
        }

        // Regardless, we should update the upper boundary
        interval[1] = min + (i + 1) * delta;
      } else {

        // We need to finalize and store the last interval
        // if necessary.
        if (haveStartedInterval) {
          haveStartedInterval = false;
          intervals.add(interval);
        }
      }
    }

    // Add last interval if there is one
    if (haveStartedInterval) {
      intervals.add(interval);
    }

    return intervals.toArray(new double[0][0]);
  }

  /**
   * Returns the natural logarithm of the density estimate at the given
   * point.
   *
   * @param value the value at which to evaluate
   * @return the natural logarithm of the density estimate at the given
   * value
   */
  public double logDensity(double value) {

    // Update boundaries if necessary
    updateBoundariesAndOrWeights();

    if (m_Boundaries == null) {
      return Math.log(Double.MIN_VALUE);
    }

    // Find the bin
    int index = Arrays.binarySearch(m_Boundaries, value);

    // Is the value outside?
    if ((index == -1) || (index == -m_Boundaries.length - 1)) {

      // Use normal density outside
      double val = 0;
      if (index == -1) { // Smaller than minimum
        val = m_TM.firstKey() - value;
      } else {
        val = value - m_TM.lastKey();
      }
      return (CONST - Math.log(m_Width) - 0.5 * (val * val / (m_Width * m_Width))) -
        Math.log(m_SumOfWeights + 2); 
    }
    
    // Is value exactly equal to right-most boundary?
    if (index == m_Boundaries.length - 1) {
      index--;
    } else {

      // Need to reverse index if necessary
      if (index < 0) {
        index = -index - 2;
      }
    }
    
    // Figure out of width
    double width = m_Boundaries[index + 1] - m_Boundaries[index];

    // Density compontent from smeared-out data point
    double densSmearedOut = 1.0 / ((m_SumOfWeights + 2) * (m_Boundaries[m_Boundaries.length - 1] -
                                                           m_Boundaries[0]));

    // Return log of density
    if (m_Weights[index] <= 0) {

      /*      System.out.println(value);
      System.out.println(this);
      System.exit(1);*/
      // Just use one smeared-out data point
      return Math.log(densSmearedOut);
    } else {
      return Math.log(densSmearedOut + m_Weights[index] / ((m_SumOfWeights + 2) * width));
    }
  }

  /**
   * Returns textual description of this estimator.
   */
  public String toString() {

    StringBuffer text = new StringBuffer();

    text.append("EqualFrequencyHistogram estimator\n\n" +
                "Bandwidth for out of range cases " + m_Width + 
                ", total weight " + m_SumOfWeights);

    if (m_Boundaries != null) {
      text.append("\nLeft boundary\tRight boundary\tWeight\n");
      for (int i = 0; i < m_Boundaries.length - 1; i++) {
        text.append(m_Boundaries[i] + "\t" + m_Boundaries[i + 1] + "\t" + m_Weights[i] + "\t" +
                    Math.exp(logDensity((m_Boundaries[i + 1] + m_Boundaries[i]) / 2)) + "\n");
      }
    }

    return text.toString();
  }

  /**
   * Main method, used for testing this class.
   */
  public static void main(String[] args) {

    // Get random number generator initialized by system
    Random r = new Random();

    // Create density estimator
    UnivariateEqualFrequencyHistogramEstimator e = new UnivariateEqualFrequencyHistogramEstimator();

    // Output the density estimator
    System.out.println(e);
    
    // Monte Carlo integration
    double sum = 0;
    for (int i = 0; i < 1000; i++) {
      sum += Math.exp(e.logDensity(r.nextDouble() * 10.0 - 5.0));
    }
    System.out.println("Approximate integral: " + 10.0 * sum / 1000);
    
    // Add Gaussian values into it
    for (int i = 0; i < 1000; i++) {
      e.addValue(0.1 * r.nextGaussian() - 3, 1);
      e.addValue(r.nextGaussian() * 0.25, 3);
    }

    // Monte Carlo integration
    sum = 0;
    int points = 10000000;
    for (int i = 0; i < points; i++) {
      double value = r.nextDouble() * 20.0 - 10.0;
      sum += Math.exp(e.logDensity(value));
    }

    // Output the density estimator
    System.out.println(e);

    System.out.println("Approximate integral: " + 20.0 * sum / points);

    // Check interval estimates
    double[][] Intervals = e.predictIntervals(0.9);
    
    System.out.println("Printing histogram intervals ---------------------");
    
    for (int k = 0; k < Intervals.length; k++) {
      System.out.println("Left: " + Intervals[k][0] + "\t Right: " + Intervals[k][1]);
    }
    
    System.out.println("Finished histogram printing intervals ---------------------");

    double Covered = 0;
    for (int i = 0; i < 1000; i++) {
      double val = -1;
      if (r.nextDouble() < 0.25) {
        val = 0.1 * r.nextGaussian() - 3.0;
      } else {
        val = r.nextGaussian() * 0.25;
      }
      for (int k = 0; k < Intervals.length; k++) {
        if (val >= Intervals[k][0] && val <= Intervals[k][1]) {
          Covered++;
          break;
        }
      }
    }
    System.out.println("Coverage at 0.9 level for histogram intervals: " + Covered / 1000);

    for (int j = 1; j < 5; j++) {
      double numTrain = Math.pow(10, j);
      System.out.println("Number of training cases: " +
                         numTrain); 

      // Compare performance to normal estimator on normally distributed data
      UnivariateEqualFrequencyHistogramEstimator eHistogram = new UnivariateEqualFrequencyHistogramEstimator();
      UnivariateNormalEstimator eNormal = new UnivariateNormalEstimator();
      
      // Add training cases
      for (int i = 0; i < numTrain; i++) {
        double val = r.nextGaussian() * 1.5 + 0.5;
        /*        if (j == 4) {
          System.err.println(val);
          }*/
        eHistogram.addValue(val, 1);
        eNormal.addValue(val, 1);
      }

      // Monte Carlo integration
      sum = 0;
      points = 10000000;
      for (int i = 0; i < points; i++) {
        double value = r.nextDouble() * 20.0 - 10.0;
        sum += Math.exp(eHistogram.logDensity(value));
      }
      System.out.println(eHistogram);
      System.out.println("Approximate integral for histogram estimator: " + 20.0 * sum / points);

      // Evaluate estimators
      double loglikelihoodHistogram = 0, loglikelihoodNormal = 0;
      for (int i = 0; i < 1000; i++) {
        double val = r.nextGaussian() * 1.5 + 0.5;
        loglikelihoodHistogram += eHistogram.logDensity(val);
        loglikelihoodNormal += eNormal.logDensity(val);
      }
      System.out.println("Loglikelihood for histogram estimator: " +
                         loglikelihoodHistogram / 1000);
      System.out.println("Loglikelihood for normal estimator: " +
                         loglikelihoodNormal / 1000);

      // Check interval estimates
      double[][] histogramIntervals = eHistogram.predictIntervals(0.95);
      double[][] normalIntervals = eNormal.predictIntervals(0.95);

      System.out.println("Printing histogram intervals ---------------------");
      
      for (int k = 0; k < histogramIntervals.length; k++) {
        System.out.println("Left: " + histogramIntervals[k][0] + "\t Right: " + histogramIntervals[k][1]);
      }

      System.out.println("Finished histogram printing intervals ---------------------");

      System.out.println("Printing normal intervals ---------------------");
      
      for (int k = 0; k < normalIntervals.length; k++) {
        System.out.println("Left: " + normalIntervals[k][0] + "\t Right: " + normalIntervals[k][1]);
      }

      System.out.println("Finished normal printing intervals ---------------------");
 
      double histogramCovered = 0;
      double normalCovered = 0;
      for (int i = 0; i < 1000; i++) {
        double val = r.nextGaussian() * 1.5 + 0.5;
        for (int k = 0; k < histogramIntervals.length; k++) {
          if (val >= histogramIntervals[k][0] && val <= histogramIntervals[k][1]) {
            histogramCovered++;
            break;
          }
        }
        for (int k = 0; k < normalIntervals.length; k++) {
          if (val >= normalIntervals[k][0] && val <= normalIntervals[k][1]) {
            normalCovered++;
            break;
          }
        }
      }
      System.out.println("Coverage at 0.95 level for histogram intervals: " + histogramCovered / 1000);
      System.out.println("Coverage at 0.95 level for normal intervals: " + normalCovered / 1000);
      
      histogramIntervals = eHistogram.predictIntervals(0.8);
      normalIntervals = eNormal.predictIntervals(0.8);
      histogramCovered = 0;
      normalCovered = 0;
      for (int i = 0; i < 1000; i++) {
        double val = r.nextGaussian() * 1.5 + 0.5;
        for (int k = 0; k < histogramIntervals.length; k++) {
          if (val >= histogramIntervals[k][0] && val <= histogramIntervals[k][1]) {
            histogramCovered++;
            break;
          }
        }
        for (int k = 0; k < normalIntervals.length; k++) {
          if (val >= normalIntervals[k][0] && val <= normalIntervals[k][1]) {
            normalCovered++;
            break;
          }
        }
      }
      System.out.println("Coverage at 0.8 level for histogram intervals: " + histogramCovered / 1000);
      System.out.println("Coverage at 0.8 level for normal intervals: " + normalCovered / 1000);
    }
  }
}