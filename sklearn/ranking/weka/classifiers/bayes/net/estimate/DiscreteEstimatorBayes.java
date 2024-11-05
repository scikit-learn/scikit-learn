
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * DiscreteEstimatorBayes.java
 * Adapted from DiscreteEstimator.java
 * 
 */
package weka.classifiers.bayes.net.estimate;

import weka.classifiers.bayes.net.search.local.Scoreable;
import weka.core.RevisionUtils;
import weka.core.Statistics;
import weka.core.Utils;
import weka.estimators.DiscreteEstimator;
import weka.estimators.Estimator;

/**
 * Symbolic probability estimator based on symbol counts and a prior.
 * 
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 1.7 $
 */
public class DiscreteEstimatorBayes extends Estimator
  implements Scoreable {

  /** for serialization */
  static final long serialVersionUID = 4215400230843212684L;
  
  /**
   * Hold the counts
   */
  protected double[] m_Counts;

  /**
   * Hold the sum of counts
   */
  protected double   m_SumOfCounts;

  /**
   * Holds number of symbols in distribution
   */
  protected int      m_nSymbols = 0;

  /**
   * Holds the prior probability
   */
  protected double   m_fPrior = 0.0;

  /**
   * Constructor
   * 
   * @param nSymbols the number of possible symbols (remember to include 0)
   * @param fPrior
   */
  public DiscreteEstimatorBayes(int nSymbols, double fPrior) {
    m_fPrior = fPrior;
    m_nSymbols = nSymbols;
    m_Counts = new double[m_nSymbols];

    for (int iSymbol = 0; iSymbol < m_nSymbols; iSymbol++) {
      m_Counts[iSymbol] = m_fPrior;
    } 

    m_SumOfCounts = m_fPrior * (double) m_nSymbols;
  }    // DiscreteEstimatorBayes

  /**
   * Add a new data value to the current estimator.
   * 
   * @param data the new data value
   * @param weight the weight assigned to the data value
   */
  public void addValue(double data, double weight) {
    m_Counts[(int) data] += weight;
    m_SumOfCounts += weight;
  } 

  /**
   * Get a probability estimate for a value
   * 
   * @param data the value to estimate the probability of
   * @return the estimated probability of the supplied value
   */
  public double getProbability(double data) {
    if (m_SumOfCounts == 0) {

      // this can only happen if numSymbols = 0 in constructor
      return 0;
    } 

    return (double) m_Counts[(int) data] / m_SumOfCounts;
  } 

  /**
   * Get a counts for a value
   * 
   * @param data the value to get the counts for
   * @return the count of the supplied value
   */
  public double getCount(double data) {
    if (m_SumOfCounts == 0) {
      // this can only happen if numSymbols = 0 in constructor
      return 0;
    } 

    return m_Counts[(int) data];
  } 
  
  /**
   * Gets the number of symbols this estimator operates with
   * 
   * @return the number of estimator symbols
   */
  public int getNumSymbols() {
    return (m_Counts == null) ? 0 : m_Counts.length;
  } 

  /**
   * Gets the log score contribution of this distribution
   * @param nType score type
   * @return the score
   */
  public double logScore(int nType, int nCardinality) {
	    double fScore = 0.0;

	    switch (nType) {

	    case (Scoreable.BAYES): {
	      for (int iSymbol = 0; iSymbol < m_nSymbols; iSymbol++) {
		fScore += Statistics.lnGamma(m_Counts[iSymbol]);
	      } 

	      fScore -= Statistics.lnGamma(m_SumOfCounts);
	      if (m_fPrior != 0.0) {
		      fScore -= m_nSymbols * Statistics.lnGamma(m_fPrior);
	    	  fScore += Statistics.lnGamma(m_nSymbols * m_fPrior);
	      }
	    } 

	      break;
		  case (Scoreable.BDeu): {
		  for (int iSymbol = 0; iSymbol < m_nSymbols; iSymbol++) {
			fScore += Statistics.lnGamma(m_Counts[iSymbol]);
		  } 

		  fScore -= Statistics.lnGamma(m_SumOfCounts);
		  //fScore -= m_nSymbols * Statistics.lnGamma(1.0);
		  //fScore += Statistics.lnGamma(m_nSymbols * 1.0);
	      fScore -= m_nSymbols * Statistics.lnGamma(1.0/(m_nSymbols * nCardinality));
	      fScore += Statistics.lnGamma(1.0/nCardinality);
		} 
		  break;

	    case (Scoreable.MDL):

	    case (Scoreable.AIC):

	    case (Scoreable.ENTROPY): {
	      for (int iSymbol = 0; iSymbol < m_nSymbols; iSymbol++) {
		double fP = getProbability(iSymbol);

		fScore += m_Counts[iSymbol] * Math.log(fP);
	      } 
	    } 

	      break;

	    default: {}
	    }

	    return fScore;
	  } 

  /**
   * Display a representation of this estimator
   * 
   * @return a string representation of the estimator
   */
  public String toString() {
    String result = "Discrete Estimator. Counts = ";

    if (m_SumOfCounts > 1) {
      for (int i = 0; i < m_Counts.length; i++) {
	result += " " + Utils.doubleToString(m_Counts[i], 2);
      } 

      result += "  (Total = " + Utils.doubleToString(m_SumOfCounts, 2) 
		+ ")\n";
    } else {
      for (int i = 0; i < m_Counts.length; i++) {
	result += " " + m_Counts[i];
      } 

      result += "  (Total = " + m_SumOfCounts + ")\n";
    } 

    return result;
  } 
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.7 $");
  }
  
  /**
   * Main method for testing this class.
   * 
   * @param argv should contain a sequence of integers which
   * will be treated as symbolic.
   */
  public static void main(String[] argv) {
    try {
      if (argv.length == 0) {
	System.out.println("Please specify a set of instances.");

	return;
      } 

      int current = Integer.parseInt(argv[0]);
      int max = current;

      for (int i = 1; i < argv.length; i++) {
	current = Integer.parseInt(argv[i]);

	if (current > max) {
	  max = current;
	} 
      } 

      DiscreteEstimator newEst = new DiscreteEstimator(max + 1, true);

      for (int i = 0; i < argv.length; i++) {
	current = Integer.parseInt(argv[i]);

	System.out.println(newEst);
	System.out.println("Prediction for " + current + " = " 
			   + newEst.getProbability(current));
	newEst.addValue(current, 1);
      } 
    } catch (Exception e) {
      System.out.println(e.getMessage());
    } 
  }    // main
 
}      // class DiscreteEstimatorBayes
