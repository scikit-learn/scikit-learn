
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
 * DiscreteEstimatorFullBayes.java
 * 
 */
package weka.classifiers.bayes.net.estimate;

import weka.core.RevisionUtils;
import weka.estimators.DiscreteEstimator;

/**
 * Symbolic probability estimator based on symbol counts and a prior.
 *  
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 1.3 $
 */
public class DiscreteEstimatorFullBayes 
  extends DiscreteEstimatorBayes {

  /** for serialization */
  static final long serialVersionUID = 6774941981423312133L;
  
  /**
   * Constructor
   * 
   * @param nSymbols the number of possible symbols (remember to include 0)
   * @param w1
   * @param w2
   * @param EmptyDist
   * @param ClassDist
   * @param fPrior
   */
  public DiscreteEstimatorFullBayes(int nSymbols, 
    double w1, double w2,
    DiscreteEstimatorBayes EmptyDist,
    DiscreteEstimatorBayes ClassDist,
    double fPrior) {
    
    super(nSymbols, fPrior);

    m_SumOfCounts = 0.0;
    for (int iSymbol = 0; iSymbol < m_nSymbols; iSymbol++) {
      double p1 = EmptyDist.getProbability(iSymbol);
      double p2 = ClassDist.getProbability(iSymbol);
      m_Counts[iSymbol] = w1 * p1 + w2 * p2;
      m_SumOfCounts += m_Counts[iSymbol];
    } 
  } // DiscreteEstimatorFullBayes
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.3 $");
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
 
}  // class DiscreteEstimatorFullBayes




