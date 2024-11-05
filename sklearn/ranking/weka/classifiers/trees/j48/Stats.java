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
 *    Stats.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Statistics;

/**
 * Class implementing a statistical routine needed by J48 to
 * compute its error estimate.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.9 $
 */
public class Stats
  implements RevisionHandler {

  /**
   * Computes estimated extra error for given total number of instances
   * and error using normal approximation to binomial distribution
   * (and continuity correction).
   *
   * @param N number of instances
   * @param e observed error
   * @param CF confidence value
   */
  public static double addErrs(double N, double e, float CF){

    // Ignore stupid values for CF
    if (CF > 0.5) {
      System.err.println("WARNING: confidence value for pruning " +
			 " too high. Error estimate not modified.");
      return 0;
    }

    // Check for extreme cases at the low end because the
    // normal approximation won't work
    if (e < 1) {

      // Base case (i.e. e == 0) from documenta Geigy Scientific
      // Tables, 6th edition, page 185
      double base = N * (1 - Math.pow(CF, 1 / N)); 
      if (e == 0) {
	return base; 
      }
    
      // Use linear interpolation between 0 and 1 like C4.5 does
      return base + e * (addErrs(N, 1, CF) - base);
    }
    
    // Use linear interpolation at the high end (i.e. between N - 0.5
    // and N) because of the continuity correction
    if (e + 0.5 >= N) {

      // Make sure that we never return anything smaller than zero
      return Math.max(N - e, 0);
    }

    // Get z-score corresponding to CF
    double z = Statistics.normalInverse(1 - CF);

    // Compute upper limit of confidence interval
    double  f = (e + 0.5) / N;
    double r = (f + (z * z) / (2 * N) +
		z * Math.sqrt((f / N) - 
			      (f * f / N) + 
			      (z * z / (4 * N * N)))) /
      (1 + (z * z) / N);

    return (r * N) - e;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.9 $");
  }
}
