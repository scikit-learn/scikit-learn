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
 *    PairedStatsCorrected.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.experiment;

import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Statistics;

/**
 * A class for storing stats on a paired comparison. This version is
 * based on the corrected resampled t-test statistic, which uses the
 * ratio of the number of test examples/the number of training examples.<p>
 *
 * For more information see:<p>
 *
 * Claude Nadeau and Yoshua Bengio, "Inference for the Generalization Error,"
 * Machine Learning, 2001.
 *
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public class PairedStatsCorrected
  extends PairedStats {

  /** The ratio used to correct the significane test */
  protected double m_testTrainRatio;

  /**
   * Creates a new PairedStatsCorrected object with the supplied
   * significance level and train/test ratio.
   *
   * @param sig the significance level for comparisons
   * @param testTrainRatio the number test examples/training examples
   */
  public PairedStatsCorrected(double sig, double testTrainRatio) {
      
    super(sig);
    m_testTrainRatio = testTrainRatio;
  }

  /**
   * Calculates the derived statistics (significance etc).
   */
  public void calculateDerived() {

    xStats.calculateDerived();
    yStats.calculateDerived();
    differencesStats.calculateDerived();

    correlation = Double.NaN;
    if (!Double.isNaN(xStats.stdDev) && !Double.isNaN(yStats.stdDev)
	&& !Utils.eq(xStats.stdDev, 0)) {
      double slope = (xySum - xStats.sum * yStats.sum / count)
	/ (xStats.sumSq - xStats.sum * xStats.mean);
      if (!Utils.eq(yStats.stdDev, 0)) {
	correlation = slope * xStats.stdDev / yStats.stdDev;
      } else {
	correlation = 1.0;
      }
    }

    if (Utils.gr(differencesStats.stdDev, 0)) {

      double tval = differencesStats.mean
	/ Math.sqrt((1 / count + m_testTrainRatio)
		    * differencesStats.stdDev * differencesStats.stdDev);
      
      if (count > 1) {
	differencesProbability = Statistics.FProbability(tval * tval, 1,
							 (int) count - 1);
      } else differencesProbability = 1;
    } else {
      if (differencesStats.sumSq == 0) {
	differencesProbability = 1.0;
      } else {
	differencesProbability = 0.0;
      }
    }
    differencesSignificance = 0;
    if (differencesProbability <= sigLevel) {
      if (xStats.mean > yStats.mean) {
	differencesSignificance = 1;
      } else {
	differencesSignificance = -1;
      }
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.5 $");
  }
}
