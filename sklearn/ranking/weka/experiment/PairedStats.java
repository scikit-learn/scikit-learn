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
 *    PairedStats.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Statistics;
import weka.core.Utils;

/**
 * A class for storing stats on a paired comparison (t-test and correlation)
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.10 $
 */
public class PairedStats
  implements RevisionHandler {
  
  /** The stats associated with the data in column 1 */
  public Stats xStats;
  
  /** The stats associated with the data in column 2 */
  public Stats yStats;
  
  /** The stats associated with the paired differences */
  public Stats differencesStats;

  /** The probability of obtaining the observed differences */
  public double differencesProbability;

  /** The correlation coefficient */
  public double correlation;

  /** The sum of the products */
  public double xySum;
  
  /** The number of data points seen */
  public double count;
  
  /**
   * A significance indicator:
   * 0 if the differences are not significant
   * > 0 if x significantly greater than y
   * < 0 if x significantly less than y
   */
  public int differencesSignificance;
  
  /** The significance level for comparisons */
  public double sigLevel;

  /** The degrees of freedom (if set programmatically) */
  protected int m_degreesOfFreedom = 0;
    
  /**
   * Creates a new PairedStats object with the supplied significance level.
   *
   * @param sig the significance level for comparisons
   */
  public PairedStats(double sig) {
      
    xStats = new Stats();
    yStats = new Stats();
    differencesStats = new Stats();
    sigLevel = sig;
  }

  /**
   * Sets the degrees of freedom (if calibration is required).
   */
  public void setDegreesOfFreedom(int d) {
   
    if (d <= 0) {
      throw new IllegalArgumentException("PairedStats: degrees of freedom must be >= 1");
    }
    m_degreesOfFreedom = d;
  }

  /**
   * Gets the degrees of freedom.
   */
  public int getDegreesOfFreedom() {

    return m_degreesOfFreedom;
  }

  /**
   * Add an observed pair of values.
   *
   * @param value1 the value from column 1
   * @param value2 the value from column 2
   */
  public void add(double value1, double value2) {

    xStats.add(value1);
    yStats.add(value2);
    differencesStats.add(value1 - value2);
    xySum += value1 * value2;
    count ++;
  }
    
  /**
   * Removes an observed pair of values.
   *
   * @param value1 the value from column 1
   * @param value2 the value from column 2
   */
  public void subtract(double value1, double value2) {

    xStats.subtract(value1);
    yStats.subtract(value2);
    differencesStats.subtract(value1 - value2);
    xySum -= value1 * value2;
    count --;
  }

    
  /**
   * Adds an array of observed pair of values.
   *
   * @param value1 the array containing values from column 1
   * @param value2 the array containing values from column 2
   */
  public void add(double value1[], double value2[]) {
    if ((value1 == null) || (value2 == null)) {
      throw new NullPointerException();
    }
    if (value1.length != value2.length) {
      throw new IllegalArgumentException("Arrays must be of the same length");
    }
    for (int i = 0; i < value1.length; i++) {
      add(value1[i], value2[i]);
    }
  }


  /**
   * Removes an array of observed pair of values.
   *
   * @param value1 the array containing values from column 1
   * @param value2 the array containing values from column 2
   */
  public void subtract(double value1[], double value2[]) {
    if ((value1 == null) || (value2 == null)) {
      throw new NullPointerException();
    }
    if (value1.length != value2.length) {
      throw new IllegalArgumentException("Arrays must be of the same length");
    }
    for (int i = 0; i < value1.length; i++) {
      subtract(value1[i], value2[i]);
    }
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
	* Math.sqrt(count)
	/ differencesStats.stdDev;

      if (m_degreesOfFreedom >= 1){
        differencesProbability = Statistics.FProbability(tval * tval, 1,
                                                         m_degreesOfFreedom);
      } else {
        if (count > 1) {
          differencesProbability = Statistics.FProbability(tval * tval, 1,
                                                           (int) count - 1);
        } else {
          differencesProbability = 1;
        }
      }
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
   * Returns statistics on the paired comparison.
   *
   * @return the t-test statistics as a string
   */
  public String toString() {

    return "Analysis for " + Utils.doubleToString(count, 0)
      + " points:\n"
      + "                "
      + "         Column 1"
      + "         Column 2"
      + "       Difference\n"
      + "Minimums        "
      + Utils.doubleToString(xStats.min, 17, 4)
      + Utils.doubleToString(yStats.min, 17, 4)
      + Utils.doubleToString(differencesStats.min, 17, 4) + '\n'
      + "Maximums        "
      + Utils.doubleToString(xStats.max, 17, 4)
      + Utils.doubleToString(yStats.max, 17, 4)
      + Utils.doubleToString(differencesStats.max, 17, 4) + '\n'
      + "Sums            "
      + Utils.doubleToString(xStats.sum, 17, 4)
      + Utils.doubleToString(yStats.sum, 17, 4)
      + Utils.doubleToString(differencesStats.sum, 17, 4) + '\n'
      + "SumSquares      "
      + Utils.doubleToString(xStats.sumSq, 17, 4)
      + Utils.doubleToString(yStats.sumSq, 17, 4)
      + Utils.doubleToString(differencesStats.sumSq, 17, 4) + '\n'
      + "Means           "
      + Utils.doubleToString(xStats.mean, 17, 4)
      + Utils.doubleToString(yStats.mean, 17, 4)
      + Utils.doubleToString(differencesStats.mean, 17, 4) + '\n'
      + "SDs             "
      + Utils.doubleToString(xStats.stdDev, 17, 4)
      + Utils.doubleToString(yStats.stdDev, 17, 4)
      + Utils.doubleToString(differencesStats.stdDev, 17, 4) + '\n'
      + "Prob(differences) "
      + Utils.doubleToString(differencesProbability, 4)
      + " (sigflag " + differencesSignificance + ")\n"
      + "Correlation       "
      + Utils.doubleToString(correlation,4) + "\n";
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.10 $");
  }

  /**
   * Tests the paired stats object from the command line.
   * reads line from stdin, expecting two values per line.
   *
   * @param args ignored.
   */
  public static void main(String [] args) {

    try {
      PairedStats ps = new PairedStats(0.05);
      java.io.LineNumberReader r = new java.io.LineNumberReader(
				   new java.io.InputStreamReader(System.in));
      String line;
      while ((line = r.readLine()) != null) {
        line = line.trim();
        if (line.equals("") || line.startsWith("@") || line.startsWith("%")) {
          continue;
        }
	java.util.StringTokenizer s 
          = new java.util.StringTokenizer(line, " ,\t\n\r\f");
	int count = 0;
	double v1 = 0, v2 = 0;
	while (s.hasMoreTokens()) {
	  double val = (new Double(s.nextToken())).doubleValue();
	  if (count == 0) {
	    v1 = val;
	  } else if (count == 1) {
	    v2 = val;
	  } else {
            System.err.println("MSG: Too many values in line \"" 
                               + line + "\", skipped.");
	    break;
	  }
	  count++;
	}
        if (count == 2) {
          ps.add(v1, v2);
        }
      }
      ps.calculateDerived();
      System.err.println(ps);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
} // PairedStats


