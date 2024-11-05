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
 *    AttributeStats.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.io.Serializable;

/**
 * A Utility class that contains summary information on an
 * the values that appear in a dataset for a particular attribute.
 *
 * @author <a href="mailto:len@reeltwo.com">Len Trigg</a>
 * @version $Revision: 5296 $
 */
public class AttributeStats
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 4434688832743939380L;
  
  //RANKING BEGIN
  public int[][] pairRankWeights;
  //RANKING END
  /** The number of int-like values */
  public int intCount = 0;
  
  /** The number of real-like values (i.e. have a fractional part) */
  public int realCount = 0;
  
  /** The number of missing values */
  public int missingCount = 0;
  
  /** The number of distinct values */
  public int distinctCount = 0;
  
  /** The number of values that only appear once */
  public int uniqueCount = 0;
  
  /** The total number of values (i.e. number of instances) */
  public int totalCount = 0;
  
  /** Stats on numeric value distributions */
  // perhaps Stats should be moved from weka.experiment to weka.core
  public weka.experiment.Stats numericStats;
  
  /** Counts of each nominal value */
  public int [] nominalCounts;
  
  /** Weight mass for each nominal value */
  public double[] nominalWeights;
    
  /**
   * Updates the counters for one more observed distinct value.
   *
   * @param value the value that has just been seen
   * @param count the number of times the value appeared
   * @param weight the weight mass of the value
   */
  protected void addDistinct(double value, int count, double weight) {
      if (count > 0) {
        if (count == 1) {
	  uniqueCount++;
	  }
        if (Utils.eq(value, (double)((int)value))) {
        	intCount += count;
        } else {
        	realCount += count;
        }
        if (nominalCounts != null) {
        	nominalCounts[(int)value] = count;
        	nominalWeights[(int)value] = weight;
        }
        if (numericStats != null) {
        	//numericStats.add(value, count);
        	numericStats.add(value, weight);
        	numericStats.calculateDerived();
        }
      }
      distinctCount++;
  }
  /**
   * Returns a human readable representation of this AttributeStats instance.
   *
   * @return a String represtinging these AttributeStats.
   */
  public String toString() {

    StringBuffer sb = new StringBuffer();
    sb.append(Utils.padLeft("Type", 4)).append(Utils.padLeft("Nom", 5));
    sb.append(Utils.padLeft("Int", 5)).append(Utils.padLeft("Real", 5));
    sb.append(Utils.padLeft("Missing", 12));
    sb.append(Utils.padLeft("Unique", 12));
    sb.append(Utils.padLeft("Dist", 6));
    if(nominalCounts != null){
      sb.append(' ');
      for (int i = 0; i < nominalCounts.length; i++) {
        sb.append(Utils.padLeft("C[" + i + "]", 5));
      }
    }
    sb.append('\n');

    long percent;
    percent = Math.round(100.0 * intCount / totalCount);
    if (nominalCounts != null) {
      sb.append(Utils.padLeft("Nom", 4)).append(' ');
      sb.append(Utils.padLeft("" + percent, 3)).append("% ");
      sb.append(Utils.padLeft("" + 0, 3)).append("% ");
    } else {
      sb.append(Utils.padLeft("Num", 4)).append(' ');
      sb.append(Utils.padLeft("" + 0, 3)).append("% ");
      sb.append(Utils.padLeft("" + percent, 3)).append("% ");
    }
    percent = Math.round(100.0 * realCount / totalCount);
    sb.append(Utils.padLeft("" + percent, 3)).append("% ");
    sb.append(Utils.padLeft("" + missingCount, 5)).append(" /");
    percent = Math.round(100.0 * missingCount / totalCount);
    sb.append(Utils.padLeft("" + percent, 3)).append("% ");
    sb.append(Utils.padLeft("" + uniqueCount, 5)).append(" /");
    percent = Math.round(100.0 * uniqueCount / totalCount);
    sb.append(Utils.padLeft("" + percent, 3)).append("% ");
    sb.append(Utils.padLeft("" + distinctCount, 5)).append(' ');
    if (nominalCounts != null) {
      for (int i = 0; i < nominalCounts.length; i++) {
        sb.append(Utils.padLeft("" + nominalCounts[i], 5));
      }
    }
    sb.append('\n');
    return sb.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5296 $");
  }
}
