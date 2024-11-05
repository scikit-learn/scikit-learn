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
 *   RemoteResult.java
 *   Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.boundaryvisualizer;

import java.io.Serializable;

/**
 * Class that encapsulates a result (and progress info) for part
 * of a distributed boundary visualization. The result of a sub-task
 * is the probabilities necessary to display one row of the final 
 * visualization.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.3 $
 * @since 1.0
 * @see Serializable
 */
public class RemoteResult
  implements Serializable {

  /** for serialization */
  private static final long serialVersionUID = 1873271280044633808L;

  /** the row number that this result corresponds to */
  private int m_rowNumber;

  /** how many pixels in a row */
  private int m_rowLength;

  /** the result - ie. the probability distributions produced by the
   * classifier for this row in the visualization */
  private double [][] m_probabilities;

  /** progress on computing this row */
  private int m_percentCompleted;

  /**
   * Creates a new <code>RemoteResult</code> instance.
   *
   * @param rowNum the row number
   * @param rowLength the number of pixels in the row
   */
  public RemoteResult(int rowNum, int rowLength) {
    m_probabilities = new double[rowLength][0];
  }
  
  /**
   * Store the classifier's distribution for a particular pixel in the
   * visualization
   *
   * @param index the pixel
   * @param distribution the probability distribution from the classifier
   */
  public void setLocationProbs(int index, double [] distribution) {
    m_probabilities[index] = distribution;
  }

  /**
   * Return the probability distributions  for this row in the visualization
   *
   * @return the probability distributions
   */
  public double [][] getProbabilities() {
    return m_probabilities;
  }

  /**
   * Set the progress for this row so far
   *
   * @param pc a percent completed value
   */
  public void setPercentCompleted(int pc) {
    m_percentCompleted = pc;
  }

  /**
   * Return the progress for this row
   *
   * @return a percent completed value
   */
  public int getPercentCompleted() {
    return m_percentCompleted;
  }
}
