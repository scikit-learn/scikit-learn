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
 * CorrelationSplitInfo.java
 * Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.m5;

import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.experiment.PairedStats;

import java.io.Serializable;

/**
 * Finds split points using correlation.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public final class CorrelationSplitInfo
  implements Cloneable, Serializable, SplitEvaluate, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 4212734895125452770L;

  /**
   * the first instance
   */
  private int    m_first;

  /**
   * the last instance
   */
  private int    m_last;
  private int    m_position;

  /**
   * the maximum impurity reduction
   */
  private double m_maxImpurity;

  /**
   * the attribute being tested
   */
  private int    m_splitAttr;

  /**
   * the best value on which to split
   */
  private double m_splitValue;

  /**
   * the number of instances
   */
  private int    m_number;

  /**
   * Constructs an object which contains the split information
   *
   * @param low the index of the first instance
   * @param high the index of the last instance
   * @param attr an attribute
   */
  public CorrelationSplitInfo(int low, int high, int attr) {
    initialize(low, high, attr);
  }

  /**
   * Makes a copy of this CorrelationSplitInfo object
   */
  public final SplitEvaluate copy() throws Exception {
    CorrelationSplitInfo s = (CorrelationSplitInfo) this.clone();

    return s;
  } 

  /**
   * Resets the object of split information
   *
   * @param low the index of the first instance
   * @param high the index of the last instance
   * @param attr the attribute
   */
  public final void initialize(int low, int high, int attr) {
    m_number = high - low + 1;
    m_first = low;
    m_last = high;
    m_position = -1;
    m_maxImpurity = -Double.MAX_VALUE;
    m_splitAttr = attr;
    m_splitValue = 0.0;
  } 

  /**
   * Finds the best splitting point for an attribute in the instances
   *
   * @param attr the splitting attribute
   * @param inst the instances
   * @exception Exception if something goes wrong
   */
  public final void attrSplit(int attr, Instances inst) throws Exception {
    int		i;
    int		len;
    int		part;
    int		low = 0;
    int		high = inst.numInstances() - 1;
    PairedStats full = new PairedStats(0.01);
    PairedStats leftSubset = new PairedStats(0.01);
    PairedStats rightSubset = new PairedStats(0.01);
    int		classIndex = inst.classIndex();
    double      leftCorr, rightCorr;
    double      leftVar, rightVar, allVar;
    double      order = 2.0;

    initialize(low, high, attr);

    if (m_number < 4) {
      return;
    } 

    len = ((high - low + 1) < 5) ? 1 : (high - low + 1) / 5;
    m_position = low;
    part = low + len - 1;

    // prime the subsets
    for (i = low; i < len; i++) {
      full.add(inst.instance(i).value(attr), 
	       inst.instance(i).value(classIndex));
      leftSubset.add(inst.instance(i).value(attr), 
		     inst.instance(i).value(classIndex));
    } 

    for (i = len; i < inst.numInstances(); i++) {
      full.add(inst.instance(i).value(attr), 
	       inst.instance(i).value(classIndex));
      rightSubset.add(inst.instance(i).value(attr), 
		      inst.instance(i).value(classIndex));
    } 

    full.calculateDerived();

    allVar = (full.yStats.stdDev * full.yStats.stdDev);
    allVar = Math.abs(allVar);
    allVar = Math.pow(allVar, (1.0 / order));

    for (i = low + len; i < high - len - 1; i++) {
      rightSubset.subtract(inst.instance(i).value(attr), 
			   inst.instance(i).value(classIndex));
      leftSubset.add(inst.instance(i).value(attr), 
		     inst.instance(i).value(classIndex));

      if (!Utils.eq(inst.instance(i + 1).value(attr), 
		    inst.instance(i).value(attr))) {
	leftSubset.calculateDerived();
	rightSubset.calculateDerived();

	leftCorr = Math.abs(leftSubset.correlation);
	rightCorr = Math.abs(rightSubset.correlation);
	leftVar = (leftSubset.yStats.stdDev * leftSubset.yStats.stdDev);
	leftVar = Math.abs(leftVar);
	leftVar = Math.pow(leftVar, (1.0 / order));
	rightVar = (rightSubset.yStats.stdDev * rightSubset.yStats.stdDev);
	rightVar = Math.abs(rightVar);
	rightVar = Math.pow(rightVar, (1.0 / order));

	double score = allVar - ((leftSubset.count / full.count) * leftVar) 
		       - ((rightSubset.count / full.count) * rightVar);

	// score /= allVar;
	leftCorr = (leftSubset.count / full.count) * leftCorr;
	rightCorr = (rightSubset.count / full.count) * rightCorr;

	double c_score = (leftCorr + rightCorr) - Math.abs(full.correlation);

	// c_score += score;
	if (!Utils.eq(score, 0.0)) {
	  if (score > m_maxImpurity) {
	    m_maxImpurity = score;
	    m_splitValue = 
	      (inst.instance(i).value(attr) + inst.instance(i + 1)
	      .value(attr)) * 0.5;
	    m_position = i;
	  } 
	} 
      } 
    } 
  } 

  /**
   * Returns the impurity of this split
   *
   * @return the impurity of this split
   */
  public double maxImpurity() {
    return m_maxImpurity;
  } 

  /**
   * Returns the attribute used in this split
   *
   * @return the attribute used in this split
   */
  public int splitAttr() {
    return m_splitAttr;
  } 

  /**
   * Returns the position of the split in the sorted values. -1 indicates that
   * a split could not be found.
   *
   * @return an <code>int</code> value
   */
  public int position() {
    return m_position;
  } 

  /**
   * Returns the split value
   *
   * @return the split value
   */
  public double splitValue() {
    return m_splitValue;
  } 
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.4 $");
  }
}
