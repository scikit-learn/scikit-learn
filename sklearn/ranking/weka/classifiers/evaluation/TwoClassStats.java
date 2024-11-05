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
 *    TwoClassStats.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.evaluation;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

/**
 * Encapsulates performance functions for two-class problems.
 *
 * @author Len Trigg (len@reeltwo.com)
 * @version $Revision: 1.9 $
 */
public class TwoClassStats
  implements RevisionHandler {

  /** The names used when converting this object to a confusion matrix */
  private static final String [] CATEGORY_NAMES = {"negative", "positive"};

  /** Pos predicted as pos */
  private double m_TruePos;

  /** Neg predicted as pos */
  private double m_FalsePos;

  /** Neg predicted as neg */
  private double m_TrueNeg;

  /** Pos predicted as neg */
  private double m_FalseNeg;

  /**
   * Creates the TwoClassStats with the given initial performance values.
   *
   * @param tp the number of correctly classified positives
   * @param fp the number of incorrectly classified negatives
   * @param tn the number of correctly classified negatives
   * @param fn the number of incorrectly classified positives
   */
  public TwoClassStats(double tp, double fp, double tn, double fn) {
      
    setTruePositive(tp); 
    setFalsePositive(fp);
    setTrueNegative(tn); 
    setFalseNegative(fn);
  }

  /** Sets the number of positive instances predicted as positive */
  public void setTruePositive(double tp) { m_TruePos = tp; }

  /** Sets the number of negative instances predicted as positive */
  public void setFalsePositive(double fp) { m_FalsePos = fp; }

  /** Sets the number of negative instances predicted as negative */
  public void setTrueNegative(double tn) { m_TrueNeg = tn; }

  /** Sets the number of positive instances predicted as negative */
  public void setFalseNegative(double fn) { m_FalseNeg = fn; }

  /** Gets the number of positive instances predicted as positive */
  public double getTruePositive() { return m_TruePos; }

  /** Gets the number of negative instances predicted as positive */
  public double getFalsePositive() { return m_FalsePos; }

  /** Gets the number of negative instances predicted as negative */
  public double getTrueNegative() { return m_TrueNeg; }

  /** Gets the number of positive instances predicted as negative */
  public double getFalseNegative() { return m_FalseNeg; }

  /**
   * Calculate the true positive rate. 
   * This is defined as<p>
   * <pre>
   * correctly classified positives
   * ------------------------------
   *       total positives
   * </pre>
   *
   * @return the true positive rate
   */
  public double getTruePositiveRate() { 
    if (0 == (m_TruePos + m_FalseNeg)) {
      return 0;
    } else {
      return m_TruePos / (m_TruePos + m_FalseNeg); 
    }
  }

  /**
   * Calculate the false positive rate. 
   * This is defined as<p>
   * <pre>
   * incorrectly classified negatives
   * --------------------------------
   *        total negatives
   * </pre>
   *
   * @return the false positive rate
   */
  public double getFalsePositiveRate() { 
    if (0 == (m_FalsePos + m_TrueNeg)) {
      return 0;
    } else {
      return m_FalsePos / (m_FalsePos + m_TrueNeg); 
    }
  }

  /**
   * Calculate the precision. 
   * This is defined as<p>
   * <pre>
   * correctly classified positives
   * ------------------------------
   *  total predicted as positive
   * </pre>
   *
   * @return the precision
   */
  public double getPrecision() { 
    if (0 == (m_TruePos + m_FalsePos)) {
      return 0;
    } else {
      return m_TruePos / (m_TruePos + m_FalsePos); 
    }
  }

  /**
   * Calculate the recall. 
   * This is defined as<p>
   * <pre>
   * correctly classified positives
   * ------------------------------
   *       total positives
   * </pre><p>
   * (Which is also the same as the truePositiveRate.)
   *
   * @return the recall
   */
  public double getRecall() { return getTruePositiveRate(); }

  /**
   * Calculate the F-Measure. 
   * This is defined as<p>
   * <pre>
   * 2 * recall * precision
   * ----------------------
   *   recall + precision
   * </pre>
   *
   * @return the F-Measure
   */
  public double getFMeasure() {

    double precision = getPrecision();
    double recall = getRecall();
    if ((precision + recall) == 0) {
      return 0;
    }
    return 2 * precision * recall / (precision + recall);
  }

  /**
   * Calculate the fallout. 
   * This is defined as<p>
   * <pre>
   * incorrectly classified negatives
   * --------------------------------
   *   total predicted as positive
   * </pre>
   *
   * @return the fallout
   */
  public double getFallout() { 
    if (0 == (m_TruePos + m_FalsePos)) {
      return 0;
    } else {
      return m_FalsePos / (m_TruePos + m_FalsePos); 
    }
  }

  /**
   * Generates a <code>ConfusionMatrix</code> representing the current
   * two-class statistics, using class names "negative" and "positive".
   *
   * @return a <code>ConfusionMatrix</code>.
   */
  public ConfusionMatrix getConfusionMatrix() {

    ConfusionMatrix cm = new ConfusionMatrix(CATEGORY_NAMES);
    cm.setElement(0, 0, m_TrueNeg);
    cm.setElement(0, 1, m_FalsePos);
    cm.setElement(1, 0, m_FalseNeg);
    cm.setElement(1, 1, m_TruePos);
    return cm;
  }

  /**
   * Returns a string containing the various performance measures
   * for the current object 
   */
  public String toString() {

    StringBuffer res = new StringBuffer();
    res.append(getTruePositive()).append(' ');
    res.append(getFalseNegative()).append(' ');
    res.append(getTrueNegative()).append(' ');
    res.append(getFalsePositive()).append(' ');
    res.append(getFalsePositiveRate()).append(' ');
    res.append(getTruePositiveRate()).append(' ');
    res.append(getPrecision()).append(' ');
    res.append(getRecall()).append(' ');
    res.append(getFMeasure()).append(' ');
    res.append(getFallout()).append(' ');
    return res.toString();
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
