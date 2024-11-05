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
 *    AssociationRule.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Abstract class for storing and manipulating an association rule.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com).
 * @version $Revision: 6485 $
 */
public abstract class AssociationRule implements Comparable<AssociationRule> {

  /**
   * Get the premise of this rule.
   * 
   * @return the premise of this rule.
   */
  public abstract Collection<Item> getPremise();

  /**
   * Get the consequence of this rule.
   * 
   * @return the consequence of this rule.
   */
  public abstract Collection<Item> getConsequence();

  /**
   * Get the name of the primary metric of this rule (e.g. confidence).
   * 
   * @return the name of the primary metric of this rule.
   */
  public abstract String getPrimaryMetricName();

  /**
   * Get the value of the metric for this rule. 
   * 
   * @return the value of the metric for this rule.
   */
  public abstract double getPrimaryMetricValue();

  /**
   * Get the value of the named metric for this rule
   * 
   * @param metricName the metric to get the value for
   * @return the value of the named metric
   * @throws Exception if the requested metric is unknown for this rule
   */
  public abstract double getNamedMetricValue(String metricName)
      throws Exception;

  /**
   * Gets the number of metrics available for this rule.
   * 
   * @return the number of metrics available for this rule
   */
  public abstract int getNumberOfMetricsForRule();

  /**
   * Return the names of the metrics available for this rule.
   * 
   * @return the names of the metrics that are available for this rule.
   */
  public abstract String[] getMetricNamesForRule();

  /**
   * Get all the available metric values for this rule. Values are
   * returned in an array with entries that correspond to the metric
   * names returned by getMetricNamesForRule().
   * 
   * @return all the available metrics for this rule.
   * @throws Exception if a metric can't be computed for some reason.
   */
  public abstract double[] getMetricValuesForRule() throws Exception;

  /**
   * Get the support for the premise.
   * 
   * @return the support for the premise.
   */
  public abstract int getPremiseSupport();

  /**
   * Get the support for the consequence.
   * 
   * @return the support for the consequence.
   */
  public abstract int getConsequenceSupport();

  /**
   * Get the total support for this rule (premise + consequence).
   * 
   * @return the total support for this rule.
   */
  public abstract int getTotalSupport();

  /**
   * Get the total number of transactions in the data.
   * 
   * @return the total number of transactions in the data.
   */
  public abstract int getTotalTransactions();
  
  /**
   * Compare this rule to the supplied rule.
   * 
   * @param other the rule to compare to.
   * @return the result of the comparison.
   */
  public int compareTo(AssociationRule other) {
    return -Double.compare(getPrimaryMetricValue(), other.getPrimaryMetricValue());
  }
  
  /**
   * Return true if this rule is equal to the supplied one.
   * 
   * @return true if this rule is the same as the supplied rule.
   */
  public boolean equals(Object other) {
    if (!(other instanceof AssociationRule)) {
      return false;
    }
    
    AssociationRule otherRule = (AssociationRule)other;
    boolean result = getPremise().equals(otherRule.getPremise()) &&
      getConsequence().equals(otherRule.getConsequence()) && 
      (getPrimaryMetricValue() == otherRule.getPrimaryMetricValue());
    
    return result;
  }
  
  public boolean containsItems(ArrayList<Item> items, boolean useOr) {
    int numItems = items.size();
    int count = 0;
    
    for (Item i : getPremise()) {
      if (items.contains(i)) {
        if (useOr) {
          return true; // can stop here
        } else {
          count++;
        }
      }
    }
    
    for (Item i : getConsequence()) {
      if (items.contains(i)) {
        if (useOr) {
          return true; // can stop here
        } else {
          count++;
        }
      }
    }
    
    if (!useOr) {
      if (count == numItems) {
        return true;
      }
    }
    
    return false;
  }
}