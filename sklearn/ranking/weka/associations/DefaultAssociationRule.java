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
 *    DefaultAssociationRule.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import java.io.Serializable;
import java.util.Collection;

import weka.core.Tag;
import weka.core.Utils;

/**
 * Class for storing and manipulating an association rule.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com).
 */
public class DefaultAssociationRule extends AssociationRule 
  implements Serializable {
  
  /** For serialization */
  private static final long serialVersionUID = -661269018702294489L;

  /** Enum for holding different metric types */
  public static enum METRIC_TYPE {
    CONFIDENCE("conf") {
      double compute(int premiseSupport, int consequenceSupport, 
          int totalSupport, int totalTransactions) {
        
        return (double)totalSupport / (double)premiseSupport;
      }
    },
    LIFT("lift") {
      double compute(int premiseSupport, int consequenceSupport, 
          int totalSupport, int totalTransactions) {
        
        double confidence = 
          METRIC_TYPE.CONFIDENCE.compute(premiseSupport, consequenceSupport, 
              totalSupport, totalTransactions);
        return confidence / ((double)consequenceSupport /
            (double)totalTransactions);
      }
    },
    LEVERAGE("lev") {
      double compute(int premiseSupport, int consequenceSupport, 
          int totalSupport, int totalTransactions) {
        
        double coverageForItemSet = (double)totalSupport /
          (double)totalTransactions;
        double expectedCoverageIfIndependent = 
          ((double)premiseSupport / (double)totalTransactions) *
          ((double)consequenceSupport / (double)totalTransactions);
        return coverageForItemSet - expectedCoverageIfIndependent;
      }
    },
    CONVICTION("conv") {
      double compute(int premiseSupport, int consequenceSupport, 
          int totalSupport, int totalTransactions) {
        
        double num = 
          (double)premiseSupport * (double)(totalTransactions - consequenceSupport) /
          (double)totalTransactions;
        double denom = premiseSupport - totalSupport + 1;
        return num / denom;
      }
    };
    
    private final String m_stringVal;
    METRIC_TYPE(String name) {
      m_stringVal = name;
    }
    
    abstract double compute(int premiseSupport, int consequenceSupport, 
        int totalSupport, int totalTransactions);
    
    public String toString() {
      return m_stringVal;
    }
    
    public String toStringMetric(int premiseSupport, int consequenceSupport,
        int totalSupport, int totalTransactions) {
      return m_stringVal + ":(" + Utils.doubleToString(compute(premiseSupport, consequenceSupport,
          totalSupport, totalTransactions), 2) + ")";
    }
    
    public String toXML(int premiseSupport, int consequenceSupport,
        int totalSupport, int totalTransactions) {
      String result = "<CRITERE name=\"" + m_stringVal + "\" value=\" " +
        Utils.doubleToString(compute(premiseSupport, consequenceSupport,
          totalSupport, totalTransactions), 2) + "\"/>";
      
      return result;
    }
  }
  
  /** Tags for display in the GUI */
  public static final Tag[] TAGS_SELECTION = {
    new Tag(METRIC_TYPE.CONFIDENCE.ordinal(), "Confidence"),
    new Tag(METRIC_TYPE.LIFT.ordinal(), "Lift"),
    new Tag(METRIC_TYPE.LEVERAGE.ordinal(), "Leverage"),
    new Tag(METRIC_TYPE.CONVICTION.ordinal(), "Conviction")
  };
  
  /** The metric type for this rule */
  protected DefaultAssociationRule.METRIC_TYPE m_metricType = METRIC_TYPE.CONFIDENCE;
  
  /** The premise of the rule */
  protected Collection<Item> m_premise;
  
  /** The consequence of the rule */
  protected Collection<Item> m_consequence;
  
  /** The support for the premise */
  protected int m_premiseSupport;
  
  /** The support for the consequence */
  protected int m_consequenceSupport;
  
  /** The total support for the item set (premise + consequence) */
  protected int m_totalSupport;
  
  /** The total number of transactions in the data */
  protected int m_totalTransactions;
  
  /**
   * Construct a new default association rule.
   * 
   * @param premise the premise of the rule
   * @param consequence the consequence of the rule
   * @param metric the metric for the rule
   * @param premiseSupport the support of the premise
   * @param consequenceSupport the support of the consequence
   * @param totalSupport the total support of the rule
   * @param totalTransactions the number of transactions in the data
   */
  public DefaultAssociationRule(Collection<Item> premise, 
      Collection<Item> consequence, METRIC_TYPE metric,
      int premiseSupport, int consequenceSupport,
      int totalSupport, int totalTransactions) {
    m_premise = premise;
    m_consequence = consequence;
    m_metricType = metric;
    m_premiseSupport = premiseSupport;
    m_consequenceSupport = consequenceSupport;
    m_totalSupport = totalSupport;
    m_totalTransactions = totalTransactions;
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getPremise()
   */
  public Collection<Item> getPremise() {
    return m_premise;
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getConsequence()
   */
  public Collection<Item> getConsequence() {
    return m_consequence;
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getPrimaryMetricName()
   */
  public String getPrimaryMetricName() {
    return TAGS_SELECTION[m_metricType.ordinal()].getReadable();
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getPrimaryMetricValue()
   */
  public double getPrimaryMetricValue() {
    return m_metricType.compute(m_premiseSupport, m_consequenceSupport, 
        m_totalSupport, m_totalTransactions);
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getNamedMetricValue(java.lang.String)
   */
  public double getNamedMetricValue(String metricName) throws Exception {
    
    DefaultAssociationRule.METRIC_TYPE requested = null;
    
    for (DefaultAssociationRule.METRIC_TYPE m : METRIC_TYPE.values()) {
      if (TAGS_SELECTION[m.ordinal()].getReadable().equals(metricName)) {
        requested = m;
      }
    }
    
    if (requested == null) {
      throw new Exception("[AssociationRule] Unknown metric: " + metricName);
    }
    
    return requested.compute(m_premiseSupport, m_consequenceSupport, 
        m_totalSupport, m_totalTransactions);
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getNumberOfMetricsForRule()
   */
  public int getNumberOfMetricsForRule() {
    return METRIC_TYPE.values().length;
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getMetricNamesForRule()
   */
  public String[] getMetricNamesForRule() {
    String[] metricNames = new String[TAGS_SELECTION.length];
    
    for (int i = 0; i < TAGS_SELECTION.length; i++) {
      metricNames[i] = TAGS_SELECTION[i].getReadable();
    }
    
    return metricNames;
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getMetricValuesForRule()
   */
  public double[] getMetricValuesForRule() throws Exception {
    double[] values = new double[TAGS_SELECTION.length];
    
    for (int i = 0; i < TAGS_SELECTION.length; i++) {
      values[i] = getNamedMetricValue(TAGS_SELECTION[i].getReadable());
    }
    
    return values;
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getPremiseSupport()
   */
  public int getPremiseSupport() {
    return m_premiseSupport;
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getConsequenceSupport()
   */
  public int getConsequenceSupport() {
    return m_consequenceSupport;
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getTotalSupport()
   */
  public int getTotalSupport() {
    return m_totalSupport;
  }
  
  /* (non-Javadoc)
   * @see weka.associations.AssociationRule#getTotalTransactions()
   */
  public int getTotalTransactions() {
    return m_totalTransactions;
  }
              
  /**
   * Get a textual description of this rule.
   * 
   * @return a textual description of this rule.
   */
  public String toString() {
    StringBuffer result = new StringBuffer();
    
    result.append(m_premise.toString() + ": " + m_premiseSupport 
        + " ==> " + m_consequence.toString() + ": " + m_totalSupport 
        + "   ");
    for (DefaultAssociationRule.METRIC_TYPE m : METRIC_TYPE.values()) {
      if (m.equals(m_metricType)) {
        result.append("<" + 
            m.toStringMetric(m_premiseSupport, m_consequenceSupport, 
                m_totalSupport, m_totalTransactions) + "> ");
      } else {
        result.append("" + 
            m.toStringMetric(m_premiseSupport, m_consequenceSupport, 
                m_totalSupport, m_totalTransactions) + " ");
      }
    }
    return result.toString();
  }
}