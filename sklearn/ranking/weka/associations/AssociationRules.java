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
 *    AssociationRules.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.associations;

import java.io.Serializable;
import java.util.List;

import weka.core.OptionHandler;
import weka.core.Utils;

/**
 * Class encapsulating a list of association rules.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6522 $
 */
public class AssociationRules implements Serializable {
  
  /** For serialization */
  private static final long serialVersionUID = 8889198755948056749L;
  
  /** The scheme that produced these rules */
  protected String m_producer = "Unknown";
  
  /** The list of rules */
  protected List<AssociationRule> m_rules;
  
  /**
   * Constructs a new AssociationRules.
   * 
   * @param rules the list of rules.
   * @param producer a string describing the scheme that produced these rules.
   */
  public AssociationRules(List<AssociationRule> rules, String producer) {
    m_rules = rules;
    m_producer = producer;
  }
  
  /**
   * Constructs a new AssociationRules.
   * 
   * @param rules the list of rules.
   * @param producer the scheme that produced the rules.
   */
  public AssociationRules(List<AssociationRule> rules, Object producer) {
    String producerString = producer.getClass().getName();
    if (producerString.startsWith("weka.associations.")) {
      producerString = producerString.substring("weka.associations.".length());
    }    
    
    if (producer instanceof OptionHandler) {
      String [] o = ((OptionHandler) producer).getOptions();
      producerString += " " + Utils.joinOptions(o);
    }
    
    m_rules = rules;
    m_producer = producerString;
  }
  
  /**
   * Constructs a new AssociationRules.
   * 
   * @param rules the list of rules.
   */
  public AssociationRules(List<AssociationRule> rules) {
    this(rules, "Unknown");
  }
  
  /**
   * Set the rules to use.
   * 
   * @param rules the rules to use.
   */
  public void setRules(List<AssociationRule> rules) {
    m_rules = rules;
  }
  
  /**
   * Get the rules.
   * 
   * @return the rules.
   */
  public List<AssociationRule> getRules() {
    return m_rules;
  }
  
  /**
   * Get the number of rules.
   * 
   * @return the number of rules.
   */
  public int getNumRules() {
    return m_rules.size();
  }
  
  
  /**
   * Set a textual description of the scheme that produced
   * these rules.
   * 
   * @param producer a textual description of the scheme that produced
   * these rules.
   */
  public void setProducer(String producer) {
    m_producer = producer;
  }
  
  /**
   * Get a string describing the scheme that produced these rules.
   * 
   * @return producer
   */
  public String getProducer() {
    return m_producer;
  }
}
