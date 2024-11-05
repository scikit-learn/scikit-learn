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
 *    FilteredAssociationRules.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import java.util.List;

import weka.filters.Filter;

/**
 * Class encapsulating a list of association rules and the preprocessing filter
 * that was applied before they were generated.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6522 $
 *
 */
public class FilteredAssociationRules extends AssociationRules {
  
  /** For serialization */
  private static final long serialVersionUID = -4223408305476916955L;

  protected Filter m_filter;
  
  protected AssociationRules m_wrappedRules;
  
  /**
   * Constructs a new FilteredAssociationRules.
   * 
   * @param producer a string describing the scheme that produced these rules.
   * @param filter the filter applied to preprocess the data used to learn the rules. 
   * @param rules the wrapped AssociationRules object.
   */
  public FilteredAssociationRules(String producer, Filter filter, AssociationRules rules) {
    super(null, producer);
    
    m_filter = filter;
    m_wrappedRules = rules;
  }
  
  /**
   * Constructs a new FilteredAssociationRules.
   * 
   * @param producer the scheme that produced the rules
   * @param filter the filter applied to preprocess the data used to learn the rules. 
   * @param rules the wrapped AssociationRules object.
   */
  public FilteredAssociationRules(Object producer, Filter filter, AssociationRules rules) {
    super(null, producer);
    
    m_filter = filter;
    m_wrappedRules = rules;
  }
  
  /**
   * Constructs a new FilteredAssociationRules.
   * 
   * @param filter the filter applied to preprocess the data used to learn the rules. 
   * @param rules the wrapped AssociationRules object.
   */
  public FilteredAssociationRules(Filter filter, AssociationRules rules) {
    super(null);
    
    m_filter = filter;
    m_wrappedRules = rules;
  }
  
  /**
   * Set the rules to use. Passes them to the wrapped AssociationRules object.
   * 
   * @param rules the rules to use.
   */
  public void setRules(List<AssociationRule> rules) {
    
    // delegate to our wrapped association rules
    m_wrappedRules.setRules(rules);
  }
  
  /**
   * Get the rules.
   * 
   * @return the rules.
   */
  public List<AssociationRule> getRules() {
    
    // delegate to our wrapped association rules
    return m_wrappedRules.getRules();
  }
  
  /**
   * Get the number of rules.
   * 
   * @return the number of rules.
   */
  public int getNumRules() {
    
    // delegate to our wrapped association rules
    return m_wrappedRules.getNumRules();
  }
  
  /**
   * Set the wrapped <code>AssociationRules</code> object to use.
   * 
   * @param rules the <code>AssociationRules</code> object to wrap.
   */
  public void setWrappedRules(AssociationRules rules) {
    m_wrappedRules = rules;
  }
  
  /**
   * Get the wrapped <code>AssociationRules</code> object.
   * 
   * @return the wrapped <code>AssociationRules</code> object.
   */
  public AssociationRules getWrappedRules() {
    return m_wrappedRules;
  }
}
