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
 *    AssociationRulesProducer.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import java.util.List;

/**
 * Interface to something that can provide a list of
 * AssociationRules.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6512 $
 */
public interface AssociationRulesProducer {
  
  /**
   * Gets the list of mined association rules.
   * 
   * @return the list of association rules discovered during mining.
   * Returns null if mining hasn't been performed yet.
   */
  AssociationRules getAssociationRules();
  
  /**
   * Gets a list of the names of the metrics output for
   * each rule. This list should be the same (in terms of
   * the names and order thereof) as that produced by
   * AssociationRule.getMetricNamesForRule().
   * 
   * @return an array of the names of the metrics available
   * for each rule learned by this producer.
   */
  String[] getRuleMetricNames();
  
  /**
   * Returns true if this AssociationRulesProducer can actually
   * produce rules. Most implementing classes will always return
   * true from this method (obviously :-)). However, an implementing
   * class that actually acts as a wrapper around things that may
   * or may not implement AssociationRulesProducer will want to
   * return false if the thing they wrap can't produce rules.
   * 
   * @return true if this producer can produce rules in its current
   * configuration
   */
  boolean canProduceRules();
}
