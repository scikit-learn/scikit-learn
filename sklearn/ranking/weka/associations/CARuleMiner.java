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
 * CARuleMiner.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import weka.core.FastVector;
import weka.core.Instances;
import weka.core.OptionHandler;

/** 
 * Interface for learning class association rules. All schemes for learning
 * class association rules implemement this interface.
 *
 * @author Stefan Mutter (mutter@cs.waikato.ac.nz)
 * @version $Revision: 1.3 $ 
 */
public interface CARuleMiner extends OptionHandler  {
    
    /**
     * Method for mining class association rules.
     * Must initialize all fields of the CARuleMiner that are not being set via options (ie. multiple calls of mineCARs
     * must always lead to the same result). Must not change the dataset in any way.
     * @param data the insatnces for which class association rules are mined
     * @throws Exception throws exception if class association rules cannot be mined
     * @return class association rules and their scoring metric in an FastVector array
     */    
 public FastVector[] mineCARs(Instances data) throws Exception;

 /**
  * Gets the instances without the class attribute
  * @return the instances withoput the class attribute
  */ 
 public Instances getInstancesNoClass();
 
 /**
  * Gets the class attribute and its values for all instances
  * @return the class attribute and its values for all instances
  */ 
 public Instances getInstancesOnlyClass();
 
 /**
  * Gets name of the scoring metric used for car mining
  * @return string containing the name of the scoring metric
  */ 
 public String metricString();
 
 /**
  * Sets the class index for the class association rule miner
  * @param index the class index
  */ 
 public void setClassIndex(int index);
 
}
