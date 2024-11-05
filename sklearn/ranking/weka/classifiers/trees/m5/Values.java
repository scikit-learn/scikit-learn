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
 *    Values.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.classifiers.trees.m5;

import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

/**
 * Stores some statistics.
 * @author Yong Wang (yongwang@cs.waikato.ac.nz)
 * @version $Revision: 1.7 $
 */
public final class Values
  implements RevisionHandler {
  
  int  numInstances;        // number of the instances
  int  missingInstances;    // number of the instances with missing values 
  int  first;               // index of the first instance
  int  last;                // index of the last instance
  int  attr;                // attribute
  double  sum;              // sum of the instances for attribute
  double  sqrSum;           // squared sum of the instances for attribute
  double  va;               // variance
  double  sd;               // standard deviation
  

  /**
   * Constructs an object which stores some statistics of the instances such 
   *      as sum, squared sum, variance, standard deviation
   * @param low the index of the first instance
   * @param high the index of the last instance
   * @param attribute the attribute
   * @param inst the instances
   */

  public Values(int low,int high,int attribute,Instances inst){
    int i,count=0;
    double value;

    numInstances = high-low+1;
    missingInstances = 0;
    first = low;
    last = high;
    attr = attribute;
    sum=0.0;
    sqrSum=0.0;
    for(i=first;i<=last;i++){
      if(inst.instance(i).isMissing(attr)==false){
	count++;
	value = inst.instance(i).value(attr);
	sum += value;
	sqrSum += value * value;
      }
      
      if(count >1){
	va = (sqrSum - sum * sum/count)/count;
	va = Math.abs(va);
	sd = Math.sqrt(va);
      }
      else {va = 0.0;  sd = 0.0;}      
    }
  }

  /**
   * Converts the stats to a string
   * @return the converted string
   */
  public final String  toString(){
    
    StringBuffer text = new StringBuffer();

    text.append("Print statistic values of instances (" + first + "-" + last + 
		"\n");
    text.append("    Number of instances:\t" + numInstances + "\n");
    text.append("    NUmber of instances with unknowns:\t" + missingInstances +
		"\n");
    text.append("    Attribute:\t\t\t:" + attr + "\n");
    text.append("    Sum:\t\t\t" + sum + "\n");
    text.append("    Squared sum:\t\t" + sqrSum + "\n");
    text.append("    Stanard Deviation:\t\t" + sd + "\n");

    return text.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.7 $");
  }
}
