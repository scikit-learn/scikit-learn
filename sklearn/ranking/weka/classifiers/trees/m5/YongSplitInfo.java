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
 *    YongSplitInfo.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.m5;

import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.Serializable;

/**
 * Stores split information.
 *
 * @author Yong Wang (yongwang@cs.waikato.ac.nz)
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public final class YongSplitInfo
  implements Cloneable, Serializable, SplitEvaluate, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 1864267581079767881L;

  private int  number;         // number of total instances
  private int  first;          // first instance index
  private int  last;           // last instance index
  private int  position;       // position of maximum impurity reduction
  private double  maxImpurity; // maximum impurity reduction
  private double  leftAve;     // left average class value
  private double rightAve;     // right average class value
  private int  splitAttr;      // spliting attribute 
  private double  splitValue;  // splitting value

  /**
   * Constructs an object which contains the split information
   * @param low the index of the first instance
   * @param high the index of the last instance
   * @param attr an attribute
   */
  public YongSplitInfo(int low,int high,int attr) {
    number = high-low+1;
    first = low;
    last = high;
    position = -1;
    maxImpurity = -1.e20;
    splitAttr = attr;      // attr < 0 is an empty object 
    splitValue = 0.0;
    Utils.SMALL = 1e-10;
  }

  /**
   * Makes a copy of this SplitInfo object
   */ 
  public final SplitEvaluate copy () throws Exception {

    YongSplitInfo s = (YongSplitInfo)this.clone();
    
    return s;
  }

  /**
   * Resets the object of split information
   * @param low the index of the first instance
   * @param high the index of the last instance
   * @param attr the attribute
   */
  public final void  initialize(int low,int high,int attr){
      
    number = high-low+1;
    first = low;
    last = high;
    position = -1;
    maxImpurity = -1.e20;
    splitAttr = attr;
    splitValue = 0.0;
  }

  /**
   * Converts the spliting information to string
   * @param inst the instances
   */
  public final String  toString(Instances inst){

    StringBuffer text =  new StringBuffer();

    text.append("Print SplitInfo:\n");
    text.append("    Instances:\t\t" + number + " (" + first + "-" + 
		position + "," + (position+1) + "-" + last + ")\n");
    text.append("    Maximum Impurity Reduction:\t" + 
		Utils.doubleToString(maxImpurity,1,4) + "\n");
    text.append("    Left average:\t" + leftAve + "\n");
    text.append("    Right average:\t" + rightAve + "\n");
    if(maxImpurity>0.0)
      text.append("    Splitting function:\t" + 
		  inst.attribute(splitAttr).name() + " = " 
		  + splitValue + "\n");
    else text.append("    Splitting function:\tnull\n");
    
    return text.toString();
  }
  
  /** 
   * Finds the best splitting point for an attribute in the instances
   * @param attr the splitting attribute
   * @param inst the instances
   * @exception Exception if something goes wrong
   */
  public final void  attrSplit(int attr,Instances inst) throws Exception {
    int i,len,count,part;
    Impurity imp;
    
    int low = 0;
    int high = inst.numInstances()-1;
    this.initialize(low,high,attr);
    if(number < 4) {
      return;
    }
    
    len = ((high-low+1)<5) ? 1 : (high-low+1) / 5; 
    
    position = low;
    
    part = low + len - 1;
    imp = new Impurity(part,attr,inst,5);
    
    count=0;
    for(i=low+len;i<=high-len-1;i++) {
      
      imp.incremental(inst.instance(i).classValue(),1);
      
      if(Utils.eq(inst.instance(i+1).value(attr),
		  inst.instance(i).value(attr)) == false) {
	count = i;
	if(imp.impurity > maxImpurity){
	  maxImpurity = imp.impurity;
	  splitValue = (inst.instance(i).value(attr) +
			inst.instance(i+1).value(attr)) * 0.5;
	  leftAve = imp.sl / imp.nl; 
	  rightAve = imp.sr / imp.nr; 
	  position=i;
	}
      } 
    }
  }

  /**
   * Returns the impurity of this split
   *
   * @return the impurity of this split
   */
  public double maxImpurity () {
    return maxImpurity;
  }

  /**
   * Returns the attribute used in this split
   *
   * @return the attribute used in this split
   */
  public int splitAttr () {
    return splitAttr;
  }

  /**
   * Returns the position of the split in the sorted values. -1 indicates that
   * a split could not be found.
   *
   * @return an <code>int</code> value
   */
  public int position () {
    return position;
  }

  /**
   * Returns the split value
   *
   * @return the split value
   */
  public double splitValue () {
    return splitValue;
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
