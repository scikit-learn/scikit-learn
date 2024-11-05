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
 * InstanceComparator.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Serializable;
import java.util.Comparator;

/**
 * A comparator for the Instance class. it can be used with or without the
 * class label. Missing values are sorted at the beginning.<br>
 * Can be used as comparator in the sorting and binary search algorithms of
 * <code>Arrays</code> and <code>Collections</code>.
 *
 * @see     Instance
 * @author  FracPete (fracpete at cs dot waikato dot ac dot nz)
 * @version $Revision: 5953 $
 * @see     java.util.Arrays
 * @see     java.util.Collections
 */
public class InstanceComparator
  implements Comparator<Instance>, Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -6589278678230949683L;
  
  /** whether to include the class in the comparison */
  protected boolean m_IncludeClass;
  
  /**
   * initializes the comparator and includes the class in the comparison 
   */
  public InstanceComparator() {
    this(true);
  }
  
  /**
   * initializes the comparator  
   */
  public InstanceComparator(boolean includeClass) {
    super();
    setIncludeClass(includeClass);
  }
  
  /**
   * sets whether the class should be included (= TRUE) in the comparison
   * 
   * @param includeClass        whether to include the class in the comparison 
   */
  public void setIncludeClass(boolean includeClass) {
    m_IncludeClass = includeClass;
  }
  
  /**
   * returns TRUE if the class is included in the comparison
   */
  public boolean getIncludeClass() {
    return m_IncludeClass;
  }

  /**
   * compares the two instances, returns -1 if o1 is smaller than o2, 0
   * if equal and +1 if greater. The method assumes that both instance objects
   * have the same attributes, they don't have to belong to the same dataset.
   * 
   * @param o1        the first instance to compare
   * @param o2        the second instance to compare
   * @return          returns -1 if o1 is smaller than o2, 0 if equal and +1 
   *                  if greater
   */
  public int compare(Instance o1, Instance o2) {
    int         result;
    Instance    inst1;
    Instance    inst2;
    int         classindex;
    int         i;
    
    inst1 = (Instance) o1;
    inst2 = (Instance) o2;
    
    // get class index
    if (inst1.classIndex() == -1)
      classindex = inst1.numAttributes() - 1;
    else
      classindex = inst1.classIndex();

    result = 0;
    for (i = 0; i < inst1.numAttributes(); i++) {
      // exclude class?
      if (!getIncludeClass() && (i == classindex))
        continue;
   
      // comparing attribute values
      // 1. special handling if missing value (NaN) is involved:
      if (inst1.isMissing(i) || inst2.isMissing(i)) {
        if (inst1.isMissing(i) && inst2.isMissing(i)) {
          continue;
        }
        else {
          if (inst1.isMissing(i))
            result = -1;
          else
            result = 1;
          break;
        }
      }
      // 2. regular values:
      else {
        if (Utils.eq(inst1.value(i), inst2.value(i))) { 
          continue;
        }
        else {
          if (inst1.value(i) < inst2.value(i))
            result = -1;
          else
            result = 1;
          break;
        }
      }
    }
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }
  
  /**
   * for testing only. takes an ARFF-filename as first argument to perform
   * some tests. 
   */
  public static void main(String[] args) throws Exception {
    Instances       inst;
    Comparator<Instance>      comp;
    
    if (args.length == 0)
      return;
    
    // read instances
    inst = new Instances(new BufferedReader(new FileReader(args[0])));
    inst.setClassIndex(inst.numAttributes() - 1);
    
    // compare incl. class
    comp = new InstanceComparator();
    System.out.println("\nIncluding the class");
    System.out.println("comparing 1. instance with 1.: " + comp.compare(inst.instance(0), inst.instance(0)));
    System.out.println("comparing 1. instance with 2.: " + comp.compare(inst.instance(0), inst.instance(1)));
    System.out.println("comparing 2. instance with 1.: " + comp.compare(inst.instance(1), inst.instance(0)));
        
    // compare excl. class
    comp = new InstanceComparator(false);
    System.out.println("\nExcluding the class");
    System.out.println("comparing 1. instance with 1.: " + comp.compare(inst.instance(0), inst.instance(0)));
    System.out.println("comparing 1. instance with 2.: " + comp.compare(inst.instance(0), inst.instance(1)));
    System.out.println("comparing 2. instance with 1.: " + comp.compare(inst.instance(1), inst.instance(0)));
  }
}
