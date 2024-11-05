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
 *    DecisionTableHashKey.java
 *    Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.classifiers.rules;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;

/**
 * Class providing hash table keys for DecisionTable
 */
public class DecisionTableHashKey 
  implements Serializable, RevisionHandler {

  /** for serialization */
  static final long serialVersionUID = 5674163500154964602L;

  /** Array of attribute values for an instance */
  private double [] attributes;

  /** True for an index if the corresponding attribute value is missing. */
  private boolean [] missing;

  /** The key */
  private int key;

  /**
   * Constructor for a hashKey
   *
   * @param t an instance from which to generate a key
   * @param numAtts the number of attributes
   * @param ignoreClass if true treat the class as a normal attribute
   * @throws Exception if something goes wrong
   */
  public DecisionTableHashKey(Instance t, int numAtts, boolean ignoreClass) throws Exception {

    int i;
    int cindex = t.classIndex();

    key = -999;
    attributes = new double [numAtts];
    missing = new boolean [numAtts];
    for (i=0;i<numAtts;i++) {
      if (i == cindex && !ignoreClass) {
        missing[i] = true;
      } else {
        if ((missing[i] = t.isMissing(i)) == false) {
          attributes[i] = t.value(i);
        }
      }
    }
  }

  /**
   * Convert a hash entry to a string
   *
   * @param t the set of instances
   * @param maxColWidth width to make the fields
   * @return string representation of the hash entry
   */
  public String toString(Instances t, int maxColWidth) {

    int i;
    int cindex = t.classIndex();
    StringBuffer text = new StringBuffer();

    for (i=0;i<attributes.length;i++) {
      if (i != cindex) {
        if (missing[i]) {
          text.append("?");
          for (int j=0;j<maxColWidth;j++) {
            text.append(" ");
          }
        } else {
          String ss = t.attribute(i).value((int)attributes[i]);
          StringBuffer sb = new StringBuffer(ss);

          for (int j=0;j < (maxColWidth-ss.length()+1); j++) {
            sb.append(" ");
          }
          text.append(sb);
        }
      }
    }
    return text.toString();
  }

  /**
   * Constructor for a hashKey
   *
   * @param t an array of feature values
   */
  public DecisionTableHashKey(double [] t) {

    int i;
    int l = t.length;

    key = -999;
    attributes = new double [l];
    missing = new boolean [l];
    for (i=0;i<l;i++) {
      if (t[i] == Double.MAX_VALUE) {
        missing[i] = true;
      } else {
        missing[i] = false;
        attributes[i] = t[i];
      }
    }
  }

  /**
   * Calculates a hash code
   *
   * @return the hash code as an integer
   */
  public int hashCode() {

    int hv = 0;

    if (key != -999)
      return key;
    for (int i=0;i<attributes.length;i++) {
      if (missing[i]) {
        hv += (i*13);
      } else {
        hv += (i * 5 * (attributes[i]+1));
      }
    }
    if (key == -999) {
      key = hv;
    }
    return hv;
  }

  /**
   * Tests if two instances are equal
   *
   * @param b a key to compare with
   * @return true if both objects are equal
   */
  public boolean equals(Object b) {

    if ((b == null) || !(b.getClass().equals(this.getClass()))) {
      return false;
    }
    boolean ok = true;
    boolean l;
    if (b instanceof DecisionTableHashKey) {
      DecisionTableHashKey n = (DecisionTableHashKey)b;
      for (int i=0;i<attributes.length;i++) {
        l = n.missing[i];
        if (missing[i] || l) {
          if ((missing[i] && !l) || (!missing[i] && l)) {
            ok = false;
            break;
          }
        } else {
          if (attributes[i] != n.attributes[i]) {
            ok = false;
            break;
          }
        }
      }
    } else {
      return false;
    }
    return ok;
  }

  /**
   * Prints the hash code
   */
  public void print_hash_code() {
    System.out.println("Hash val: "+hashCode());
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5928 $");
  }
}
