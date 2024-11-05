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
 * Test.java
 * Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.datagenerators;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.Serializable;

/** 
 * Class to represent a test. <br/>
 * <br/>
 * The string representation of the test can be supplied in standard notation
 * or for a subset of types of attributes  in Prolog notation.<br/>
 *
 * Following examples for all possible tests that can be represented by
 * this class, given in standard notation.<br/>
 * <br/>
 * Examples of tests for numeric attributes:<br/>
 * B &gt;= 2.333<br/>        B &lt; 4.56<br/>
 * <br/>
 * Examples of tests for nominal attributes with more then 2 values:<br/>
 * A = rain <br/>            A != rain<br/>
 * <br/>
 * Examples of tests for nominal attribute with exactly 2 values:<br/>
 * A = false <br/>            A = true<br/>
 * <br/>
 * <br/>
 * The Prolog notation is only supplied for numeric attributes and
 * for nominal attributes that have the values "true" and "false".<br/>
 * <br/>
 * Following examples for the Prolog notation provided.<br/>
 * <br/>
 * Examples of tests for numeric attributes:<br/>
 * The same as for standard notation above.<br/>
 * <br/>
 * Examples of tests for nominal attributes with values "true"and "false":<br/>
 * A<br/>
 * not(A)<br/>
 * <br/>
 * (Other nominal attributes are not supported by the Prolog notation.)<br/>
 * <br/>
 *
 * @author Gabi Schmidberger (gabi@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $ 
 **/

public class Test 
  implements Serializable, RevisionHandler {

  /** for serialization */
  static final long serialVersionUID = -8890645875887157782L;
  
  /** the attribute index */
  int m_AttIndex;
  
  /** the split */
  double m_Split;
  
  /** whether to negate the test */
  boolean m_Not;
  
  /** the dataset */
  Instances m_Dataset;
	
  /** 
   * Constructor 
   * 
   * @param i the attribute index
   * @param s the split
   * @param dataset the dataset
   */
  public Test(int i, double s, Instances dataset) { 
    m_AttIndex = i; 
    m_Split = s;
    m_Dataset = dataset;
   
    m_Not = false;
  }

  /** 
   * Constructor 
   * 
   * @param i the attribute index
   * @param s the split
   * @param dataset the dataset
   * @param n whether to negate the test
   */
  public Test(int i, double s, Instances dataset, boolean n) {
    m_AttIndex = i;
    m_Split = s;
    m_Dataset = dataset;
    m_Not = n;
  }

  /** 
   * Negates the test.
   *
   * @return the test itself negated
   */
  public Test getNot() { // returns a modified copy
    return new Test(m_AttIndex, m_Split, m_Dataset, m_Not ? false : true);
    }

  /**
   * Determines whether an instance passes the test.
   *
   * @param inst the instance
   * @return true if the instance satisfies the test, false otherwise
   * @throws Exception if something goes wrong
   */   
  public boolean passesTest(Instance inst) throws Exception {
    if (inst.isMissing(m_AttIndex)) return false; // missing values fail
	
    boolean isNominal = inst.attribute(m_AttIndex).isNominal() || inst.attribute(m_AttIndex).isRanking();
    double attribVal = inst.value(m_AttIndex);
    if (!m_Not) {
      if (isNominal) {
        if (((int) attribVal) != ((int) m_Split)) return false;
      }
      else if (attribVal >= m_Split) return false;
    } else {
      if (isNominal) {
	if (((int) attribVal) == ((int) m_Split)) return false;
      }
      else if (attribVal < m_Split) return false;
    }
    return true;
  }

  /**
   * Returns the test represented by a string.
   *
   * @return a string representing the test
   */   
  public String toString() {
    return (m_Dataset.attribute(m_AttIndex).name() + " " +
    testComparisonString());
  }

  /**
   * Returns the test represented by a string in Prolog notation.
   *
   * @return a string representing the test in Prolog notation
   */   
  public String toPrologString() {
    Attribute att = m_Dataset.attribute(m_AttIndex);
    StringBuffer str = new StringBuffer();
    String attName = m_Dataset.attribute(m_AttIndex).name();
    if (att.isNumeric()) {
      str = str.append(attName + " ");
      if (m_Not) str = str.append(">= " + Utils.doubleToString(m_Split, 3));
      else str = str.append("< " + Utils.doubleToString(m_Split, 3));
    } else {
      String value = att.value((int)m_Split);
    
      if (value == "false") { str = str.append("not(" + attName + ")"); }      
      else { str = str.append(attName); }
    }
  return str.toString();
  }

  /**
   * Gives a string representation of the test, starting from the comparison
   * symbol.
   *
   * @return a string representing the test
   */   
  private String testComparisonString() {
    Attribute att = m_Dataset.attribute(m_AttIndex);
    if (att.isNumeric()) {
      return ((m_Not ? ">= " : "< ") + Utils.doubleToString(m_Split,3));
    }
    else {
      if (att.numValues() != 2) 
        return ((m_Not ? "!= " : "= ") + att.value((int)m_Split));
      else return ("= " 
                   + (m_Not ?
      att.value((int)m_Split == 0 ? 1 : 0) : att.value((int)m_Split)));
    }
  }

  /**
   * Gives a string representation of the test in Prolog notation, starting
   * from the comparison symbol.
   *
   * @return a string representing the test in Prolog notation
   */   
  private String testPrologComparisonString() {
    Attribute att = m_Dataset.attribute(m_AttIndex);
    if (att.isNumeric()) {
      return ((m_Not ? ">= " : "< ") + Utils.doubleToString(m_Split,3));
    }
    else {
      if (att.numValues() != 2) 
        return ((m_Not ? "!= " : "= ") + att.value((int)m_Split));
      else return ("= " 
                   + (m_Not ? att.value((int)m_Split == 0 ? 1 : 0) 
                          : att.value((int)m_Split)));
    }
  }

  /**
   * Compares the test with the test that is given as parameter.
   *
   * @param t the test the object is compared to
   * @return true if the two Tests are equal 
   */   
  public boolean equalTo(Test t) {
    return (m_AttIndex == t.m_AttIndex && m_Split == t.m_Split && m_Not == t.m_Not);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.5 $");
  }
}

