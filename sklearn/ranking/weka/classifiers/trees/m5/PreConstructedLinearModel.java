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
 *    RuleNode.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.m5;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.Serializable;

/**
 * This class encapsulates a linear regression function. It is a classifier
 * but does not learn the function itself, instead it is constructed with
 * coefficients and intercept obtained elsewhere. The buildClassifier method
 * must still be called however as this stores a copy of the training data's
 * header for use in printing the model to the console.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5928 $
 */
public class PreConstructedLinearModel 
  extends AbstractClassifier 
  implements Serializable {
  
  /** for serialization */
  static final long serialVersionUID = 2030974097051713247L;
  
  /** The coefficients */
  private double [] m_coefficients;

  /** The intercept */
  private double m_intercept;

  /** Holds the instances header for printing the model */
  private Instances m_instancesHeader;

  /** number of coefficients in the model */
  private int m_numParameters;

  /**
   * Constructor
   *
   * @param coeffs an array of coefficients
   * @param intercept the intercept
   */
  public PreConstructedLinearModel(double [] coeffs, double intercept) {
    m_coefficients = coeffs;
    m_intercept = intercept;
    int count = 0;
    for (int i = 0; i < coeffs.length; i++) {
      if (coeffs[i] != 0) {
	count++;
      }
    }
    m_numParameters = count;
  }

  /**
   * Builds the classifier. In this case all that is done is that a
   * copy of the training instances header is saved.
   *
   * @param instances an <code>Instances</code> value
   * @exception Exception if an error occurs
   */
  public void buildClassifier(Instances instances) throws Exception {
    m_instancesHeader = new Instances(instances, 0);
  }

  /**
   * Predicts the class of the supplied instance using the linear model.
   *
   * @param inst the instance to make a prediction for
   * @return the prediction
   * @exception Exception if an error occurs
   */
  public double classifyInstance(Instance inst) throws Exception {
    double result = 0;

    //    System.out.println(inst);
    for (int i = 0; i < m_coefficients.length; i++) {
      if (i != inst.classIndex() && !inst.isMissing(i)) {
	//	System.out.println(inst.value(i)+" "+m_coefficients[i]);
	result += m_coefficients[i] * inst.value(i);
      }
    }
    
    result += m_intercept;
    return result;
  }
  
  /**
   * Return the number of parameters (coefficients) in the linear model
   *
   * @return the number of parameters
   */
  public int numParameters() {
    return m_numParameters;
  }
  
  /**
   * Return the array of coefficients
   *
   * @return the coefficients
   */
  public double [] coefficients() {
    return m_coefficients;
  }

  /**
   * Return the intercept
   *
   * @return the intercept
   */
  public double intercept() {
    return m_intercept;
  }

  /**
   * Returns a textual description of this linear model
   *
   * @return String containing a description of this linear model
   */
  public String toString() {
    StringBuffer b = new StringBuffer();
    b.append("\n"+m_instancesHeader.classAttribute().name() + " = ");
    boolean first = true;
    for (int i = 0; i < m_coefficients.length; i++) {
      if (m_coefficients[i] != 0.0) {
	double c = m_coefficients[i];
	if (first) {
	  b.append("\n\t" + Utils.doubleToString(c, 12, 4).trim() + " * " 
		   + m_instancesHeader.attribute(i).name() + " ");
	  first = false;
	} else {
	  b.append("\n\t" + ((m_coefficients[i] < 0) ? 
			   "- " + Utils.doubleToString(Math.abs(c), 12, 4).trim() : "+ "
		   + Utils.doubleToString(Math.abs(c), 12, 4).trim()) + " * "
		   + m_instancesHeader.attribute(i).name() + " ");
	}
      }
    }
    
    b.append("\n\t" + ((m_intercept < 0) ? "- " : "+ ")
	     + Utils.doubleToString(Math.abs(m_intercept), 12, 4).trim());
    return b.toString();
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
