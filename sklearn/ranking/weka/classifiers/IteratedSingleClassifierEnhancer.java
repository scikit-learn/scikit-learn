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
 *    IteratedSingleClassifierEnhancer.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers;

import weka.core.Instances;
import weka.core.Option;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/**
 * Abstract utility class for handling settings common to
 * meta classifiers that build an ensemble from a single base learner.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6041 $
 */
public abstract class IteratedSingleClassifierEnhancer
  extends SingleClassifierEnhancer {

  /** for serialization */
  private static final long serialVersionUID = -6217979135443319724L;

  /** Array for storing the generated base classifiers. */
  protected Classifier[] m_Classifiers;

  /** The number of iterations. */
  protected int m_NumIterations = 10;

  /**
   * Stump method for building the classifiers.
   *
   * @param data the training data to be used for generating the
   * bagged classifier.
   * @exception Exception if the classifier could not be built successfully
   */
  public void buildClassifier(Instances data) throws Exception {

    if (m_Classifier == null) {
      throw new Exception("A base classifier has not been specified!");
    }
    m_Classifiers = AbstractClassifier.makeCopies(m_Classifier, m_NumIterations);
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(2);

    newVector.addElement(new Option(
          "\tNumber of iterations.\n"
          + "\t(default 10)",
          "I", 1, "-I <num>"));

    Enumeration enu = super.listOptions();
    while (enu.hasMoreElements()) {
      newVector.addElement(enu.nextElement());
    }
    return newVector.elements();
  }

  /**
   * Parses a given list of options. Valid options are:<p>
   *
   * -W classname <br>
   * Specify the full class name of the base learner.<p>
   *
   * -I num <br>
   * Set the number of iterations (default 10). <p>
   *
   * Options after -- are passed to the designated classifier.<p>
   *
   * @param options the list of options as an array of strings
   * @exception Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    String iterations = Utils.getOption('I', options);
    if (iterations.length() != 0) {
      setNumIterations(Integer.parseInt(iterations));
    } else {
      setNumIterations(10);
    }

    super.setOptions(options);
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] superOptions = super.getOptions();
    String [] options = new String [superOptions.length + 2];

    int current = 0;
    options[current++] = "-I";
    options[current++] = "" + getNumIterations();

    System.arraycopy(superOptions, 0, options, current,
        superOptions.length);

    return options;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numIterationsTipText() {
    return "The number of iterations to be performed.";
  }

  /**
   * Sets the number of bagging iterations
   */
  public void setNumIterations(int numIterations) {

    m_NumIterations = numIterations;
  }

  /**
   * Gets the number of bagging iterations
   *
   * @return the maximum number of bagging iterations
   */
  public int getNumIterations() {

    return m_NumIterations;
  }
}
