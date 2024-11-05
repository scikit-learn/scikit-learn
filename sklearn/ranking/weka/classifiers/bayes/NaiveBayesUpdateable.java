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
 *    NaiveBayesUpdateable.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.bayes;

import weka.classifiers.UpdateableClassifier;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;

/**
 <!-- globalinfo-start -->
 * Class for a Naive Bayes classifier using estimator classes. This is the updateable version of NaiveBayes.<br/>
 * This classifier will use a default precision of 0.1 for numeric attributes when buildClassifier is called with zero training instances.<br/>
 * <br/>
 * For more information on Naive Bayes classifiers, see<br/>
 * <br/>
 * George H. John, Pat Langley: Estimating Continuous Distributions in Bayesian Classifiers. In: Eleventh Conference on Uncertainty in Artificial Intelligence, San Mateo, 338-345, 1995.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{John1995,
 *    address = {San Mateo},
 *    author = {George H. John and Pat Langley},
 *    booktitle = {Eleventh Conference on Uncertainty in Artificial Intelligence},
 *    pages = {338-345},
 *    publisher = {Morgan Kaufmann},
 *    title = {Estimating Continuous Distributions in Bayesian Classifiers},
 *    year = {1995}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -K
 *  Use kernel density estimator rather than normal
 *  distribution for numeric attributes</pre>
 * 
 * <pre> -D
 *  Use supervised discretization to process numeric attributes
 * </pre>
 * 
 * <pre> -O
 *  Display model in old format (good when there are many classes)
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.11 $
 */
public class NaiveBayesUpdateable extends NaiveBayes 
  implements UpdateableClassifier {
  
  /** for serialization */
  static final long serialVersionUID = -5354015843807192221L;
 
  /**
   * Returns a string describing this classifier
   * @return a description of the classifier suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Class for a Naive Bayes classifier using estimator classes. This is the "
      +"updateable version of NaiveBayes.\n"
      +"This classifier will use a default precision of 0.1 for numeric attributes "
      +"when buildClassifier is called with zero training instances.\n\n"
      +"For more information on Naive Bayes classifiers, see\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    return super.getTechnicalInformation();
  }

  /**
   * Set whether supervised discretization is to be used.
   *
   * @param newblah true if supervised discretization is to be used.
   */
  public void setUseSupervisedDiscretization(boolean newblah) {

    if (newblah) {
      throw new IllegalArgumentException("Can't use discretization " + 
					 "in NaiveBayesUpdateable!");
    }
    m_UseDiscretization = false;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.11 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv the options
   */
  public static void main(String [] argv) {
    runClassifier(new NaiveBayesUpdateable(), argv);
  }
}

