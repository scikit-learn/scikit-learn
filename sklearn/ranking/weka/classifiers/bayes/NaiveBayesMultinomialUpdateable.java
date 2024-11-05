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
 *    NaiveBayesMultinomialUpdateable.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *    Copyright (C) 2007 Jiang Su (incremental version)
 */

package weka.classifiers.bayes;

import weka.classifiers.UpdateableClassifier;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;

/**
 <!-- globalinfo-start -->
 * Class for building and using a multinomial Naive Bayes classifier. For more information see,<br/>
 * <br/>
 * Andrew Mccallum, Kamal Nigam: A Comparison of Event Models for Naive Bayes Text Classification. In: AAAI-98 Workshop on 'Learning for Text Categorization', 1998.<br/>
 * <br/>
 * The core equation for this classifier:<br/>
 * <br/>
 * P[Ci|D] = (P[D|Ci] x P[Ci]) / P[D] (Bayes rule)<br/>
 * <br/>
 * where Ci is class i and D is a document.<br/>
 * <br/>
 * Incremental version of the algorithm.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Mccallum1998,
 *    author = {Andrew Mccallum and Kamal Nigam},
 *    booktitle = {AAAI-98 Workshop on 'Learning for Text Categorization'},
 *    title = {A Comparison of Event Models for Naive Bayes Text Classification},
 *    year = {1998}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 <!-- options-end -->
 *
 * @author Andrew Golightly (acg4@cs.waikato.ac.nz)
 * @author Bernhard Pfahringer (bernhard@cs.waikato.ac.nz)
 * @author Jiang Su
 * @version $Revision: 1.3 $
 */
public class NaiveBayesMultinomialUpdateable
  extends NaiveBayesMultinomial
  implements UpdateableClassifier {

  /** for serialization */
  private static final long serialVersionUID = -7204398796974263186L;
  
  /** the word count per class */
  protected double[] m_wordsPerClass;
  
  /**
   * Returns a string describing this classifier
   * 
   * @return 		a description of the classifier suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return
        super.globalInfo() + "\n\n"
      + "Incremental version of the algorithm.";
  }

  /**
   * Generates the classifier.
   *
   * @param instances 	set of instances serving as training data
   * @throws Exception 	if the classifier has not been generated successfully
   */
  public void buildClassifier(Instances instances) throws Exception {
    // can classifier handle the data?
    getCapabilities().testWithFail(instances);

    // remove instances with missing class
    instances = new Instances(instances);
    instances.deleteWithMissingClass();

    m_headerInfo = new Instances(instances, 0);
    m_numClasses = instances.numClasses();
    m_numAttributes = instances.numAttributes();
    m_probOfWordGivenClass = new double[m_numClasses][];
    m_wordsPerClass = new double[m_numClasses];
    m_probOfClass = new double[m_numClasses];

    // initialising the matrix of word counts
    // NOTE: Laplace estimator introduced in case a word that does not 
    // appear for a class in the training set does so for the test set
    double laplace = 1;
    for (int c = 0; c < m_numClasses; c++) {
      m_probOfWordGivenClass[c] = new double[m_numAttributes];
      m_probOfClass[c]   = laplace;
      m_wordsPerClass[c] = laplace * m_numAttributes;
      for(int att = 0; att<m_numAttributes; att++) {
	m_probOfWordGivenClass[c][att] = laplace;
      }
    }

    for (int i = 0; i < instances.numInstances(); i++)
      updateClassifier(instances.instance(i));
  }

  /**
   * Updates the classifier with the given instance.
   *
   * @param instance 	the new training instance to include in the model
   * @throws Exception 	if the instance could not be incorporated in
   * 			the model.
   */
  public void updateClassifier(Instance instance) throws Exception {
    int classIndex = (int) instance.value(instance.classIndex());
    m_probOfClass[classIndex] += instance.weight();

    for (int a = 0; a < instance.numValues(); a++) {
      if (instance.index(a) == instance.classIndex() ||
	  instance.isMissing(a))
	continue;

      double numOccurences = instance.valueSparse(a) * instance.weight();
      if (numOccurences < 0)
	throw new Exception(
	    "Numeric attribute values must all be greater or equal to zero.");
      m_wordsPerClass[classIndex] += numOccurences;
      m_probOfWordGivenClass[classIndex][instance.index(a)] += numOccurences;
    }
  }

  /**
   * Calculates the class membership probabilities for the given test
   * instance.
   *
   * @param instance 	the instance to be classified
   * @return 		predicted class probability distribution
   * @throws Exception 	if there is a problem generating the prediction
   */
  public double[] distributionForInstance(Instance instance) throws Exception {
    double[] probOfClassGivenDoc = new double[m_numClasses];

    // calculate the array of log(Pr[D|C])
    double[] logDocGivenClass = new double[m_numClasses];
    for (int c = 0; c < m_numClasses; c++) {
      logDocGivenClass[c] += Math.log(m_probOfClass[c]);
      int allWords = 0;
      for (int i = 0; i < instance.numValues(); i++) {
	if (instance.index(i) == instance.classIndex())
	  continue;
	double frequencies = instance.valueSparse(i);
	allWords += frequencies;
	logDocGivenClass[c] += frequencies *
	Math.log(m_probOfWordGivenClass[c][instance.index(i)]);
      }
      logDocGivenClass[c] -= allWords * Math.log(m_wordsPerClass[c]);
    }

    double max = logDocGivenClass[Utils.maxIndex(logDocGivenClass)];
    for (int i = 0; i < m_numClasses; i++)
      probOfClassGivenDoc[i] = Math.exp(logDocGivenClass[i] - max);

    Utils.normalize(probOfClassGivenDoc);

    return probOfClassGivenDoc;
  }

  /**
   * Returns a string representation of the classifier.
   *
   * @return 		a string representation of the classifier
   */
  public String toString() {
    StringBuffer result = new StringBuffer();

    result.append("The independent probability of a class\n");
    result.append("--------------------------------------\n");

    for (int c = 0; c < m_numClasses; c++)
      result.append(m_headerInfo.classAttribute().value(c)).append("\t").
      append(Double.toString(m_probOfClass[c])).append("\n");

    result.append("\nThe probability of a word given the class\n");
    result.append("-----------------------------------------\n\t");

    for (int c = 0; c < m_numClasses; c++)
      result.append(m_headerInfo.classAttribute().value(c)).append("\t");

    result.append("\n");

    for (int w = 0; w < m_numAttributes; w++) {
      result.append(m_headerInfo.attribute(w).name()).append("\t");
      for (int c = 0; c < m_numClasses; c++)
	result.append(
	    Double.toString(Math.exp(m_probOfWordGivenClass[c][w]))).append("\t");
      result.append("\n");
    }

    return result.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.3 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param args 	the options
   */
  public static void main(String[] args) {
    runClassifier(new NaiveBayesMultinomialUpdateable(), args);
  }
}
