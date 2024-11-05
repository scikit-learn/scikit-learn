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
 *    EvaluationUtils.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.evaluation;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.util.Random;

/**
 * Contains utility functions for generating lists of predictions in 
 * various manners.
 *
 * @author Len Trigg (len@reeltwo.com)
 * @version $Revision: 5928 $
 */
public class EvaluationUtils
  implements RevisionHandler {

  /** Seed used to randomize data in cross-validation */
  private int m_Seed = 1;

  /** Sets the seed for randomization during cross-validation */
  public void setSeed(int seed) { m_Seed = seed; }

  /** Gets the seed for randomization during cross-validation */
  public int getSeed() { return m_Seed; }
  
  /**
   * Generate a bunch of predictions ready for processing, by performing a
   * cross-validation on the supplied dataset.
   *
   * @param classifier the Classifier to evaluate
   * @param data the dataset
   * @param numFolds the number of folds in the cross-validation.
   * @exception Exception if an error occurs
   */
  public FastVector getCVPredictions(Classifier classifier, 
                                     Instances data, 
                                     int numFolds) 
    throws Exception {

    FastVector predictions = new FastVector();
    Instances runInstances = new Instances(data);
    Random random = new Random(m_Seed);
    runInstances.randomize(random);
    if ((runInstances.classAttribute().isNominal() || runInstances.classAttribute().isRanking()) && (numFolds > 1)) {
      runInstances.stratify(numFolds);
    }
    int inst = 0;
    for (int fold = 0; fold < numFolds; fold++) {
      Instances train = runInstances.trainCV(numFolds, fold, random);
      Instances test = runInstances.testCV(numFolds, fold);
      FastVector foldPred = getTrainTestPredictions(classifier, train, test);
      predictions.appendElements(foldPred);
    } 
    return predictions;
  }

  /**
   * Generate a bunch of predictions ready for processing, by performing a
   * evaluation on a test set after training on the given training set.
   *
   * @param classifier the Classifier to evaluate
   * @param train the training dataset
   * @param test the test dataset
   * @exception Exception if an error occurs
   */
  public FastVector getTrainTestPredictions(Classifier classifier, 
                                            Instances train, Instances test) 
    throws Exception {
    
    classifier.buildClassifier(train);
    return getTestPredictions(classifier, test);
  }

  /**
   * Generate a bunch of predictions ready for processing, by performing a
   * evaluation on a test set assuming the classifier is already trained.
   *
   * @param classifier the pre-trained Classifier to evaluate
   * @param test the test dataset
   * @exception Exception if an error occurs
   */
  public FastVector getTestPredictions(Classifier classifier, 
                                       Instances test) 
    throws Exception {
    
    FastVector predictions = new FastVector();
    for (int i = 0; i < test.numInstances(); i++) {
      if (!test.instance(i).classIsMissing()) {
        predictions.addElement(getPrediction(classifier, test.instance(i)));
      }
    }
    return predictions;
  }

  
  /**
   * Generate a single prediction for a test instance given the pre-trained
   * classifier.
   *
   * @param classifier the pre-trained Classifier to evaluate
   * @param test the test instance
   * @exception Exception if an error occurs
   */
  public Prediction getPrediction(Classifier classifier,
                                  Instance test)
    throws Exception {
   
    double actual = test.classValue();
    double [] dist = classifier.distributionForInstance(test);
    if (test.classAttribute().isNominal() || test.classAttribute().isRanking()) {
      return new NominalPrediction(actual, dist, test.weight());
    } else {
      return new NumericPrediction(actual, dist[0], test.weight());
    }
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

