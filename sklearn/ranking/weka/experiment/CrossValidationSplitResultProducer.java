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
 *    CrossValidationSplitResultProducer.java
 *    Copyright (C) 1999, 2009 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.AdditionalMeasureProducer;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.File;
import java.util.Calendar;
import java.util.Enumeration;
import java.util.Random;
import java.util.TimeZone;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Carries out one split of a repeated k-fold cross-validation, using the set SplitEvaluator to generate some results. Note that the run number is actually the nth split of a repeated k-fold cross-validation, i.e. if k=10, run number 100 is the 10th fold of the 10th cross-validation run. This producer's sole purpose is to allow more fine-grained distribution of cross-validation experiments. If the class attribute is nominal, the dataset is stratified.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -X &lt;number of folds&gt;
 *  The number of folds to use for the cross-validation.
 *  (default 10)</pre>
 * 
 * <pre> -D
 * Save raw split evaluator output.</pre>
 * 
 * <pre> -O &lt;file/directory name/path&gt;
 *  The filename where raw output will be stored.
 *  If a directory name is specified then then individual
 *  outputs will be gzipped, otherwise all output will be
 *  zipped to the named file. Use in conjuction with -D. (default splitEvalutorOut.zip)</pre>
 * 
 * <pre> -W &lt;class name&gt;
 *  The full class name of a SplitEvaluator.
 *  eg: weka.experiment.ClassifierSplitEvaluator</pre>
 * 
 * <pre> 
 * Options specific to split evaluator weka.experiment.ClassifierSplitEvaluator:
 * </pre>
 * 
 * <pre> -W &lt;class name&gt;
 *  The full class name of the classifier.
 *  eg: weka.classifiers.bayes.NaiveBayes</pre>
 * 
 * <pre> -C &lt;index&gt;
 *  The index of the class for which IR statistics
 *  are to be output. (default 1)</pre>
 * 
 * <pre> -I &lt;index&gt;
 *  The index of an attribute to output in the
 *  results. This attribute should identify an
 *  instance in order to know which instances are
 *  in the test set of a cross validation. if 0
 *  no output (default 0).</pre>
 * 
 * <pre> -P
 *  Add target and prediction columns to the result
 *  for each fold.</pre>
 * 
 * <pre> 
 * Options specific to classifier weka.classifiers.rules.ZeroR:
 * </pre>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 <!-- options-end -->
 * 
 * All options after -- will be passed to the split evaluator.
 *
 * @author Len Trigg 
 * @author Eibe Frank
 * @version $Revision: 5828 $
 */
public class CrossValidationSplitResultProducer 
  extends CrossValidationResultProducer {
  
  /** for serialization */
  static final long serialVersionUID = 1403798164046795073L;
  
  /**
   * Returns a string describing this result producer
   * @return a description of the result producer suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Carries out one split of a repeated k-fold cross-validation, "
      + "using the set SplitEvaluator to generate some results. "
      + "Note that the run number is actually the nth split of a repeated "
      + "k-fold cross-validation, i.e. if k=10, run number 100 is the 10th "
      + "fold of the 10th cross-validation run. This producer's sole purpose "
      + "is to allow more fine-grained distribution of cross-validation "
      + "experiments. If the class attribute is nominal, the dataset is stratified.";
  }
  
  /**
   * Gets the keys for a specified run number. Different run
   * numbers correspond to different randomizations of the data. Keys
   * produced should be sent to the current ResultListener
   *
   * @param run the run number to get keys for.
   * @throws Exception if a problem occurs while getting the keys
   */
  public void doRunKeys(int run) throws Exception {
    if (m_Instances == null) {
      throw new Exception("No Instances set");
    }

    // Add in some fields to the key like run and fold number, dataset name
    Object [] seKey = m_SplitEvaluator.getKey();
    Object [] key = new Object [seKey.length + 3];
    key[0] = Utils.backQuoteChars(m_Instances.relationName());
    key[2] = "" + (((run - 1) % m_NumFolds) + 1);
    key[1] = "" + (((run - 1) / m_NumFolds) + 1);
    System.arraycopy(seKey, 0, key, 3, seKey.length);
    if (m_ResultListener.isResultRequired(this, key)) {
      try {
        m_ResultListener.acceptResult(this, key, null);
      } catch (Exception ex) {
        // Save the train and test datasets for debugging purposes?
        throw ex;
      }
    }
  }

  /**
   * Gets the results for a specified run number. Different run
   * numbers correspond to different randomizations of the data. Results
   * produced should be sent to the current ResultListener
   *
   * @param run the run number to get results for.
   * @throws Exception if a problem occurs while getting the results
   */
  public void doRun(int run) throws Exception {

    if (getRawOutput()) {
      if (m_ZipDest == null) {
	m_ZipDest = new OutputZipper(m_OutputFile);
      }
    }

    if (m_Instances == null) {
      throw new Exception("No Instances set");
    }

    // Compute run and fold number from given run
    int fold = (run - 1) % m_NumFolds;
    run = ((run - 1) / m_NumFolds) + 1; 
    

    // Randomize on a copy of the original dataset
    Instances runInstances = new Instances(m_Instances);
    Random random = new Random(run);
    runInstances.randomize(random);
    if (runInstances.classAttribute().isNominal() || runInstances.classAttribute().isRanking()) {
      runInstances.stratify(m_NumFolds);
    }

    // Add in some fields to the key like run and fold number, dataset name
    Object [] seKey = m_SplitEvaluator.getKey();
    Object [] key = new Object [seKey.length + 3];
    key[0] =  Utils.backQuoteChars(m_Instances.relationName());
    key[1] = "" + run;
    key[2] = "" + (fold + 1);
    System.arraycopy(seKey, 0, key, 3, seKey.length);
    if (m_ResultListener.isResultRequired(this, key)) {
      Instances train = runInstances.trainCV(m_NumFolds, fold, random);
      Instances test = runInstances.testCV(m_NumFolds, fold);
      try {
        Object [] seResults = m_SplitEvaluator.getResult(train, test);
        Object [] results = new Object [seResults.length + 1];
        results[0] = getTimestamp();
        System.arraycopy(seResults, 0, results, 1,
                         seResults.length);
        if (m_debugOutput) {
          String resultName = (""+run+"."+(fold+1)+"."
                               + Utils.backQuoteChars(runInstances.relationName())
                               +"."
                               +m_SplitEvaluator.toString()).replace(' ','_');
          resultName = Utils.removeSubstring(resultName, 
                                             "weka.classifiers.");
          resultName = Utils.removeSubstring(resultName, 
                                             "weka.filters.");
          resultName = Utils.removeSubstring(resultName, 
                                             "weka.attributeSelection.");
          m_ZipDest.zipit(m_SplitEvaluator.getRawResultOutput(), resultName);
        }
        m_ResultListener.acceptResult(this, key, results);
      } catch (Exception ex) {
        // Save the train and test datasets for debugging purposes?
        throw ex;
      }
    }
  }

  /**
   * Gets a text descrption of the result producer.
   *
   * @return a text description of the result producer.
   */
  public String toString() {

    String result = "CrossValidationSplitResultProducer: ";
    result += getCompatibilityState();
    if (m_Instances == null) {
      result += ": <null Instances>";
    } else {
      result += ": " + Utils.backQuoteChars(m_Instances.relationName());
    }
    return result;
  }

  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5828 $");
  }
} // CrossValidationSplitResultProducer

