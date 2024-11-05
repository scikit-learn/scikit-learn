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
 *    CostCurve.java
 *    Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.evaluation;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

/**
 * Generates points illustrating probablity cost tradeoffs that can be 
 * obtained by varying the threshold value between classes. For example, 
 * the typical threshold value of 0.5 means the predicted probability of 
 * "positive" must be higher than 0.5 for the instance to be predicted as 
 * "positive".
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5987 $
 */

public class CostCurve 
  implements RevisionHandler {

  /** The name of the relation used in cost curve datasets */
  public static final String RELATION_NAME = "CostCurve";

  /** attribute name: Probability Cost Function */
  public static final String PROB_COST_FUNC_NAME = "Probability Cost Function";
  /** attribute name: Normalized Expected Cost */
  public static final String NORM_EXPECTED_COST_NAME = "Normalized Expected Cost";
  /** attribute name: Threshold */
  public static final String THRESHOLD_NAME = "Threshold";

  /**
   * Calculates the performance stats for the default class and return 
   * results as a set of Instances. The
   * structure of these Instances is as follows:<p> <ul> 
   * <li> <b>Probability Cost Function </b>
   * <li> <b>Normalized Expected Cost</b>
   * <li> <b>Threshold</b> contains the probability threshold that gives
   * rise to the previous performance values. 
   * </ul> <p>
   *
   * @see TwoClassStats
   * @param predictions the predictions to base the curve on
   * @return datapoints as a set of instances, null if no predictions
   * have been made.
   */
  public Instances getCurve(FastVector predictions) {

    if (predictions.size() == 0) {
      return null;
    }
    return getCurve(predictions, 
                    ((NominalPrediction)predictions.elementAt(0))
                    .distribution().length - 1);
  }

  /**
   * Calculates the performance stats for the desired class and return 
   * results as a set of Instances.
   *
   * @param predictions the predictions to base the curve on
   * @param classIndex index of the class of interest.
   * @return datapoints as a set of instances.
   */
  public Instances getCurve(FastVector predictions, int classIndex) {

    if ((predictions.size() == 0) ||
        (((NominalPrediction)predictions.elementAt(0))
         .distribution().length <= classIndex)) {
      return null;
    }
    
    ThresholdCurve tc = new ThresholdCurve();
    Instances threshInst = tc.getCurve(predictions, classIndex);

    Instances insts = makeHeader();
    int fpind = threshInst.attribute(ThresholdCurve.FP_RATE_NAME).index();
    int tpind = threshInst.attribute(ThresholdCurve.TP_RATE_NAME).index();
    int threshind = threshInst.attribute(ThresholdCurve.THRESHOLD_NAME).index();
    
    double [] vals;
    double fpval, tpval, thresh;
    for (int i = 0; i< threshInst.numInstances(); i++) {
      fpval = threshInst.instance(i).value(fpind);
      tpval = threshInst.instance(i).value(tpind);
      thresh = threshInst.instance(i).value(threshind);
      vals = new double [3];
      vals[0] = 0; vals[1] = fpval; vals[2] = thresh;
      insts.add(new DenseInstance(1.0, vals));
      vals = new double [3];
      vals[0] = 1; vals[1] = 1.0 - tpval; vals[2] = thresh;
      insts.add(new DenseInstance(1.0, vals));
    }
    
    return insts;
  }

  /**
   * generates the header
   * 
   * @return the header
   */
  private Instances makeHeader() {

    FastVector fv = new FastVector();
    fv.addElement(new Attribute(PROB_COST_FUNC_NAME));
    fv.addElement(new Attribute(NORM_EXPECTED_COST_NAME));
    fv.addElement(new Attribute(THRESHOLD_NAME));
    return new Instances(RELATION_NAME, fv, 100);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }

  /**
   * Tests the CostCurve generation from the command line.
   * The classifier is currently hardcoded. Pipe in an arff file.
   *
   * @param args currently ignored
   */
  public static void main(String [] args) {

    try {
      
      Instances inst = new Instances(new java.io.InputStreamReader(System.in));
      
      inst.setClassIndex(inst.numAttributes() - 1);
      CostCurve cc = new CostCurve();
      EvaluationUtils eu = new EvaluationUtils();
      Classifier classifier = new weka.classifiers.functions.Logistic();
      FastVector predictions = new FastVector();
      for (int i = 0; i < 2; i++) { // Do two runs.
	eu.setSeed(i);
	predictions.appendElements(eu.getCVPredictions(classifier, inst, 10));
	//System.out.println("\n\n\n");
      }
      Instances result = cc.getCurve(predictions);
      System.out.println(result);
      
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
