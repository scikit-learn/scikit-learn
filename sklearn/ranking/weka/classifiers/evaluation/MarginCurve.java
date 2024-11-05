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
 *    MarginCurve.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.evaluation;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

/**
 * Generates points illustrating the prediction margin. The margin is defined
 * as the difference between the probability predicted for the actual class and
 * the highest probability predicted for the other classes. One hypothesis
 * as to the good performance of boosting algorithms is that they increaes the
 * margins on the training data and this gives better performance on test data.
 *
 * @author Len Trigg (len@reeltwo.com)
 * @version $Revision: 5987 $
 */
public class MarginCurve
  implements RevisionHandler {

  /**
   * Calculates the cumulative margin distribution for the set of
   * predictions, returning the result as a set of Instances. The
   * structure of these Instances is as follows:<p> <ul> 
   * <li> <b>Margin</b> contains the margin value (which should be plotted
   * as an x-coordinate) 
   * <li> <b>Current</b> contains the count of instances with the current 
   * margin (plot as y axis)
   * <li> <b>Cumulative</b> contains the count of instances with margin
   * less than or equal to the current margin (plot as y axis)
   * </ul> <p>
   *
   * @return datapoints as a set of instances, null if no predictions
   * have been made.  
   */
  public Instances getCurve(FastVector predictions) {

    if (predictions.size() == 0) {
      return null;
    }

    Instances insts = makeHeader();
    double [] margins = getMargins(predictions);
    int [] sorted = Utils.sort(margins);
    int binMargin = 0;
    int totalMargin = 0;
    insts.add(makeInstance(-1, binMargin, totalMargin));
    for (int i = 0; i < sorted.length; i++) {
      double current = margins[sorted[i]];
      double weight = ((NominalPrediction)predictions.elementAt(sorted[i]))
        .weight();
      totalMargin += weight;
      binMargin += weight;
      if (true) {
        insts.add(makeInstance(current, binMargin, totalMargin));
        binMargin = 0;
      }
    }
    return insts;
  }

  /**
   * Pulls all the margin values out of a vector of NominalPredictions.
   *
   * @param predictions a FastVector containing NominalPredictions
   * @return an array of margin values.
   */
  private double [] getMargins(FastVector predictions) {

    // sort by predicted probability of the desired class.
    double [] margins = new double [predictions.size()];
    for (int i = 0; i < margins.length; i++) {
      NominalPrediction pred = (NominalPrediction)predictions.elementAt(i);
      margins[i] = pred.margin();
    }
    return margins;
  }

  /**
   * Creates an Instances object with the attributes we will be calculating.
   *
   * @return the Instances structure.
   */
  private Instances makeHeader() {

    FastVector fv = new FastVector();
    fv.addElement(new Attribute("Margin"));
    fv.addElement(new Attribute("Current"));
    fv.addElement(new Attribute("Cumulative"));
    return new Instances("MarginCurve", fv, 100);
  }
  
  /**
   * Creates an Instance object with the attributes calculated.
   *
   * @param margin the margin for this data point.
   * @param current the number of instances with this margin.
   * @param cumulative the number of instances with margin less than or equal
   * to this margin.
   * @return the Instance object.
   */
  private Instance makeInstance(double margin, int current, int cumulative) {

    int count = 0;
    double [] vals = new double[3];
    vals[count++] = margin;
    vals[count++] = current;
    vals[count++] = cumulative;
    return new DenseInstance(1.0, vals);
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
   * Tests the MarginCurve generation from the command line.
   * The classifier is currently hardcoded. Pipe in an arff file.
   *
   * @param args currently ignored
   */
  public static void main(String [] args) {

    try {
      Utils.SMALL = 0;
      Instances inst = new Instances(new java.io.InputStreamReader(System.in));
      inst.setClassIndex(inst.numAttributes() - 1);
      MarginCurve tc = new MarginCurve();
      EvaluationUtils eu = new EvaluationUtils();
      weka.classifiers.meta.LogitBoost classifier 
        = new weka.classifiers.meta.LogitBoost();
      classifier.setNumIterations(20);
      FastVector predictions 
        = eu.getTrainTestPredictions(classifier, inst, inst);
      Instances result = tc.getCurve(predictions);
      System.out.println(result);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
