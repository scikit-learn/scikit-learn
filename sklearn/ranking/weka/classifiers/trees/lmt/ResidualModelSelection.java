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
 *    ResidualModelSelection.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.lmt;

import weka.classifiers.trees.j48.ClassifierSplitModel;
import weka.classifiers.trees.j48.Distribution;
import weka.classifiers.trees.j48.ModelSelection;
import weka.classifiers.trees.j48.NoSplit;
import weka.core.Instances;
import weka.core.RevisionUtils;

/**
 * Helper class for logistic model trees (weka.classifiers.trees.lmt.LMT) to implement the 
 * splitting criterion based on residuals.
 * 
 * @author Niels Landwehr
 * @version $Revision: 1.4 $
 */
public class ResidualModelSelection
  extends ModelSelection {

  /** for serialization */
  private static final long serialVersionUID = -293098783159385148L;

  /** Minimum number of instances for leaves*/
  protected int m_minNumInstances;

  /** Minimum information gain for split*/
  protected double m_minInfoGain;    

  /**
   * Constructor to create ResidualModelSelection object. 
   * @param minNumInstances minimum number of instances for leaves
   */
  public ResidualModelSelection(int minNumInstances) {
    m_minNumInstances = minNumInstances;
    m_minInfoGain = 1.0E-4;
  }

  /**Method not in use*/
  public void cleanup() {
    //method not in use
  }

  /**
   * Selects split based on residuals for the given dataset.
   */
  public final ClassifierSplitModel selectModel(Instances data, 
      double[][] dataZs, double[][] dataWs) throws Exception{

    int numAttributes = data.numAttributes();

    if (numAttributes < 2) throw new Exception("Can't select Model without non-class attribute");
    if (data.numInstances() < m_minNumInstances) return new NoSplit(new Distribution(data));


    double bestGain = -Double.MAX_VALUE;
    int bestAttribute = -1;

    //try split on every attribute
    for (int i = 0; i < numAttributes; i++) {
      if (i != data.classIndex()) {

	//build split
	ResidualSplit split = new ResidualSplit(i);	    
	split.buildClassifier(data, dataZs, dataWs);

	if (split.checkModel(m_minNumInstances)){

	  //evaluate split 
	  double gain = split.entropyGain();	
	  if (gain > bestGain) {
	    bestGain = gain;
	    bestAttribute = i;
	  }
	}
      }    	    
    }     

    if (bestGain >= m_minInfoGain){
      //return best split
      ResidualSplit split = new ResidualSplit(bestAttribute);
      split.buildClassifier(data, dataZs, dataWs);	
      return split;	    
    } else {	    
      //could not find any split with enough information gain
      return new NoSplit(new Distribution(data));	    
    }
  }

  /**Method not in use*/
  public final ClassifierSplitModel selectModel(Instances train) {
    //method not in use
    return null;
  }

  /**Method not in use*/
  public final ClassifierSplitModel selectModel(Instances train, Instances test) {
    //method not in use
    return null;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.4 $");
  }
}
