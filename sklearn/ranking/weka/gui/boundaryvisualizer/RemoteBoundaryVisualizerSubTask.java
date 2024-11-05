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
 *   RemoteBoundaryVisualizerSubTask.java
 *   Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.boundaryvisualizer;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Utils;
import weka.experiment.Task;
import weka.experiment.TaskStatusInfo;

import java.util.Random;

/**
 * Class that encapsulates a sub task for distributed boundary
 * visualization. Produces probability distributions for each pixel
 * in one row of the visualization.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 5987 $
 * @since 1.0
 * @see Task
 */
public class RemoteBoundaryVisualizerSubTask implements Task {

  // status information for this sub task
  private TaskStatusInfo m_status = new TaskStatusInfo();

  // the result of this sub task
  private RemoteResult m_result;

  // which row are we doing
  private int m_rowNumber;

  // width and height of the visualization
  private int m_panelHeight;
  private int m_panelWidth;

  // the classifier to use
  private Classifier m_classifier;

  // the kernel density estimator
  private DataGenerator m_dataGenerator;

  // the training data
  private Instances m_trainingData;

  // attributes for visualizing on (fixed dimensions)
  private int m_xAttribute;
  private int m_yAttribute;

  // pixel width and height in terms of attribute values
  private double m_pixHeight;
  private double m_pixWidth;

  // min, max of these attributes
  private double m_minX;
  private double m_minY;
  private double m_maxX;
  private double m_maxY;

  // number of samples to take from each region in the fixed dimensions
  private int m_numOfSamplesPerRegion = 2;

  // number of samples per kernel = base ^ (# non-fixed dimensions)
  private int m_numOfSamplesPerGenerator;
  private double m_samplesBase = 2.0;

  // A random number generator 
  private Random m_random;

  private double [] m_weightingAttsValues;
  private boolean [] m_attsToWeightOn;
  private double [] m_vals;
  private double [] m_dist;
  private Instance m_predInst;
  
  /**
   * Set the row number for this sub task
   *
   * @param rn the row number
   */
  public void setRowNumber(int rn) {
    m_rowNumber = rn;
  }

  /**
   * Set the width of the visualization
   *
   * @param pw the width
   */
  public void setPanelWidth(int pw) {
    m_panelWidth = pw;
  }

  /**
   * Set the height of the visualization
   *
   * @param ph the height
   */
  public void setPanelHeight(int ph) {
    m_panelHeight = ph;
  }

  /**
   * Set the height of a pixel
   *
   * @param ph the height of a pixel
   */
  public void setPixHeight(double ph) {
    m_pixHeight = ph;
  }

  /**
   * Set the width of a pixel
   *
   * @param pw the width of a pixel
   */
  public void setPixWidth(double pw) {
    m_pixWidth = pw;
  }

  /**
   * Set the classifier to use
   *
   * @param dc the classifier
   */
  public void setClassifier(Classifier dc) {
    m_classifier = dc;
  }

  /**
   * Set the density estimator to use
   *
   * @param dg the density estimator
   */
  public void setDataGenerator(DataGenerator dg) {
    m_dataGenerator = dg;
  }

  /**
   * Set the training data
   *
   * @param i the training data
   */
  public void setInstances(Instances i) {
    m_trainingData = i;
  }

  /**
   * Set the minimum and maximum values of the x axis fixed dimension
   *
   * @param minx a <code>double</code> value
   * @param maxx a <code>double</code> value
   */
  public void setMinMaxX(double minx, double maxx) {
    m_minX = minx; m_maxX = maxx;
  }

  /**
   * Set the minimum and maximum values of the y axis fixed dimension
   *
   * @param miny a <code>double</code> value
   * @param maxy a <code>double</code> value
   */
  public void setMinMaxY(double miny, double maxy) {
    m_minY = miny; m_maxY = maxy;
  }

  /**
   * Set the x axis fixed dimension
   *
   * @param xatt an <code>int</code> value
   */
  public void setXAttribute(int xatt) {
    m_xAttribute = xatt;
  }

  /**
   * Set the y axis fixed dimension
   *
   * @param yatt an <code>int</code> value
   */
  public void setYAttribute(int yatt) {
    m_yAttribute = yatt;
  }

  /**
   * Set the number of points to uniformly sample from a region (fixed
   * dimensions).
   *
   * @param num an <code>int</code> value
   */
  public void setNumSamplesPerRegion(int num) {
    m_numOfSamplesPerRegion = num;
  }

  /**
   * Set the base for computing the number of samples to obtain from each
   * generator. number of samples = base ^ (# non fixed dimensions)
   *
   * @param ksb a <code>double</code> value
   */
  public void setGeneratorSamplesBase(double ksb) {
    m_samplesBase = ksb;
  }

  /**
   * Perform the sub task
   */
  public void execute() {

    m_random = new Random(m_rowNumber * 11);
    m_dataGenerator.setSeed(m_rowNumber * 11);
    m_result = new RemoteResult(m_rowNumber, m_panelWidth);
    m_status.setTaskResult(m_result);
    m_status.setExecutionStatus(TaskStatusInfo.PROCESSING);

    try {
      m_numOfSamplesPerGenerator = 
	(int)Math.pow(m_samplesBase, m_trainingData.numAttributes()-3);
      if (m_trainingData == null) {
	throw new Exception("No training data set (BoundaryPanel)");
      }
      if (m_classifier == null) {
	throw new Exception("No classifier set (BoundaryPanel)");
      }
      if (m_dataGenerator == null) {
	throw new Exception("No data generator set (BoundaryPanel)");
      }
      if (m_trainingData.attribute(m_xAttribute).isNominal() || m_trainingData.attribute(m_xAttribute).isRanking() || 
	m_trainingData.attribute(m_yAttribute).isNominal() || m_trainingData.attribute(m_yAttribute).isRanking()) {
	throw new Exception("Visualization dimensions must be numeric "
			    +"(RemoteBoundaryVisualizerSubTask)");
      }
      
      m_attsToWeightOn = new boolean[m_trainingData.numAttributes()];
      m_attsToWeightOn[m_xAttribute] = true;
      m_attsToWeightOn[m_yAttribute] = true;
      
      // generate samples
      m_weightingAttsValues = new double [m_attsToWeightOn.length];
      m_vals = new double[m_trainingData.numAttributes()];
      m_predInst = new DenseInstance(1.0, m_vals);
      m_predInst.setDataset(m_trainingData);

      System.err.println("Executing row number "+m_rowNumber);
      for (int j = 0; j < m_panelWidth; j++) {
	double [] preds = calculateRegionProbs(j, m_rowNumber);
	m_result.setLocationProbs(j, preds);
	m_result.
	  setPercentCompleted((int)(100 * ((double)j / (double)m_panelWidth)));
      }
    } catch (Exception ex) {
      m_status.setExecutionStatus(TaskStatusInfo.FAILED);
      m_status.setStatusMessage("Row "+m_rowNumber+" failed.");
      System.err.print(ex);
      return;
    }

    // finished
    m_status.setExecutionStatus(TaskStatusInfo.FINISHED);
    m_status.setStatusMessage("Row "+m_rowNumber+" completed successfully.");
  }


  private double [] calculateRegionProbs(int j, int i) throws Exception {
    double [] sumOfProbsForRegion = 
      new double [m_trainingData.classAttribute().numValues()];

    for (int u = 0; u < m_numOfSamplesPerRegion; u++) {
      
      double [] sumOfProbsForLocation = 
	new double [m_trainingData.classAttribute().numValues()];
      
      m_weightingAttsValues[m_xAttribute] = getRandomX(j);
      m_weightingAttsValues[m_yAttribute] = getRandomY(m_panelHeight-i-1);
      
      m_dataGenerator.setWeightingValues(m_weightingAttsValues);
      
      double [] weights = m_dataGenerator.getWeights();
      double sumOfWeights = Utils.sum(weights);
      int [] indices = Utils.sort(weights);
      
      // Prune 1% of weight mass
      int [] newIndices = new int[indices.length];
      double sumSoFar = 0; 
      double criticalMass = 0.99 * sumOfWeights;
      int index = weights.length - 1; int counter = 0;
      for (int z = weights.length - 1; z >= 0; z--) {
	newIndices[index--] = indices[z];
	sumSoFar += weights[indices[z]];
	counter++;
	if (sumSoFar > criticalMass) {
	  break;
	}
      }
      indices = new int[counter];
      System.arraycopy(newIndices, index + 1, indices, 0, counter);
      
      for (int z = 0; z < m_numOfSamplesPerGenerator; z++) {
        
	m_dataGenerator.setWeightingValues(m_weightingAttsValues);
	double [][] values = m_dataGenerator.generateInstances(indices);
        
	for (int q = 0; q < values.length; q++) {
	  if (values[q] != null) {
	    System.arraycopy(values[q], 0, m_vals, 0, m_vals.length);
	    m_vals[m_xAttribute] = m_weightingAttsValues[m_xAttribute];
	    m_vals[m_yAttribute] = m_weightingAttsValues[m_yAttribute];
            
	    // classify the instance
	    m_dist = m_classifier.distributionForInstance(m_predInst);

	    for (int k = 0; k < sumOfProbsForLocation.length; k++) {
	      sumOfProbsForLocation[k] += (m_dist[k] * weights[q]); 
	    }
	  }
	}
      }
      
      for (int k = 0; k < sumOfProbsForRegion.length; k++) {
	sumOfProbsForRegion[k] += (sumOfProbsForLocation[k] * sumOfWeights); 
      }
    }
    
    // average
    Utils.normalize(sumOfProbsForRegion);

    // cache
    double [] tempDist = new double[sumOfProbsForRegion.length];
    System.arraycopy(sumOfProbsForRegion, 0, tempDist, 
		     0, sumOfProbsForRegion.length);
		
    return tempDist;
  }

  /**
   * Return a random x attribute value contained within
   * the pix'th horizontal pixel
   *
   * @param pix the horizontal pixel number
   * @return a value in attribute space
   */
  private double getRandomX(int pix) {

    double minPix =  m_minX + (pix * m_pixWidth);

    return minPix + m_random.nextDouble() * m_pixWidth;
  }

  /**
   * Return a random y attribute value contained within
   * the pix'th vertical pixel
   *
   * @param pix the vertical pixel number
   * @return a value in attribute space
   */
  private double getRandomY(int pix) {
    
    double minPix = m_minY + (pix * m_pixHeight);
    
    return minPix +  m_random.nextDouble() * m_pixHeight;
  }
  
  /**
   * Return status information for this sub task
   *
   * @return a <code>TaskStatusInfo</code> value
   */
  public TaskStatusInfo getTaskStatus() {
    return m_status;
  }
}
