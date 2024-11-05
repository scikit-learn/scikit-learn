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
 *    PredictionAppender.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.clusterers.DensityBasedClusterer;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;

import java.awt.BorderLayout;
import java.beans.EventSetDescriptor;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

import javax.swing.JPanel;

/**
 * Bean that can can accept batch or incremental classifier events
 * and produce dataset or instance events which contain instances with
 * predictions appended.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 6804 $
 */
public class PredictionAppender
  extends JPanel
  implements DataSource, TrainingSetProducer, TestSetProducer, Visible, BeanCommon,
	     EventConstraints, BatchClassifierListener,
	     IncrementalClassifierListener, BatchClustererListener, Serializable {

  /** for serialization */
  private static final long serialVersionUID = -2987740065058976673L;

  /**
   * Objects listenening for dataset events
   */
  protected Vector m_dataSourceListeners = new Vector();

  /**
   * Objects listening for instances events
   */
  protected Vector m_instanceListeners = new Vector();
  
  /**
   * Objects listening for training set events
   */
  protected Vector m_trainingSetListeners = new Vector();;
  
  /**
   * Objects listening for test set events
   */
  protected Vector m_testSetListeners = new Vector();

  /**
   * Non null if this object is a target for any events.
   */
  protected Object m_listenee = null;

  /**
   * Format of instances to be produced.
   */
  protected Instances m_format;

  protected BeanVisual m_visual = 
    new BeanVisual("PredictionAppender", 
		   BeanVisual.ICON_PATH+"PredictionAppender.gif",
		   BeanVisual.ICON_PATH+"PredictionAppender_animated.gif");

  /**
   * Append classifier's predicted probabilities (if the class is discrete
   * and the classifier is a distribution classifier)
   */
  protected boolean m_appendProbabilities;

  protected transient weka.gui.Logger m_logger;

  /**
   * Global description of this bean
   *
   * @return a <code>String</code> value
   */
  public String globalInfo() {
    return "Accepts batch or incremental classifier events and "
      +"produces a new data set with classifier predictions appended.";
  }

  /**
   * Creates a new <code>PredictionAppender</code> instance.
   */
  public PredictionAppender() {
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
  }

  /**
   * Set a custom (descriptive) name for this bean
   * 
   * @param name the name to use
   */
  public void setCustomName(String name) {
    m_visual.setText(name);
  }

  /**
   * Get the custom (descriptive) name for this bean (if one has been set)
   * 
   * @return the custom name (or the default name)
   */
  public String getCustomName() {
    return m_visual.getText();
  }

  /**
   * Return a tip text suitable for displaying in a GUI
   *
   * @return a <code>String</code> value
   */
  public String appendPredictedProbabilitiesTipText() {
    return "append probabilities rather than labels for discrete class "
      +"predictions";
  }

  /**
   * Return true if predicted probabilities are to be appended rather
   * than class value
   *
   * @return a <code>boolean</code> value
   */
  public boolean getAppendPredictedProbabilities() {
    return m_appendProbabilities;
  }

  /**
   * Set whether to append predicted probabilities rather than
   * class value (for discrete class data sets)
   *
   * @param ap a <code>boolean</code> value
   */
  public void setAppendPredictedProbabilities(boolean ap) {
    m_appendProbabilities = ap;
  }

  /**
   * Add a training set listener
   *
   * @param tsl a <code>TrainingSetListener</code> value
   */
  public void addTrainingSetListener(TrainingSetListener tsl) {
    // TODO Auto-generated method stub
    m_trainingSetListeners.addElement(tsl);
    // pass on any format that we might have determined so far
    if (m_format != null) {
      TrainingSetEvent e = new TrainingSetEvent(this, m_format);
      tsl.acceptTrainingSet(e);
    }
  }

  /**
   * Remove a training set listener
   *
   * @param tsl a <code>TrainingSetListener</code> value
   */
  public void removeTrainingSetListener(TrainingSetListener tsl) {   
    m_trainingSetListeners.removeElement(tsl);
  }

  /**
   * Add a test set listener
   *
   * @param tsl a <code>TestSetListener</code> value
   */
  public void addTestSetListener(TestSetListener tsl) {
    m_testSetListeners.addElement(tsl);
//  pass on any format that we might have determined so far
    if (m_format != null) {
      TestSetEvent e = new TestSetEvent(this, m_format);
      tsl.acceptTestSet(e);
    }
  }

  /**
   * Remove a test set listener
   *
   * @param tsl a <code>TestSetListener</code> value
   */
  public void removeTestSetListener(TestSetListener tsl) {
    m_testSetListeners.removeElement(tsl);
  }
  
  /**
   * Add a datasource listener
   *
   * @param dsl a <code>DataSourceListener</code> value
   */
  public synchronized void addDataSourceListener(DataSourceListener dsl) {
    m_dataSourceListeners.addElement(dsl);
    // pass on any format that we might have determined so far
    if (m_format != null) {
      DataSetEvent e = new DataSetEvent(this, m_format);
      dsl.acceptDataSet(e);
    }
  }
  
  /**
   * Remove a datasource listener
   *
   * @param dsl a <code>DataSourceListener</code> value
   */
  public synchronized void removeDataSourceListener(DataSourceListener dsl) {
    m_dataSourceListeners.remove(dsl);
  }

  /**
   * Add an instance listener
   *
   * @param dsl a <code>InstanceListener</code> value
   */
  public synchronized void addInstanceListener(InstanceListener dsl) {
    m_instanceListeners.addElement(dsl);
    // pass on any format that we might have determined so far
    if (m_format != null) {
      InstanceEvent e = new InstanceEvent(this, m_format);
      dsl.acceptInstance(e);
    }
  }
  
  /**
   * Remove an instance listener
   *
   * @param dsl a <code>InstanceListener</code> value
   */
  public synchronized void removeInstanceListener(InstanceListener dsl) {
    m_instanceListeners.remove(dsl);
  }

  /**
   * Set the visual for this data source
   *
   * @param newVisual a <code>BeanVisual</code> value
   */
  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  /**
   * Get the visual being used by this data source.
   *
   */
  public BeanVisual getVisual() {
    return m_visual;
  }

  /**
   * Use the default images for a data source
   *
   */
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"PredictionAppender.gif",
		       BeanVisual.ICON_PATH+"PredictionAppender_animated.gif");
  }

  protected InstanceEvent m_instanceEvent;

  
  /**
   * Accept and process an incremental classifier event
   *
   * @param e an <code>IncrementalClassifierEvent</code> value
   */
  public void acceptClassifier(IncrementalClassifierEvent e) {
    weka.classifiers.Classifier classifier = e.getClassifier();
    Instance currentI = e.getCurrentInstance();
    int status = e.getStatus();
    int oldNumAtts = 0;
    if (status == IncrementalClassifierEvent.NEW_BATCH) {
      oldNumAtts = e.getStructure().numAttributes();
    } else {
      oldNumAtts = currentI.dataset().numAttributes();
    }
    if (status == IncrementalClassifierEvent.NEW_BATCH) {
      m_instanceEvent = new InstanceEvent(this, null, 0);
      // create new header structure
      Instances oldStructure = new Instances(e.getStructure(), 0);
      //String relationNameModifier = oldStructure.relationName()
	//+"_with predictions";
      String relationNameModifier = "_with predictions";
	//+"_with predictions";
       if (!m_appendProbabilities 
	   || oldStructure.classAttribute().isNumeric()) {
	 try {
	   m_format = makeDataSetClass(oldStructure, oldStructure, classifier,
						     relationNameModifier);
	 } catch (Exception ex) {
	   ex.printStackTrace();
	   return;
	 }
       } else if (m_appendProbabilities) {
	 try {
	   m_format = 
	     makeDataSetProbabilities(oldStructure, oldStructure, classifier,
				      relationNameModifier);

	 } catch (Exception ex) {
	   ex.printStackTrace();
	   return;
	 }
       }
       // Pass on the structure
       m_instanceEvent.setStructure(m_format);
       notifyInstanceAvailable(m_instanceEvent);
       return;
    }

    double[] instanceVals = new double [m_format.numAttributes()];
    Instance newInst = null;
    try {
      // process the actual instance
      for (int i = 0; i < oldNumAtts; i++) {
	instanceVals[i] = currentI.value(i);
      }
      if (!m_appendProbabilities 
	  || currentI.dataset().classAttribute().isNumeric()) {
	double predClass = 
	  classifier.classifyInstance(currentI);
	instanceVals[instanceVals.length - 1] = predClass;
      } else if (m_appendProbabilities) {
	double [] preds = classifier.distributionForInstance(currentI);
	for (int i = oldNumAtts; i < instanceVals.length; i++) {
	  instanceVals[i] = preds[i-oldNumAtts];
	}      
      }      
    } catch (Exception ex) {
      ex.printStackTrace();
      return;
    } finally {
      newInst = new DenseInstance(currentI.weight(), instanceVals);
      newInst.setDataset(m_format);
      m_instanceEvent.setInstance(newInst);
      m_instanceEvent.setStatus(status);
      // notify listeners
      notifyInstanceAvailable(m_instanceEvent);
    }

    if (status == IncrementalClassifierEvent.BATCH_FINISHED) {
      // clean up
      //      m_incrementalStructure = null;
      m_instanceEvent = null;
    }
  }

  /**
   * Accept and process a batch classifier event
   *
   * @param e a <code>BatchClassifierEvent</code> value
   */
  public void acceptClassifier(BatchClassifierEvent e) {
    if (m_dataSourceListeners.size() > 0 
	|| m_trainingSetListeners.size() > 0
	|| m_testSetListeners.size() > 0) {

      if (e.getTestSet() == null) {
        // can't append predictions
        return;
      }

      Instances testSet = e.getTestSet().getDataSet();
      Instances trainSet = e.getTrainSet().getDataSet();
      int setNum = e.getSetNumber();
      int maxNum = e.getMaxSetNumber();

      weka.classifiers.Classifier classifier = e.getClassifier();
      String relationNameModifier = "_set_"+e.getSetNumber()+"_of_"
	+e.getMaxSetNumber();
      if (!m_appendProbabilities || testSet.classAttribute().isNumeric()) {
	try {
	  Instances newTestSetInstances = makeDataSetClass(testSet, trainSet, 
	      classifier, relationNameModifier);
	  Instances newTrainingSetInstances = makeDataSetClass(trainSet, trainSet, 
	      classifier, relationNameModifier);
	  
	  if (m_trainingSetListeners.size() > 0) {
	    TrainingSetEvent tse = new TrainingSetEvent(this,
		new Instances(newTrainingSetInstances, 0));
	    tse.m_setNumber = setNum;
	    tse.m_maxSetNumber = maxNum;
	    notifyTrainingSetAvailable(tse);
	    // fill in predicted values
            for (int i = 0; i < trainSet.numInstances(); i++) {
              double predClass = 
        	classifier.classifyInstance(trainSet.instance(i));
              newTrainingSetInstances.instance(i).setValue(newTrainingSetInstances.numAttributes()-1,
        	  predClass);
            }
            tse = new TrainingSetEvent(this,
        	newTrainingSetInstances);
            tse.m_setNumber = setNum;
            tse.m_maxSetNumber = maxNum;
            notifyTrainingSetAvailable(tse);
	  }
	  
	  if (m_testSetListeners.size() > 0) {
	    TestSetEvent tse = new TestSetEvent(this,
		new Instances(newTestSetInstances, 0));
	    tse.m_setNumber = setNum;
	    tse.m_maxSetNumber = maxNum;
	    notifyTestSetAvailable(tse);
	  }
	  if (m_dataSourceListeners.size() > 0) {
	    notifyDataSetAvailable(new DataSetEvent(this, new Instances(newTestSetInstances,0)));
	  }
          if (e.getTestSet().isStructureOnly()) {
	    m_format = newTestSetInstances;
	  }
          if (m_dataSourceListeners.size() > 0 || m_testSetListeners.size() > 0) {
            // fill in predicted values
            for (int i = 0; i < testSet.numInstances(); i++) {
              Instance tempInst = testSet.instance(i);
              
              // if the class value is missing, then copy the instance
              // and set the data set to the training data. This is
              // just in case this test data was loaded from a CSV file
              // with all missing values for a nominal class (in this
              // case we have no information on the legal class values
              // in the test data)
              if (tempInst.isMissing(tempInst.classIndex()) && 
                  !(classifier instanceof weka.classifiers.misc.InputMappedClassifier)) {
                tempInst = (Instance)testSet.instance(i).copy();
                tempInst.setDataset(trainSet);
              }
              double predClass = 
        	classifier.classifyInstance(tempInst);
              newTestSetInstances.instance(i).setValue(newTestSetInstances.numAttributes()-1,
        	  predClass);
            }
          }
	  // notify listeners
          if (m_testSetListeners.size() > 0) {
            TestSetEvent tse = new TestSetEvent(this, newTestSetInstances);
            tse.m_setNumber = setNum;
            tse.m_maxSetNumber = maxNum;
            notifyTestSetAvailable(tse);
          }
          if (m_dataSourceListeners.size() > 0) {
            notifyDataSetAvailable(new DataSetEvent(this, newTestSetInstances));            
          }
	  return;
	} catch (Exception ex) {
	  ex.printStackTrace();
	}
      }
      if (m_appendProbabilities) {
	try {
	  Instances newTestSetInstances = 
	    makeDataSetProbabilities(testSet, trainSet,
				     classifier,relationNameModifier);
	  Instances newTrainingSetInstances = 
	    makeDataSetProbabilities(trainSet, trainSet,
				     classifier,relationNameModifier);
	  if (m_trainingSetListeners.size() > 0) {
	    TrainingSetEvent tse = new TrainingSetEvent(this,
		new Instances(newTrainingSetInstances, 0));
	    tse.m_setNumber = setNum;
	    tse.m_maxSetNumber = maxNum;
	    notifyTrainingSetAvailable(tse);
//	    fill in predicted probabilities
	    for (int i = 0; i < trainSet.numInstances(); i++) {
	      double [] preds = classifier.
	      distributionForInstance(trainSet.instance(i));
	      for (int j = 0; j < trainSet.classAttribute().numValues(); j++) {
		newTrainingSetInstances.instance(i).setValue(trainSet.numAttributes()+j,
		    preds[j]);
	      }
	    }
	    tse = new TrainingSetEvent(this,
        	newTrainingSetInstances);
            tse.m_setNumber = setNum;
            tse.m_maxSetNumber = maxNum;
            notifyTrainingSetAvailable(tse);
	  }
	  if (m_testSetListeners.size() > 0) {
	    TestSetEvent tse = new TestSetEvent(this,
		new Instances(newTestSetInstances, 0));
	    tse.m_setNumber = setNum;
	    tse.m_maxSetNumber = maxNum;
	    notifyTestSetAvailable(tse);
	  }
	  if (m_dataSourceListeners.size() > 0) {
	    notifyDataSetAvailable(new DataSetEvent(this, new Instances(newTestSetInstances,0)));
	  }
          if (e.getTestSet().isStructureOnly()) {
	    m_format = newTestSetInstances;
	  }
          if (m_dataSourceListeners.size() > 0 || m_testSetListeners.size() > 0) {
            // fill in predicted probabilities
            for (int i = 0; i < testSet.numInstances(); i++) {
              Instance tempInst = testSet.instance(i);
              
              // if the class value is missing, then copy the instance
              // and set the data set to the training data. This is
              // just in case this test data was loaded from a CSV file
              // with all missing values for a nominal class (in this
              // case we have no information on the legal class values
              // in the test data)
              if (tempInst.isMissing(tempInst.classIndex()) && 
                  !(classifier instanceof weka.classifiers.misc.InputMappedClassifier)) {
                tempInst = (Instance)testSet.instance(i).copy();
                tempInst.setDataset(trainSet);
              }
              
              double [] preds = classifier.
              distributionForInstance(tempInst);
              for (int j = 0; j < tempInst.classAttribute().numValues(); j++) {
        	newTestSetInstances.instance(i).setValue(testSet.numAttributes()+j,
        	    preds[j]);
              }
            }
          }
          
          // notify listeners
          if (m_testSetListeners.size() > 0) {
            TestSetEvent tse = new TestSetEvent(this, newTestSetInstances);
            tse.m_setNumber = setNum;
            tse.m_maxSetNumber = maxNum;
            notifyTestSetAvailable(tse);
          }
          if (m_dataSourceListeners.size() > 0) {
            notifyDataSetAvailable(new DataSetEvent(this, newTestSetInstances));
          }
	} catch (Exception ex) {
	  ex.printStackTrace();
	}
      }
    }
  }
  
  
  /**
   * Accept and process a batch clusterer event
   *
   * @param e a <code>BatchClassifierEvent</code> value
   */
  public void acceptClusterer(BatchClustererEvent e) {
    if (m_dataSourceListeners.size() > 0
        || m_trainingSetListeners.size() > 0
	|| m_testSetListeners.size() > 0) {

      if(e.getTestSet().isStructureOnly()) {
        return;
      }
      Instances testSet = e.getTestSet().getDataSet();

      weka.clusterers.Clusterer clusterer = e.getClusterer();
      String test;
      if(e.getTestOrTrain() == 0) {
        test = "test";
      } else {
        test = "training";
      }
      String relationNameModifier = "_"+test+"_"+e.getSetNumber()+"_of_"
	+e.getMaxSetNumber();
      if (!m_appendProbabilities || !(clusterer instanceof DensityBasedClusterer)) {
	if(m_appendProbabilities && !(clusterer instanceof DensityBasedClusterer)){
          System.err.println("Only density based clusterers can append probabilities. Instead cluster will be assigned for each instance.");
          if (m_logger != null) {
            m_logger.logMessage("[PredictionAppender] "
                + statusMessagePrefix() + " Only density based clusterers can "
                +"append probabilities. Instead cluster will be assigned for each "
                +"instance.");
            m_logger.statusMessage(statusMessagePrefix()
                +"WARNING: Only density based clusterers can append probabilities. "
                +"Instead cluster will be assigned for each instance.");
          }
        }
        try {
	  Instances newInstances = makeClusterDataSetClass(testSet, clusterer,
                                                           relationNameModifier);

          // data source listeners get both train and test sets
          if (m_dataSourceListeners.size() > 0) {
            notifyDataSetAvailable(new DataSetEvent(this, new Instances(newInstances,0)));
          }

          if (m_trainingSetListeners.size() > 0 && e.getTestOrTrain() > 0) {
             TrainingSetEvent tse = 
               new TrainingSetEvent(this, new Instances(newInstances, 0));
             tse.m_setNumber = e.getSetNumber();
             tse.m_maxSetNumber = e.getMaxSetNumber();
	    notifyTrainingSetAvailable(tse);
          }

          if (m_testSetListeners.size() > 0 && e.getTestOrTrain() == 0) {
             TestSetEvent tse = 
               new TestSetEvent(this, new Instances(newInstances, 0));
             tse.m_setNumber = e.getSetNumber();
             tse.m_maxSetNumber = e.getMaxSetNumber();
	    notifyTestSetAvailable(tse);
          }
          
	  // fill in predicted values
	  for (int i = 0; i < testSet.numInstances(); i++) {
	    double predCluster = 
	      clusterer.clusterInstance(testSet.instance(i));
	    newInstances.instance(i).setValue(newInstances.numAttributes()-1,
					      predCluster);
	  }
	  // notify listeners
          if (m_dataSourceListeners.size() > 0) {
            notifyDataSetAvailable(new DataSetEvent(this, newInstances));
          }
          if (m_trainingSetListeners.size() > 0 && e.getTestOrTrain() > 0) {
             TrainingSetEvent tse = 
               new TrainingSetEvent(this, newInstances);
             tse.m_setNumber = e.getSetNumber();
             tse.m_maxSetNumber = e.getMaxSetNumber();
	    notifyTrainingSetAvailable(tse);
          }
          if (m_testSetListeners.size() > 0 && e.getTestOrTrain() == 0) {
             TestSetEvent tse = 
               new TestSetEvent(this, newInstances);
             tse.m_setNumber = e.getSetNumber();
             tse.m_maxSetNumber = e.getMaxSetNumber();
	    notifyTestSetAvailable(tse);
          }

	  return;
	} catch (Exception ex) {
	  ex.printStackTrace();
	}
      }
      else{
	try {
	  Instances newInstances = 
	    makeClusterDataSetProbabilities(testSet,
                                            clusterer,relationNameModifier);
	  notifyDataSetAvailable(new DataSetEvent(this, new Instances(newInstances,0)));
          
	  // fill in predicted probabilities
	  for (int i = 0; i < testSet.numInstances(); i++) {
	    double [] probs = clusterer.
	      distributionForInstance(testSet.instance(i));
	    for (int j = 0; j < clusterer.numberOfClusters(); j++) {
	      newInstances.instance(i).setValue(testSet.numAttributes()+j,
						probs[j]);
	    }
	  }
	  // notify listeners
	  notifyDataSetAvailable(new DataSetEvent(this, newInstances));
	} catch (Exception ex) {
	  ex.printStackTrace();
	}
      }
    }
  }

  private Instances 
    makeDataSetProbabilities(Instances insts, Instances format,
			     weka.classifiers.Classifier classifier,
			     String relationNameModifier) 
  throws Exception {
    
    // adjust structure for InputMappedClassifier (if necessary)
    if (classifier instanceof weka.classifiers.misc.InputMappedClassifier) {
      format = 
        ((weka.classifiers.misc.InputMappedClassifier)classifier).
        getModelHeader(new Instances(format, 0));
    }
    
    String classifierName = classifier.getClass().getName();
    classifierName = classifierName.
      substring(classifierName.lastIndexOf('.')+1, classifierName.length());
    int numOrigAtts = insts.numAttributes();
    Instances newInstances = new Instances(insts);
    for (int i = 0; i < format.classAttribute().numValues(); i++) {
      weka.filters.unsupervised.attribute.Add addF = new
	weka.filters.unsupervised.attribute.Add();
      addF.setAttributeIndex("last");
      addF.setAttributeName(classifierName+"_prob_"+format.classAttribute().value(i));
      addF.setInputFormat(newInstances);
      newInstances = weka.filters.Filter.useFilter(newInstances, addF);
    }
    newInstances.setRelationName(insts.relationName()+relationNameModifier);
    return newInstances;
  }

  private Instances makeDataSetClass(Instances insts, Instances structure,
				     weka.classifiers.Classifier classifier,
				     String relationNameModifier) 
  throws Exception {
    
    // adjust structure for InputMappedClassifier (if necessary)
    if (classifier instanceof weka.classifiers.misc.InputMappedClassifier) {
      structure = 
        ((weka.classifiers.misc.InputMappedClassifier)classifier).
        getModelHeader(new Instances(structure, 0));
    }
    
    weka.filters.unsupervised.attribute.Add addF = new
      weka.filters.unsupervised.attribute.Add();
    addF.setAttributeIndex("last");
    String classifierName = classifier.getClass().getName();
    classifierName = classifierName.
      substring(classifierName.lastIndexOf('.')+1, classifierName.length());
    addF.setAttributeName("class_predicted_by: "+classifierName);
    if (structure.classAttribute().isNominal() || structure.classAttribute().isRanking()) {
      String classLabels = "";
      Enumeration enu = structure.classAttribute().enumerateValues();
      classLabels += (String)enu.nextElement();
      while (enu.hasMoreElements()) {
	classLabels += ","+(String)enu.nextElement();
      }
      addF.setNominalLabels(classLabels);
    }
    addF.setInputFormat(insts);


    Instances newInstances = 
      weka.filters.Filter.useFilter(insts, addF);
    newInstances.setRelationName(insts.relationName()+relationNameModifier);
    return newInstances;
  }
  
  private Instances 
    makeClusterDataSetProbabilities(Instances format,
			     weka.clusterers.Clusterer clusterer,
			     String relationNameModifier) 
  throws Exception {
    int numOrigAtts = format.numAttributes();
    Instances newInstances = new Instances(format);
    for (int i = 0; i < clusterer.numberOfClusters(); i++) {
      weka.filters.unsupervised.attribute.Add addF = new
	weka.filters.unsupervised.attribute.Add();
      addF.setAttributeIndex("last");
      addF.setAttributeName("prob_cluster"+i);
      addF.setInputFormat(newInstances);
      newInstances = weka.filters.Filter.useFilter(newInstances, addF);
    }
    newInstances.setRelationName(format.relationName()+relationNameModifier);
    return newInstances;
  }

  private Instances makeClusterDataSetClass(Instances format,
				     weka.clusterers.Clusterer clusterer,
				     String relationNameModifier) 
  throws Exception {
    
    weka.filters.unsupervised.attribute.Add addF = new
      weka.filters.unsupervised.attribute.Add();
    addF.setAttributeIndex("last");
    String clustererName = clusterer.getClass().getName();
    clustererName = clustererName.
      substring(clustererName.lastIndexOf('.')+1, clustererName.length());
    addF.setAttributeName("assigned_cluster: "+clustererName);
    //if (format.classAttribute().isNominal()) {
    String clusterLabels = "0";
      /*Enumeration enu = format.classAttribute().enumerateValues();
      clusterLabels += (String)enu.nextElement();
      while (enu.hasMoreElements()) {
	clusterLabels += ","+(String)enu.nextElement();
      }*/
    for(int i = 1; i <= clusterer.numberOfClusters()-1; i++)
        clusterLabels += ","+i;
    addF.setNominalLabels(clusterLabels);
    //}
    addF.setInputFormat(format);


    Instances newInstances = 
      weka.filters.Filter.useFilter(format, addF);
    newInstances.setRelationName(format.relationName()+relationNameModifier);
    return newInstances;
  }

  /**
   * Notify all instance listeners that an instance is available
   *
   * @param e an <code>InstanceEvent</code> value
   */
  protected void notifyInstanceAvailable(InstanceEvent e) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_instanceListeners.clone();
    }
    
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((InstanceListener)l.elementAt(i)).acceptInstance(e);
      }
    }
  }

  /**
   * Notify all Data source listeners that a data set is available
   *
   * @param e a <code>DataSetEvent</code> value
   */
  protected void notifyDataSetAvailable(DataSetEvent e) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_dataSourceListeners.clone();
    }
    
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((DataSourceListener)l.elementAt(i)).acceptDataSet(e);
      }
    }
  }
  
  /**
   * Notify all test set listeners that a test set is available
   *
   * @param e a <code>TestSetEvent</code> value
   */
  protected void notifyTestSetAvailable(TestSetEvent e) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_testSetListeners.clone();
    }
    
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((TestSetListener)l.elementAt(i)).acceptTestSet(e);
      }
    }
  }
  
  /**
   * Notify all test set listeners that a test set is available
   *
   * @param e a <code>TestSetEvent</code> value
   */
  protected void notifyTrainingSetAvailable(TrainingSetEvent e) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_trainingSetListeners.clone();
    }
    
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((TrainingSetListener)l.elementAt(i)).acceptTrainingSet(e);
      }
    }
  }

  /**
   * Set a logger
   *
   * @param logger a <code>weka.gui.Logger</code> value
   */
  public void setLog(weka.gui.Logger logger) {
    m_logger = logger;
  }

  public void stop() {
    // tell the listenee (upstream bean) to stop
    if (m_listenee instanceof BeanCommon) {
      ((BeanCommon)m_listenee).stop();
    }
  }
  
  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return false;
  }

  /**
   * Returns true if, at this time, 
   * the object will accept a connection according to the supplied
   * event name
   *
   * @param eventName the event
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {
    return (m_listenee == null);
  }

  /**
   * Returns true if, at this time, 
   * the object will accept a connection according to the supplied
   * EventSetDescriptor
   *
   * @param esd the EventSetDescriptor
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(EventSetDescriptor esd) {
    return connectionAllowed(esd.getName());
  }

  /**
   * Notify this object that it has been registered as a listener with
   * a source with respect to the supplied event name
   *
   * @param eventName
   * @param source the source with which this object has been registered as
   * a listener
   */
  public synchronized void connectionNotification(String eventName,
						  Object source) {
    if (connectionAllowed(eventName)) {
      m_listenee = source;
    }
  }

  /**
   * Notify this object that it has been deregistered as a listener with
   * a source with respect to the supplied event name
   *
   * @param eventName the event name
   * @param source the source with which this object has been registered as
   * a listener
   */
  public synchronized void disconnectionNotification(String eventName,
						     Object source) {
    if (m_listenee == source) {
      m_listenee = null;
      m_format = null; // assume any calculated instance format if now invalid
    }
  }

  /**
   * Returns true, if at the current time, the named event could
   * be generated. Assumes that supplied event names are names of
   * events that could be generated by this bean.
   *
   * @param eventName the name of the event in question
   * @return true if the named event could be generated at this point in
   * time
   */
  public boolean eventGeneratable(String eventName) {
    if (m_listenee == null) {
      return false;
    }

    if (m_listenee instanceof EventConstraints) {
      if (eventName.equals("instance")) {
	if (!((EventConstraints)m_listenee).
	    eventGeneratable("incrementalClassifier")) {
	  return false;
	}
      }
      if (eventName.equals("dataSet") 
	  || eventName.equals("trainingSet") 
	  || eventName.equals("testSet")) {
	if (((EventConstraints)m_listenee).
	    eventGeneratable("batchClassifier")) {
	  return true;
	}
	if (((EventConstraints)m_listenee).eventGeneratable("batchClusterer")) {
	  return true;
	}
	return false;
      }
    }
    return true;
  }
  
  private String statusMessagePrefix() {
    return getCustomName() + "$" + hashCode() + "|";
  }
}
