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
 *    DecisionTable.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.rules;

import weka.attributeSelection.ASSearch;
import weka.attributeSelection.BestFirst;
import weka.attributeSelection.SubsetEvaluator;
import weka.attributeSelection.ASEvaluation;
import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Evaluation;
import weka.classifiers.lazy.IBk;
import weka.core.AdditionalMeasureProducer;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class for building and using a simple decision table majority classifier.<br/>
 * <br/>
 * For more information see: <br/>
 * <br/>
 * Ron Kohavi: The Power of Decision Tables. In: 8th European Conference on Machine Learning, 174-189, 1995.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Kohavi1995,
 *    author = {Ron Kohavi},
 *    booktitle = {8th European Conference on Machine Learning},
 *    pages = {174-189},
 *    publisher = {Springer},
 *    title = {The Power of Decision Tables},
 *    year = {1995}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -S &lt;search method specification&gt;
 *  Full class name of search method, followed
 *  by its options.
 *  eg: "weka.attributeSelection.BestFirst -D 1"
 *  (default weka.attributeSelection.BestFirst)</pre>
 * 
 * <pre> -X &lt;number of folds&gt;
 *  Use cross validation to evaluate features.
 *  Use number of folds = 1 for leave one out CV.
 *  (Default = leave one out CV)</pre>
 * 
 * <pre> -E &lt;acc | rmse | mae | auc&gt;
 *  Performance evaluation measure to use for selecting attributes.
 *  (Default = accuracy for discrete class and rmse for numeric class)</pre>
 * 
 * <pre> -I
 *  Use nearest neighbour instead of global table majority.</pre>
 * 
 * <pre> -R
 *  Display decision table rules.
 * </pre>
 * 
 * <pre> 
 * Options specific to search method weka.attributeSelection.BestFirst:
 * </pre>
 * 
 * <pre> -P &lt;start set&gt;
 *  Specify a starting set of attributes.
 *  Eg. 1,3,5-7.</pre>
 * 
 * <pre> -D &lt;0 = backward | 1 = forward | 2 = bi-directional&gt;
 *  Direction of search. (default = 1).</pre>
 * 
 * <pre> -N &lt;num&gt;
 *  Number of non-improving nodes to
 *  consider before terminating search.</pre>
 * 
 * <pre> -S &lt;num&gt;
 *  Size of lookup cache for evaluated subsets.
 *  Expressed as a multiple of the number of
 *  attributes in the data set. (default = 1)</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5987 $ 
 */
public class DecisionTable 
  extends AbstractClassifier 
  implements OptionHandler, WeightedInstancesHandler, 
             AdditionalMeasureProducer, TechnicalInformationHandler {

  /** for serialization */
  static final long serialVersionUID = 2888557078165701326L;

  /** The hashtable used to hold training instances */
  protected Hashtable m_entries;

  /** The class priors to use when there is no match in the table */
  protected double [] m_classPriorCounts;
  protected double [] m_classPriors;

  /** Holds the final feature set */
  protected int [] m_decisionFeatures;

  /** Discretization filter */
  protected Filter m_disTransform;

  /** Filter used to remove columns discarded by feature selection */
  protected Remove m_delTransform;

  /** IB1 used to classify non matching instances rather than majority class */
  protected IBk m_ibk;

  /** Holds the original training instances */
  protected Instances m_theInstances;

  /** Holds the final feature selected set of instances */
  protected Instances m_dtInstances;

  /** The number of attributes in the dataset */
  protected int m_numAttributes;

  /** The number of instances in the dataset */
  private int m_numInstances;

  /** Class is nominal */
  protected boolean m_classIsNominal;

  /** Use the IBk classifier rather than majority class */
  protected boolean m_useIBk;

  /** Display Rules */
  protected boolean m_displayRules;

  /** Number of folds for cross validating feature sets */
  private int m_CVFolds;

  /** Random numbers for use in cross validation */
  private Random m_rr;

  /** Holds the majority class */
  protected double m_majority;

  /** The search method to use */
  protected ASSearch m_search = new BestFirst();

  /** Our own internal evaluator */
  protected ASEvaluation m_evaluator;

  /** The evaluation object used to evaluate subsets */
  protected Evaluation m_evaluation;

  /** default is accuracy for discrete class and RMSE for numeric class */
  public static final int EVAL_DEFAULT = 1;
  public static final int EVAL_ACCURACY = 2;
  public static final int EVAL_RMSE = 3;
  public static final int EVAL_MAE = 4;
  public static final int EVAL_AUC = 5;

  public static final Tag [] TAGS_EVALUATION = {
    new Tag(EVAL_DEFAULT, "Default: accuracy (discrete class); RMSE (numeric class)"),
    new Tag(EVAL_ACCURACY, "Accuracy (discrete class only"),
    new Tag(EVAL_RMSE, "RMSE (of the class probabilities for discrete class)"),
    new Tag(EVAL_MAE, "MAE (of the class probabilities for discrete class)"),
    new Tag(EVAL_AUC, "AUC (area under the ROC curve - discrete class only)")
  };

  protected int m_evaluationMeasure = EVAL_DEFAULT;

  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return  
    "Class for building and using a simple decision table majority "
    + "classifier.\n\n"
    + "For more information see: \n\n"
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
    TechnicalInformation 	result;

    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Ron Kohavi");
    result.setValue(Field.TITLE, "The Power of Decision Tables");
    result.setValue(Field.BOOKTITLE, "8th European Conference on Machine Learning");
    result.setValue(Field.YEAR, "1995");
    result.setValue(Field.PAGES, "174-189");
    result.setValue(Field.PUBLISHER, "Springer");

    return result;
  }
  
  /**
   * Inserts an instance into the hash table
   *
   * @param inst instance to be inserted
   * @param instA to create the hash key from
   * @throws Exception if the instance can't be inserted
   */
  private void insertIntoTable(Instance inst, double [] instA)
  throws Exception {

    double [] tempClassDist2;
    double [] newDist;
    DecisionTableHashKey thekey;

    if (instA != null) {
      thekey = new DecisionTableHashKey(instA);
    } else {
      thekey = new DecisionTableHashKey(inst, inst.numAttributes(), false);
    }

    // see if this one is already in the table
    tempClassDist2 = (double []) m_entries.get(thekey);
    if (tempClassDist2 == null) {
      if (m_classIsNominal) {
	newDist = new double [m_theInstances.classAttribute().numValues()];
	
	//Leplace estimation
	for (int i = 0; i < m_theInstances.classAttribute().numValues(); i++) {
	  newDist[i] = 1.0;
	}
	
	newDist[(int)inst.classValue()] = inst.weight();

	// add to the table
	m_entries.put(thekey, newDist);
      } else {
	newDist = new double [2];
	newDist[0] = inst.classValue() * inst.weight();
	newDist[1] = inst.weight();

	// add to the table
	m_entries.put(thekey, newDist);
      }
    } else { 

      // update the distribution for this instance
      if (m_classIsNominal) {
	tempClassDist2[(int)inst.classValue()]+=inst.weight();

	// update the table
	m_entries.put(thekey, tempClassDist2);
      } else  {
	tempClassDist2[0] += (inst.classValue() * inst.weight());
	tempClassDist2[1] += inst.weight();

	// update the table
	m_entries.put(thekey, tempClassDist2);
      }
    }
  }

  /**
   * Classifies an instance for internal leave one out cross validation
   * of feature sets
   *
   * @param instance instance to be "left out" and classified
   * @param instA feature values of the selected features for the instance
   * @return the classification of the instance
   * @throws Exception if something goes wrong
   */
  double evaluateInstanceLeaveOneOut(Instance instance, double [] instA)
  throws Exception {

    DecisionTableHashKey thekey;
    double [] tempDist;
    double [] normDist;

    thekey = new DecisionTableHashKey(instA);
    if (m_classIsNominal) {

      // if this one is not in the table
      if ((tempDist = (double [])m_entries.get(thekey)) == null) {
	throw new Error("This should never happen!");
      } else {
	normDist = new double [tempDist.length];
	System.arraycopy(tempDist,0,normDist,0,tempDist.length);
	normDist[(int)instance.classValue()] -= instance.weight();

	// update the table
	// first check to see if the class counts are all zero now
	boolean ok = false;
	for (int i=0;i<normDist.length;i++) {
	  if (Utils.gr(normDist[i],1.0)) {
	    ok = true;
	    break;
	  }
	}

//	downdate the class prior counts
	m_classPriorCounts[(int)instance.classValue()] -= 
	  instance.weight();
	double [] classPriors = m_classPriorCounts.clone();
	Utils.normalize(classPriors);
	if (!ok) { // majority class
	  normDist = classPriors;
	}

	m_classPriorCounts[(int)instance.classValue()] += 
	  instance.weight();

	//if (ok) {
	Utils.normalize(normDist);
	if (m_evaluationMeasure == EVAL_AUC) {
	  m_evaluation.evaluateModelOnceAndRecordPrediction(normDist, instance);						
	} else {
	  m_evaluation.evaluateModelOnce(normDist, instance);
	}
	return Utils.maxIndex(normDist);
	/*} else {
	  normDist = new double [normDist.length];
	  normDist[(int)m_majority] = 1.0;
	  if (m_evaluationMeasure == EVAL_AUC) {
	    m_evaluation.evaluateModelOnceAndRecordPrediction(normDist, instance);						
	  } else {
	    m_evaluation.evaluateModelOnce(normDist, instance);
	  }
	  return m_majority;
	} */
      }
      //      return Utils.maxIndex(tempDist);
    } else {

      // see if this one is already in the table
      if ((tempDist = (double[])m_entries.get(thekey)) != null) {
	normDist = new double [tempDist.length];
	System.arraycopy(tempDist,0,normDist,0,tempDist.length);
	normDist[0] -= (instance.classValue() * instance.weight());
	normDist[1] -= instance.weight();
	if (Utils.eq(normDist[1],0.0)) {
	  double [] temp = new double[1];
	  temp[0] = m_majority;
	  m_evaluation.evaluateModelOnce(temp, instance);
	  return m_majority;
	} else {
	  double [] temp = new double[1];
	  temp[0] = normDist[0] / normDist[1];
	  m_evaluation.evaluateModelOnce(temp, instance);
	  return temp[0];
	}
      } else {
	throw new Error("This should never happen!");
      }
    }

    // shouldn't get here 
    // return 0.0;
  }

  /**
   * Calculates the accuracy on a test fold for internal cross validation
   * of feature sets
   *
   * @param fold set of instances to be "left out" and classified
   * @param fs currently selected feature set
   * @return the accuracy for the fold
   * @throws Exception if something goes wrong
   */
  double evaluateFoldCV(Instances fold, int [] fs) throws Exception {

    int i;
    int ruleCount = 0;
    int numFold = fold.numInstances();
    int numCl = m_theInstances.classAttribute().numValues();
    double [][] class_distribs = new double [numFold][numCl];
    double [] instA = new double [fs.length];
    double [] normDist;
    DecisionTableHashKey thekey;
    double acc = 0.0;
    int classI = m_theInstances.classIndex();
    Instance inst;

    if (m_classIsNominal) {
      normDist = new double [numCl];
    } else {
      normDist = new double [2];
    }

    // first *remove* instances
    for (i=0;i<numFold;i++) {
      inst = fold.instance(i);
      for (int j=0;j<fs.length;j++) {
	if (fs[j] == classI) {
	  instA[j] = Double.MAX_VALUE; // missing for the class
	} else if (inst.isMissing(fs[j])) {
	  instA[j] = Double.MAX_VALUE;
	} else{
	  instA[j] = inst.value(fs[j]);
	}
      }
      thekey = new DecisionTableHashKey(instA);
      if ((class_distribs[i] = (double [])m_entries.get(thekey)) == null) {
	throw new Error("This should never happen!");
      } else {
	if (m_classIsNominal) {
	  class_distribs[i][(int)inst.classValue()] -= inst.weight();
	} else {
	  class_distribs[i][0] -= (inst.classValue() * inst.weight());
	  class_distribs[i][1] -= inst.weight();
	}
	ruleCount++;
      }
      m_classPriorCounts[(int)inst.classValue()] -= 
	inst.weight();	
    }
    double [] classPriors = m_classPriorCounts.clone();
    Utils.normalize(classPriors);

    // now classify instances
    for (i=0;i<numFold;i++) {
      inst = fold.instance(i);
      System.arraycopy(class_distribs[i],0,normDist,0,normDist.length);
      if (m_classIsNominal) {
	boolean ok = false;
	for (int j=0;j<normDist.length;j++) {
	  if (Utils.gr(normDist[j],1.0)) {
	    ok = true;
	    break;
	  }
	}

	if (!ok) { // majority class
	  normDist = classPriors.clone();
	}

//	if (ok) {
	Utils.normalize(normDist);
	if (m_evaluationMeasure == EVAL_AUC) {
	  m_evaluation.evaluateModelOnceAndRecordPrediction(normDist, inst);						
	} else {
	  m_evaluation.evaluateModelOnce(normDist, inst);
	}
	/*	} else {					
	  normDist[(int)m_majority] = 1.0;
	  if (m_evaluationMeasure == EVAL_AUC) {
	    m_evaluation.evaluateModelOnceAndRecordPrediction(normDist, inst);						
	  } else {
	    m_evaluation.evaluateModelOnce(normDist, inst);					
	  }
	} */
      } else {
	if (Utils.eq(normDist[1],0.0)) {
	  double [] temp = new double[1];
	  temp[0] = m_majority;
	  m_evaluation.evaluateModelOnce(temp, inst);
	} else {
	  double [] temp = new double[1];
	  temp[0] = normDist[0] / normDist[1];
	  m_evaluation.evaluateModelOnce(temp, inst);
	}
      }
    }

    // now re-insert instances
    for (i=0;i<numFold;i++) {
      inst = fold.instance(i);

      m_classPriorCounts[(int)inst.classValue()] += 
	inst.weight();

      if (m_classIsNominal) {
	class_distribs[i][(int)inst.classValue()] += inst.weight();
      } else {
	class_distribs[i][0] += (inst.classValue() * inst.weight());
	class_distribs[i][1] += inst.weight();
      }
    }
    return acc;
  }


  /**
   * Evaluates a feature subset by cross validation
   *
   * @param feature_set the subset to be evaluated
   * @param num_atts the number of attributes in the subset
   * @return the estimated accuracy
   * @throws Exception if subset can't be evaluated
   */
  protected double estimatePerformance(BitSet feature_set, int num_atts)
  throws Exception {

    m_evaluation = new Evaluation(m_theInstances);
    int i;
    int [] fs = new int [num_atts];

    double [] instA = new double [num_atts];
    int classI = m_theInstances.classIndex();

    int index = 0;
    for (i=0;i<m_numAttributes;i++) {
      if (feature_set.get(i)) {
	fs[index++] = i;
      }
    }

    // create new hash table
    m_entries = new Hashtable((int)(m_theInstances.numInstances() * 1.5));

    // insert instances into the hash table
    for (i=0;i<m_numInstances;i++) {

      Instance inst = m_theInstances.instance(i);
      for (int j=0;j<fs.length;j++) {
	if (fs[j] == classI) {
	  instA[j] = Double.MAX_VALUE; // missing for the class
	} else if (inst.isMissing(fs[j])) {
	  instA[j] = Double.MAX_VALUE;
	} else {
	  instA[j] = inst.value(fs[j]);
	}
      }
      insertIntoTable(inst, instA);
    }


    if (m_CVFolds == 1) {

      // calculate leave one out error
      for (i=0;i<m_numInstances;i++) {
	Instance inst = m_theInstances.instance(i);
	for (int j=0;j<fs.length;j++) {
	  if (fs[j] == classI) {
	    instA[j] = Double.MAX_VALUE; // missing for the class
	  } else if (inst.isMissing(fs[j])) {
	    instA[j] = Double.MAX_VALUE;
	  } else {
	    instA[j] = inst.value(fs[j]);
	  }
	}
	evaluateInstanceLeaveOneOut(inst, instA);				
      }
    } else {
      m_theInstances.randomize(m_rr);
      m_theInstances.stratify(m_CVFolds);

      // calculate 10 fold cross validation error
      for (i=0;i<m_CVFolds;i++) {
	Instances insts = m_theInstances.testCV(m_CVFolds,i);
	evaluateFoldCV(insts, fs);
      }
    }

    switch (m_evaluationMeasure) {
    case EVAL_DEFAULT:
      if (m_classIsNominal) {
	return m_evaluation.pctCorrect();
      }
      return -m_evaluation.rootMeanSquaredError();
    case EVAL_ACCURACY:
      return m_evaluation.pctCorrect();
    case EVAL_RMSE:
      return -m_evaluation.rootMeanSquaredError();
    case EVAL_MAE:
      return -m_evaluation.meanAbsoluteError();
    case EVAL_AUC:
      double [] classPriors = m_evaluation.getClassPriors();
      Utils.normalize(classPriors);
      double weightedAUC = 0;
      for (i = 0; i < m_theInstances.classAttribute().numValues(); i++) {
	double tempAUC = m_evaluation.areaUnderROC(i);
	if (!Utils.isMissingValue(tempAUC)) {
	  weightedAUC += (classPriors[i] * tempAUC);
	} else {
	  System.err.println("Undefined AUC!!");
	}
      }
      return weightedAUC;
    }
    // shouldn't get here
    return 0.0;
  }

  /**
   * Returns a String representation of a feature subset
   *
   * @param sub BitSet representation of a subset
   * @return String containing subset
   */
  private String printSub(BitSet sub) {

    String s="";
    for (int jj=0;jj<m_numAttributes;jj++) {
      if (sub.get(jj)) {
	s += " "+(jj+1);
      }
    }
    return s;
  }

  /**
   * Resets the options.
   */
  protected void resetOptions()  {

    m_entries = null;
    m_decisionFeatures = null;
    m_useIBk = false;
    m_CVFolds = 1;
    m_displayRules = false;
    m_evaluationMeasure = EVAL_DEFAULT;
  }

  /**
   * Constructor for a DecisionTable
   */
  public DecisionTable() {

    resetOptions();
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(7);

    newVector.addElement(new Option(
	"\tFull class name of search method, followed\n"
	+ "\tby its options.\n"
	+ "\teg: \"weka.attributeSelection.BestFirst -D 1\"\n"
	+ "\t(default weka.attributeSelection.BestFirst)",
	"S", 1, "-S <search method specification>"));

    newVector.addElement(new Option(
	"\tUse cross validation to evaluate features.\n" +
	"\tUse number of folds = 1 for leave one out CV.\n" +
	"\t(Default = leave one out CV)",
	"X", 1, "-X <number of folds>"));

    newVector.addElement(new Option(
	"\tPerformance evaluation measure to use for selecting attributes.\n" +
	"\t(Default = accuracy for discrete class and rmse for numeric class)",
	"E", 1, "-E <acc | rmse | mae | auc>"));

    newVector.addElement(new Option(
	"\tUse nearest neighbour instead of global table majority.",
	"I", 0, "-I"));

    newVector.addElement(new Option(
	"\tDisplay decision table rules.\n",
	"R", 0, "-R")); 

    newVector.addElement(new Option(
	"",
	"", 0, "\nOptions specific to search method "
	+ m_search.getClass().getName() + ":"));
    Enumeration enu = ((OptionHandler)m_search).listOptions();
    while (enu.hasMoreElements()) {
      newVector.addElement(enu.nextElement());
    }
    return newVector.elements();
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String crossValTipText() {
    return "Sets the number of folds for cross validation (1 = leave one out).";
  }

  /**
   * Sets the number of folds for cross validation (1 = leave one out)
   *
   * @param folds the number of folds
   */
  public void setCrossVal(int folds) {

    m_CVFolds = folds;
  }

  /**
   * Gets the number of folds for cross validation
   *
   * @return the number of cross validation folds
   */
  public int getCrossVal() {

    return m_CVFolds;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String useIBkTipText() {
    return "Sets whether IBk should be used instead of the majority class.";
  }

  /**
   * Sets whether IBk should be used instead of the majority class
   *
   * @param ibk true if IBk is to be used
   */
  public void setUseIBk(boolean ibk) {

    m_useIBk = ibk;
  }

  /**
   * Gets whether IBk is being used instead of the majority class
   *
   * @return true if IBk is being used
   */
  public boolean getUseIBk() {

    return m_useIBk;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String displayRulesTipText() {
    return "Sets whether rules are to be printed.";
  }

  /**
   * Sets whether rules are to be printed
   *
   * @param rules true if rules are to be printed
   */
  public void setDisplayRules(boolean rules) {

    m_displayRules = rules;
  }

  /**
   * Gets whether rules are being printed
   *
   * @return true if rules are being printed
   */
  public boolean getDisplayRules() {

    return m_displayRules;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String searchTipText() {
    return "The search method used to find good attribute combinations for the "
    + "decision table.";
  }
  /**
   * Sets the search method to use
   * 
   * @param search
   */
  public void setSearch(ASSearch search) {
    m_search = search;
  }

  /**
   * Gets the current search method
   * 
   * @return the search method used
   */
  public ASSearch getSearch() {
    return m_search;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String evaluationMeasureTipText() {
    return "The measure used to evaluate the performance of attribute combinations "
    + "used in the decision table.";
  }
  /**
   * Gets the currently set performance evaluation measure used for selecting
   * attributes for the decision table
   * 
   * @return the performance evaluation measure
   */
  public SelectedTag getEvaluationMeasure() {
    return new SelectedTag(m_evaluationMeasure, TAGS_EVALUATION);
  }

  /**
   * Sets the performance evaluation measure to use for selecting attributes
   * for the decision table
   * 
   * @param newMethod the new performance evaluation metric to use
   */
  public void setEvaluationMeasure(SelectedTag newMethod) {
    if (newMethod.getTags() == TAGS_EVALUATION) {
      m_evaluationMeasure = newMethod.getSelectedTag().getID();
    }
  }

  /**
   * Parses the options for this object. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -S &lt;search method specification&gt;
   *  Full class name of search method, followed
   *  by its options.
   *  eg: "weka.attributeSelection.BestFirst -D 1"
   *  (default weka.attributeSelection.BestFirst)</pre>
   * 
   * <pre> -X &lt;number of folds&gt;
   *  Use cross validation to evaluate features.
   *  Use number of folds = 1 for leave one out CV.
   *  (Default = leave one out CV)</pre>
   * 
   * <pre> -E &lt;acc | rmse | mae | auc&gt;
   *  Performance evaluation measure to use for selecting attributes.
   *  (Default = accuracy for discrete class and rmse for numeric class)</pre>
   * 
   * <pre> -I
   *  Use nearest neighbour instead of global table majority.</pre>
   * 
   * <pre> -R
   *  Display decision table rules.
   * </pre>
   * 
   * <pre> 
   * Options specific to search method weka.attributeSelection.BestFirst:
   * </pre>
   * 
   * <pre> -P &lt;start set&gt;
   *  Specify a starting set of attributes.
   *  Eg. 1,3,5-7.</pre>
   * 
   * <pre> -D &lt;0 = backward | 1 = forward | 2 = bi-directional&gt;
   *  Direction of search. (default = 1).</pre>
   * 
   * <pre> -N &lt;num&gt;
   *  Number of non-improving nodes to
   *  consider before terminating search.</pre>
   * 
   * <pre> -S &lt;num&gt;
   *  Size of lookup cache for evaluated subsets.
   *  Expressed as a multiple of the number of
   *  attributes in the data set. (default = 1)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    String optionString;

    resetOptions();

    optionString = Utils.getOption('X',options);
    if (optionString.length() != 0) {
      m_CVFolds = Integer.parseInt(optionString);
    }

    m_useIBk = Utils.getFlag('I',options);

    m_displayRules = Utils.getFlag('R',options);

    optionString = Utils.getOption('E', options);
    if (optionString.length() != 0) {
      if (optionString.equals("acc")) {
	setEvaluationMeasure(new SelectedTag(EVAL_ACCURACY, TAGS_EVALUATION));
      } else if (optionString.equals("rmse")) {
	setEvaluationMeasure(new SelectedTag(EVAL_RMSE, TAGS_EVALUATION));
      } else if (optionString.equals("mae")) {
	setEvaluationMeasure(new SelectedTag(EVAL_MAE, TAGS_EVALUATION));
      } else if (optionString.equals("auc")) {
	setEvaluationMeasure(new SelectedTag(EVAL_AUC, TAGS_EVALUATION));
      } else {
	throw new IllegalArgumentException("Invalid evaluation measure");
      }
    }

    String searchString = Utils.getOption('S', options);
    if (searchString.length() == 0)
      searchString = weka.attributeSelection.BestFirst.class.getName();
    String [] searchSpec = Utils.splitOptions(searchString);
    if (searchSpec.length == 0) {
      throw new IllegalArgumentException("Invalid search specification string");
    }
    String searchName = searchSpec[0];
    searchSpec[0] = "";
    setSearch(ASSearch.forName(searchName, searchSpec));
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [9];
    int current = 0;

    options[current++] = "-X"; options[current++] = "" + m_CVFolds;

    if (m_evaluationMeasure != EVAL_DEFAULT) {
      options[current++] = "-E";
      switch (m_evaluationMeasure) {
      case EVAL_ACCURACY:
	options[current++] = "acc";
	break;
      case EVAL_RMSE:
	options[current++] = "rmse";
	break;
      case EVAL_MAE:
	options[current++] = "mae";
	break;
      case EVAL_AUC:
	options[current++] = "auc";
	break;
      }
    }
    if (m_useIBk) {
      options[current++] = "-I";
    }
    if (m_displayRules) {
      options[current++] = "-R";
    }

    options[current++] = "-S";
    options[current++] = "" + getSearchSpec();

    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Gets the search specification string, which contains the class name of
   * the search method and any options to it
   *
   * @return the search string.
   */
  protected String getSearchSpec() {

    ASSearch s = getSearch();
    if (s instanceof OptionHandler) {
      return s.getClass().getName() + " "
      + Utils.joinOptions(((OptionHandler)s).getOptions());
    }
    return s.getClass().getName();
  }

  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.DATE_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enable(Capability.NOMINAL_CLASS);
    if (m_evaluationMeasure != EVAL_ACCURACY && m_evaluationMeasure != EVAL_AUC) {
      result.enable(Capability.NUMERIC_CLASS);
      result.enable(Capability.DATE_CLASS);
    }
    
    result.enable(Capability.MISSING_CLASS_VALUES);

    return result;
  }
  
  private class DummySubsetEvaluator extends ASEvaluation implements SubsetEvaluator {
    /** for serialization */
    private static final long serialVersionUID = 3927442457704974150L;
      
    public void buildEvaluator(Instances data) throws Exception {
    }

    public double evaluateSubset(BitSet subset) throws Exception {

      int fc = 0;
      for (int jj = 0;jj < m_numAttributes; jj++) {
        if (subset.get(jj)) {
          fc++;
        }
      }

      return estimatePerformance(subset, fc);
    }
  }

  /**
   * Sets up a dummy subset evaluator that basically just delegates
   * evaluation to the estimatePerformance method in DecisionTable
   */
  protected void setUpEvaluator() throws Exception {
    m_evaluator = new DummySubsetEvaluator();
  }

  protected boolean m_saveMemory = true;
  /**
   * Generates the classifier.
   *
   * @param data set of instances serving as training data 
   * @throws Exception if the classifier has not been generated successfully
   */
  public void buildClassifier(Instances data) throws Exception {

    // can classifier handle the data?
    getCapabilities().testWithFail(data);

    // remove instances with missing class
    m_theInstances = new Instances(data);
    m_theInstances.deleteWithMissingClass();

    m_rr = new Random(1);

    if (m_theInstances.classAttribute().isNominal() || m_theInstances.classAttribute().isRanking())  {//	 Set up class priors
      m_classPriorCounts = 
	new double [data.classAttribute().numValues()];
      Arrays.fill(m_classPriorCounts, 1.0);
      for (int i = 0; i <data.numInstances(); i++) {
	Instance curr = data.instance(i);
	m_classPriorCounts[(int)curr.classValue()] += 
	  curr.weight();
      }
      m_classPriors = m_classPriorCounts.clone();
      Utils.normalize(m_classPriors);
    }

    setUpEvaluator();

    if (m_theInstances.classAttribute().isNumeric()) {
      m_disTransform = new weka.filters.unsupervised.attribute.Discretize();
      m_classIsNominal = false;

      // use binned discretisation if the class is numeric
      ((weka.filters.unsupervised.attribute.Discretize)m_disTransform).
      setBins(10);
      ((weka.filters.unsupervised.attribute.Discretize)m_disTransform).
      setInvertSelection(true);

      // Discretize all attributes EXCEPT the class 
      String rangeList = "";
      rangeList+=(m_theInstances.classIndex()+1);
      //System.out.println("The class col: "+m_theInstances.classIndex());

      ((weka.filters.unsupervised.attribute.Discretize)m_disTransform).
      setAttributeIndices(rangeList);
    } else {
      m_disTransform = new weka.filters.supervised.attribute.Discretize();
      ((weka.filters.supervised.attribute.Discretize)m_disTransform).setUseBetterEncoding(true);
      m_classIsNominal = true;
    }

    m_disTransform.setInputFormat(m_theInstances);
    m_theInstances = Filter.useFilter(m_theInstances, m_disTransform);

    m_numAttributes = m_theInstances.numAttributes();
    m_numInstances = m_theInstances.numInstances();
    m_majority = m_theInstances.meanOrMode(m_theInstances.classAttribute());

    // Perform the search
    int [] selected = m_search.search(m_evaluator, m_theInstances);

    m_decisionFeatures = new int [selected.length+1];
    System.arraycopy(selected, 0, m_decisionFeatures, 0, selected.length);
    m_decisionFeatures[m_decisionFeatures.length-1] = m_theInstances.classIndex();

    // reduce instances to selected features
    m_delTransform = new Remove();
    m_delTransform.setInvertSelection(true);

    // set features to keep
    m_delTransform.setAttributeIndicesArray(m_decisionFeatures); 
    m_delTransform.setInputFormat(m_theInstances);
    m_dtInstances = Filter.useFilter(m_theInstances, m_delTransform);

    // reset the number of attributes
    m_numAttributes = m_dtInstances.numAttributes();

    // create hash table
    m_entries = new Hashtable((int)(m_dtInstances.numInstances() * 1.5));

    // insert instances into the hash table
    for (int i = 0; i < m_numInstances; i++) {
      Instance inst = m_dtInstances.instance(i);
      insertIntoTable(inst, null);
    }

    // Replace the global table majority with nearest neighbour?
    if (m_useIBk) {
      m_ibk = new IBk();
      m_ibk.buildClassifier(m_theInstances);
    }

    // Save memory
    if (m_saveMemory) {
      m_theInstances = new Instances(m_theInstances, 0);
      m_dtInstances = new Instances(m_dtInstances, 0);
    }
    m_evaluation = null;
  }

  /**
   * Calculates the class membership probabilities for the given 
   * test instance.
   *
   * @param instance the instance to be classified
   * @return predicted class probability distribution
   * @throws Exception if distribution can't be computed
   */
  public double [] distributionForInstance(Instance instance)
  throws Exception {

    DecisionTableHashKey thekey;
    double [] tempDist;
    double [] normDist;

    m_disTransform.input(instance);
    m_disTransform.batchFinished();
    instance = m_disTransform.output();

    m_delTransform.input(instance);
    m_delTransform.batchFinished();
    instance = m_delTransform.output();

    thekey = new DecisionTableHashKey(instance, instance.numAttributes(), false);

    // if this one is not in the table
    if ((tempDist = (double [])m_entries.get(thekey)) == null) {
      if (m_useIBk) {
	tempDist = m_ibk.distributionForInstance(instance);
      } else {
	if (!m_classIsNominal) {
	  tempDist = new double[1];
	  tempDist[0] = m_majority;
	} else {
	  tempDist = m_classPriors.clone();
	  /*tempDist = new double [m_theInstances.classAttribute().numValues()];
	  tempDist[(int)m_majority] = 1.0; */
	}
      }
    } else {
      if (!m_classIsNominal) {
	normDist = new double[1];
	normDist[0] = (tempDist[0] / tempDist[1]);
	tempDist = normDist;
      } else {

	// normalise distribution
	normDist = new double [tempDist.length];
	System.arraycopy(tempDist,0,normDist,0,tempDist.length);
	Utils.normalize(normDist);
	tempDist = normDist;
      }
    }
    return tempDist;
  }

  /**
   * Returns a string description of the features selected
   *
   * @return a string of features
   */
  public String printFeatures() {

    int i;
    String s = "";

    for (i=0;i<m_decisionFeatures.length;i++) {
      if (i==0) {
	s = ""+(m_decisionFeatures[i]+1);
      } else {
	s += ","+(m_decisionFeatures[i]+1);
      }
    }
    return s;
  }

  /**
   * Returns the number of rules
   * @return the number of rules
   */
  public double measureNumRules() {
    return m_entries.size();
  }

  /**
   * Returns an enumeration of the additional measure names
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector newVector = new Vector(1);
    newVector.addElement("measureNumRules");
    return newVector.elements();
  }

  /**
   * Returns the value of the named measure
   * @param additionalMeasureName the name of the measure to query for its value
   * @return the value of the named measure
   * @throws IllegalArgumentException if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) {
    if (additionalMeasureName.compareToIgnoreCase("measureNumRules") == 0) {
      return measureNumRules();
    } else {
      throw new IllegalArgumentException(additionalMeasureName 
	  + " not supported (DecisionTable)");
    }
  }

  /**
   * Returns a description of the classifier.
   *
   * @return a description of the classifier as a string.
   */
  public String toString() {

    if (m_entries == null) {
      return "Decision Table: No model built yet.";
    } else {
      StringBuffer text = new StringBuffer();

      text.append("Decision Table:"+
	  "\n\nNumber of training instances: "+m_numInstances+
	  "\nNumber of Rules : "+m_entries.size()+"\n");

      if (m_useIBk) {
	text.append("Non matches covered by IB1.\n");
      } else {
	text.append("Non matches covered by Majority class.\n");
      }

      text.append(m_search.toString());
      /*text.append("Best first search for feature set,\nterminated after "+
					m_maxStale+" non improving subsets.\n"); */

      text.append("Evaluation (for feature selection): CV ");
      if (m_CVFolds > 1) {
	text.append("("+m_CVFolds+" fold) ");
      } else {
	text.append("(leave one out) ");
      }
      text.append("\nFeature set: "+printFeatures());

      if (m_displayRules) {

	// find out the max column width
	int maxColWidth = 0;
	for (int i=0;i<m_dtInstances.numAttributes();i++) {
	  if (m_dtInstances.attribute(i).name().length() > maxColWidth) {
	    maxColWidth = m_dtInstances.attribute(i).name().length();
	  }

	  if (m_classIsNominal || (i != m_dtInstances.classIndex())) {
	    Enumeration e = m_dtInstances.attribute(i).enumerateValues();
	    while (e.hasMoreElements()) {
	      String ss = (String)e.nextElement();
	      if (ss.length() > maxColWidth) {
		maxColWidth = ss.length();
	      }
	    }
	  }
	}

	text.append("\n\nRules:\n");
	StringBuffer tm = new StringBuffer();
	for (int i=0;i<m_dtInstances.numAttributes();i++) {
	  if (m_dtInstances.classIndex() != i) {
	    int d = maxColWidth - m_dtInstances.attribute(i).name().length();
	    tm.append(m_dtInstances.attribute(i).name());
	    for (int j=0;j<d+1;j++) {
	      tm.append(" ");
	    }
	  }
	}
	tm.append(m_dtInstances.attribute(m_dtInstances.classIndex()).name()+"  ");

	for (int i=0;i<tm.length()+10;i++) {
	  text.append("=");
	}
	text.append("\n");
	text.append(tm);
	text.append("\n");
	for (int i=0;i<tm.length()+10;i++) {
	  text.append("=");
	}
	text.append("\n");

	Enumeration e = m_entries.keys();
	while (e.hasMoreElements()) {
	  DecisionTableHashKey tt = (DecisionTableHashKey)e.nextElement();
	  text.append(tt.toString(m_dtInstances,maxColWidth));
	  double [] ClassDist = (double []) m_entries.get(tt);

	  if (m_classIsNominal) {
	    int m = Utils.maxIndex(ClassDist);
	    try {
	      text.append(m_dtInstances.classAttribute().value(m)+"\n");
	    } catch (Exception ee) {
	      System.out.println(ee.getMessage());
	    }
	  } else {
	    text.append((ClassDist[0] / ClassDist[1])+"\n");
	  }
	}

	for (int i=0;i<tm.length()+10;i++) {
	  text.append("=");
	}
	text.append("\n");
	text.append("\n");
      }
      return text.toString();
    }
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
   * Main method for testing this class.
   *
   * @param argv the command-line options
   */
  public static void main(String [] argv) {
    runClassifier(new DecisionTable(), argv);
  }
}
