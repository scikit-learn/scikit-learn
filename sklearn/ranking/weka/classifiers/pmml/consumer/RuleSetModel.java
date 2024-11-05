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
 *    RuleSetModel.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.pmml.consumer;

import java.io.Serializable;
import java.util.ArrayList;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.classifiers.pmml.consumer.TreeModel.MiningFunction;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.pmml.MiningSchema;

/**
 * Class implementing import of PMML RuleSetModel. Can be used as a Weka
 * classifier for prediction only (buildClassifier() raises an Exception).
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 5987 $
 */
public class RuleSetModel extends PMMLClassifier {
  
  /** For serialization */
  private static final long serialVersionUID = 1993161168811020547L;

  /**
   * Abstract inner base class for Rules
   */
  static abstract class Rule implements Serializable {
    
    /** For serialization */
    private static final long serialVersionUID = 6236231263477446102L;
    
    /** The predicate for this rule */
    protected TreeModel.Predicate m_predicate;
    
    public Rule(Element ruleE, MiningSchema miningSchema) throws Exception {
      // Set up the predicate
      m_predicate = TreeModel.Predicate.getPredicate(ruleE, miningSchema);
    }
    
    /**
     * Collect the rule(s) that fire for the supplied incoming instance
     * 
     * @param input a vector of independent and derived independent variables
     * @param ruleCollection the array list to add any firing rules into
     */
    public abstract void fires(double[] input, ArrayList<SimpleRule> ruleCollection);
    
    /**
     * Get a textual description of this Rule
     * 
     * @param prefix prefix string (typically some number of spaces) to prepend
     * @param indent the number of additional spaces to add to the prefix
     * @return a description of this Rule as a String
     */
    public abstract String toString(String prefix, int indent);
    
  }
  
  /**
   * Inner class for representing simple rules
   */
  static class SimpleRule extends Rule {
    
    /** For serialization */
    private static final long serialVersionUID = -2612893679476049682L;

    /** The ID for the rule (optional) */
    protected String m_ID;
    
    /** The predicted value when the rule fires (required) */
    protected String m_scoreString;
    
    /** 
     * The predicted value as a number (regression) or index (classification)
     * when the rule fires (required)
     */
    protected double m_score = Utils.missingValue();
    
    /** The number of training/test instances on which the rule fired (optional) */
    protected double m_recordCount = Utils.missingValue();
    
    /** 
     * The number of training/test instances on which the rule fired and the
     * prediction was correct (optional)
     */
    protected double m_nbCorrect = Utils.missingValue();
    
    /** The confidence of the rule (optional) */
    protected double m_confidence = Utils.missingValue();
    
    /** The score distributions for this rule (if any) */
    protected ArrayList<TreeModel.ScoreDistribution> m_scoreDistributions = 
      new ArrayList<TreeModel.ScoreDistribution>();
    
    /**
     *  The relative importance of the rule. May or may not be equal to the
     * confidence (optional).
     */
    protected double m_weight = Utils.missingValue();
    
    public String toString(String prefix, int indent) {
      StringBuffer temp = new StringBuffer();
      
      for (int i = 0; i < indent; i++) {
        prefix += " ";
      }
      
      temp.append(prefix + "Simple rule: " + m_predicate + "\n");
      temp.append(prefix + " => " + m_scoreString + "\n");
      if (!Utils.isMissingValue(m_recordCount)) {
        temp.append(prefix + " recordCount: " + m_recordCount + "\n");
      }
      if (!Utils.isMissingValue(m_nbCorrect)) {
        temp.append(prefix + "   nbCorrect: " + m_nbCorrect + "\n");
      }
      if (!Utils.isMissingValue(m_confidence)) {
        temp.append(prefix + "  confidence: " + m_confidence + "\n");
      }
      if (!Utils.isMissingValue(m_weight)) {
        temp.append(prefix + "      weight: " + m_weight + "\n");
      }
      
      return temp.toString();
    }
    
    public String toString() {
      return toString("", 0);
    }
        
    /**
     * Constructor for a simple rule
     * 
     * @param ruleE the XML element holding the simple rule
     * @param miningSchema the mining schema to use
     * @throws Exception if something goes wrong
     */
    public SimpleRule(Element ruleE, MiningSchema miningSchema) throws Exception {
      super(ruleE, miningSchema);
      
      String id = ruleE.getAttribute("id");
      if (id != null && id.length() > 0) {
        m_ID = id;
      }
      
      m_scoreString = ruleE.getAttribute("score");
      Attribute classAtt = miningSchema.getFieldsAsInstances().classAttribute(); 
      if (classAtt.isNumeric()) {
        m_score = Double.parseDouble(m_scoreString);
      } else {
        if (classAtt.indexOfValue(m_scoreString) < 0) {
          throw new Exception("[SimpleRule] class value " + m_scoreString + 
              "does not exist in class attribute " + classAtt.name());
        }
        m_score = classAtt.indexOfValue(m_scoreString);
      }
      
      String recordCount = ruleE.getAttribute("recordCount");
      if (recordCount != null && recordCount.length() > 0) {
        m_recordCount = Double.parseDouble(recordCount);
      }
      
      String nbCorrect = ruleE.getAttribute("nbCorrect");
      if (nbCorrect != null && nbCorrect.length() > 0) {
        m_nbCorrect = Double.parseDouble(nbCorrect);
      }
      
      String confidence = ruleE.getAttribute("confidence");
      if (confidence != null && confidence.length() > 0) {
        m_confidence = Double.parseDouble(confidence);
      }
      
      String weight = ruleE.getAttribute("weight");
      if (weight != null && weight.length() > 0) {
        m_weight = Double.parseDouble(weight);
      }
      
      // get the ScoreDistributions (if any)
      if (miningSchema.getFieldsAsInstances().classAttribute().isNominal() || miningSchema.getFieldsAsInstances().classAttribute().isRanking()) {
        // see if we have any ScoreDistribution entries
        NodeList scoreChildren = ruleE.getChildNodes();
                
        for (int i = 0; i < scoreChildren.getLength(); i++) {
          Node child = scoreChildren.item(i);
          if (child.getNodeType() == Node.ELEMENT_NODE) {
            String tagName = ((Element)child).getTagName();
            if (tagName.equals("ScoreDistribution")) {
              TreeModel.ScoreDistribution newDist = 
                new TreeModel.ScoreDistribution((Element)child, 
                  miningSchema, m_recordCount);
              m_scoreDistributions.add(newDist);
            }
          }
        }
        
        // check that we have as many score distribution elements as there
        // are class labels in the data
        if (m_scoreDistributions.size() > 0 && 
            m_scoreDistributions.size() != 
              miningSchema.getFieldsAsInstances().classAttribute().numValues()) {
          throw new Exception("[SimpleRule] Number of score distribution elements is "
              + " different than the number of class labels!");
        }
        
        //backfit the confidence values (if necessary)
        if (Utils.isMissingValue(m_recordCount)) {
          double baseCount = 0;
          for (TreeModel.ScoreDistribution s : m_scoreDistributions) {
            baseCount += s.getRecordCount();
          }
          
          for (TreeModel.ScoreDistribution s : m_scoreDistributions) {
            s.deriveConfidenceValue(baseCount);
          }
        }
      }
    }
    
    /**
     * Collect the rule(s) that fire for the supplied incoming instance
     * 
     * @param input a vector of independent and derived independent variables
     * @param ruleCollection the array list to add any firing rules into
     */
    public void fires(double[] input, ArrayList<SimpleRule> ruleCollection) {
      if (m_predicate.evaluate(input) == TreeModel.Predicate.Eval.TRUE) {
        ruleCollection.add(this);
      }      
    }
    
    /**
     * Score the incoming instance
     * 
     * @param instance a vector containing the incoming independent and
     * derived independent variables
     * @param classAtt the class attribute
     * @param rsm the rule selection method (ignored by simple rules)
     * @return a probability distribution over the class labels or
     * the predicted value (in element zero of the array if the class is numeric)
     * @throws Exception if something goes wrong
     */
    public double[] score(double[] instance, Attribute classAtt) 
      throws Exception {
      
      double[] preds;
      if (classAtt.isNumeric()) {
        preds = new double[1];
        preds[0] = m_score;
      } else {
        preds = new double[classAtt.numValues()];
        if (m_scoreDistributions.size() > 0) {
          for (TreeModel.ScoreDistribution s : m_scoreDistributions) {
            preds[s.getClassLabelIndex()] = s.getConfidence();
          }
        } else if (!Utils.isMissingValue(m_confidence)) {
          preds[classAtt.indexOfValue(m_scoreString)] = m_confidence;
        } else {
          preds[classAtt.indexOfValue(m_scoreString)] = 1.0;
        }
      }      
      
      return preds;
    }
    
    /**
     * Get the weight of the rule
     * 
     * @return the weight of the rule
     */
    public double getWeight() {
      return m_weight;
    }
    
    /**
     * Get the ID of the rule
     * 
     * @return the ID of the rule
     */
    public String getID() {
      return m_ID;
    }
    
    /**
     * Get the predicted value of this rule (either a number
     * for regression problems or an index of a class label for
     * classification problems)
     * 
     * @return the predicted value of this rule
     */
    public double getScore() {
      return m_score;
    }
  }
  
  /**
   * Inner class representing a compound rule
   */
  static class CompoundRule extends Rule {
    
    /** For serialization */
    private static final long serialVersionUID = -2853658811459970718L;
    
    /** The child rules of this compound rule */
    ArrayList<Rule> m_childRules = new ArrayList<Rule>();
    
    public String toString(String prefix, int indent) {
      StringBuffer temp = new StringBuffer();

      for (int i = 0; i < indent; i++) {
        prefix += " ";
      }
      
      temp.append(prefix + "Compound rule: " + m_predicate + "\n");
      
      for (Rule r : m_childRules) {
        temp.append(r.toString(prefix, indent + 1));
      }

      return temp.toString();
    }
    
    public String toString() {
      return toString("", 0);
    }
    
    /**
     * Constructor.
     * 
     * @param ruleE XML node holding the rule
     * @param miningSchema the mining schema to use
     * @throws Exception if something goes wrong
     */
    public CompoundRule(Element ruleE, MiningSchema miningSchema) throws Exception {
      
      // get the Predicate
      super(ruleE, miningSchema);
      
      // get the nested rules
      NodeList ruleChildren = ruleE.getChildNodes();
      for (int i = 0; i < ruleChildren.getLength(); i++) {
        Node child = ruleChildren.item(i);
        if (child.getNodeType() == Node.ELEMENT_NODE) {
          String tagName = ((Element)child).getTagName();
          if (tagName.equals("SimpleRule")) {
            Rule childRule = new SimpleRule(((Element)child), miningSchema);
            m_childRules.add(childRule);
          } else if (tagName.equals("CompoundRule")) {
            Rule childRule = new CompoundRule(((Element)child), miningSchema);
            m_childRules.add(childRule);
          }
        }
      }
    }
    
    /**
     * Collect the rule(s) that fire for the supplied incoming instance
     * 
     * @param input a vector of independent and derived independent variables
     * @param ruleCollection the array list to add any firing rules into
     */
    public void fires(double[] input, ArrayList<SimpleRule> ruleCollection) {
      
      // evaluate our predicate first
      if (m_predicate.evaluate(input) == TreeModel.Predicate.Eval.TRUE) {
        // now check the child rules
        for (Rule r : m_childRules) {
          r.fires(input, ruleCollection);
        }
      }      
    }
  }
  
  /**
   * Inner class representing a set of rules 
   */
  static class RuleSet implements Serializable {
    
    /** For serialization */
    private static final long serialVersionUID = -8718126887943074376L;

    enum RuleSelectionMethod {
      WEIGHTEDSUM("weightedSum"),
      WEIGHTEDMAX("weightedMax"),
      FIRSTHIT("firstHit");
      
      private final String m_stringVal;
      
      RuleSelectionMethod(String name) {
        m_stringVal = name;
      }
      
      public String toString() {
        return m_stringVal;
      }
    }
    
    /** 
     * The number of training/test cases to which the ruleset was
     * applied to generate support and confidence measures for individual
     * rules (optional)
     */
    private double m_recordCount = Utils.missingValue();
    
    /** 
     * The number of training/test cases for which the default
     * score is correct (optional)
     */
    private double m_nbCorrect = Utils.missingValue();
    
    /** 
     * The default value to predict when no rule in the
     * ruleset fires (as a String; optional) 
     * */
    private String m_defaultScore;
    
    /** 
     * The default value to predict (either a real value or an
     * index) 
     * */
    private double m_defaultPrediction = Utils.missingValue();
    
    /** 
     * The default distribution to predict when no rule in the
     * ruleset fires (nominal class only, optional)
     */
    private ArrayList<TreeModel.ScoreDistribution> m_scoreDistributions =
      new ArrayList<TreeModel.ScoreDistribution>();
    
    /**
     * The default confidence value to return along with a score
     * when no rules in the set fire (optional)
     */
    private double m_defaultConfidence = Utils.missingValue();
    
    /** The active rule selection method */
    private RuleSelectionMethod m_currentMethod;
    
    /** The selection of rule selection methods allowed */
    private ArrayList<RuleSelectionMethod> m_availableRuleSelectionMethods = 
      new ArrayList<RuleSelectionMethod>();
    
    /** The rules contained in the rule set */
    private ArrayList<Rule> m_rules = new ArrayList<Rule>();
    
    /* (non-Javadoc)
     * @see java.lang.Object#toString()
     */
    public String toString() {
      StringBuffer temp = new StringBuffer();
      
      temp.append("Rule selection method: " + m_currentMethod + "\n");
      if (m_defaultScore != null) {
        temp.append("Default prediction: " + m_defaultScore + "\n");
        
        if (!Utils.isMissingValue(m_recordCount)) {
          temp.append("       recordCount: " + m_recordCount + "\n");
        }
        if (!Utils.isMissingValue(m_nbCorrect)) {
          temp.append("         nbCorrect: " + m_nbCorrect + "\n");
        }
        if (!Utils.isMissingValue(m_defaultConfidence)) {
          temp.append(" defaultConfidence: " + m_defaultConfidence + "\n");
        }
        
        temp.append("\n");
      }
      
      for (Rule r : m_rules) {
        temp.append(r + "\n");
      }
      
      return temp.toString();
    }
    
    /**
     * Constructor for a RuleSet.
     * 
     * @param ruleSetNode the XML node holding the RuleSet
     * @param miningSchema the mining schema to use
     * @throws Exception if something goes wrong
     */
    public RuleSet(Element ruleSetNode, MiningSchema miningSchema) 
      throws Exception {
      
      String recordCount = ruleSetNode.getAttribute("recordCount");
      if (recordCount != null && recordCount.length() > 0) {
        m_recordCount = Double.parseDouble(recordCount);
      }
      
      String nbCorrect = ruleSetNode.getAttribute("nbCorrect");
      if (nbCorrect != null & nbCorrect.length() > 0) {
        m_nbCorrect = Double.parseDouble(nbCorrect);
      }
      
      String defaultScore = ruleSetNode.getAttribute("defaultScore");
      if (defaultScore != null && defaultScore.length() > 0) {
        m_defaultScore = defaultScore;
        
        Attribute classAtt = miningSchema.getFieldsAsInstances().classAttribute();
        if (classAtt == null) {
          throw new Exception("[RuleSet] class attribute not set!");
        }
        
        if (classAtt.isNumeric()) {
          m_defaultPrediction = Double.parseDouble(defaultScore);
        } else {
          if (classAtt.indexOfValue(defaultScore) < 0) {
            throw new Exception("[RuleSet] class value " + defaultScore + 
                " not found!");
          }
          m_defaultPrediction = classAtt.indexOfValue(defaultScore);
        }
      }
      
      String defaultConfidence = ruleSetNode.getAttribute("defaultConfidence");
      if (defaultConfidence != null && defaultConfidence.length() > 0) {
        m_defaultConfidence = Double.parseDouble(defaultConfidence);
      }
      
      // get the rule selection methods
      NodeList selectionNL = ruleSetNode.getElementsByTagName("RuleSelectionMethod");
      for (int i = 0; i < selectionNL.getLength(); i++) {
        Node selectN = selectionNL.item(i);
        if (selectN.getNodeType() == Node.ELEMENT_NODE) {
          Element sN = (Element)selectN;
          String criterion = sN.getAttribute("criterion");
          for (RuleSelectionMethod m : RuleSelectionMethod.values()) {
            if (m.toString().equals(criterion)) {
              m_availableRuleSelectionMethods.add(m);
              if (i == 0) {
                // set the default (first specified one)
                m_currentMethod = m;
              }
            }
          }
        }
      }
      
      if (miningSchema.getFieldsAsInstances().classAttribute().isNominal() || miningSchema.getFieldsAsInstances().classAttribute().isRanking()) {
        // see if we have any ScoreDistribution entries
        NodeList scoreChildren = ruleSetNode.getChildNodes();
        for (int i = 0; i < scoreChildren.getLength(); i++) {
          Node child = scoreChildren.item(i);
          if (child.getNodeType() == Node.ELEMENT_NODE) {
            String tagName = ((Element)child).getTagName();
            if (tagName.equals("ScoreDistribution")) {
              TreeModel.ScoreDistribution newDist = 
                new TreeModel.ScoreDistribution((Element)child, 
                  miningSchema, m_recordCount);
              m_scoreDistributions.add(newDist);
            }
          }
        }
        
        //backfit the confidence values (if necessary)
        if (Utils.isMissingValue(m_recordCount)) {
          double baseCount = 0;
          for (TreeModel.ScoreDistribution s : m_scoreDistributions) {
            baseCount += s.getRecordCount();
          }
          
          for (TreeModel.ScoreDistribution s : m_scoreDistributions) {
            s.deriveConfidenceValue(baseCount);
          }
        }
      }
      
      // Get the rules in this rule set
      NodeList ruleChildren = ruleSetNode.getChildNodes();
      for (int i = 0; i < ruleChildren.getLength(); i++) {
        Node child = ruleChildren.item(i);
        if (child.getNodeType() == Node.ELEMENT_NODE) {
          String tagName = ((Element)child).getTagName();
          if (tagName.equals("SimpleRule")) {
            Rule tempRule = new SimpleRule(((Element)child), miningSchema);
            m_rules.add(tempRule);
          } else if (tagName.equals("CompoundRule")) {
            Rule tempRule = new CompoundRule(((Element)child), miningSchema);
            m_rules.add(tempRule);
          }
        }
      }
    }
    
    /**
     * Score an incoming instance by collecting all rules that fire.
     * 
     * @param instance a vector of incoming attribte and derived field values
     * @param classAtt the class attribute
     * @return a predicted probability distribution
     * @throws Exception is something goes wrong
     */
    protected double[] score(double[] instance, Attribute classAtt)
      throws Exception {
      
      double[] preds = null;
      if (classAtt.isNumeric()) {
        preds = new double[1];
      } else {
        preds = new double[classAtt.numValues()];
      }
      
      // holds the rules that fire for this test case
      ArrayList<SimpleRule> firingRules = new ArrayList<SimpleRule>();
      
      for (Rule r : m_rules) {
        r.fires(instance, firingRules);
      }
      
      if (firingRules.size() > 0) {
        if (m_currentMethod == RuleSelectionMethod.FIRSTHIT) {
          preds = firingRules.get(0).score(instance, classAtt);
        } else if (m_currentMethod == RuleSelectionMethod.WEIGHTEDMAX) {
          double wMax = Double.NEGATIVE_INFINITY;
          SimpleRule best = null;
          for (SimpleRule s : firingRules) {
            if (Utils.isMissingValue(s.getWeight())) {
              throw new Exception("[RuleSet] Scoring criterion is WEIGHTEDMAX, but " +
              		"rule " + s.getID() + " does not have a weight defined!");
            }
            if (s.getWeight() > wMax) {
              wMax = s.getWeight();
              best = s;
            }
          }
          if (best == null) {
            throw new Exception("[RuleSet] Unable to determine the best rule under " +
            		"the WEIGHTEDMAX criterion!");
          }
          preds = best.score(instance, classAtt);          
        } else if (m_currentMethod == RuleSelectionMethod.WEIGHTEDSUM) {
          double sumOfWeights = 0;
          for (SimpleRule s : firingRules) {
            if (Utils.isMissingValue(s.getWeight())) {
              throw new Exception("[RuleSet] Scoring criterion is WEIGHTEDSUM, but " +
                        "rule " + s.getID() + " does not have a weight defined!");
            }            
            if (classAtt.isNumeric()) {
              sumOfWeights += s.getWeight();
              preds[0] += (s.getScore() * s.getWeight());
            } else {
              preds[(int)s.getScore()] += s.getWeight();
            }
          }
          if (classAtt.isNumeric()) {
            if (sumOfWeights == 0) {
              throw new Exception("[RuleSet] Sum of weights is zero!");
            }
            preds[0] /= sumOfWeights;
          } else {
            // array gets normalized in the distributionForInstance() method
          }
        }
      } else {
        // default prediction
        if (classAtt.isNumeric()) {
          preds[0] = m_defaultPrediction;
        } else {
          if (m_scoreDistributions.size() > 0) {
            for (TreeModel.ScoreDistribution s : m_scoreDistributions) {
              preds[s.getClassLabelIndex()] = s.getConfidence();
            }
          } else if (!Utils.isMissingValue(m_defaultConfidence)) {
            preds[(int)m_defaultPrediction] = m_defaultConfidence;
          } else {
            preds[(int)m_defaultPrediction] = 1.0;
          }
        }
      }
      
      return preds;
    }
  }
  
  /** The mining function */
  protected MiningFunction m_functionType = MiningFunction.CLASSIFICATION;
  
  /** The model name (if defined) */
  protected String m_modelName;
  
  /** The algorithm name (if defined) */
  protected String m_algorithmName;
  
  /** The set of rules */
  protected RuleSet m_ruleSet;

  /**
   * Constructor for a RuleSetModel
   * 
   * @param model the XML element encapsulating the RuleSetModel
   * @param dataDictionary the data dictionary to use
   * @param miningSchema the mining schema to use
   * @throws Exception if something goes wrong
   */
  public RuleSetModel(Element model, Instances dataDictionary,
      MiningSchema miningSchema) throws Exception {
    
    super(dataDictionary, miningSchema);
    
    if (!getPMMLVersion().equals("3.2")) {
      // TODO: might have to throw an exception and only support 3.2
    }
    
    String fn = model.getAttribute("functionName");
    if (fn.equals("regression")) {
      m_functionType = MiningFunction.REGRESSION;
    }
    
    String modelName = model.getAttribute("modelName");
    if (modelName != null && modelName.length() > 0) {
      m_modelName = modelName;
    }
    
    String algoName = model.getAttribute("algorithmName");
    if (algoName != null && algoName.length() > 0) {
      m_algorithmName = algoName;
    }    
    
    NodeList ruleset = model.getElementsByTagName("RuleSet");
    if (ruleset.getLength() == 1) {
      Node ruleSetNode = ruleset.item(0);
      if (ruleSetNode.getNodeType() == Node.ELEMENT_NODE) {
        m_ruleSet = new RuleSet((Element)ruleSetNode, miningSchema);
      }
    } else {
      throw new Exception ("[RuleSetModel] Should only have a single RuleSet!");
    }
  }
  
  /**                                                                                                             
   * Classifies the given test instance. The instance has to belong to a                                          
   * dataset when it's being classified.                                                          
   *                                                                                                              
   * @param inst the instance to be classified                                                                
   * @return the predicted most likely class for the instance or                                                  
   * Utils.missingValue() if no prediction is made                                                             
   * @exception Exception if an error occurred during the prediction                                              
   */
  public double[] distributionForInstance(Instance inst) throws Exception {
    if (!m_initialized) {
      mapToMiningSchema(inst.dataset());
    }
    double[] preds = null;
    
    if (m_miningSchema.getFieldsAsInstances().classAttribute().isNumeric()) {
      preds = new double[1];
    } else {
      preds = new double[m_miningSchema.getFieldsAsInstances().classAttribute().numValues()];
    }
    
    double[] incoming = m_fieldsMap.instanceToSchema(inst, m_miningSchema);
    
    preds = m_ruleSet.score(incoming, 
        m_miningSchema.getFieldsAsInstances().classAttribute());
    
    if (m_miningSchema.getFieldsAsInstances().classAttribute().isNominal() || m_miningSchema.getFieldsAsInstances().classAttribute().isRanking()) {
      Utils.normalize(preds);
    }
    
    return preds;
  }
  
  /**
   * Return a textual description of this model.
   * 
   * @return a textual description of this model
   */
  public String toString() {
    StringBuffer temp = new StringBuffer();
    
    temp.append("PMML version " + getPMMLVersion());
    if (!getCreatorApplication().equals("?")) {
      temp.append("\nApplication: " + getCreatorApplication());
    }
    temp.append("\nPMML Model: RuleSetModel");
    temp.append("\n\n");
    temp.append(m_miningSchema);
    
    if (m_algorithmName != null) {
      temp.append("\nAlgorithm: " + m_algorithmName + "\n");
    }
    
    temp.append(m_ruleSet);
  
    return temp.toString();
  }

  /**
   * Get the revision string for this class
   * 
   * @return the revision string
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }
  
}
