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
 *    TreeModel.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.pmml.consumer;

import java.io.Serializable;
import java.util.ArrayList;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.core.Attribute;
import weka.core.Drawable;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.pmml.*;

/**
 * Class implementing import of PMML TreeModel. Can be used as a Weka
 * classifier for prediction (buildClassifier() raises and Exception).
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 5987 $;
 */
public class TreeModel extends PMMLClassifier implements Drawable {
  
  /**
   * For serialization
   */
  private static final long serialVersionUID = -2065158088298753129L;

  /**
   * Inner class representing the ScoreDistribution element
   */
  static class ScoreDistribution implements Serializable {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = -123506262094299933L;

    /** The class label for this distribution element */
    private String m_classLabel;
    
    /** The index of the class label */
    private int m_classLabelIndex = -1;
    
    /** The count for this label */
    private double m_recordCount;
    
    /** The optional confidence value */
    private double m_confidence = Utils.missingValue();
    
    /**
     * Construct a ScoreDistribution entry
     * 
     * @param scoreE the node containing the distribution
     * @param miningSchema the mining schema
     * @param baseCount the number of records at the node that owns this 
     * distribution entry
     * @throws Exception if something goes wrong
     */
    protected ScoreDistribution(Element scoreE, MiningSchema miningSchema, double baseCount) 
      throws Exception {
      // get the label
      m_classLabel = scoreE.getAttribute("value");
      Attribute classAtt = miningSchema.getFieldsAsInstances().classAttribute();
      if (classAtt == null || classAtt.indexOfValue(m_classLabel) < 0) {
        throw new Exception("[ScoreDistribution] class attribute not set or class value " +
            m_classLabel + " not found!");
      }
      
      m_classLabelIndex = classAtt.indexOfValue(m_classLabel);
      
      // get the frequency
      String recordC = scoreE.getAttribute("recordCount");
      m_recordCount = Double.parseDouble(recordC);
      
      // get the optional confidence
      String confidence = scoreE.getAttribute("confidence");
      if (confidence != null && confidence.length() > 0) {
        m_confidence = Double.parseDouble(confidence);        
      } else if (!Utils.isMissingValue(baseCount) && baseCount > 0) {
        m_confidence = m_recordCount / baseCount;
      }
    }
    
    /**
     * Backfit confidence value (does nothing if the confidence
     * value is already set).
     * 
     * @param baseCount the total number of records (supplied either
     * explicitly from the node that owns this distribution entry
     * or most likely computed from summing the recordCounts of all
     * the distribution entries in the distribution that owns this
     * entry).
     */
    void deriveConfidenceValue(double baseCount) {
      if (Utils.isMissingValue(m_confidence) && 
          !Utils.isMissingValue(baseCount) && 
          baseCount > 0) {
        m_confidence = m_recordCount / baseCount;
      }
    }
    
    String getClassLabel() {
      return m_classLabel;
    }
    
    int getClassLabelIndex() {
      return m_classLabelIndex;
    }
    
    double getRecordCount() {
      return m_recordCount;
    }
    
    double getConfidence() {
      return m_confidence;
    }
    
    public String toString() {
      return m_classLabel + ": " + m_recordCount 
        + " (" + Utils.doubleToString(m_confidence, 2) + ") ";
    }
  }
  
  /**
   * Base class for Predicates
   */
  static abstract class Predicate implements Serializable {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = 1035344165452733887L;

    enum Eval {
      TRUE,
      FALSE,
      UNKNOWN;
    }
    
    /**
     * Evaluate this predicate.
     * 
     * @param input the input vector of attribute and derived field values.
     * 
     * @return the evaluation status of this predicate.
     */
    abstract Eval evaluate(double[] input);
    
    protected String toString(int level, boolean cr) {
      return toString(level);
    }
    
    protected String toString(int level) {
      StringBuffer text = new StringBuffer();
      for (int j = 0; j < level; j++) {
        text.append("|   ");
      }
      
      return text.append(toString()).toString();
    }
    
    static Eval booleanToEval(boolean missing, boolean result) {
      if (missing) {
        return Eval.UNKNOWN;
      } else if (result) {
        return Eval.TRUE;
      } else {
        return Eval.FALSE;
      }
    }
    
    /**
     * Factory method to return the appropriate predicate for
     * a given node in the tree.
     * 
     * @param nodeE the XML node encapsulating the tree node.
     * @param miningSchema the mining schema in use
     * @return a Predicate
     * @throws Exception of something goes wrong.
     */
    static Predicate getPredicate(Element nodeE, 
        MiningSchema miningSchema) throws Exception {
      
      Predicate result = null;
      NodeList children = nodeE.getChildNodes();
      for (int i = 0; i < children.getLength(); i++) {
        Node child = children.item(i);
        if (child.getNodeType() == Node.ELEMENT_NODE) {
          String tagName = ((Element)child).getTagName();
          if (tagName.equals("True")) {
            result = new True();
            break;
          } else if (tagName.equals("False")) {
            result = new False();
            break;
          } else if (tagName.equals("SimplePredicate")) {
            result = new SimplePredicate((Element)child, miningSchema);
            break;
          } else if (tagName.equals("CompoundPredicate")) {
            result = new CompoundPredicate((Element)child, miningSchema);
            break;
          } else if (tagName.equals("SimpleSetPredicate")) {
           result = new SimpleSetPredicate((Element)child, miningSchema);
           break;
          }
        }
      }
      
      if (result == null) {
        throw new Exception("[Predicate] unknown or missing predicate type in node");
      }
      
      return result;
    }
  }
  
  /**
   * Simple True Predicate
   */
  static class True extends Predicate {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = 1817942234610531627L;

    public Predicate.Eval evaluate(double[] input) {
      return Predicate.Eval.TRUE;
    }
    
    public String toString() {
      return "True: ";
    }
  }
  
  /**
   * Simple False Predicate
   */
  static class False extends Predicate {
    
    /**
     * For serialization 
     */
    private static final long serialVersionUID = -3647261386442860365L;

    public Predicate.Eval evaluate(double[] input) {
      return Predicate.Eval.FALSE;
    }
    
    public String toString() {
      return "False: ";
    }
  }
  
  /**
   * Class representing the SimplePredicate
   */
  static class SimplePredicate extends Predicate {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = -6156684285069327400L;

    enum Operator {
      EQUAL("equal") {
        Predicate.Eval evaluate(double[] input, double value, int fieldIndex) {
          return Predicate.booleanToEval(Utils.isMissingValue(input[fieldIndex]), 
              weka.core.Utils.eq(input[fieldIndex], value));
        }
        
        String shortName() {
          return "==";
        }
      },
      NOTEQUAL("notEqual")
       {
        Predicate.Eval evaluate(double[] input, double value, int fieldIndex) {
          return Predicate.booleanToEval(Utils.isMissingValue(input[fieldIndex]), 
              (input[fieldIndex] != value));
        }
        
        String shortName() {
          return "!=";
        }
      },
      LESSTHAN("lessThan")  {
        Predicate.Eval evaluate(double[] input, double value, int fieldIndex) {
          return Predicate.booleanToEval(Utils.isMissingValue(input[fieldIndex]),
              (input[fieldIndex] < value));
        }
        
        String shortName() {
          return "<";
        }
      },
      LESSOREQUAL("lessOrEqual") {
        Predicate.Eval evaluate(double[] input, double value, int fieldIndex) {
          return Predicate.booleanToEval(Utils.isMissingValue(input[fieldIndex]),
              (input[fieldIndex] <= value));
        }
        
        String shortName() {
          return "<=";
        }
      },
      GREATERTHAN("greaterThan") {
        Predicate.Eval evaluate(double[] input, double value, int fieldIndex) {
          return Predicate.booleanToEval(Utils.isMissingValue(input[fieldIndex]),
              (input[fieldIndex] > value));
        }
        
        String shortName() {
          return ">";
        }
      },
      GREATEROREQUAL("greaterOrEqual") {
        Predicate.Eval evaluate(double[] input, double value, int fieldIndex) {
          return Predicate.booleanToEval(Utils.isMissingValue(input[fieldIndex]),
              (input[fieldIndex] >= value));
        }
        
        String shortName() {
          return ">=";
        }
      },
      ISMISSING("isMissing") {
        Predicate.Eval evaluate(double[] input, double value, int fieldIndex) {
          return Predicate.booleanToEval(false,
              Utils.isMissingValue(input[fieldIndex]));
        }
        
        String shortName() {
          return toString();
        }
      },
      ISNOTMISSING("isNotMissing") {
        Predicate.Eval evaluate(double[] input, double value, int fieldIndex) {
          return Predicate.booleanToEval(false, !Utils.isMissingValue(input[fieldIndex]));
        }
        
        String shortName() {
          return toString();
        }
      };
      
      abstract Predicate.Eval evaluate(double[] input, double value, int fieldIndex);
      abstract String shortName();
      
      private final String m_stringVal;
      
      Operator(String name) {
        m_stringVal = name;
      }
            
      public String toString() {
        return m_stringVal;
      }
    }
    
    /** the field that we are comparing against */
    int m_fieldIndex = -1;
    
    /** the name of the field */
    String m_fieldName;
    
    /** true if the field is nominal */
    boolean m_isNominal;
    
    /** the value as a string (if nominal) */
    String m_nominalValue;
    
    /** the value to compare against (if nominal it holds the index of the value) */
    double m_value;
    
    /** the operator to use */
    Operator m_operator;
        
    public SimplePredicate(Element simpleP, 
        MiningSchema miningSchema) throws Exception {
      Instances totalStructure = miningSchema.getFieldsAsInstances();
      
      // get the field name and set up the index
      String fieldS = simpleP.getAttribute("field");
      Attribute att = totalStructure.attribute(fieldS);
      if (att == null) {
        throw new Exception("[SimplePredicate] unable to find field " + fieldS
            + " in the incoming instance structure!");
      }
      
      // find the index
      int index = -1;
      for (int i = 0; i < totalStructure.numAttributes(); i++) {
        if (totalStructure.attribute(i).name().equals(fieldS)) {
          index = i;
          m_fieldName = totalStructure.attribute(i).name();
          break;
        }
      }
      m_fieldIndex = index;
      if (att.isNominal() || att.isRanking()) {
        m_isNominal = true;
      }
      
      // get the operator
      String oppS = simpleP.getAttribute("operator");
      for (Operator o : Operator.values()) {
        if (o.toString().equals(oppS)) {
          m_operator = o;
          break;
        }
      }
      
      if (m_operator != Operator.ISMISSING && m_operator != Operator.ISNOTMISSING) {
        String valueS = simpleP.getAttribute("value");
        if (att.isNumeric()) {
          m_value = Double.parseDouble(valueS);
        } else {
          m_nominalValue = valueS;
          m_value = att.indexOfValue(valueS);
          if (m_value < 0) {
            throw new Exception("[SimplePredicate] can't find value " + valueS + " in nominal " +
                "attribute " + att.name());
          }
        }
      }
    }
    
    public Predicate.Eval evaluate(double[] input) {
      return m_operator.evaluate(input, m_value, m_fieldIndex);
    }
        
    public String toString() {
      StringBuffer temp = new StringBuffer();
      
      temp.append(m_fieldName + " " + m_operator.shortName());
      if (m_operator != Operator.ISMISSING && m_operator != Operator.ISNOTMISSING) {
        temp.append(" " + ((m_isNominal) ? m_nominalValue : "" + m_value));
      }
      
      return temp.toString();
    }
  }
  
  /**
   * Class representing the CompoundPredicate
   */
  static class CompoundPredicate extends Predicate {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = -3332091529764559077L;

    enum BooleanOperator {
      OR("or") {
        Predicate.Eval evaluate(ArrayList<Predicate> constituents, double[] input) {
          Predicate.Eval currentStatus = Predicate.Eval.FALSE;
          for (Predicate p : constituents) {
            Predicate.Eval temp = p.evaluate(input);
            if (temp == Predicate.Eval.TRUE) {
              currentStatus = temp;
              break;
            } else if (temp == Predicate.Eval.UNKNOWN) {
              currentStatus = temp;
            }            
          }
          return currentStatus;
        }
      },
      AND("and") {
        Predicate.Eval evaluate(ArrayList<Predicate> constituents, double[] input) {
          Predicate.Eval currentStatus = Predicate.Eval.TRUE;
          for (Predicate p : constituents) {
            Predicate.Eval temp = p.evaluate(input);
            if (temp == Predicate.Eval.FALSE) {
              currentStatus = temp;
              break;
            } else if (temp == Predicate.Eval.UNKNOWN) {
              currentStatus = temp;
            }
          }          
          return currentStatus;
        }
      },
      XOR("xor") {
        Predicate.Eval evaluate(ArrayList<Predicate> constituents, double[] input) {
          Predicate.Eval currentStatus = constituents.get(0).evaluate(input);
          if (currentStatus != Predicate.Eval.UNKNOWN) {
            for (int i = 1; i < constituents.size(); i++) {
              Predicate.Eval temp = constituents.get(i).evaluate(input);
              if (temp == Predicate.Eval.UNKNOWN) {
                currentStatus = temp;
                break;
              } else {
                if (currentStatus != temp) {
                  currentStatus = Predicate.Eval.TRUE;
                } else {
                  currentStatus = Predicate.Eval.FALSE;
                }
              }
            }
          }
          return currentStatus;
        }
      },
      SURROGATE("surrogate") {
        Predicate.Eval evaluate(ArrayList<Predicate> constituents, double[] input) {
          Predicate.Eval currentStatus = constituents.get(0).evaluate(input);
          
          int i = 1;
          while (currentStatus == Predicate.Eval.UNKNOWN) {
            currentStatus = constituents.get(i).evaluate(input);            
          }
          
          // return false if all our surrogates evaluate to unknown.
          if (currentStatus == Predicate.Eval.UNKNOWN) {
            currentStatus = Predicate.Eval.FALSE;
          }
          
          return currentStatus;
        }
      };
      
      abstract Predicate.Eval evaluate(ArrayList<Predicate> constituents, double[] input);
      
      private final String m_stringVal;
      
      BooleanOperator(String name) {
        m_stringVal = name;
      }
      
      public String toString() {
        return m_stringVal;
      }
    }
    
    /** the constituent Predicates */
    ArrayList<Predicate> m_components = new ArrayList<Predicate>();
    
    /** the boolean operator */
    BooleanOperator m_booleanOperator;
        
    public CompoundPredicate(Element compoundP, 
        MiningSchema miningSchema) throws Exception {
//      Instances totalStructure = miningSchema.getFieldsAsInstances();
      
      String booleanOpp = compoundP.getAttribute("booleanOperator");
      for (BooleanOperator b : BooleanOperator.values()) {
        if (b.toString().equals(booleanOpp)) {
          m_booleanOperator = b;
        }
      }
      
      // now get all the encapsulated operators
      NodeList children = compoundP.getChildNodes();
      for (int i = 0; i < children.getLength(); i++) {
        Node child = children.item(i);
        if (child.getNodeType() == Node.ELEMENT_NODE) {
          String tagName = ((Element)child).getTagName();
          if (tagName.equals("True")) {
            m_components.add(new True());
          } else if (tagName.equals("False")) {
            m_components.add(new False());
          } else if (tagName.equals("SimplePredicate")) {
            m_components.add(new SimplePredicate((Element)child, miningSchema));
          } else if (tagName.equals("CompoundPredicate")) {
            m_components.add(new CompoundPredicate((Element)child, miningSchema));
          } else {
            m_components.add(new SimpleSetPredicate((Element)child, miningSchema));
          }
        }
      }
    }
    
    public Predicate.Eval evaluate(double[] input) {
      return m_booleanOperator.evaluate(m_components, input);
    }
    
    public String toString() {
      return toString(0, false);
    }
    
    public String toString(int level, boolean cr) {
      StringBuffer text = new StringBuffer();
      for (int j = 0; j < level; j++) {
        text.append("|   ");
      }
      
      text.append("Compound [" + m_booleanOperator.toString() + "]");
      if (cr) {
        text.append("\\n");
      } else {
        text.append("\n");
      }
      for (int i = 0; i < m_components.size(); i++) {
        text.append(m_components.get(i).toString(level, cr).replace(":", ""));
        if (i != m_components.size()-1) {
          if (cr) {
            text.append("\\n");
          } else {
            text.append("\n");
          }
        }
      }
      
      return text.toString();
    }
  }
  
  /**
   * Class representing the SimpleSetPredicate
   */
  static class SimpleSetPredicate extends Predicate {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = -2711995401345708486L;

    enum BooleanOperator {
        IS_IN("isIn") {
          Predicate.Eval evaluate(double[] input, int fieldIndex, 
              Array set, Attribute nominalLookup) {            
            if (set.getType() == Array.ArrayType.STRING) {
              String value = "";
              if (!Utils.isMissingValue(input[fieldIndex])) {
                value = nominalLookup.value((int)input[fieldIndex]);
              }
              return Predicate.booleanToEval(Utils.isMissingValue(input[fieldIndex]), 
                  set.contains(value));
            } else if (set.getType() == Array.ArrayType.NUM ||
                set.getType() == Array.ArrayType.REAL) {
              return Predicate.booleanToEval(Utils.isMissingValue(input[fieldIndex]), 
                set.contains(input[fieldIndex]));
            }
            return Predicate.booleanToEval(Utils.isMissingValue(input[fieldIndex]), 
                set.contains((int)input[fieldIndex]));
          }
        },
        IS_NOT_IN("isNotIn") {
          Predicate.Eval evaluate(double[] input, int fieldIndex,
              Array set, Attribute nominalLookup) {
            Predicate.Eval result = IS_IN.evaluate(input, fieldIndex, set, nominalLookup);
            if (result == Predicate.Eval.FALSE) {
              result = Predicate.Eval.TRUE;
            } else if (result == Predicate.Eval.TRUE) {
              result = Predicate.Eval.FALSE;
            }
            
            return result;
          }
        };
        
        abstract Predicate.Eval evaluate(double[] input, int fieldIndex, 
            Array set, Attribute nominalLookup);
        
        private final String m_stringVal;
        
        BooleanOperator(String name) {
          m_stringVal = name;
        }
        
        public String toString() {
          return m_stringVal;
        }
    }
    
    /** the field to reference */
    int m_fieldIndex = -1;
    
    /** the name of the field */
    String m_fieldName;
    
    /** is the referenced field nominal? */
    boolean m_isNominal = false;
    
    /** the attribute to lookup nominal values from */
    Attribute m_nominalLookup;
    
    /** the boolean operator */
    BooleanOperator m_operator = BooleanOperator.IS_IN;
    
    /** the array holding the set of values */
    Array m_set;
        
    public SimpleSetPredicate(Element setP, 
        MiningSchema miningSchema) throws Exception {
      Instances totalStructure = miningSchema.getFieldsAsInstances();
      
      // get the field name and set up the index
      String fieldS = setP.getAttribute("field");
      Attribute att = totalStructure.attribute(fieldS);
      if (att == null) {
        throw new Exception("[SimplePredicate] unable to find field " + fieldS
            + " in the incoming instance structure!");
      }
      
      // find the index
      int index = -1;
      for (int i = 0; i < totalStructure.numAttributes(); i++) {
        if (totalStructure.attribute(i).name().equals(fieldS)) {
          index = i;
          m_fieldName = totalStructure.attribute(i).name();
          break;
        }
      }
      m_fieldIndex = index;
      if (att.isNominal() || att.isRanking()) {
        m_isNominal = true;
        m_nominalLookup = att;
      }
  
      // need to scan the children looking for an array type
      NodeList children = setP.getChildNodes();
      for (int i = 0; i < children.getLength(); i++) {
        Node child = children.item(i);
        if (child.getNodeType() == Node.ELEMENT_NODE) {
          if (Array.isArray((Element)child)) {
            // found the array
            m_set = Array.create((Element)child);
            break;
          }
        }
      }

      if (m_set == null) {
        throw new Exception("[SimpleSetPredictate] couldn't find an " +
        "array containing the set values!");
      }
      
      // check array type against field type
      if (m_set.getType() == Array.ArrayType.STRING &&
          !m_isNominal) {
        throw new Exception("[SimpleSetPredicate] referenced field " +
            totalStructure.attribute(m_fieldIndex).name() + 
            " is numeric but array type is string!");
      } else if (m_set.getType() != Array.ArrayType.STRING && 
          m_isNominal) {
        throw new Exception("[SimpleSetPredicate] referenced field " +
            totalStructure.attribute(m_fieldIndex).name() +
            " is nominal but array type is numeric!");
      }      
    }
    
    public Predicate.Eval evaluate(double[] input) {
      return m_operator.evaluate(input, m_fieldIndex, m_set, m_nominalLookup);
    }
    
    public String toString() {
      StringBuffer temp = new StringBuffer();
      
      temp.append(m_fieldName + " " + m_operator.toString() + " ");
      temp.append(m_set.toString());
      
      return temp.toString();
    }
  }
  
  /**
   * Class for handling a Node in the tree
   */
  class TreeNode implements Serializable {
    // TODO: perhaps implement a class called Statistics that contains Partitions?
        
    /**
     * For serialization
     */
    private static final long serialVersionUID = 3011062274167063699L;

    /** ID for this node */
    private String m_ID = "" + this.hashCode();
    
    /** The score as a string */
    private String m_scoreString;
    
    /** The index of this predicted value (if class is nominal) */
    private int m_scoreIndex = -1;
    
    /** The score as a number (if target is numeric) */
    private double m_scoreNumeric = Utils.missingValue();
    
    /** The record count at this node (if defined) */
    private double m_recordCount = Utils.missingValue();
    
    /** The ID of the default child (if applicable) */
    private String m_defaultChildID;
    
    /** Holds the node of the default child (if defined) */
    private TreeNode m_defaultChild;
    
    /** The distribution for labels (classification) */
    private ArrayList<ScoreDistribution> m_scoreDistributions = 
      new ArrayList<ScoreDistribution>();
    
    /** The predicate for this node */
    private Predicate m_predicate;
    
    /** The children of this node */
    private ArrayList<TreeNode> m_childNodes = new ArrayList<TreeNode>();
    
    
    protected TreeNode(Element nodeE, MiningSchema miningSchema) throws Exception {
      Attribute classAtt = miningSchema.getFieldsAsInstances().classAttribute();
      
      // get the ID
      String id = nodeE.getAttribute("id");
      if (id != null && id.length() > 0) {
        m_ID = id;
      }
      
      // get the score for this node
      String scoreS = nodeE.getAttribute("score");
      if (scoreS != null && scoreS.length() > 0) {
        m_scoreString = scoreS;
        
        // try to parse as a number in case we 
        // are part of a regression tree
        if (classAtt.isNumeric()) {
          try {
            m_scoreNumeric = Double.parseDouble(scoreS);
          } catch (NumberFormatException ex) {
            throw new Exception("[TreeNode] class is numeric but unable to parse score " 
                + m_scoreString + " as a number!");
          }
        } else {
          // store the index of this class value
          m_scoreIndex = classAtt.indexOfValue(m_scoreString);
          
          if (m_scoreIndex < 0) {
            throw new Exception("[TreeNode] can't find match for predicted value " 
                + m_scoreString + " in class attribute!");
          }
        }
      }
      
      // get the record count if defined
      String recordC = nodeE.getAttribute("recordCount");
      if (recordC != null && recordC.length() > 0) {
        m_recordCount = Double.parseDouble(recordC);
      }
      
      // get the default child (if applicable)
      String defaultC = nodeE.getAttribute("defaultChild");
      if (defaultC != null && defaultC.length() > 0) {
        m_defaultChildID = defaultC;
      }
      
      //TODO: Embedded model (once we support model composition)
      
      // Now get the ScoreDistributions (if any and mining function 
      // is classification) at this level
      if (m_functionType == MiningFunction.CLASSIFICATION) {
        getScoreDistributions(nodeE, miningSchema);
      }
      
      // Now get the Predicate
      m_predicate = Predicate.getPredicate(nodeE, miningSchema);
      
      // Now get the child Node(s)
      getChildNodes(nodeE, miningSchema);
      
      // If we have a default child specified, find it now
      if (m_defaultChildID != null) {
        for (TreeNode t : m_childNodes) {
          if (t.getID().equals(m_defaultChildID)) {
            m_defaultChild = t;
            break;
          }
        }
      }
    }
    
    private void getChildNodes(Element nodeE, MiningSchema miningSchema) throws Exception {
      NodeList children = nodeE.getChildNodes();
      
      for (int i = 0; i < children.getLength(); i++) {
        Node child = children.item(i);
        if (child.getNodeType() == Node.ELEMENT_NODE) {
          String tagName = ((Element)child).getTagName();
          if (tagName.equals("Node")) {
            TreeNode tempN = new TreeNode((Element)child, miningSchema);
            m_childNodes.add(tempN);
          }
        }
      }
    }
    
    private void getScoreDistributions(Element nodeE, 
        MiningSchema miningSchema) throws Exception {
      
      NodeList scoreChildren = nodeE.getChildNodes();
      for (int i = 0; i < scoreChildren.getLength(); i++) {
        Node child = scoreChildren.item(i);
        if (child.getNodeType() == Node.ELEMENT_NODE) {
          String tagName = ((Element)child).getTagName();
          if (tagName.equals("ScoreDistribution")) {
            ScoreDistribution newDist = new ScoreDistribution((Element)child, 
                miningSchema, m_recordCount);
            m_scoreDistributions.add(newDist);
          }
        }
      }
      
      // backfit the confidence values
      if (Utils.isMissingValue(m_recordCount)) {
        double baseCount = 0;
        for (ScoreDistribution s : m_scoreDistributions) {
          baseCount += s.getRecordCount();
        }
        
        for (ScoreDistribution s : m_scoreDistributions) {
          s.deriveConfidenceValue(baseCount);
        }
      }
    }
        
    /**
     * Get the score value as a string.
     * 
     * @return the score value as a String.
     */
    protected String getScore() {
      return m_scoreString;
    }
    
    /**
     * Get the score value as a number (regression trees only).
     * 
     * @return the score as a number
     */
    protected double getScoreNumeric() {
      return m_scoreNumeric;
    }
    
    /**
     * Get the ID of this node.
     * 
     * @return the ID of this node.
     */
    protected String getID() {
      return m_ID;
    }
    
    /**
     * Get the Predicate at this node.
     * 
     * @return the predicate at this node.
     */
    protected Predicate getPredicate() {
      return m_predicate;
    }
    
    /**
     * Get the record count at this node.
     * 
     * @return the record count at this node.
     */
    protected double getRecordCount() {
      return m_recordCount;
    }
    
    protected void dumpGraph(StringBuffer text) throws Exception {
      text.append("N" + m_ID + " ");
      if (m_scoreString != null) {
        text.append("[label=\"score=" + m_scoreString);
      }
      
      if (m_scoreDistributions.size() > 0 && m_childNodes.size() == 0) {
        text.append("\\n");
        for (ScoreDistribution s : m_scoreDistributions) {
          text.append(s + "\\n");
        }
      }
      
      text.append("\"");
      
      if (m_childNodes.size() == 0) {
        text.append(" shape=box style=filled");
        
      }
      
      text.append("]\n");
      
      for (TreeNode c : m_childNodes) {
        text.append("N" + m_ID +"->" + "N" + c.getID());
        text.append(" [label=\"" + c.getPredicate().toString(0, true));
        text.append("\"]\n");
        c.dumpGraph(text);
      }
    }
    
    public String toString() {
      StringBuffer text = new StringBuffer();
      
      // print out the root
      dumpTree(0, text);

      return text.toString();
    }
    
    protected void dumpTree(int level, StringBuffer text) {
      if (m_childNodes.size() > 0) {

        for (int i = 0; i < m_childNodes.size(); i++) {
          text.append("\n");
          
/*          for (int j = 0; j < level; j++) {
            text.append("|   ");
          } */
          
          // output the predicate for this child node
          TreeNode child = m_childNodes.get(i);
          text.append(child.getPredicate().toString(level, false));
          
          // process recursively
          child.dumpTree(level + 1 , text);          
        }
      } else {
        // leaf
        text.append(": ");
        if (!Utils.isMissingValue(m_scoreNumeric)) {
          text.append(m_scoreNumeric);
        } else {
          text.append(m_scoreString + " ");
          if (m_scoreDistributions.size() > 0) {
            text.append("[");
            for (ScoreDistribution s : m_scoreDistributions) {
              text.append(s);
            }
            text.append("]");
          } else {
            text.append(m_scoreString);
          }
        }
      }
    }
    
    /**
     * Score an incoming instance. Invokes a missing value handling strategy.
     * 
     * @param instance a vector of incoming attribute and derived field values.
     * @param classAtt the class attribute
     * @return a predicted probability distribution.
     * @throws Exception if something goes wrong.
     */
    protected double[] score(double[] instance, Attribute classAtt) throws Exception {
      double[] preds = null;
      
      if (classAtt.isNumeric()) {
        preds = new double[1];
      } else {
        preds = new double[classAtt.numValues()];
      }
      
      // leaf?
      if (m_childNodes.size() == 0) {
        doLeaf(classAtt, preds);
      } else {
        // process the children
        switch (TreeModel.this.m_missingValueStrategy) {
        case NONE:
          preds = missingValueStrategyNone(instance, classAtt);
          break;
        case LASTPREDICTION:
          preds = missingValueStrategyLastPrediction(instance, classAtt);
          break;
        case DEFAULTCHILD:
          preds = missingValueStrategyDefaultChild(instance, classAtt);
          break;
        default:
          throw new Exception("[TreeModel] not implemented!");
        }
      }
      
      return preds;
    }
    
    /**
     * Compute the predictions for a leaf.
     * 
     * @param classAtt the class attribute
     * @param preds an array to hold the predicted probabilities.
     * @throws Exception if something goes wrong.
     */
    protected void doLeaf(Attribute classAtt, double[] preds) throws Exception {
      if (classAtt.isNumeric()) {
        preds[0] = m_scoreNumeric;
      } else {
        if (m_scoreDistributions.size() == 0) {
          preds[m_scoreIndex] = 1.0;
        } else {
          // collect confidences from the score distributions
          for (ScoreDistribution s : m_scoreDistributions) {
            preds[s.getClassLabelIndex()] = s.getConfidence();
          }
        }
      }
    }
    
    /**
     * Evaluate on the basis of the no true child strategy.
     * 
     * @param classAtt the class attribute.
     * @param preds an array to hold the predicted probabilities.
     * @throws Exception if something goes wrong.
     */
    protected void doNoTrueChild(Attribute classAtt, double[] preds) 
      throws Exception {
      if (TreeModel.this.m_noTrueChildStrategy == 
        NoTrueChildStrategy.RETURNNULLPREDICTION) {
        for (int i = 0; i < classAtt.numValues(); i++) {
          preds[i] = Utils.missingValue();
        }
      } else {
        // return the predictions at this node
        doLeaf(classAtt, preds);
      }
    }
    
    /**
     * Compute predictions and optionally invoke the weighted confidence
     * missing value handling strategy.
     * 
     * @param instance the incoming vector of attribute and derived field values.
     * @param classAtt the class attribute.
     * @return the predicted probability distribution.
     * @throws Exception if something goes wrong.
     */
    protected double[] missingValueStrategyWeightedConfidence(double[] instance,
        Attribute classAtt) throws Exception {
      
      if (classAtt.isNumeric()) {
        throw new Exception("[TreeNode] missing value strategy weighted confidence, "
            + "but class is numeric!");
      }
      
      double[] preds = null;
      TreeNode trueNode = null;
      boolean strategyInvoked = false;
      int nodeCount = 0;
      
      // look at the evaluation of the child predicates
      for (TreeNode c : m_childNodes) {
        if (c.getPredicate().evaluate(instance) == Predicate.Eval.TRUE) {
          // note the first child to evaluate to true
          if (trueNode == null) {
            trueNode = c;
          }
          nodeCount++;
        } else if (c.getPredicate().evaluate(instance) == Predicate.Eval.UNKNOWN) {
          strategyInvoked = true;
          nodeCount++;
        }
      }
      
      if (strategyInvoked) {
        // we expect to combine nodeCount distributions
        double[][] dists = new double[nodeCount][];
        double[] weights = new double[nodeCount];
        
        // collect the distributions and weights
        int count = 0;
        for (TreeNode c : m_childNodes) {
          if (c.getPredicate().evaluate(instance) == Predicate.Eval.TRUE ||
              c.getPredicate().evaluate(instance) == Predicate.Eval.UNKNOWN) {
            
            weights[count] = c.getRecordCount();
            if (Utils.isMissingValue(weights[count])) {
              throw new Exception("[TreeNode] weighted confidence missing value " +
              		"strategy invoked, but no record count defined for node " +
              		c.getID());
            }            
            dists[count++] = c.score(instance, classAtt);
          }
        }
        
        // do the combination
        preds = new double[classAtt.numValues()];
        for (int i = 0; i < classAtt.numValues(); i++) {
          for (int j = 0; j < nodeCount; j++) {
            preds[i] += ((weights[j] / m_recordCount) * dists[j][i]); 
          }
        }
      } else {
        if (trueNode != null) {
          preds = trueNode.score(instance, classAtt);
        } else {
          doNoTrueChild(classAtt, preds);
        }
      }
      
      return preds;
    }
    
    protected double[] freqCountsForAggNodesStrategy(double[] instance,
        Attribute classAtt) throws Exception {
    
      double[] counts = new double[classAtt.numValues()];
      
      if (m_childNodes.size() > 0) {
        // collect the counts
        for (TreeNode c : m_childNodes) {
          if (c.getPredicate().evaluate(instance) == Predicate.Eval.TRUE ||
              c.getPredicate().evaluate(instance) == Predicate.Eval.UNKNOWN) {

            double[] temp = c.freqCountsForAggNodesStrategy(instance, classAtt);
            for (int i = 0; i < classAtt.numValues(); i++) {
              counts[i] += temp[i];
            }
          }
        }
      } else {
        // process the score distributions
        if (m_scoreDistributions.size() == 0) {
          throw new Exception("[TreeModel] missing value strategy aggregate nodes:" +
          		" no score distributions at leaf " + m_ID);
        }
        for (ScoreDistribution s : m_scoreDistributions) {
          counts[s.getClassLabelIndex()] = s.getRecordCount();
        }
      }
            
      return counts;
    }
    
    /**
     * Compute predictions and optionally invoke the aggregate nodes
     * missing value handling strategy.
     * 
     * @param instance the incoming vector of attribute and derived field values.
     * @param classAtt the class attribute.
     * @return the predicted probability distribution.
     * @throws Exception if something goes wrong.
     */
    protected double[] missingValueStrategyAggregateNodes(double[] instance,
        Attribute classAtt) throws Exception {
      
      if (classAtt.isNumeric()) {
        throw new Exception("[TreeNode] missing value strategy aggregate nodes, "
            + "but class is numeric!");
      }

      double[] preds = null;
      TreeNode trueNode = null;
      boolean strategyInvoked = false;
      int nodeCount = 0;
      
      // look at the evaluation of the child predicates
      for (TreeNode c : m_childNodes) {
        if (c.getPredicate().evaluate(instance) == Predicate.Eval.TRUE) {
          // note the first child to evaluate to true
          if (trueNode == null) {
            trueNode = c;
          }
          nodeCount++;
        } else if (c.getPredicate().evaluate(instance) == Predicate.Eval.UNKNOWN) {
          strategyInvoked = true;
          nodeCount++;
        }
      }
      
      if (strategyInvoked) {
        double[] aggregatedCounts = 
          freqCountsForAggNodesStrategy(instance, classAtt);
        
        // normalize
        Utils.normalize(aggregatedCounts);
        preds = aggregatedCounts;
      } else {
        if (trueNode != null) {
          preds = trueNode.score(instance, classAtt);
        } else {
          doNoTrueChild(classAtt, preds);
        }
      }
      
      return preds;             
    }
    
    /**
     * Compute predictions and optionally invoke the default child
     * missing value handling strategy.
     * 
     * @param instance the incoming vector of attribute and derived field values.
     * @param classAtt the class attribute.
     * @return the predicted probability distribution.
     * @throws Exception if something goes wrong.
     */
    protected double[] missingValueStrategyDefaultChild(double[] instance, 
        Attribute classAtt) throws Exception {
      
      double[] preds = null;
      boolean strategyInvoked = false;
      
      // look for a child whose predicate evaluates to TRUE
      for (TreeNode c : m_childNodes) {
        if (c.getPredicate().evaluate(instance) == Predicate.Eval.TRUE) {
          preds = c.score(instance, classAtt);
          break;
        } else if (c.getPredicate().evaluate(instance) == Predicate.Eval.UNKNOWN) {
          strategyInvoked = true;
        }
      }
      
      // no true child found
      if (preds == null) {
        if (!strategyInvoked) {
          doNoTrueChild(classAtt, preds);
        } else {
          // do the strategy
          
          // NOTE: we don't actually implement the missing value penalty since
          // we always return a full probability distribution.
          if (m_defaultChild != null) {
            preds = m_defaultChild.score(instance, classAtt);
          } else {
            throw new Exception("[TreeNode] missing value strategy is defaultChild, but " +
            		"no default child has been specified in node " + m_ID);
          }
        }
      }
                  
      return preds;
    }
    
    /**
     * Compute predictions and optionally invoke the last prediction
     * missing value handling strategy.
     * 
     * @param instance the incoming vector of attribute and derived field values.
     * @param classAtt the class attribute.
     * @return the predicted probability distribution.
     * @throws Exception if something goes wrong.
     */
    protected double[] missingValueStrategyLastPrediction(double[] instance, 
        Attribute classAtt) throws Exception {
      
      double[] preds = null;
      boolean strategyInvoked = false;
      
      // look for a child whose predicate evaluates to TRUE
      for (TreeNode c : m_childNodes) {
        if (c.getPredicate().evaluate(instance) == Predicate.Eval.TRUE) {
          preds = c.score(instance, classAtt);
          break;
        } else if (c.getPredicate().evaluate(instance) == Predicate.Eval.UNKNOWN) {
          strategyInvoked = true;
        }
      }
      
      // no true child found
      if (preds == null) {
        preds = new double[classAtt.numValues()];
        if (!strategyInvoked) {
          // no true child
          doNoTrueChild(classAtt, preds);
        } else {
          // do the strategy
          doLeaf(classAtt, preds);
        }
      }
      
      return preds;
    }
    
    /**
     * Compute predictions and optionally invoke the null prediction
     * missing value handling strategy.
     * 
     * @param instance the incoming vector of attribute and derived field values.
     * @param classAtt the class attribute.
     * @return the predicted probability distribution.
     * @throws Exception if something goes wrong.
     */
    protected double[] missingValueStrategyNullPrediction(double[] instance,
        Attribute classAtt) throws Exception {
      
      double[] preds = null;
      boolean strategyInvoked = false;
      
      // look for a child whose predicate evaluates to TRUE
      for (TreeNode c : m_childNodes) {
        if (c.getPredicate().evaluate(instance) == Predicate.Eval.TRUE) {
          preds = c.score(instance, classAtt);
          break;
        } else if (c.getPredicate().evaluate(instance) == Predicate.Eval.UNKNOWN) {
          strategyInvoked = true;
        }
      }
      
      // no true child found
      if (preds == null) {
        preds = new double[classAtt.numValues()];
        if (!strategyInvoked) {
          doNoTrueChild(classAtt, preds);
        } else {
          // do the strategy
          for (int i = 0; i < classAtt.numValues(); i++) {
            preds[i] = Utils.missingValue();
          }
        }
      }
      
      return preds;
    }
    
    /**
     * Compute predictions and optionally invoke the "none"
     * missing value handling strategy (invokes no true child).
     * 
     * @param instance the incoming vector of attribute and derived field values.
     * @param classAtt the class attribute.
     * @return the predicted probability distribution.
     * @throws Exception if something goes wrong.
     */
    protected double[] missingValueStrategyNone(double[] instance, Attribute classAtt)
      throws Exception {
      
      double[] preds = null;
      
      // look for a child whose predicate evaluates to TRUE
      for (TreeNode c : m_childNodes) {
        if (c.getPredicate().evaluate(instance) == Predicate.Eval.TRUE) {
          preds = c.score(instance, classAtt);
          break;
        }
      }
      
      if (preds == null) {
        preds = new double[classAtt.numValues()];
        
        // no true child strategy
        doNoTrueChild(classAtt, preds);
      }
      
      return preds;
    }
  }
  
  /**
   * Enumerated type for the mining function
   */
  enum MiningFunction {
    CLASSIFICATION,
    REGRESSION;
  }
  
  enum MissingValueStrategy {
    LASTPREDICTION("lastPrediction"),
    NULLPREDICTION("nullPrediction"),
    DEFAULTCHILD("defaultChild"),
    WEIGHTEDCONFIDENCE("weightedConfidence"),
    AGGREGATENODES("aggregateNodes"),
    NONE("none");
    
    private final String m_stringVal;
    
    MissingValueStrategy(String name) {
      m_stringVal = name;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }
  
  enum NoTrueChildStrategy {
    RETURNNULLPREDICTION("returnNullPrediction"),
    RETURNLASTPREDICTION("returnLastPrediction");
    
    private final String m_stringVal;
    
    NoTrueChildStrategy(String name) {
      m_stringVal = name;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }
  
  enum SplitCharacteristic {
    BINARYSPLIT("binarySplit"),
    MULTISPLIT("multiSplit");
  
    private final String m_stringVal;
    
    SplitCharacteristic(String name) {
      m_stringVal = name;
    }
    
    public String toString() {
      return m_stringVal;
    }  
  }
  
  /** The mining function */
  protected MiningFunction m_functionType = MiningFunction.CLASSIFICATION;
  
  /** The missing value strategy */
  protected MissingValueStrategy m_missingValueStrategy = MissingValueStrategy.NONE;
  
  /** 
   * The missing value penalty (if defined). 
   * We don't actually make use of this since we always return 
   * full probability distributions.
   */
  protected double m_missingValuePenalty = Utils.missingValue();
  
  /** The no true child strategy to use */
  protected NoTrueChildStrategy m_noTrueChildStrategy = NoTrueChildStrategy.RETURNNULLPREDICTION;
  
  /** The splitting type */
  protected SplitCharacteristic m_splitCharacteristic = SplitCharacteristic.MULTISPLIT;
  
  /** The root of the tree */
  protected TreeNode m_root;
  
  public TreeModel(Element model, Instances dataDictionary, 
      MiningSchema miningSchema) throws Exception {
    
    super(dataDictionary, miningSchema);
    
    if (!getPMMLVersion().equals("3.2")) {
      // TODO: might have to throw an exception and only support 3.2
    }
    
    String fn = model.getAttribute("functionName");
    if (fn.equals("regression")) {
      m_functionType = MiningFunction.REGRESSION;
    }
    
    // get the missing value strategy (if any)
    String missingVS = model.getAttribute("missingValueStrategy");
    if (missingVS != null && missingVS.length() > 0) {
      for (MissingValueStrategy m : MissingValueStrategy.values()) {
        if (m.toString().equals(missingVS)) {
          m_missingValueStrategy = m;
          break;
        }
      }
    }

    // get the missing value penalty (if any)
    String missingP = model.getAttribute("missingValuePenalty");
    if (missingP != null && missingP.length() > 0) {
      // try to parse as a number
      try {
        m_missingValuePenalty = Double.parseDouble(missingP);
      } catch (NumberFormatException ex) {
        System.err.println("[TreeModel] WARNING: " +
          "couldn't parse supplied missingValuePenalty as a number");
      }
    }

    String splitC = model.getAttribute("splitCharacteristic");

    if (splitC != null && splitC.length() > 0) {
      for (SplitCharacteristic s : SplitCharacteristic.values()) {
        if (s.toString().equals(splitC)) {
          m_splitCharacteristic = s;
          break;
        }
      }
    }
    
    // find the root node of the tree
    NodeList children = model.getChildNodes();
    for (int i = 0; i < children.getLength(); i++) {
      Node child = children.item(i);
      if (child.getNodeType() == Node.ELEMENT_NODE) {
        String tagName = ((Element)child).getTagName();
        if (tagName.equals("Node")) {
          m_root = new TreeNode((Element)child, miningSchema);          
          break;
        }
      }
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
    
    preds = m_root.score(incoming, m_miningSchema.getFieldsAsInstances().classAttribute());
    
   return preds; 
  }
  
  public String toString() {
    StringBuffer temp = new StringBuffer();

    temp.append("PMML version " + getPMMLVersion());
    if (!getCreatorApplication().equals("?")) {
      temp.append("\nApplication: " + getCreatorApplication());
    }
    temp.append("\nPMML Model: TreeModel");
    temp.append("\n\n");
    temp.append(m_miningSchema);
    
    temp.append("Split-type: " + m_splitCharacteristic + "\n");
    temp.append("No true child strategy: " + m_noTrueChildStrategy + "\n");
    temp.append("Missing value strategy: " + m_missingValueStrategy + "\n");
    
    temp.append(m_root.toString());
    
    return temp.toString();
  }
  
  public String graph() throws Exception {
    StringBuffer text = new StringBuffer();
    text.append("digraph PMMTree {\n");
    
    m_root.dumpGraph(text);
    
    text.append("}\n");
    
    return text.toString();
  }

  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }

  public int graphType() {
    return Drawable.TREE;
  }
}
