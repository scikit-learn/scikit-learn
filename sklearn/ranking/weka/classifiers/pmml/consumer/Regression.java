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
 *    Regression.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.pmml.consumer;

import java.io.Serializable;
import java.util.ArrayList;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.pmml.*;

/**
 * Class implementing import of PMML Regression model. Can be
 * used as a Weka classifier for prediction (buildClassifier()
 * raises an Exception).
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com
 * @version $Revision: 6018 $
 */
public class Regression extends PMMLClassifier
  implements Serializable {

  /** For serialization */
  private static final long serialVersionUID = -5551125528409488634L;

  /**
   * Inner class for encapsulating a regression table
   */
  static class RegressionTable implements Serializable {

    /** For serialization */
    private static final long serialVersionUID = -5259866093996338995L;

    /**
     * Abstract inner base class for different predictor types.
     */
    abstract static class Predictor implements Serializable {

      /** For serialization */
      private static final long serialVersionUID = 7043831847273383618L;
      
      /** Name of this predictor */
      protected String m_name;
      
      /** 
       * Index of the attribute in the mining schema that corresponds to this
       * predictor
       */
      protected int m_miningSchemaAttIndex = -1;
      
      /** Coefficient for this predictor */
      protected double m_coefficient = 1.0;
      
      /**
       * Constructs a new Predictor.
       * 
       * @param predictor the <code>Element</code> encapsulating this predictor
       * @param miningSchema the mining schema as an Instances object
       * @throws Exception if there is a problem constructing this Predictor
       */
      protected Predictor(Element predictor, Instances miningSchema) throws Exception {
        m_name = predictor.getAttribute("name");
        for (int i = 0; i < miningSchema.numAttributes(); i++) {
          Attribute temp = miningSchema.attribute(i);
          if (temp.name().equals(m_name)) {
            m_miningSchemaAttIndex = i;
          }
        }
        
        if (m_miningSchemaAttIndex == -1) {
          throw new Exception("[Predictor] unable to find matching attribute for "
                              + "predictor " + m_name);
        }

        String coeff = predictor.getAttribute("coefficient");
        if (coeff.length() > 0) {
          m_coefficient = Double.parseDouble(coeff);
        }
      }

      /**
       * Returns a textual description of this predictor applicable
       * to all sub classes.
       */
      public String toString() {
        return Utils.doubleToString(m_coefficient, 12, 4) + " * ";
      }

      /**
       * Abstract add method. Adds this predictor into the sum for the
       * current prediction.
       * 
       * @param preds the prediction computed so far. For regression, it is a
       * single element array; for classification it is a multi-element array
       * @param input the input instance's values
       */
      public abstract void add(double[] preds, double[] input);
    }

    /**
     * Inner class for a numeric predictor
     */
    protected class NumericPredictor extends Predictor {
      /**
       * For serialization
       */
      private static final long serialVersionUID = -4335075205696648273L;
      
      /** The exponent*/
      protected double m_exponent = 1.0;

      /**
       * Constructs a NumericPredictor.
       * 
       * @param predictor the <code>Element</code> holding the predictor
       * @param miningSchema the mining schema as an Instances object
       * @throws Exception if something goes wrong while constructing this
       * predictor
       */
      protected NumericPredictor(Element predictor, 
                              Instances miningSchema) throws Exception {
        super(predictor, miningSchema);
        
        String exponent = predictor.getAttribute("exponent");
        if (exponent.length() > 0) {
          m_exponent = Double.parseDouble(exponent);
        }
      }

      /**
       * Return a textual description of this predictor.
       */
      public String toString() {
        String output = super.toString();
        output += m_name;
        if (m_exponent > 1.0 || m_exponent < 1.0) {
          output += "^" + Utils.doubleToString(m_exponent, 4);
        }
        return output;
      }

      /**
       * Adds this predictor into the sum for the
       * current prediction.
       * 
       * @param preds the prediction computed so far. For regression, it is a
       * single element array; for classification it is a multi-element array
       * @param input the input instance's values
       */
      public void add(double[] preds, double[] input) {
        if (m_targetCategory == -1) {
          preds[0] += m_coefficient * Math.pow(input[m_miningSchemaAttIndex], m_exponent);
        } else {
          preds[m_targetCategory] += 
            m_coefficient * Math.pow(input[m_miningSchemaAttIndex], m_exponent);
        }
      }
    }

    /**
     * Inner class encapsulating a categorical predictor.
     */
    protected class CategoricalPredictor extends Predictor {
      
      /**For serialization */
      private static final long serialVersionUID = 3077920125549906819L;
      
      /** The attribute value for this predictor */
      protected String m_valueName;
      
      /** The index of the attribute value for this predictor */
      protected int m_valueIndex = -1;

      /**
       * Constructs a CategoricalPredictor.
       * 
       * @param predictor the <code>Element</code> containing the predictor
       * @param miningSchema the mining schema as an Instances object
       * @throws Exception if something goes wrong while constructing
       * this predictor
       */
      protected CategoricalPredictor(Element predictor,
                                  Instances miningSchema) throws Exception {
        super(predictor, miningSchema);
        
        String valName = predictor.getAttribute("value");
        if (valName.length() == 0) {
          throw new Exception("[CategoricalPredictor] attribute value not specified!");
        }
        
        m_valueName = valName;

        Attribute att = miningSchema.attribute(m_miningSchemaAttIndex);
        if (att.isString()) {
          // means that there were no Value elements defined in the 
          // data dictionary (and hence the mining schema).
          // We add our value here.
          att.addStringValue(m_valueName);
        }
        m_valueIndex = att.indexOfValue(m_valueName);
        /*        for (int i = 0; i < att.numValues(); i++) {
          if (att.value(i).equals(m_valueName)) {
            m_valueIndex = i;
          }
          }*/

        if (m_valueIndex == -1) {
          throw new Exception("[CategoricalPredictor] unable to find value "
                              + m_valueName + " in mining schema attribute "
                              + att.name());
        }
      }

      /**
       * Return a textual description of this predictor.
       */
      public String toString() {
        String output = super.toString();
        output += m_name + "=" + m_valueName;
        return output;
      }

      /**
       * Adds this predictor into the sum for the
       * current prediction.
       * 
       * @param preds the prediction computed so far. For regression, it is a
       * single element array; for classification it is a multi-element array
       * @param input the input instance's values
       */
      public void add(double[] preds, double[] input) {
        
        // if the value is equal to the one in the input then add the coefficient
        if (m_valueIndex == (int)input[m_miningSchemaAttIndex]) {
          if (m_targetCategory == -1) {
            preds[0] += m_coefficient;
          } else {
            preds[m_targetCategory] += m_coefficient;
          }
        }
      }
    }

    /**
     * Inner class to handle PredictorTerms.
     */
    protected class PredictorTerm implements Serializable {

      /** For serialization */
      private static final long serialVersionUID = 5493100145890252757L;

      /** The coefficient for this predictor term */
      protected double m_coefficient = 1.0;

      /** the indexes of the terms to be multiplied */
      protected int[] m_indexes;

      /** The names of the terms (attributes) to be multiplied */
      protected String[] m_fieldNames;

      /**
       * Construct a new PredictorTerm.
       * 
       * @param predictorTerm the <code>Element</code> describing the predictor term
       * @param miningSchema the mining schema as an Instances object
       * @throws Exception if something goes wrong while constructing this
       * predictor term
       */
      protected PredictorTerm(Element predictorTerm, 
                              Instances miningSchema) throws Exception {

        String coeff = predictorTerm.getAttribute("coefficient");
        if (coeff != null && coeff.length() > 0) {
          try {
            m_coefficient = Double.parseDouble(coeff);
          } catch (IllegalArgumentException ex) {
            throw new Exception("[PredictorTerm] unable to parse coefficient");
          }
        }
        
        NodeList fields = predictorTerm.getElementsByTagName("FieldRef");
        if (fields.getLength() > 0) {
          m_indexes = new int[fields.getLength()];
          m_fieldNames = new String[fields.getLength()];

          for (int i = 0; i < fields.getLength(); i++) {
            Node fieldRef = fields.item(i);
            if (fieldRef.getNodeType() == Node.ELEMENT_NODE) {
              String fieldName = ((Element)fieldRef).getAttribute("field");
              if (fieldName != null && fieldName.length() > 0) {
                boolean found = false;
                // look for this field in the mining schema
                for (int j = 0; j < miningSchema.numAttributes(); j++) {
                  if (miningSchema.attribute(j).name().equals(fieldName)) {
                    
                    // all referenced fields MUST be numeric
                    if (!miningSchema.attribute(j).isNumeric()) {
                      throw new Exception("[PredictorTerm] field is not continuous: "
                                          + fieldName);
                    }
                    found = true;
                    m_indexes[i] = j;
                    m_fieldNames[i] = fieldName;
                    break;
                  }
                }
                if (!found) {
                  throw new Exception("[PredictorTerm] Unable to find field "
                                      + fieldName + " in mining schema!");
                }
              }
            }
          }
        }
      }

      /**
       * Return a textual description of this predictor term.
       */
      public String toString() {
        StringBuffer result = new StringBuffer();
        result.append("(" + Utils.doubleToString(m_coefficient, 12, 4));
        for (int i = 0; i < m_fieldNames.length; i++) {
          result.append(" * " + m_fieldNames[i]);
        }
        result.append(")");
        return result.toString();
      }

      /**
       * Adds this predictor term into the sum for the
       * current prediction.
       * 
       * @param preds the prediction computed so far. For regression, it is a
       * single element array; for classification it is a multi-element array
       * @param input the input instance's values
       */
      public void add(double[] preds, double[] input) {
        int indx = 0;
        if (m_targetCategory != -1) {
          indx = m_targetCategory;
        }

        double result = m_coefficient;
        for (int i = 0; i < m_indexes.length; i++) {
          result *= input[m_indexes[i]];
        }
        preds[indx] += result;
      }
    }
    
    /** Constant for regression model type */
    public static final int REGRESSION = 0;
    
    /** Constant for classification model type */
    public static final int CLASSIFICATION = 1;

    /** The type of function - regression or classification */
    protected int m_functionType = REGRESSION;
    
    /** The mining schema */
    protected MiningSchema m_miningSchema;
        
    /** The intercept */
    protected double m_intercept = 0.0;
    
    /** classification only */
    protected int m_targetCategory = -1;

    /** Numeric and categorical predictors */
    protected ArrayList<Predictor> m_predictors = 
      new ArrayList<Predictor>();

    /** Interaction terms */
    protected ArrayList<PredictorTerm> m_predictorTerms =
      new ArrayList<PredictorTerm>();

    /**
     * Return a textual description of this RegressionTable.
     */
    public String toString() {
      Instances miningSchema = m_miningSchema.getFieldsAsInstances();
      StringBuffer temp = new StringBuffer();
      temp.append("Regression table:\n");
      temp.append(miningSchema.classAttribute().name());
      if (m_functionType == CLASSIFICATION) {
        temp.append("=" + miningSchema.
                    classAttribute().value(m_targetCategory));
      }

      temp.append(" =\n\n");
      
      // do the predictors
      for (int i = 0; i < m_predictors.size(); i++) {
        temp.append(m_predictors.get(i).toString() + " +\n");
      }
      
      // do the predictor terms
      for (int i = 0; i < m_predictorTerms.size(); i++) {
        temp.append(m_predictorTerms.get(i).toString() + " +\n");
      }

      temp.append(Utils.doubleToString(m_intercept, 12, 4));
      temp.append("\n\n");

      return temp.toString();
    }

    /**
     * Construct a regression table from an <code>Element</code>
     *
     * @param table the table to encapsulate
     * @param functionType the type of function 
     * (regression or classification)
     * to use
     * @param mSchema the mining schema
     * @throws Exception if there is a problem while constructing
     * this regression table
     */
    protected RegressionTable(Element table, 
                           int functionType,
                           MiningSchema mSchema) throws Exception {

      m_miningSchema = mSchema;
      m_functionType = functionType;

      Instances miningSchema = m_miningSchema.getFieldsAsInstances();

      // get the intercept
      String intercept = table.getAttribute("intercept");
      if (intercept.length() > 0) {
        m_intercept = Double.parseDouble(intercept);
      }

      // get the target category (if classification)
      if (m_functionType == CLASSIFICATION) {
        // target category MUST be defined
        String targetCat = table.getAttribute("targetCategory");
        if (targetCat.length() > 0) {
          Attribute classA = miningSchema.classAttribute();
          for (int i = 0; i < classA.numValues(); i++) {
            if (classA.value(i).equals(targetCat)) {
              m_targetCategory = i;
            }
          }
        } 
        if (m_targetCategory == -1) {
          throw new Exception("[RegressionTable] No target categories defined for classification");
        }
      }

      // read all the numeric predictors
      NodeList numericPs = table.getElementsByTagName("NumericPredictor");
      for (int i = 0; i < numericPs.getLength(); i++) {
        Node nP = numericPs.item(i);
        if (nP.getNodeType() == Node.ELEMENT_NODE) {
          NumericPredictor numP = new NumericPredictor((Element)nP, miningSchema);
          m_predictors.add(numP);
        }
      }

      // read all the categorical predictors
      NodeList categoricalPs = table.getElementsByTagName("CategoricalPredictor");
      for (int i = 0; i < categoricalPs.getLength(); i++) {
        Node cP = categoricalPs.item(i);
        if (cP.getNodeType() == Node.ELEMENT_NODE) {
          CategoricalPredictor catP = new CategoricalPredictor((Element)cP, miningSchema);
          m_predictors.add(catP);
        }
      }

      // read all the PredictorTerms
      NodeList predictorTerms = table.getElementsByTagName("PredictorTerm");
      for (int i = 0; i < predictorTerms.getLength(); i++) {
        Node pT = predictorTerms.item(i);
        PredictorTerm predT = new PredictorTerm((Element)pT, miningSchema);
        m_predictorTerms.add(predT);
      }
    }

    public void predict(double[] preds, double[] input) {
      if (m_targetCategory == -1) {
        preds[0] = m_intercept;
      } else {
        preds[m_targetCategory] = m_intercept;
      }
      
      // add the predictors
      for (int i = 0; i < m_predictors.size(); i++) {
        Predictor p = m_predictors.get(i);
        p.add(preds, input);
      }

      // add the PredictorTerms
      for (int i = 0; i < m_predictorTerms.size(); i++) {
        PredictorTerm pt = m_predictorTerms.get(i);
        pt.add(preds, input);
      }
    }
  }

  /** Description of the algorithm */
  protected String m_algorithmName;

  /** The regression tables for this regression */
  protected RegressionTable[] m_regressionTables;

  /**
   * Enum for the normalization methods.
   */
  enum Normalization {
    NONE, SIMPLEMAX, SOFTMAX, LOGIT, PROBIT, CLOGLOG,
      EXP, LOGLOG, CAUCHIT}

  /** The normalization to use */
  protected Normalization m_normalizationMethod = Normalization.NONE;

  /**
   * Constructs a new PMML Regression.
   * 
   * @param model the <code>Element</code> containing the regression model
   * @param dataDictionary the data dictionary as an Instances object
   * @param miningSchema the mining schema
   * @throws Exception if there is a problem constructing this Regression
   */
  public Regression(Element model, Instances dataDictionary,
                    MiningSchema miningSchema) throws Exception {
    super(dataDictionary, miningSchema);
    
    int functionType = RegressionTable.REGRESSION;

    // determine function name first
    String fName = model.getAttribute("functionName");
    
    if (fName.equals("regression")) {
      functionType = RegressionTable.REGRESSION;
    } else if (fName.equals("classification")) {
      functionType = RegressionTable.CLASSIFICATION;
    } else {
      throw new Exception("[PMML Regression] Function name not defined in pmml!");
    }

    // do we have an algorithm name?
    String algName = model.getAttribute("algorithmName");
    if (algName != null && algName.length() > 0) {
      m_algorithmName = algName;
    }

    // determine normalization method (if any)
    m_normalizationMethod = determineNormalization(model);

    setUpRegressionTables(model, functionType);

    // convert any string attributes in the mining schema
    //miningSchema.convertStringAttsToNominal();
  }

  /**
   * Create all the RegressionTables for this model.
   * 
   * @param model the <code>Element</code> holding this regression model
   * @param functionType the type of function (regression or
   * classification)
   * @throws Exception if there is a problem setting up the regression
   * tables
   */
  private void setUpRegressionTables(Element model,
                                     int functionType) throws Exception {
    NodeList tableList = model.getElementsByTagName("RegressionTable");
    
    if (tableList.getLength() == 0) {
      throw new Exception("[Regression] no regression tables defined!");
    }

    m_regressionTables = new RegressionTable[tableList.getLength()];
    
    for (int i = 0; i < tableList.getLength(); i++) {
      Node table = tableList.item(i);
      if (table.getNodeType() == Node.ELEMENT_NODE) {
        RegressionTable tempRTable = 
          new RegressionTable((Element)table, 
                              functionType, 
                              m_miningSchema);
        m_regressionTables[i] = tempRTable;
      }
    }
  }

  /**
   * Return the type of normalization used for this regression
   * 
   * @param model the <code>Element</code> holding the model
   * @return the normalization used in this regression
   */
  private static Normalization determineNormalization(Element model) {
    
    Normalization normMethod = Normalization.NONE;

    String normName = model.getAttribute("normalizationMethod");
    if (normName.equals("simplemax")) {
      normMethod = Normalization.SIMPLEMAX;
    } else if (normName.equals("softmax")) {
      normMethod = Normalization.SOFTMAX;
    } else if (normName.equals("logit")) {
      normMethod = Normalization.LOGIT;
    } else if (normName.equals("probit")) {
      normMethod = Normalization.PROBIT;
    } else if (normName.equals("cloglog")) {
      normMethod = Normalization.CLOGLOG;
    } else if (normName.equals("exp")) {
      normMethod = Normalization.EXP;
    } else if (normName.equals("loglog")) {
      normMethod = Normalization.LOGLOG;
    } else if (normName.equals("cauchit")) {
      normMethod = Normalization.CAUCHIT;
    } 
    return normMethod;
  }

  /**
   * Return a textual description of this Regression model.
   */
  public String toString() {
    StringBuffer temp = new StringBuffer();
    temp.append("PMML version " + getPMMLVersion());
    if (!getCreatorApplication().equals("?")) {
      temp.append("\nApplication: " + getCreatorApplication());
    }
    if (m_algorithmName != null) {
      temp.append("\nPMML Model: " + m_algorithmName);
    }
    temp.append("\n\n");
    temp.append(m_miningSchema);

    for (RegressionTable table : m_regressionTables) {
      temp.append(table);
    }
    
    if (m_normalizationMethod != Normalization.NONE) {
      temp.append("Normalization: " + m_normalizationMethod);
    }
    temp.append("\n");

    return temp.toString();
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

    // create an array of doubles that holds values from the incoming
    // instance; in order of the fields in the mining schema. We will
    // also handle missing values and outliers here.
    //    System.err.println(inst);
    double[] incoming = m_fieldsMap.instanceToSchema(inst, m_miningSchema);

    // scan for missing values. If there are still missing values after instanceToSchema(),
    // then missing value handling has been deferred to the PMML scheme. The specification
    // (Regression PMML 3.2) seems to contradict itself with regards to classification and categorical
    // variables. In one place it states that if a categorical variable is missing then 
    // variable_name=value is 0 for any value. Further down in the document it states: "if
    // one or more of the y_j cannot be evaluated because the value in one of the referenced
    // fields is missing, then the following formulas (for computing p_j) do not apply. In
    // that case the predictions are defined by the priorProbability values in the Target
    // element".

    // In this implementation we will default to information in the Target element (default
    // value for numeric prediction and prior probabilities for classification). If there is
    // no Target element defined, then an Exception is thrown.

    boolean hasMissing = false;
    for (int i = 0; i < incoming.length; i++) {
      if (i != m_miningSchema.getFieldsAsInstances().classIndex() && 
          Utils.isMissingValue(incoming[i])) {
        hasMissing = true;
        break;
      }
    }

    if (hasMissing) {
      if (!m_miningSchema.hasTargetMetaData()) {
        String message = "[Regression] WARNING: Instance to predict has missing value(s) but "
          + "there is no missing value handling meta data and no "
          + "prior probabilities/default value to fall back to. No "
          + "prediction will be made (" 
          + ((m_miningSchema.getFieldsAsInstances().classAttribute().isNominal() ||
        	  m_miningSchema.getFieldsAsInstances().classAttribute().isRanking()||
              m_miningSchema.getFieldsAsInstances().classAttribute().isString())
              ? "zero probabilities output)."
              : "NaN output).");
        if (m_log == null) {
          System.err.println(message);
        } else {
          m_log.logMessage(message);
        }
        if (m_miningSchema.getFieldsAsInstances().classAttribute().isNumeric()) {
          preds[0] = Utils.missingValue();
        }
        return preds;
      } else {
        // use prior probablilities/default value
        TargetMetaInfo targetData = m_miningSchema.getTargetMetaData();
        if (m_miningSchema.getFieldsAsInstances().classAttribute().isNumeric()) {
          preds[0] = targetData.getDefaultValue();
        } else {
          Instances miningSchemaI = m_miningSchema.getFieldsAsInstances();
          for (int i = 0; i < miningSchemaI.classAttribute().numValues(); i++) {
            preds[i] = targetData.getPriorProbability(miningSchemaI.classAttribute().value(i));
          }
        }
        return preds;
      }
    } else {
      // loop through the RegressionTables
      for (int i = 0; i < m_regressionTables.length; i++) {
        m_regressionTables[i].predict(preds, incoming);
      }
 
      // Now apply the normalization
      switch (m_normalizationMethod) {
      case NONE:
        // nothing to be done
        break;
      case SIMPLEMAX:
        Utils.normalize(preds);
        break;
      case SOFTMAX:
        for (int i = 0; i < preds.length; i++) {
          preds[i] = Math.exp(preds[i]);
        }
        if (preds.length == 1) {
          // hack for those models that do binary logistic regression as
          // a numeric prediction model
          preds[0] = preds[0] / (preds[0] + 1.0);
        } else {
          Utils.normalize(preds);
        }
        break;
      case LOGIT:
        for (int i = 0; i < preds.length; i++) {
          preds[i] = 1.0 / (1.0 + Math.exp(-preds[i]));
        }
        Utils.normalize(preds);
        break;
      case PROBIT:
        for (int i = 0; i < preds.length; i++) {
          preds[i] = weka.core.matrix.Maths.pnorm(preds[i]);
        }
        Utils.normalize(preds);
        break;
      case CLOGLOG:
        // note this is supposed to be illegal for regression
        for (int i = 0; i < preds.length; i++) {
          preds[i] = 1.0 - Math.exp(-Math.exp(-preds[i]));
        }
        Utils.normalize(preds);
        break;
      case EXP:
        for (int i = 0; i < preds.length; i++) {
          preds[i] = Math.exp(preds[i]);
        }
        Utils.normalize(preds);
        break;
      case LOGLOG:
        // note this is supposed to be illegal for regression
        for (int i = 0; i < preds.length; i++) {
          preds[i] = Math.exp(-Math.exp(-preds[i]));
        }
        Utils.normalize(preds);
        break;
      case CAUCHIT:
        for (int i = 0; i < preds.length; i++) {
          preds[i] = 0.5 + (1.0 / Math.PI) * Math.atan(preds[i]);
        }
        Utils.normalize(preds);
        break;
      default:
          throw new Exception("[Regression] unknown normalization method");
      }

      // If there is a Target defined, and this is a numeric prediction problem,
      // then apply any min, max, rescaling etc.
      if (m_miningSchema.getFieldsAsInstances().classAttribute().isNumeric()
          && m_miningSchema.hasTargetMetaData()) {
        TargetMetaInfo targetData = m_miningSchema.getTargetMetaData();
        preds[0] = targetData.applyMinMaxRescaleCast(preds[0]);
      }
    }
    
    return preds;
  }

  /* (non-Javadoc)
   * @see weka.core.RevisionHandler#getRevision()
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6018 $");
  }
}
