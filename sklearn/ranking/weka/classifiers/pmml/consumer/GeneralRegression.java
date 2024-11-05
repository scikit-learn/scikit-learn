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
 *    GeneralRegression.java
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
 * Class implementing import of PMML General Regression model. Can be
 * used as a Weka classifier for prediction (buildClassifier()
 * raises an Exception).
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 5987 $
 */
public class GeneralRegression extends PMMLClassifier
  implements Serializable {

  /**
   * For serialization
   */
  private static final long serialVersionUID = 2583880411828388959L;

  /**
   * Enumerated type for the model type.
   */
  enum ModelType {

    // same type of model
    REGRESSION ("regression"), 
      GENERALLINEAR ("generalLinear"), 
      MULTINOMIALLOGISTIC ("multinomialLogistic"),
      ORDINALMULTINOMIAL ("ordinalMultinomial"), 
      GENERALIZEDLINEAR ("generalizedLinear");

    private final String m_stringVal;
    ModelType(String name) {
      m_stringVal = name;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }
  
  // the model type
  protected ModelType m_modelType = ModelType.REGRESSION;

  // the model name (if defined)
  protected String m_modelName;
    
  // the algorithm name (if defined)
  protected String m_algorithmName;

  // the function type (regression or classification)
  protected int m_functionType = Regression.RegressionTable.REGRESSION;

  /**
   * Enumerated type for the cumulative link function
   * (ordinal multinomial model type only).
   */
  enum CumulativeLinkFunction {
    NONE ("none") {
      double eval(double value, double offset) {
        return Double.NaN; // no evaluation defined in this case!
      }
    },
    LOGIT ("logit") {
      double eval(double value, double offset) {
        return 1.0 / (1.0 + Math.exp(-(value + offset)));
      }
    },
    PROBIT ("probit") {
      double eval(double value, double offset) {
        return weka.core.matrix.Maths.pnorm(value + offset); 
      }
    },
    CLOGLOG ("cloglog") {
      double eval(double value, double offset) {
        return 1.0 - Math.exp(-Math.exp(value + offset));
      }
    },
    LOGLOG ("loglog") {
      double eval(double value, double offset) {
        return Math.exp(-Math.exp(-(value + offset))); 
      }
    },
    CAUCHIT ("cauchit") {
      double eval(double value, double offset) {
        return 0.5 + (1.0 / Math.PI) * Math.atan(value + offset);
      }
    };

    /**
     * Evaluation function.
     * 
     * @param value the raw response value
     * @param offset the offset to add to the raw value 
     * @return the result of the link function
     */
    abstract double eval(double value, double offset);
    
    private final String m_stringVal;
    
    /**
     * Constructor
     * 
     * @param name textual name for this enum
     */
    CumulativeLinkFunction(String name) {
      m_stringVal = name;
    }
    
    /* (non-Javadoc)
     * @see java.lang.Enum#toString()
     */
    public String toString() {
      return m_stringVal;
    }
  }
  
  // cumulative link function (ordinal multinomial only)
  protected CumulativeLinkFunction m_cumulativeLinkFunction 
    = CumulativeLinkFunction.NONE;


  /**
   * Enumerated type for the link function (general linear and
   * generalized linear model types only).
   */
  enum LinkFunction {
    NONE ("none") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return Double.NaN; // no evaluation defined in this case!
      }
    },
    CLOGLOG ("cloglog") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return (1.0 - Math.exp(-Math.exp(value + offset))) * trials;
      }
    },
    IDENTITY ("identity") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return (value + offset) * trials;
      }
    },
    LOG ("log") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return Math.exp(value + offset) * trials;
      }
    },
    LOGC ("logc") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return (1.0 - Math.exp(value + offset)) * trials;
      }
    },
    LOGIT ("logit") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return (1.0 / (1.0 + Math.exp(-(value + offset)))) * trials;
      }
    },
    LOGLOG ("loglog") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return Math.exp(-Math.exp(-(value + offset))) * trials;
      }
    },
    NEGBIN ("negbin") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return (1.0 / (distParam * (Math.exp(-(value + offset)) - 1.0))) * trials;
      }
    },
    ODDSPOWER ("oddspower") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return (linkParam < 0.0 || linkParam > 0.0)
        ? (1.0 / (1.0 + Math.pow(1.0 + linkParam * (value + offset), (-1.0 / linkParam)))) * trials
        : (1.0 / (1.0 + Math.exp(-(value + offset)))) * trials;
      }
    },
    POWER ("power") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return (linkParam < 0.0 || linkParam > 0.0)
        ? Math.pow(value + offset, (1.0 / linkParam)) * trials
            : Math.exp(value + offset) * trials;
      }
    },
    PROBIT ("probit") {
      double eval(double value, double offset, double trials,
                  double distParam, double linkParam) {
        return weka.core.matrix.Maths.pnorm(value + offset) * trials;
      }
    };

    /**
     * Evaluation function.
     * 
     * @param value the raw response value
     * @param offset the offset to add to the raw value
     * @param trials the trials value to multiply the result by
     * @param distParam the distribution parameter (negbin only)
     * @param linkParam the link parameter (power and oddspower only) 
     * @return the result of the link function
     */
    abstract double eval(double value, double offset, double trials, 
                         double distParam, double linkParam);
    
    private final String m_stringVal;
    
    /**
     * Constructor.
     * 
     * @param name the textual name of this link function
     */
    LinkFunction(String name) {
      m_stringVal = name;
    }

    /* (non-Javadoc)
     * @see java.lang.Enum#toString()
     */
    public String toString() {
      return m_stringVal;
    }
  }
  
  // link function (generalLinear model type only)
  protected LinkFunction m_linkFunction = LinkFunction.NONE;
  protected double m_linkParameter = Double.NaN;
  protected String m_trialsVariable;
  protected double m_trialsValue = Double.NaN;

  /**
   * Enumerated type for the distribution (general linear
   * and generalized linear model types only).
   */
  enum Distribution {
    NONE ("none"),
    NORMAL ("normal"),
    BINOMIAL ("binomial"),
    GAMMA ("gamma"),
    INVGAUSSIAN ("igauss"),
    NEGBINOMIAL ("negbin"),
    POISSON ("poisson");

    private final String m_stringVal;
    Distribution(String name) {
      m_stringVal = name;
    }

    /* (non-Javadoc)
     * @see java.lang.Enum#toString()
     */
    public String toString() {
      return m_stringVal;
    }
  }
  
  // generalLinear and generalizedLinear model type only
  protected Distribution m_distribution = Distribution.NORMAL;

  // ancillary parameter value for the negative binomial distribution
  protected double m_distParameter = Double.NaN;

  // if present, this variable is used during scoring generalizedLinear/generalLinear or
  // ordinalMultinomial models
  protected String m_offsetVariable;

  // if present, this variable is used during scoring generalizedLinear/generalLinear or
  // ordinalMultinomial models. It works like a user-specified intercept.
  // At most, only one of offsetVariable or offsetValue may be specified.
  protected double m_offsetValue = Double.NaN;

  /**
   * Small inner class to hold the name of a parameter plus
   * its optional descriptive label
   */
  static class Parameter implements Serializable {
    // ESCA-JAVA0096:
    /** For serialization */
    // CHECK ME WITH serialver
    private static final long serialVersionUID = 6502780192411755341L;

    protected String m_name = null;
    protected String m_label = null;
  }

  // List of model parameters
  protected ArrayList<Parameter> m_parameterList = new ArrayList<Parameter>();

  /**
   * Small inner class to hold the name of a factor or covariate,
   * plus the index of the attribute it corresponds to in the
   * mining schema.
   */
  static class Predictor implements Serializable {
    /** For serialization */
    // CHECK ME WITH serialver
    private static final long serialVersionUID = 6502780192411755341L;

    protected String m_name = null;
    protected int m_miningSchemaIndex = -1;
    
    public String toString() {
      return m_name;
    }
  }
  
  // FactorList
  protected ArrayList<Predictor> m_factorList = new ArrayList<Predictor>();

  // CovariateList
  protected ArrayList<Predictor> m_covariateList = new ArrayList<Predictor>();

  /**
   * Small inner class to hold details on a predictor-to-parameter
   * correlation.
   */
  static class PPCell implements Serializable {
    /** For serialization */
    // CHECK ME WITH serialver
    private static final long serialVersionUID = 6502780192411755341L;
    
    protected String m_predictorName = null;
    protected String m_parameterName = null;

    // either the exponent of a numeric attribute or the index of
    // a discrete value
    protected double m_value = 0;

    // optional. The default is for all target categories to
    // share the same PPMatrix.
    // TO-DO: implement multiple PPMatrixes 
    protected String m_targetCategory = null;
    
  }
  
  // PPMatrix (predictor-to-parameter matrix)
  // rows = parameters, columns = predictors (attributes)
  protected PPCell[][] m_ppMatrix;

  /**
   * Small inner class to hold a single entry in the 
   * ParamMatrix (parameter matrix).
   */
  static class PCell implements Serializable {
    
    /** For serialization */
    // CHECK ME WITH serialver
    private static final long serialVersionUID = 6502780192411755341L;

    // may be null for numeric target. May also be null if this coefficent
    // applies to all target categories.
    protected String m_targetCategory = null;
    protected String m_parameterName = null;
    // coefficient
    protected double m_beta = 0.0;
    // optional degrees of freedom
    protected int m_df = -1;
  }
  
  // ParamMatrix. rows = target categories (only one if target is numeric),
  // columns = parameters (in order that they occur in the parameter list).
  protected PCell[][] m_paramMatrix;

  /**
   * Constructs a GeneralRegression classifier.
   * 
   * @param model the Element that holds the model definition
   * @param dataDictionary the data dictionary as a set of Instances
   * @param miningSchema the mining schema
   * @throws Exception if there is a problem constructing the general regression
   * object from the PMML.
   */
  public GeneralRegression(Element model, Instances dataDictionary,
                           MiningSchema miningSchema) throws Exception {

    super(dataDictionary, miningSchema);
 
    // get the model type
    String mType = model.getAttribute("modelType");
    boolean found = false;
    for (ModelType m : ModelType.values()) {
      if (m.toString().equals(mType)) {
        m_modelType = m;
        found = true;
        break;
      }      
    }
    if (!found) {
      throw new Exception("[GeneralRegression] unknown model type: " + mType);
    }

    if (m_modelType == ModelType.ORDINALMULTINOMIAL) {
      // get the cumulative link function
      String cLink = model.getAttribute("cumulativeLink");
      found = false;
      for (CumulativeLinkFunction c : CumulativeLinkFunction.values()) {
        if (c.toString().equals(cLink)) {
          m_cumulativeLinkFunction = c;
          found = true;
          break;
        }
      }
      if (!found) {
        throw new Exception("[GeneralRegression] cumulative link function " + cLink);
      }
    } else if (m_modelType == ModelType.GENERALIZEDLINEAR || 
                m_modelType == ModelType.GENERALLINEAR) {
      // get the link function
      String link = model.getAttribute("linkFunction");
      found = false;
      for (LinkFunction l : LinkFunction.values()) {
        if (l.toString().equals(link)) {
          m_linkFunction = l;
          found = true;
          break;
        }
      }
      if (!found) {
        throw new Exception("[GeneralRegression] unknown link function " + link);
      }

      // get the link parameter
      String linkP = model.getAttribute("linkParameter");
      if (linkP != null && linkP.length() > 0) {
        try {
          m_linkParameter = Double.parseDouble(linkP);
        } catch (IllegalArgumentException ex) {
          throw new Exception("[GeneralRegression] unable to parse the link parameter");
        }
      }

      // get the trials variable
      String trials = model.getAttribute("trialsVariable");
      if (trials != null && trials.length() > 0) {
        m_trialsVariable = trials;
      }

      // get the trials value
      String trialsV = model.getAttribute("trialsValue");
      if (trialsV != null && trialsV.length() > 0) {
        try {
          m_trialsValue = Double.parseDouble(trialsV);
        } catch (IllegalArgumentException ex) {
          throw new Exception("[GeneralRegression] unable to parse the trials value"); 
        }
      }
    }
  
    String mName = model.getAttribute("modelName");
    if (mName != null && mName.length() > 0) {
      m_modelName = mName;
    }

    String fName = model.getAttribute("functionName");
    if (fName.equals("classification")) {
      m_functionType = Regression.RegressionTable.CLASSIFICATION;
    }

    String algName = model.getAttribute("algorithmName");
    if (algName != null && algName.length() > 0) {
      m_algorithmName = algName;
    }

    String distribution = model.getAttribute("distribution");
    if (distribution != null && distribution.length() > 0) {
      found = false;
      for (Distribution d : Distribution.values()) {
        if (d.toString().equals(distribution)) {
          m_distribution = d;
          found = true;
          break;
        }
      }
      if (!found) {
        throw new Exception("[GeneralRegression] unknown distribution type " + distribution);
      }
    }

    String distP = model.getAttribute("distParameter");
    if (distP != null && distP.length() > 0) {
      try {
        m_distParameter = Double.parseDouble(distP);
      } catch (IllegalArgumentException ex) {
        throw new Exception("[GeneralRegression] unable to parse the distribution parameter");
      }
    }

    String offsetV = model.getAttribute("offsetVariable");
    if (offsetV != null && offsetV.length() > 0) {
       m_offsetVariable = offsetV;
    }

    String offsetVal = model.getAttribute("offsetValue");
    if (offsetVal != null && offsetVal.length() > 0) {
      try {
        m_offsetValue = Double.parseDouble(offsetVal);
      } catch (IllegalArgumentException ex) {
        throw new Exception("[GeneralRegression] unable to parse the offset value");
      }
    }

    // get the parameter list
    readParameterList(model);
    
    // get the factors and covariates
    readFactorsAndCovariates(model, "FactorList");
    readFactorsAndCovariates(model, "CovariateList");

    // read the PPMatrix
    readPPMatrix(model);

    // read the parameter estimates
    readParamMatrix(model);
  }

  /**
   * Read the list of parameters.
   *
   * @param model the Element that contains the model
   * @throws Exception if there is some problem with extracting the
   * parameters.
   */
  protected void readParameterList(Element model) throws Exception {
    NodeList paramL = model.getElementsByTagName("ParameterList");

    // should be just one parameter list
    if (paramL.getLength() == 1) {
      Node paramN = paramL.item(0);
      if (paramN.getNodeType() == Node.ELEMENT_NODE) {
        NodeList parameterList = ((Element)paramN).getElementsByTagName("Parameter");
        for (int i = 0; i < parameterList.getLength(); i++) {
          Node parameter = parameterList.item(i);
          if (parameter.getNodeType() == Node.ELEMENT_NODE) {
            Parameter p = new Parameter();
            p.m_name = ((Element)parameter).getAttribute("name");
            String label = ((Element)parameter).getAttribute("label");
            if (label != null && label.length() > 0) {
              p.m_label = label;
            }
            m_parameterList.add(p);
          }
        }
      }
    } else {
      throw new Exception("[GeneralRegression] more than one parameter list!");
    }
  }

  /**
   * Read the lists of factors and covariates.
   *
   * @param model the Element that contains the model
   * @param factorOrCovariate holds the String "FactorList" or
   * "CovariateList"
   * @throws Exception if there is a factor or covariate listed
   * that isn't in the mining schema
   */
  protected void readFactorsAndCovariates(Element model, 
                                          String factorOrCovariate) 
    throws Exception {
    Instances miningSchemaI = m_miningSchema.getFieldsAsInstances();

    NodeList factorL = model.getElementsByTagName(factorOrCovariate);
    if (factorL.getLength() == 1) { // should be 0 or 1 FactorList element
      Node factor = factorL.item(0);
      if (factor.getNodeType() == Node.ELEMENT_NODE) {
        NodeList predL = ((Element)factor).getElementsByTagName("Predictor");
        for (int i = 0; i < predL.getLength(); i++) {
          Node pred = predL.item(i);
          if (pred.getNodeType() == Node.ELEMENT_NODE) {
            Predictor p = new Predictor();
            p.m_name = ((Element)pred).getAttribute("name");
            // find the index of this predictor in the mining schema
            boolean found = false;
            for (int j = 0; j < miningSchemaI.numAttributes(); j++) {
              if (miningSchemaI.attribute(j).name().equals(p.m_name)) {
                found = true;
                p.m_miningSchemaIndex = j;
                break;
              }
            }
            if (found) {
              if (factorOrCovariate.equals("FactorList")) {
                m_factorList.add(p);
              } else {
                m_covariateList.add(p);
              }
            } else {
              throw new Exception("[GeneralRegression] reading factors and covariates - "
                                  + "unable to find predictor " +
                                  p.m_name + " in the mining schema");
            }
          }
        }
      }
    } else if (factorL.getLength() > 1){
      throw new Exception("[GeneralRegression] more than one " + factorOrCovariate
                          + "! ");
    }
  }

  /**
   * Read the PPMatrix from the xml. Does not handle multiple PPMatrixes yet.
   *
   * @param model the Element that contains the model
   * @throws Exception if there is a problem parsing cell values.
   */
  protected void readPPMatrix(Element model) throws Exception {
    Instances miningSchemaI = m_miningSchema.getFieldsAsInstances();
    
    NodeList matrixL = model.getElementsByTagName("PPMatrix");

    // should be exactly one PPMatrix
    if (matrixL.getLength() == 1) {
      // allocate space for the matrix
      // column that corresponds to the class will be empty (and will be missed out
      // when printing the model).
      m_ppMatrix = new PPCell[m_parameterList.size()][miningSchemaI.numAttributes()];

      Node ppM = matrixL.item(0);
      if (ppM.getNodeType() == Node.ELEMENT_NODE) {
        NodeList cellL = ((Element)ppM).getElementsByTagName("PPCell");
        for (int i = 0; i < cellL.getLength(); i++) {
          Node cell = cellL.item(i);
          if (cell.getNodeType() == Node.ELEMENT_NODE) {
            String predictorName = ((Element)cell).getAttribute("predictorName");
            String parameterName = ((Element)cell).getAttribute("parameterName");
            String value = ((Element)cell).getAttribute("value");
            double expOrIndex = -1;
            int predictorIndex = -1;
            int parameterIndex = -1;
            for (int j = 0; j < m_parameterList.size(); j++) {
              if (m_parameterList.get(j).m_name.equals(parameterName)) {
                parameterIndex = j;
                break;
              }
            }
            if (parameterIndex == -1) {
              throw new Exception("[GeneralRegression] unable to find parameter name "
                                  + parameterName + " in parameter list");
            }

            Predictor p = getCovariate(predictorName);
            if (p != null) {
              try {
                expOrIndex = Double.parseDouble(value);
                predictorIndex = p.m_miningSchemaIndex;
              } catch (IllegalArgumentException ex) {
                throw new Exception("[GeneralRegression] unable to parse PPCell value: "
                                    + value);
              }
            } else {
              // try as a factor
              p = getFactor(predictorName);
              if (p != null) {
                // An example pmml file from DMG seems to suggest that it
                // is possible for a continuous variable in the mining schema
                // to be treated as a factor, so we have to check for this
                if (miningSchemaI.attribute(p.m_miningSchemaIndex).isNumeric()) {
                  // parse this value as a double. It will be treated as a value
                  // to match rather than an exponent since we are dealing with
                  // a factor here
                  try {
                    expOrIndex = Double.parseDouble(value);
                  } catch (IllegalArgumentException ex) {
                    throw new Exception("[GeneralRegresion] unable to parse PPCell value: "
                                        + value);
                  }
                } else {
                  // it is a nominal attribute in the mining schema so find
                  // the index that correponds to this value
                  Attribute att = miningSchemaI.attribute(p.m_miningSchemaIndex); 
                  expOrIndex = att.indexOfValue(value);
                  if (expOrIndex == -1) {
                    throw new Exception("[GeneralRegression] unable to find PPCell value "
                                        + value + " in mining schema attribute "
                                        + att.name());
                  }
                }
              } else {
                throw new Exception("[GeneralRegression] cant find predictor "
                                    + predictorName + "in either the factors list "
                                    + "or the covariates list");
              }
              predictorIndex = p.m_miningSchemaIndex;
            }

            // fill in cell value
            PPCell ppc = new PPCell();
            ppc.m_predictorName = predictorName; ppc.m_parameterName = parameterName;
            ppc.m_value = expOrIndex;

            // TO-DO: ppc.m_targetCategory (when handling for multiple PPMatrixes is implemented)
            m_ppMatrix[parameterIndex][predictorIndex] = ppc;
          }
        }
      }
    } else {
      throw new Exception("[GeneralRegression] more than one PPMatrix!");
    }
  }

  private Predictor getCovariate(String predictorName) {
    for (int i = 0; i < m_covariateList.size(); i++) {
      if (predictorName.equals(m_covariateList.get(i).m_name)) {
        return m_covariateList.get(i);
      }
    }
    return null;
  }

  private Predictor getFactor(String predictorName) {
    for (int i = 0; i < m_factorList.size(); i++) {
      if (predictorName.equals(m_factorList.get(i).m_name)) {
        return m_factorList.get(i);
      }
    }
    return null;
  }

  /**
   * Read the parameter matrix from the xml.
   * 
   * @param model Element that holds the model
   * @throws Exception if a problem is encountered during extraction of
   * the parameter matrix
   */
  private void readParamMatrix(Element model) throws Exception {

    Instances miningSchemaI = m_miningSchema.getFieldsAsInstances();
    Attribute classAtt = miningSchemaI.classAttribute();
    // used when function type is classification but class attribute is numeric
    // in the mining schema. We will assume that there is a Target specified in
    // the pmml that defines the legal values for this class.
    ArrayList<String> targetVals = null;

    NodeList matrixL = model.getElementsByTagName("ParamMatrix");
    if (matrixL.getLength() != 1) {
      throw new Exception("[GeneralRegression] more than one ParamMatrix!");
    }
    Element matrix = (Element)matrixL.item(0);


    // check for the case where the class in the mining schema is numeric,
    // but this attribute is treated as discrete
    if (m_functionType == Regression.RegressionTable.CLASSIFICATION &&
        classAtt.isNumeric()) {
      // try and convert the class attribute to nominal. For this to succeed
      // there has to be a Target element defined in the PMML.
      if (!m_miningSchema.hasTargetMetaData()) {
        throw new Exception("[GeneralRegression] function type is classification and "
                            + "class attribute in mining schema is numeric, however, "
                            + "there is no Target element "
                            + "specifying legal discrete values for the target!");

      }

      if (m_miningSchema.getTargetMetaData().getOptype() 
          != TargetMetaInfo.Optype.CATEGORICAL) {
        throw new Exception("[GeneralRegression] function type is classification and "
                            + "class attribute in mining schema is numeric, however "
                            + "Target element in PMML does not have optype categorical!");
      }

      // OK now get legal values
      targetVals = m_miningSchema.getTargetMetaData().getValues();
      if (targetVals.size() == 0) {
        throw new Exception("[GeneralRegression] function type is classification and "
                            + "class attribute in mining schema is numeric, however "
                            + "Target element in PMML does not have any discrete values "
                            + "defined!");
      }

      // Finally, convert the class in the mining schema to nominal
      m_miningSchema.convertNumericAttToNominal(miningSchemaI.classIndex(), targetVals);
    }
    
    // allocate space for the matrix 
    m_paramMatrix = 
        new PCell[(classAtt.isNumeric())
                  ? 1
                  : classAtt.numValues()][m_parameterList.size()];

    NodeList pcellL = matrix.getElementsByTagName("PCell");
    for (int i = 0; i < pcellL.getLength(); i++) {
      // indicates that that this beta applies to all target categories
      // or target is numeric
      int targetCategoryIndex = -1;
      int parameterIndex = -1;
      Node pcell = pcellL.item(i);
      if (pcell.getNodeType() == Node.ELEMENT_NODE) {
        String paramName = ((Element)pcell).getAttribute("parameterName");
        String targetCatName = ((Element)pcell).getAttribute("targetCategory");
        String coefficient = ((Element)pcell).getAttribute("beta");
        String df = ((Element)pcell).getAttribute("df");

        for (int j = 0; j < m_parameterList.size(); j++) {
          if (m_parameterList.get(j).m_name.equals(paramName)) {
            parameterIndex = j;
            // use the label if defined
            if (m_parameterList.get(j).m_label != null) {
              paramName = m_parameterList.get(j).m_label;
            }
            break;
          }
        }
        if (parameterIndex == -1) {
          throw new Exception("[GeneralRegression] unable to find parameter name "
                              + paramName + " in parameter list");
        }

        if (targetCatName != null && targetCatName.length() > 0) {
          if (classAtt.isNominal() || classAtt.isRanking() || classAtt.isString()) {
            targetCategoryIndex = classAtt.indexOfValue(targetCatName);
          } else {
            throw new Exception("[GeneralRegression] found a PCell with a named "
                                + "target category: " + targetCatName
                                + " but class attribute is numeric in "
                                + "mining schema");
          }
        }

        PCell p = new PCell();
        if (targetCategoryIndex != -1) {
          p.m_targetCategory = targetCatName;
        }
        p.m_parameterName = paramName;
        try {
          p.m_beta = Double.parseDouble(coefficient);
        } catch (IllegalArgumentException ex) {
          throw new Exception("[GeneralRegression] unable to parse beta value "
                              + coefficient + " as a double from PCell");
        }
        if (df != null && df.length() > 0) {
          try {
            p.m_df = Integer.parseInt(df);
          } catch (IllegalArgumentException ex) {
            throw new Exception("[GeneralRegression] unable to parse df value "
                              + df + " as an int from PCell");
          }
        }
        
        if (targetCategoryIndex != -1) {
          m_paramMatrix[targetCategoryIndex][parameterIndex] = p;
        } else {
          // this PCell to all target categories (covers numeric class, in
          // which case there will be only one row in the matrix anyway)
          for (int j = 0; j < m_paramMatrix.length; j++) {
            m_paramMatrix[j][parameterIndex] = p;
          }
        }
      }
    }
  }

  /**
   * Return a textual description of this general regression.
   * 
   * @return a description of this general regression
   */
  public String toString() {
    StringBuffer temp = new StringBuffer();
    temp.append("PMML version " + getPMMLVersion());
    if (!getCreatorApplication().equals("?")) {
      temp.append("\nApplication: " + getCreatorApplication());
    }
    temp.append("\nPMML Model: " + m_modelType);
    temp.append("\n\n");
    temp.append(m_miningSchema);

    if (m_factorList.size() > 0) {
      temp.append("Factors:\n");
      for (Predictor p : m_factorList) {
        temp.append("\t" + p + "\n");
      }
    }
    temp.append("\n");
    if (m_covariateList.size() > 0) {
      temp.append("Covariates:\n");
      for (Predictor p : m_covariateList) {
        temp.append("\t" + p + "\n");
      }
    }
    temp.append("\n");
    
    printPPMatrix(temp);
    temp.append("\n");
    printParameterMatrix(temp);
    
    // do the link function stuff
    temp.append("\n");
    
    if (m_linkFunction != LinkFunction.NONE) {
      temp.append("Link function: " + m_linkFunction);
      if (m_offsetVariable != null) {
        temp.append("\n\tOffset variable " + m_offsetVariable);
      } else if (!Double.isNaN(m_offsetValue)) {
        temp.append("\n\tOffset value " + m_offsetValue);
      }
      
      if (m_trialsVariable != null) {
        temp.append("\n\tTrials variable " + m_trialsVariable);
      } else if (!Double.isNaN(m_trialsValue)) {
        temp.append("\n\tTrials value " + m_trialsValue);
      }
      
      if (m_distribution != Distribution.NONE) {
        temp.append("\nDistribution: " + m_distribution);
      }
      
      if (m_linkFunction == LinkFunction.NEGBIN &&
          m_distribution == Distribution.NEGBINOMIAL &&
          !Double.isNaN(m_distParameter)) {
        temp.append("\n\tDistribution parameter " + m_distParameter);
      }
      
      if (m_linkFunction == LinkFunction.POWER ||
          m_linkFunction == LinkFunction.ODDSPOWER) {
        if (!Double.isNaN(m_linkParameter)) {
          temp.append("\n\nLink parameter " + m_linkParameter);
        }
      }
    }
    
    if (m_cumulativeLinkFunction != CumulativeLinkFunction.NONE) {
      temp.append("Cumulative link function: " + m_cumulativeLinkFunction);
      
      if (m_offsetVariable != null) {
        temp.append("\n\tOffset variable " + m_offsetVariable);
      } else if (!Double.isNaN(m_offsetValue)) {
        temp.append("\n\tOffset value " + m_offsetValue);
      }
    }
    temp.append("\n");
    
    return temp.toString();
  }
  
  /**
   * Format and print the PPMatrix to the supplied StringBuffer.
   * 
   * @param buff the StringBuffer to append to
   */
  protected void printPPMatrix(StringBuffer buff) {
    Instances miningSchemaI = m_miningSchema.getFieldsAsInstances();
    int maxAttWidth = 0;
    for (int i = 0; i < miningSchemaI.numAttributes(); i++) {
      Attribute a = miningSchemaI.attribute(i);
      if (a.name().length() > maxAttWidth) {
        maxAttWidth = a.name().length();
      }
    }

    // check the width of the values
    for (int i = 0; i < m_parameterList.size(); i++) {
      for (int j = 0; j < miningSchemaI.numAttributes(); j++) {
        if (m_ppMatrix[i][j] != null) {
          double width = Math.log(Math.abs(m_ppMatrix[i][j].m_value)) /
            Math.log(10.0);
          if (width < 0) {
            width = 1;
          }
          // decimal + # decimal places + 1
          width += 2.0;
          if ((int)width > maxAttWidth) {
            maxAttWidth = (int)width;
          }
          if (miningSchemaI.attribute(j).isNominal() || miningSchemaI.attribute(j).isRanking()|| 
              miningSchemaI.attribute(j).isString()) {
            // check the width of this value
            String val = miningSchemaI.attribute(j).value((int)m_ppMatrix[i][j].m_value) + " ";
            if (val.length() > maxAttWidth) {
              maxAttWidth = val.length();
            }
          }
        }
      }
    }

    // get the max parameter width
    int maxParamWidth = "Parameter  ".length();
    for (Parameter p : m_parameterList) {
      String temp = (p.m_label != null)
        ? p.m_label + " "
        : p.m_name + " ";

      if (temp.length() > maxParamWidth) {
        maxParamWidth = temp.length();
      }
    }

    buff.append("Predictor-to-Parameter matrix:\n");
    buff.append(PMMLUtils.pad("Predictor", " ", (maxParamWidth + (maxAttWidth * 2 + 2))
                              - "Predictor".length(), true));
    buff.append("\n" + PMMLUtils.pad("Parameter", " ", maxParamWidth - "Parameter".length(), false));
    // attribute names
    for (int i = 0; i < miningSchemaI.numAttributes(); i++) {
      if (i != miningSchemaI.classIndex()) {
        String attName = miningSchemaI.attribute(i).name();
        buff.append(PMMLUtils.pad(attName, " ", maxAttWidth + 1 - attName.length(), true));
      }
    }
    buff.append("\n");

    for (int i = 0; i < m_parameterList.size(); i++) {
      Parameter param = m_parameterList.get(i);
      String paramS = (param.m_label != null)
        ? param.m_label
        : param.m_name;
      buff.append(PMMLUtils.pad(paramS, " ", 
                                maxParamWidth - paramS.length(), false));
      for (int j = 0; j < miningSchemaI.numAttributes(); j++) {
        if (j != miningSchemaI.classIndex()) {
          PPCell p = m_ppMatrix[i][j];
          String val = " ";
          if (p != null) {
            if (miningSchemaI.attribute(j).isNominal() || miningSchemaI.attribute(j).isRanking()||
                miningSchemaI.attribute(j).isString()) {
              val = miningSchemaI.attribute(j).value((int)p.m_value);
            } else {
              val = "" + Utils.doubleToString(p.m_value, maxAttWidth, 4).trim();
            }
          }
          buff.append(PMMLUtils.pad(val, " ", maxAttWidth + 1 - val.length(), true));
        }
      }
      buff.append("\n");
    }
  }

  /**
   * Format and print the parameter matrix to the supplied StringBuffer.
   * 
   * @param buff the StringBuffer to append to
   */
  protected void printParameterMatrix(StringBuffer buff) {
    Instances miningSchemaI = m_miningSchema.getFieldsAsInstances();

    // get the maximum class value width (nominal)
    int maxClassWidth = miningSchemaI.classAttribute().name().length();
    if (miningSchemaI.classAttribute().isNominal() || miningSchemaI.classAttribute().isRanking() 
        || miningSchemaI.classAttribute().isString()) {
      for (int i = 0; i < miningSchemaI.classAttribute().numValues(); i++) {
        if (miningSchemaI.classAttribute().value(i).length() > maxClassWidth) {
          maxClassWidth = miningSchemaI.classAttribute().value(i).length();
        }
      }
    }

    // get the maximum parameter name/label width
    int maxParamWidth = 0;
    for (int i = 0; i < m_parameterList.size(); i++) {
      Parameter p = m_parameterList.get(i);
      String val = (p.m_label != null)
        ? p.m_label + " "
        : p.m_name + " ";
      if (val.length() > maxParamWidth) {
        maxParamWidth = val.length();
      }
    }

    // get the max beta value width
    int maxBetaWidth = "Coeff.".length();
    for (int i = 0; i < m_paramMatrix.length; i++) {
      for (int j = 0; j < m_parameterList.size(); j++) {
        PCell p = m_paramMatrix[i][j];
        if (p != null) {
          double width = Math.log(Math.abs(p.m_beta)) / Math.log(10);
          if (width < 0) {
            width = 1;
          }
          // decimal + # decimal places + 1
          width += 7.0;
          if ((int)width > maxBetaWidth) {
            maxBetaWidth = (int)width;
          }
        }
      }
    }

    buff.append("Parameter estimates:\n");
    buff.append(PMMLUtils.pad(miningSchemaI.classAttribute().name(), " ", 
                              maxClassWidth + maxParamWidth + 2 - 
                              miningSchemaI.classAttribute().name().length(), false));
    buff.append(PMMLUtils.pad("Coeff.", " ", maxBetaWidth + 1 - "Coeff.".length(), true));
    buff.append(PMMLUtils.pad("df", " ", maxBetaWidth - "df".length(), true));
    buff.append("\n");
    for (int i = 0; i < m_paramMatrix.length; i++) {
      // scan for non-null entry for this class value
      boolean ok = false;
      for (int j = 0; j < m_parameterList.size(); j++) {
        if (m_paramMatrix[i][j] != null) {
          ok = true;
        }
      }
      if (!ok) {
        continue;
      }
      // first the class value (if nominal)
      String cVal = (miningSchemaI.classAttribute().isNominal() || miningSchemaI.classAttribute().isRanking() ||
          miningSchemaI.classAttribute().isString())
        ? miningSchemaI.classAttribute().value(i)
        : " ";
      buff.append(PMMLUtils.pad(cVal, " ", maxClassWidth - cVal.length(), false));     
      buff.append("\n");
      for (int j = 0; j < m_parameterList.size(); j++) {
        PCell p = m_paramMatrix[i][j];
        if (p != null) {
          String label = p.m_parameterName;
          buff.append(PMMLUtils.pad(label, " ", maxClassWidth + maxParamWidth + 2 -
                                    label.length(), true));
          String betaS = Utils.doubleToString(p.m_beta, maxBetaWidth, 4).trim();
          buff.append(PMMLUtils.pad(betaS, " ", maxBetaWidth + 1 - betaS.length(), true));
          String dfS = Utils.doubleToString(p.m_df, maxBetaWidth, 4).trim();
          buff.append(PMMLUtils.pad(dfS, " ", maxBetaWidth - dfS.length(), true));
          buff.append("\n");
        }
      }
    }
  }
  
  /**
   * Construct the incoming parameter vector based on the values
   * in the incoming test instance.
   * 
   * @param incomingInst the values of the incoming test instance
   * @return the populated parameter vector ready to be multiplied against
   * the vector of coefficients.
   * @throws Exception if there is some problem whilst constructing the
   * parameter vector
   */
  private double[] incomingParamVector(double[] incomingInst) throws Exception {
    Instances miningSchemaI = m_miningSchema.getFieldsAsInstances();
    double[] incomingPV = new double[m_parameterList.size()];
    
    for (int i = 0; i < m_parameterList.size(); i++) {
      //
      // default is that this row represents the intercept.
      // this will be the case if there are all null entries in this row
      incomingPV[i] = 1.0;

      // loop over the attributes (predictors)
      for (int j = 0; j < miningSchemaI.numAttributes(); j++) {        
        PPCell cellEntry = m_ppMatrix[i][j];
        Predictor p = null;
        if (cellEntry != null) {
          if ((p = getFactor(cellEntry.m_predictorName)) != null) {
            if ((int)incomingInst[p.m_miningSchemaIndex] == (int)cellEntry.m_value) {
              incomingPV[i] *= 1.0; // we have a match
            } else {
              incomingPV[i] *= 0.0;
            }
          } else if ((p = getCovariate(cellEntry.m_predictorName)) != null) {
              incomingPV[i] *= Math.pow(incomingInst[p.m_miningSchemaIndex], cellEntry.m_value);
          } else {
            throw new Exception("[GeneralRegression] can't find predictor "
                + cellEntry.m_predictorName + " in either the list of factors or covariates");
          }
        }
      }
    }
    
    return incomingPV;
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
    double[] incoming = m_fieldsMap.instanceToSchema(inst, m_miningSchema);
    
    // In this implementation we will default to information in the Target element (default
    // value for numeric prediction and prior probabilities for classification). If there is
    // no Target element defined, then an Exception is thrown.

    boolean hasMissing = false;
    for (int i = 0; i < incoming.length; i++) {
      if (i != m_miningSchema.getFieldsAsInstances().classIndex() && 
          Double.isNaN(incoming[i])) {
        hasMissing = true;
        break;
      }
    }
    
    if (hasMissing) {
      if (!m_miningSchema.hasTargetMetaData()) {
        String message = "[GeneralRegression] WARNING: Instance to predict has missing value(s) but "
          + "there is no missing value handling meta data and no "
          + "prior probabilities/default value to fall back to. No "
          + "prediction will be made (" 
          + ((m_miningSchema.getFieldsAsInstances().classAttribute().isNominal() || m_miningSchema.getFieldsAsInstances().classAttribute().isRanking()
              || m_miningSchema.getFieldsAsInstances().classAttribute().isString())
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
      // construct input parameter vector here
      double[] inputParamVector = incomingParamVector(incoming);
      computeResponses(incoming, inputParamVector, preds);
    }
    
    return preds;
  }
  
  /**
   * Compute the responses for the function given the parameter values corresponding
   * to the current incoming instance.
   * 
   * @param incomingInst raw incoming instance values (after missing value
   * replacement and outlier treatment)
   * @param incomingParamVector incoming instance values mapped to parameters
   * @param responses will contain the responses computed by the function
   * @throws Exception if something goes wrong
   */
  private void computeResponses(double[] incomingInst, 
                                double[] incomingParamVector,
                                double[] responses) throws Exception {
    for (int i = 0; i < responses.length; i++) {
      for (int j = 0; j < m_parameterList.size(); j++) {
        // a row of the parameter matrix should have all non-null entries
        // except for the last class (in the case of classification) which
        // should have just an intercept of 0. Need to handle the case where
        // no intercept has been defined in the pmml file for the last class
        PCell p = m_paramMatrix[i][j];
        if (p == null) {
          responses[i] += 0.0 * incomingParamVector[j];
        } else {
          responses[i] += incomingParamVector[j] * p.m_beta;
        }
      }
    }
    
    switch(m_modelType) {
    case MULTINOMIALLOGISTIC:
      computeProbabilitiesMultinomialLogistic(responses);
      break;
    case REGRESSION:
      // nothing to be done
      break;
    case GENERALLINEAR:
    case GENERALIZEDLINEAR:
      if (m_linkFunction != LinkFunction.NONE) {
        computeResponseGeneralizedLinear(incomingInst, responses);
      } else {
        throw new Exception("[GeneralRegression] no link function specified!");
      }
      break;
    case ORDINALMULTINOMIAL:
      if (m_cumulativeLinkFunction != CumulativeLinkFunction.NONE) {
        computeResponseOrdinalMultinomial(incomingInst, responses);
      } else {
        throw new Exception("[GeneralRegression] no cumulative link function specified!");
      }
      break;
      default:
        throw new Exception("[GeneralRegression] unknown model type");
    }
  }
  
  /**
   * Computes probabilities for the multinomial logistic model type.
   * 
   * @param responses will hold the responses computed by the function.
   */
  private static void computeProbabilitiesMultinomialLogistic(double[] responses) {
    double[] r = responses.clone();
    for (int j = 0; j < r.length; j++) {
      double sum = 0;
      boolean overflow = false;
      for (int k = 0; k < r.length; k++) {
        if (r[k] - r[j] > 700) {
          overflow = true;
          break;
        }
        sum += Math.exp(r[k] - r[j]);
      }
      if (overflow) {
        responses[j] = 0.0;
      } else {
        responses[j] = 1.0 / sum;
      }
    }
  }
  
  /**
   * Computes responses for the general linear and generalized linear model
   * types.
   * 
   * @param incomingInst the raw incoming instance values (after missing value
   * replacement and outlier treatment etc).
   * @param responses will hold the responses computed by the function
   * @throws Exception if a problem occurs. 
   */
  private void computeResponseGeneralizedLinear(double[] incomingInst, 
                                                double[] responses) 
    throws Exception {
    double[] r = responses.clone();
    
    double offset = 0;
    if (m_offsetVariable != null) {
      Attribute offsetAtt = 
        m_miningSchema.getFieldsAsInstances().attribute(m_offsetVariable);
      if (offsetAtt == null) {
        throw new Exception("[GeneralRegression] unable to find offset variable "
            + m_offsetVariable + " in the mining schema!");
      }
      offset = incomingInst[offsetAtt.index()];
    } else if (!Double.isNaN(m_offsetValue)) {
      offset = m_offsetValue;
    }
    
    double trials = 1;
    if (m_trialsVariable != null) {
      Attribute trialsAtt = m_miningSchema.getFieldsAsInstances().attribute(m_trialsVariable);
      if (trialsAtt == null) {
        throw new Exception("[GeneralRegression] unable to find trials variable "
            + m_trialsVariable + " in the mining schema!");
      }
      trials = incomingInst[trialsAtt.index()];
    } else if (!Double.isNaN(m_trialsValue)) {
      trials = m_trialsValue;
    }
    
    double distParam = 0;
    if (m_linkFunction == LinkFunction.NEGBIN && 
        m_distribution == Distribution.NEGBINOMIAL) {
      if (Double.isNaN(m_distParameter)) {
        throw new Exception("[GeneralRegression] no distribution parameter defined!");
      }
      distParam = m_distParameter;
    }
    
    double linkParam = 0;
    if (m_linkFunction == LinkFunction.POWER || 
        m_linkFunction == LinkFunction.ODDSPOWER) {
      if (Double.isNaN(m_linkParameter)) {
        throw new Exception("[GeneralRegression] no link parameter defined!");
      }
      linkParam = m_linkParameter;
    }
   
    for (int i = 0; i < r.length; i++) {
      responses[i] = m_linkFunction.eval(r[i], offset, trials, distParam, linkParam);
    }
  }
    
  /**
   * Computes responses for the ordinal multinomial model type.
   * 
   * @param incomingInst the raw incoming instance values (after missing value
   * replacement and outlier treatment etc).
   * @param responses will hold the responses computed by the function
   * @throws Exception if a problem occurs. 
   */
  private void computeResponseOrdinalMultinomial(double[] incomingInst, 
                                                  double[] responses) throws Exception {
    
    double[] r = responses.clone();
    
    double offset = 0;
    if (m_offsetVariable != null) {
      Attribute offsetAtt = 
        m_miningSchema.getFieldsAsInstances().attribute(m_offsetVariable);
      if (offsetAtt == null) {
        throw new Exception("[GeneralRegression] unable to find offset variable "
            + m_offsetVariable + " in the mining schema!");
      }
      offset = incomingInst[offsetAtt.index()];
    } else if (!Double.isNaN(m_offsetValue)) {
      offset = m_offsetValue;
    }
    
    for (int i = 0; i < r.length; i++) {
      if (i == 0) {
        responses[i] = m_cumulativeLinkFunction.eval(r[i], offset);
   
      } else if (i == (r.length - 1)) {
        responses[i] = 1.0 - responses[i - 1];
      } else {
        responses[i] = m_cumulativeLinkFunction.eval(r[i], offset) - responses[i - 1];
      }
    }
  }

  /* (non-Javadoc)
   * @see weka.core.RevisionHandler#getRevision()
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }
}
