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
 *    SupportVectorMachineModel.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.pmml.consumer;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.classifiers.pmml.consumer.NeuralNetwork.MiningFunction;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.pmml.MiningSchema;
import weka.core.pmml.TargetMetaInfo;
import weka.core.pmml.VectorDictionary;
import weka.core.pmml.VectorInstance;
import weka.gui.Logger;

/**
 * Implements a PMML SupportVectorMachineModel
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6528 $
 */
public class SupportVectorMachineModel extends PMMLClassifier 
  implements Serializable {
  
  /** For serialization */
  private static final long serialVersionUID = 6225095165118374296L;

  /**
   * Abstract base class for kernels
   */
  static abstract class Kernel implements Serializable {
    
    /** The log object to use */
    protected Logger m_log = null;
    
    protected Kernel(Logger log) {
      m_log = log;
    }
    
    /** For serialization */
    private static final long serialVersionUID = -6696443459968934767L;

    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector instance
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public abstract double evaluate(VectorInstance x, VectorInstance y) 
      throws Exception;
    
    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector (as an array of values)
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public abstract double evaluate(VectorInstance x, double[] y) 
      throws Exception;
    
    /**
     * Factory method returning a new Kernel constructed from 
     * the supplied XML element
     * 
     * @param svmMachineModelElement the XML element containing the kernel
     * @param log the logging object to use
     * @return a new Kernel
     * @throws Exception if something goes wrong
     */
    public static Kernel getKernel(Element svmMachineModelElement,
        Logger log) throws Exception {

      
      NodeList kList = svmMachineModelElement.getElementsByTagName("LinearKernelType");
      if (kList.getLength() > 0) {
        return new LinearKernel(log);
      }
      
      kList = svmMachineModelElement.getElementsByTagName("PolynomialKernelType");
      if (kList.getLength() > 0) {
        return new PolynomialKernel((Element)kList.item(0), log);
      }
      
      kList = svmMachineModelElement.getElementsByTagName("RadialBasisKernelType");
      if (kList.getLength() > 0) {
        return new RadialBasisKernel((Element)kList.item(0), log);
      }
      
      kList = svmMachineModelElement.getElementsByTagName("SigmoidKernelType");
      if (kList.getLength() > 0) {
        return new SigmoidKernel((Element)kList.item(0), log);
      }
      
      throw new Exception("[Kernel] Can't find a kernel that I recognize!");
    }
  }
  
  /**
   * Subclass of Kernel implementing a simple linear (dot product) kernel
   */
  static class LinearKernel extends Kernel implements Serializable {
    
    public LinearKernel(Logger log) {
      super(log);
    }
    
    public LinearKernel() {
      super(null);
    }
    
    /** For serialization */
    private static final long serialVersionUID = 8991716708484953837L;

    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector instance
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public double evaluate(VectorInstance x, VectorInstance y) 
      throws Exception {
      return x.dotProduct(y);
    }
    
    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector (as an array of values)
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public double evaluate(VectorInstance x, double[] y)
      throws Exception {
      return x.dotProduct(y);
    }
    
    /**
     * Return a textual description of this kernel
     * 
     * @return a string describing this kernel
     */
    public String toString() {
      return "Linear kernel: K(x,y) = <x,y>";
    }
  }
  
  /**
   * Subclass of Kernel implementing a polynomial kernel
   */
  static class PolynomialKernel extends Kernel implements Serializable {
    
    /** For serialization */
    private static final long serialVersionUID = -616176630397865281L;
    
    protected double m_gamma = 1;
    protected double m_coef0 = 1;
    protected double m_degree = 1;
    
    public PolynomialKernel(Element polyNode) {
      this(polyNode, null);
    }
    
    public PolynomialKernel(Element polyNode, Logger log) {
      super(log);
      
      String gammaString = polyNode.getAttribute("gamma");
      if (gammaString != null && gammaString.length() > 0) {
        try {
          m_gamma = Double.parseDouble(gammaString);
        } catch (NumberFormatException e) {
          String message = "[PolynomialKernel] : WARNING, can't parse "
            + "gamma attribute. Using default value of 1.";
          if (m_log == null) {
            System.err.println(message);
          } else {
            m_log.logMessage(message);
          }
        }
      }
      
      String coefString = polyNode.getAttribute("coef0");
      if (coefString != null && coefString.length() > 0) {
        try {
          m_coef0 = Double.parseDouble(coefString);
        } catch (NumberFormatException e) {
          String message = "[PolynomialKernel] : WARNING, can't parse "
            + "coef0 attribute. Using default value of 1.";
          if (m_log == null) {
            System.err.println(message);
          } else {
            m_log.logMessage(message);
          }
        }
      }
      
      String degreeString = polyNode.getAttribute("degree");
      if (degreeString != null && degreeString.length() > 0) {
        try {
          m_degree = Double.parseDouble(degreeString);
        } catch (NumberFormatException e) {
          String message = "[PolynomialKernel] : WARNING, can't parse "
            + "degree attribute. Using default value of 1.";
          if (m_log == null) {
            System.err.println(message);
          } else {
            m_log.logMessage(message);
          }
        }
      }
    }
    
    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector instance
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public double evaluate(VectorInstance x, VectorInstance y)
      throws Exception {
      double dotProd = x.dotProduct(y);
      return Math.pow(m_gamma * dotProd + m_coef0, m_degree);      
    }
    
    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector (as an array of values)
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public double evaluate(VectorInstance x, double[] y)
      throws Exception {
      double dotProd = x.dotProduct(y);
      return Math.pow(m_gamma * dotProd + m_coef0, m_degree);      
    }
    
    /**
     * Return a textual description of this kernel
     * 
     * @return a string describing this kernel
     */
    public String toString() {
      return "Polynomial kernel: K(x,y) = (" + m_gamma + " * <x,y> + " 
        + m_coef0 +")^" + m_degree;
    }
  }
  
  /**
   * Subclass of Kernel implementing a radial basis function kernel
   */
  static class RadialBasisKernel extends Kernel implements Serializable {
    
    /** For serialization */
    private static final long serialVersionUID = -3834238621822239042L;
    
    protected double m_gamma = 1;
    
    public RadialBasisKernel(Element radialElement) {
      this(radialElement, null);
    }
    
    public RadialBasisKernel(Element radialElement, Logger log) {
      super(log);
      
      String gammaString = radialElement.getAttribute("gamma");
      if (gammaString != null && gammaString.length() > 0) {
        try {
          m_gamma = Double.parseDouble(gammaString);
        } catch (NumberFormatException e) {
          String message = "[RadialBasisKernel] : WARNING, can't parse "
            + "gamma attribute. Using default value of 1.";
          if (m_log == null) {
            System.err.println(message);
          } else {
            m_log.logMessage(message);
          }
        }
      }
    }
   
    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector instance
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public double evaluate(VectorInstance x, VectorInstance y)
    throws Exception {
      VectorInstance diff = x.subtract(y);
      double result = -m_gamma * diff.dotProduct(diff);

      return Math.exp(result);    
    }

    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector (as an array of values)
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public double evaluate(VectorInstance x, double[] y)
      throws Exception {
      VectorInstance diff = x.subtract(y);
//      System.err.println("diff: " + diff.getValues());
      double result = -m_gamma * diff.dotProduct(diff);
//      System.err.println("Result: " + result);
      return Math.exp(result);
    }
    
    /**
     * Return a textual description of this kernel
     * 
     * @return a string describing this kernel
     */
    public String toString() {
      return "Radial kernel: K(x,y) = exp(-" + m_gamma + " * ||x - y||^2)"; 
    }    
  }
  
  /**
   * Subclass of Kernel implementing a sigmoid function
   */
  static class SigmoidKernel extends Kernel implements Serializable {
    
    /** For serialization */
    private static final long serialVersionUID = 8713475894705750117L;

    protected double m_gamma = 1;
    
    protected double m_coef0 = 1;
    
    public SigmoidKernel(Element sigElement) {
      this(sigElement, null);
    }
    
    public SigmoidKernel(Element sigElement, Logger log) {
      super(log);
      
      String gammaString = sigElement.getAttribute("gamma");
      if (gammaString != null && gammaString.length() > 0) {
        try {
          m_gamma = Double.parseDouble(gammaString);
        } catch (NumberFormatException e) {
          String message = "[SigmoidKernel] : WARNING, can't parse "
            + "gamma attribute. Using default value of 1.";
          if (m_log == null) {
            System.err.println(message);
          } else {
            m_log.logMessage(message);
          }
        }
      }
      
      String coefString = sigElement.getAttribute("coef0");
      if (coefString != null && coefString.length() > 0) {
        try {
          m_coef0 = Double.parseDouble(coefString);
        } catch (NumberFormatException e) {
          String message = "[SigmoidKernel] : WARNING, can't parse "
            + "coef0 attribute. Using default value of 1.";
          if (m_log == null) {
            System.err.println(message);
          } else {
            m_log.logMessage(message);
          }
        }
      }
    }
    
    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector instance
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public double evaluate(VectorInstance x, VectorInstance y)
      throws Exception {
      
      double dotProd = x.dotProduct(y);
      double z = m_gamma * dotProd + m_coef0;
      double a = Math.exp(z);
      double b = Math.exp(-z);
      return ((a - b) / (a + b));      
    }

    /**
     * Compute the result of the kernel evaluation on the supplied vectors
     * 
     * @param x the first vector instance
     * @param y the second vector (as an array of values)
     * @return the result of the kernel evaluation
     * @throws Exception if something goes wrong
     */
    public double evaluate(VectorInstance x, double[] y)
      throws Exception {
      
      double dotProd = x.dotProduct(y);
      double z = m_gamma * dotProd + m_coef0;
      double a = Math.exp(z);
      double b = Math.exp(-z);
      return ((a - b) / (a + b));
    }
    
    /**
     * Return a textual description of this kernel
     * 
     * @return a string describing this kernel
     */
    public String toString() {
      return "Sigmoid kernel: K(x,y) = tanh(" + m_gamma + " * <x,y> + " 
        + m_coef0 +")";
    }
  }
  
  /**
   * Inner class implementing a single binary (classification) SVM
   */
  static class SupportVectorMachine implements Serializable {
    
    /** For serialization */
    private static final long serialVersionUID = -7650496802836815608L;

    /** The target label for classification problems */
    protected String m_targetCategory;
    
    /** The index of the global alternate target label for classification */
    protected int m_globalAlternateTargetCategoryIndex = -1;
    
    /** The index of the target label for classification problems */
    protected int m_targetCategoryIndex = -1;
        
    /** PMML 4.0 - index of the alternate category for one-vs-one */
    protected int m_localAlternateTargetCategoryIndex = -1;
    
    /** PMML 4.0 - local threshold (overrides the global one, if set) */
    protected double m_localThreshold = Double.MAX_VALUE; // default - not set
        
    /** The mining schema */
    protected MiningSchema m_miningSchema;
    
    /** The log object to use (if supplied) */
    protected Logger m_log;
    
    /** 
     * True if this is a linear machine expressed in terms of the original
     * attributes
     */
    protected boolean m_coeffsOnly = false;
    
    /** The support vectors used by this machine (if not linear with coeffs only */
    protected List<VectorInstance> m_supportVectors = 
      new ArrayList<VectorInstance>();
    
    /** The constant term - b */
    protected double m_intercept = 0;
    
    /** The coefficients for the vectors */
    protected double[] m_coefficients;
    
    /**
     * Computes the prediction from this support vector machine. For classification,
     * it fills in the appropriate element in the preds array with either a 0 or a 1.
     * If the output of the machine is < 0, then a 1 is entered into the array
     * element corresponding to this machine's targetCategory, otherwise a 0 is entered.
     * Note that this is different to the scoring procedure in the 3.2 spec (see the comments
     * in the source code of this method for more information about this).
     * 
     * @param input the test instance to predict
     * @param kernel the kernel to use
     * @param vecDict the vector dictionary (null if the machine is linear and expressed
     * in terms of attribute weights)
     * @param preds the prediction array to fill in
     * @param cMethod the classification method to use (for classification problems)
     * @param globalThreshold the global threshold (used if there is no local threshold
     * for this machine)
     * @throws Exception if something goes wrong
     */
    public void distributionForInstance(double[] input, Kernel kernel, 
        VectorDictionary vecDict, double[] preds, classificationMethod cMethod,
        double globalThreshold) 
      throws Exception {
      int targetIndex = 0;
      
      if (!m_coeffsOnly) {
        // get an array that only holds the values that are referenced
        // by the support vectors
        input = vecDict.incomingInstanceToVectorFieldVals(input);
      }
      
      if (m_miningSchema.getFieldsAsInstances().classAttribute().isNominal() || m_miningSchema.getFieldsAsInstances().classAttribute().isRanking()) {
        targetIndex = m_targetCategoryIndex;
      }
      
      double result = 0;
      for (int i = 0; i < m_coefficients.length; i++) {
   //     System.err.println("X : " + m_supportVectors.get(i).getValues());
        double val = 0;
        if (!m_coeffsOnly) {
          val = kernel.evaluate(m_supportVectors.get(i), input);
        } else {
          val = input[i];
        }
        val *= m_coefficients[i];
//        System.err.println("val " + val);
        result += val;
      }
      result += m_intercept;
      
/*      if (result < 0 && m_miningSchema.getFieldsAsInstances().classAttribute().isNominal()) {
          System.err.println("[SupportVectorMachine] result (" 
            + result + ") is less than zero" +
        		" for a nominal class value!");
        result = 0;
      } else if (result > 1 && m_miningSchema.getFieldsAsInstances().classAttribute().isNominal()) {
        System.err.println("[SupportVectorMachine] result (" 
            + result + ") is greater than one" +
                        " for a nominal class value!");
        result = 1;
      }
      System.err.println("result " + result); */
      
      // TODO revisit this when I actually find out what is going on with
      // the Zementis model (the only easily available SVM model out there
      // at this time).
      
      if (cMethod == classificationMethod.NONE || 
          m_miningSchema.getFieldsAsInstances().classAttribute().isNumeric()) {
        // ----------------------------------------------------------------------
        // The PMML 3.2 spec states that the output of the machine should lie
        // between 0 and 1 for a binary classification case with 1 corresponding
        // to the machine's targetCategory and 0 corresponding to the
        // alternateBinaryTargetCategory (implying a threshold of 0.5). This seems kind of
        // non standard, and indeed, in the 4.0 spec, it has changed to output < threshold
        // corresponding to the machine's targetCategory and >= threshold corresponding
        // to the alternateBinaryTargetCategory. What has been implemented here is
        // the later with a default threshold of 0 (since the 3.2 spec doesn't have
        // a way to specify a threshold). The example SVM PMML model from Zementis,
        // which is PMML version 3.2, produces output between -1 and 1 (give or take).
        // Implementing the 3.2 scoring mechanism as described in the spec (and truncating 
        // output at 0 and 1) results in all the predicted labels getting flipped on the data
        // used to construct the Zementis model!
        //
        // April 2010 - the Zementis guys have emailed me to say that their model
        // has been prepared for PMML 4.0, even though it states it is 3.2 in the
        // XML.
        // ----------------------------------------------------------------------

        if (m_miningSchema.getFieldsAsInstances().classAttribute().isNominal() || m_miningSchema.getFieldsAsInstances().classAttribute().isRanking()) {
          if (result < 0) {
            preds[targetIndex] = 1;
          } else {
            preds[targetIndex] = 0;
          }
        } else {
          preds[targetIndex] = result;
        }
      } else {
        // PMML 4.0
        if (cMethod == classificationMethod.ONE_AGAINST_ALL) {
          // smallest value output by a machine is the predicted class
          preds[targetIndex] = result;
        } else {
          // one-vs-one
          double threshold = (m_localThreshold < Double.MAX_VALUE)
            ? m_localThreshold
            : globalThreshold;
          
          // vote
          if (result < threshold) {
            preds[targetIndex]++;
          } else {
            int altCat = (m_localAlternateTargetCategoryIndex != -1)
              ? m_localAlternateTargetCategoryIndex
              : m_globalAlternateTargetCategoryIndex;
            
            preds[altCat]++;
          }
        }
      }
//      preds[targetIndex] = result;
    }
        
    /**
     * Constructs a new SupportVectorMachine from the supplied XML element
     * 
     * @param machineElement the XML element containing the SVM
     * @param miningSchema the mining schema for the PMML model
     * @param dictionary the VectorDictionary from which to look up the support
     * vectors used by this machine (may be null if the machine is linear and
     * expressed in terms of attribute weights)
     * @param svmRep the representation of this SVM (uses support vectors or is linear
     * an uses attribute weights)
     * @param altCategoryInd the index of the global alternateBinaryTarget (if classification)
     * @param log the log object to use
     * @throws Exception if something goes wrong
     */
    public SupportVectorMachine(Element machineElement, 
        MiningSchema miningSchema, VectorDictionary dictionary, SVM_representation svmRep, 
        int altCategoryInd, Logger log) throws Exception {

      m_miningSchema = miningSchema;
      m_log = log;
      
      String targetCat = machineElement.getAttribute("targetCategory");
      if (targetCat != null && targetCat.length() > 0) {
        m_targetCategory = targetCat;
        Attribute classAtt = m_miningSchema.getFieldsAsInstances().classAttribute(); 
        if (classAtt.isNominal() || classAtt.isRanking()) {
          int index = classAtt.indexOfValue(m_targetCategory);
          
          if (index < 0) {
            throw new Exception("[SupportVectorMachine] : can't find target category: "
                + m_targetCategory + " in the class attribute!");
          }
          
          m_targetCategoryIndex = index;
          
          // now check for the PMML 4.0 alternateTargetCategory
          String altTargetCat = machineElement.getAttribute("alternateTargetCategory");
          if (altTargetCat != null && altTargetCat.length() > 0) {
            index = classAtt.indexOfValue(altTargetCat);
            if (index < 0) {
              throw new Exception("[SupportVectorMachine] : can't find alternate target category: "
                  + altTargetCat + " in the class attribute!");
            }            
            m_localAlternateTargetCategoryIndex = index;
          } else {
            // set the global one
            m_globalAlternateTargetCategoryIndex = altCategoryInd;
          }
          
        } else {
          throw new Exception("[SupportVectorMachine] : target category supplied " +
          		"but class attribute is numeric!");
        }
      } else {
        if (m_miningSchema.getFieldsAsInstances().classAttribute().isNominal() || m_miningSchema.getFieldsAsInstances().classAttribute().isRanking()) {
          m_targetCategoryIndex = (altCategoryInd == 0)
            ? 1
            : 0;
          
          m_globalAlternateTargetCategoryIndex = altCategoryInd;
          System.err.println("Setting target index for machine to " + m_targetCategoryIndex);
        }
      }
      
      if (svmRep == SVM_representation.SUPPORT_VECTORS) {
        // get the vectors
        NodeList vectorsL = 
          machineElement.getElementsByTagName("SupportVectors");
        if (vectorsL.getLength() > 0) {
          Element vectors = (Element)vectorsL.item(0);
          NodeList allTheVectorsL = 
            vectors.getElementsByTagName("SupportVector");
          for (int i = 0; i < allTheVectorsL.getLength(); i++) {
            Node vec = allTheVectorsL.item(i);
            String vecId = ((Element)vec).getAttribute("vectorId");
            VectorInstance suppV = dictionary.getVector(vecId);
            if (suppV == null) {
              throw new Exception("[SupportVectorMachine] : can't find " +
                  "vector with ID: " + vecId + " in the " +
              "vector dictionary!");
            }
            m_supportVectors.add(suppV);
          }
        }
      } else {
        m_coeffsOnly = true;
      }
      
      // get the coefficients
      NodeList coefficientsL = machineElement.getElementsByTagName("Coefficients");
      // should be just one list of coefficients
      if (coefficientsL.getLength() != 1) {
        throw new Exception("[SupportVectorMachine] Should be just one list of " +
        		"coefficients per binary SVM!");
      }
      Element cL = (Element)coefficientsL.item(0);
      String intercept = cL.getAttribute("absoluteValue");
      if (intercept != null && intercept.length() > 0) {
        m_intercept = Double.parseDouble(intercept);
      }
            
      // now get the individual coefficient elements
      NodeList coeffL = cL.getElementsByTagName("Coefficient");
      if (coeffL.getLength() == 0) {
        throw new Exception("[SupportVectorMachine] No coefficients defined!");
      }
      
      m_coefficients = new double[coeffL.getLength()];
      
      for (int i = 0; i < coeffL.getLength(); i++) {
        Element coeff = ((Element) coeffL.item(i));
        String val = coeff.getAttribute("value");
        m_coefficients[i] = Double.parseDouble(val);
      }
    }
    
    /**
     * Get a textual description of this binary SVM
     * 
     * @return a description of this SVM as a string.
     */
    public String toString() {
      StringBuffer temp = new StringBuffer();
      
      temp.append("Binary SVM");
      if (m_miningSchema.getFieldsAsInstances().classAttribute().isNominal() || m_miningSchema.getFieldsAsInstances().classAttribute().isRanking()) {
        temp.append(" (target category = " + m_targetCategory + ")");
        if (m_localAlternateTargetCategoryIndex != -1) {
          temp.append("\n (alternate category = " 
              + m_miningSchema.getFieldsAsInstances().classAttribute().
                value(m_localAlternateTargetCategoryIndex) + ")");
        }
      }
      temp.append("\n\n");
      
      for (int i = 0; i < m_supportVectors.size(); i++) {
        temp.append("\n" + m_coefficients[i] + " * [" 
            + m_supportVectors.get(i).getValues() + " * X]");
      }
      
      if (m_intercept >= 0) {
        temp.append("\n +" + m_intercept);
      } else {
        temp.append("\n " + m_intercept);
      }
      return temp.toString();
    }
  }
  
  static enum SVM_representation {
    SUPPORT_VECTORS,
    COEFFICIENTS; // for the inputs if machine is linear and expressed in terms of the attributes
  }
  
  static enum classificationMethod {
    NONE, // PMML 3.x
    ONE_AGAINST_ALL, // PMML 4.0 default
    ONE_AGAINST_ONE;
  }
  
  /** The mining function **/
  protected MiningFunction m_functionType = MiningFunction.CLASSIFICATION;
  
  /** The classification method (PMML 4.0) */
  protected classificationMethod m_classificationMethod = 
    classificationMethod.NONE; // PMML 3.x (only handles binary problems)
  
  /** The model name (if defined) */
  protected String m_modelName;
  
  /** The algorithm name (if defined) */
  protected String m_algorithmName;
  
  /** The dictionary of support vectors */
  protected VectorDictionary m_vectorDictionary;
  
  /** The kernel function to use */
  protected Kernel m_kernel;
  
  /** The individual binary SVMs */
  protected List<SupportVectorMachine> m_machines = 
    new ArrayList<SupportVectorMachine>();
  
  /** The other class index (in the case of a single binary SVM - PMML 3.2). */
  protected int m_alternateBinaryTargetCategory = -1;
  
  /** Do we have support vectors, or just attribute coefficients for a linear machine? */
  protected SVM_representation m_svmRepresentation = SVM_representation.SUPPORT_VECTORS;
  
  /** PMML 4.0 threshold value */
  protected double m_threshold = 0;
  
  /**
   * Construct a new SupportVectorMachineModel encapsulating the information provided
   * in the PMML document.
   * 
   * @param model the SVM element from the PMML document
   * @param dataDictionary the data dictionary
   * @param miningSchema the mining schema
   * @throws Exception if the model can't be constructed from the PMML
   */
  public SupportVectorMachineModel(Element model, Instances dataDictionary,
      MiningSchema miningSchema) throws Exception {
    
    super(dataDictionary, miningSchema);
    
    if (!getPMMLVersion().equals("3.2")) {
      //TODO might have to throw an exception and only support 3.2
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
    
    String svmRep = model.getAttribute("svmRepresentation");
    if (svmRep != null && svmRep.length() > 0) {
      if (svmRep.equals("Coefficients")) {
        m_svmRepresentation = SVM_representation.COEFFICIENTS;
      }
    }
    
    String altTargetCat = model.getAttribute("alternateBinaryTargetCategory");
    if (altTargetCat != null && altTargetCat.length() > 0) {
      int altTargetInd = 
        m_miningSchema.getFieldsAsInstances().classAttribute().indexOfValue(altTargetCat);
      
      if (altTargetInd < 0) {
        throw new Exception("[SupportVectorMachineModel] can't find alternate " +
        		"target value " + altTargetCat);
      }
      m_alternateBinaryTargetCategory = altTargetInd;
    }
    
    // PMML 4.0
    String thresholdS = model.getAttribute("threshold");
    if (thresholdS != null && thresholdS.length() > 0) {
      m_threshold = Double.parseDouble(thresholdS);
    }
    
    // PMML 4.0
    if (getPMMLVersion().startsWith("4.")) {
      m_classificationMethod = classificationMethod.ONE_AGAINST_ALL; // default for PMML 4.0
    }
    
    String classificationMethodS = model.getAttribute("classificationMethod");
    if (classificationMethodS != null && classificationMethodS.length()> 0) {
      if (classificationMethodS.equals("OneAgainstOne")) {
        m_classificationMethod = classificationMethod.ONE_AGAINST_ONE;
      }
    }
    
    if (m_svmRepresentation == SVM_representation.SUPPORT_VECTORS) {
      m_vectorDictionary = VectorDictionary.getVectorDictionary(model, miningSchema);
    }
    
    m_kernel = Kernel.getKernel(model, m_log);
    if (m_svmRepresentation == SVM_representation.COEFFICIENTS && 
        !(m_kernel instanceof LinearKernel)) {
      throw new Exception("[SupportVectorMachineModel] representation is " +
      		"coefficients, but kernel is not linear!");
    }
    
    // Get the individual machines
    NodeList machineL = model.getElementsByTagName("SupportVectorMachine");
    if (machineL.getLength() == 0) {
      throw new Exception("[SupportVectorMachineModel] No binary SVMs" +
      		" defined in model file!");
    }
    for (int i = 0; i < machineL.getLength(); i++) {
      Node machine = machineL.item(i);
      SupportVectorMachine newMach = 
        new SupportVectorMachine((Element)machine, 
            m_miningSchema, m_vectorDictionary, m_svmRepresentation,
            m_alternateBinaryTargetCategory, m_log);
      m_machines.add(newMach);
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
      for (int i = 0; i < preds.length; i++) {
        preds[i] = -1; // mark all entries as not calculated to begin with
      }
    }
    
    double[] incoming = m_fieldsMap.instanceToSchema(inst, m_miningSchema);
    
    boolean hasMissing = false;
    for (int i = 0; i < incoming.length; i++) {
      if (i != m_miningSchema.getFieldsAsInstances().classIndex() && 
          Double.isNaN(incoming[i])) {
        hasMissing = true;
        //System.err.println("Missing value for att : " + m_miningSchema.getFieldsAsInstances().attribute(i).name());
        break;
      }
    }
    
    if (hasMissing) {
      if (!m_miningSchema.hasTargetMetaData()) {
        String message = "[SupportVectorMachineModel] WARNING: Instance to predict has missing value(s) but "
          + "there is no missing value handling meta data and no "
          + "prior probabilities/default value to fall back to. No "
          + "prediction will be made (" 
          + ((m_miningSchema.getFieldsAsInstances().classAttribute().isNominal()
        	  ||m_miningSchema.getFieldsAsInstances().classAttribute().isRanking()
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
      for (SupportVectorMachine m : m_machines) {
        m.distributionForInstance(incoming, m_kernel, m_vectorDictionary, 
            preds, m_classificationMethod, m_threshold);
      }
    }
    
    if (m_classificationMethod != classificationMethod.NONE &&
        (m_miningSchema.getFieldsAsInstances().classAttribute().isNominal()||m_miningSchema.getFieldsAsInstances().classAttribute().isRanking())) {
      // PMML 4.0
      if (m_classificationMethod == classificationMethod.ONE_AGAINST_ALL) {
        // find the minimum value
        int minI = Utils.minIndex(preds);
        preds = new double[preds.length];
        preds[minI] = 1.0;
      } else {
        // nothing to do for one-against-one - just normalize the
        // votes
      }
    }
    
    if (m_machines.size() == preds.length - 1) {
      double total = 0;
      int unset = -1;
      for (int i = 0; i < preds.length; i++) {
        if (preds[i] != -1) {
          total += preds[i];
        } else {
          unset = i;
        }
      }
      
      if (total > 1.0) {
        throw new Exception("[SupportVectorMachineModel] total of probabilities" +
        		" is greater than 1!");
      }
      
      preds[unset] = 1.0 - total;
    }
    
    if (preds.length > 1) {
      Utils.normalize(preds);
    }
    
    return preds;
  }
  
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6528 $");
  }
  
  /**
   * Get a textual description of this SupportVectorMachineModel
   * 
   * @return a description of this SupportVectorMachineModel as a string
   */
  public String toString() {
    StringBuffer temp = new StringBuffer();
    
    temp.append("PMML version " + getPMMLVersion());
    if (!getCreatorApplication().equals("?")) {
      temp.append("\nApplication: " + getCreatorApplication());
    }
    temp.append("\nPMML Model: Support Vector Machine Model");
    
    temp.append("\n\n");
    temp.append(m_miningSchema);
    
    temp.append("Kernel: \n\t");
    temp.append(m_kernel);
    temp.append("\n");
    
    if (m_classificationMethod != classificationMethod.NONE) {
      temp.append("Multi-class classifcation using ");
      if (m_classificationMethod == classificationMethod.ONE_AGAINST_ALL) {
        temp.append("one-against-all");
      } else {
        temp.append("one-against-one");
      }
      temp.append("\n\n");
    }
    
    for (SupportVectorMachine v : m_machines) {
      temp.append("\n" + v);
    }
    
    return temp.toString();
  }
}
