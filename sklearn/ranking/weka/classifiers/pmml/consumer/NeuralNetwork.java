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
 *    NeuralNetwork.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.pmml.consumer;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

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
 * Class implementing import of PMML Neural Network model. Can be used as a Weka
 * classifier for prediction (buildClassifier() raises an Exception).
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision 1.0 $
 */
public class NeuralNetwork extends PMMLClassifier {
  
  /**
   * For serialization
   */
  private static final long serialVersionUID = -4545904813133921249L;

  /**
   * Small inner class for a NeuralInput (essentially just
   * wraps a DerivedField and adds an ID)
   */
  static class NeuralInput implements Serializable {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = -1902233762824835563L;
    
    /** Field that this input refers to */
    private DerivedFieldMetaInfo m_field;
    
    /** ID string */
    private String m_ID = null;
    
    private String getID() {
      return m_ID;
    }
    
    protected NeuralInput(Element input, MiningSchema miningSchema) throws Exception {
      m_ID = input.getAttribute("id");
      
      NodeList fL = input.getElementsByTagName("DerivedField");
      if (fL.getLength() != 1) {
        throw new Exception("[NeuralInput] expecting just one derived field!");
      }
      
      Element dF = (Element)fL.item(0);
      Instances allFields = miningSchema.getFieldsAsInstances();
      ArrayList<Attribute> fieldDefs = new ArrayList<Attribute>();
      for (int i = 0; i < allFields.numAttributes(); i++) {
        fieldDefs.add(allFields.attribute(i));
      }
      m_field = new DerivedFieldMetaInfo(dF, fieldDefs, miningSchema.getTransformationDictionary());
    }
    
    protected double getValue(double[] incoming) throws Exception {
      return m_field.getDerivedValue(incoming);
    }
    
    public String toString() {
      StringBuffer temp = new StringBuffer();
      
      temp.append("Nueral input (" + getID() + ")\n");
      temp.append(m_field);
      
      return temp.toString();
    }
  }
  
  /**
   * Inner class representing a layer in the network.
   */
  class NeuralLayer implements Serializable {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = -8386042001675763922L;

    /** The number of neurons in this layer */
    private int m_numNeurons = 0;
    
    /** Activation function (if defined, overrides one in NeuralNetwork) */
    private ActivationFunction m_layerActivationFunction = null;
    
    /** Threshold (if defined overrides one in NeuralNetwork) */
    private double m_layerThreshold = Double.NaN; 
    
    /** Width (if defined overrides one in NeuralNetwork) */
    private double m_layerWidth = Double.NaN;
    
    /** Altitude (if defined overrides one in NeuralNetwork) */
    private double m_layerAltitude = Double.NaN;
    
    /** Normalization (if defined overrides one in NeuralNetwork) */
    private Normalization m_layerNormalization = null;
    
    /** The neurons at this hidden layer */
    private Neuron[] m_layerNeurons = null;
    
    /** Stores the output of this layer (for given inputs) */
    private HashMap<String, Double> m_layerOutput = new HashMap<String, Double>();
    
    protected NeuralLayer(Element layerE) {
      
      String activationFunction = layerE.getAttribute("activationFunction");
      if (activationFunction != null && activationFunction.length() > 0) {
        for (ActivationFunction a : ActivationFunction.values()) {
          if (a.toString().equals(activationFunction)) {
            m_layerActivationFunction = a;
            break;
          }
        }
      } else {
        // use the network-level activation function
        m_layerActivationFunction = m_activationFunction;
      }
      
      String threshold = layerE.getAttribute("threshold");
      if (threshold != null && threshold.length() > 0) {
        m_layerThreshold = Double.parseDouble(threshold);
      } else {
        // use network-level threshold
        m_layerThreshold = m_threshold;
      }
      
      String width = layerE.getAttribute("width");
      if (width != null && width.length() > 0) {
        m_layerWidth = Double.parseDouble(width);
      } else {
        // use network-level width
        m_layerWidth = m_width;
      }
      
      String altitude = layerE.getAttribute("altitude");
      if (altitude != null && altitude.length() > 0) {
        m_layerAltitude = Double.parseDouble(altitude);
      } else {
        // use network-level altitude
        m_layerAltitude = m_altitude;
      }
      
      String normMethod = layerE.getAttribute("normalizationMethod");
      if (normMethod != null && normMethod.length() > 0) {
        for (Normalization n : Normalization.values()) {
          if (n.toString().equals(normMethod)) {
            m_layerNormalization = n;
            break;
          }
        }
      } else {
        // use network-level normalization method
        m_layerNormalization = m_normalizationMethod;
      }
      
      NodeList neuronL = layerE.getElementsByTagName("Neuron");
      m_numNeurons = neuronL.getLength();
      m_layerNeurons = new Neuron[m_numNeurons];
      for (int i = 0; i < neuronL.getLength(); i++) {
        Node neuronN = neuronL.item(i);
        if (neuronN.getNodeType() == Node.ELEMENT_NODE) {
          m_layerNeurons[i] = new Neuron((Element)neuronN, this);
        }
      }
    }
    
    protected ActivationFunction getActivationFunction() {
      return m_layerActivationFunction;
    }
    
    protected double getThreshold() {
      return m_layerThreshold;
    }
    
    protected double getWidth() {
      return m_layerWidth;
    }
    
    protected double getAltitude() {
      return m_layerAltitude;
    }
    
    protected Normalization getNormalization() {
      return m_layerNormalization;
    }
    
    /**
     * Compute the output values for this layer.
     * 
     * @param incoming the incoming values
     * @return the output values for this layer
     * @throws Exception if there is a problem computing the outputs
     */
    protected HashMap<String, Double> computeOutput(HashMap<String, Double> incoming) 
      throws Exception {
      
      m_layerOutput.clear();
      
      double normSum = 0;
      for (int i = 0; i < m_layerNeurons.length; i++) {
        double neuronOut = m_layerNeurons[i].getValue(incoming);
        String neuronID = m_layerNeurons[i].getID();

        if (m_layerNormalization == Normalization.SOFTMAX) {
          normSum += Math.exp(neuronOut);
        } else if (m_layerNormalization == Normalization.SIMPLEMAX) {
          normSum += neuronOut;
        }
        //System.err.println("Inserting ID " + neuronID + " " + neuronOut);
        m_layerOutput.put(neuronID, neuronOut);
      }
      
      // apply the normalization (if necessary)
      if (m_layerNormalization != Normalization.NONE) {
        for (int i = 0; i < m_layerNeurons.length; i++) {
          double val = m_layerOutput.get(m_layerNeurons[i].getID());
//          System.err.println("Normalizing ID " + m_layerNeurons[i].getID() + " " + val);
          if (m_layerNormalization == Normalization.SOFTMAX) {
            val = Math.exp(val) / normSum;
          } else {
            val = (val / normSum);
          }
          m_layerOutput.put(m_layerNeurons[i].getID(), val);
        }
      }
      return m_layerOutput;
    }
    
    public String toString() {
      StringBuffer temp = new StringBuffer();
      
      temp.append("activation: " + getActivationFunction() + "\n");
      if (!Double.isNaN(getThreshold())) {
        temp.append("threshold: " + getThreshold() + "\n");
      }
      if (!Double.isNaN(getWidth())) {
        temp.append("width: " + getWidth() + "\n");
      }
      if (!Double.isNaN(getAltitude())) {
        temp.append("altitude: " + getAltitude() + "\n");
      }
      temp.append("normalization: " + m_layerNormalization + "\n");
      for (int i = 0; i < m_numNeurons; i++) {
        temp.append(m_layerNeurons[i] + "\n");
      }

      return temp.toString();
    }
  }
  
  /**
   * Inner class encapsulating a Neuron
   */
  static class Neuron implements Serializable {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = -3817434025682603443L;

    /** ID string */
    private String m_ID = null;
    
    /** The layer we belong to (for accessing activation function, threshold etc.) */
    private NeuralLayer m_layer;
    
    /** The bias */
    private double m_bias = 0.0;
    
    /** The width (if defined overrides the one in NeuralLayer or NeuralNetwork) */
    private double m_neuronWidth = Double.NaN;
    
    /** The altitude (if defined overrides the one in NeuralLayer or NeuralNetwork) */
    private double m_neuronAltitude = Double.NaN;
    
    /** The IDs of the neurons/neural inputs that we are connected to */
    private String[] m_connectionIDs = null;
    
    /** The weights corresponding to the connections */
    private double[] m_weights = null;
    
    protected Neuron(Element neuronE, NeuralLayer layer) {
      m_layer = layer;
      
      m_ID = neuronE.getAttribute("id");
      
      String bias = neuronE.getAttribute("bias");
      if (bias != null && bias.length() > 0) {
        m_bias = Double.parseDouble(bias);
      }
      
      String width = neuronE.getAttribute("width");
      if (width != null && width.length() > 0) {
        m_neuronWidth = Double.parseDouble(width);
      }
      
      String altitude = neuronE.getAttribute("altitude");
      if (altitude != null && altitude.length() > 0) {
        m_neuronAltitude = Double.parseDouble(altitude);
      }
      
      // get the connection details
      NodeList conL = neuronE.getElementsByTagName("Con");
      m_connectionIDs = new String[conL.getLength()];
      m_weights = new double[conL.getLength()];
      for (int i = 0; i < conL.getLength(); i++) {
        Node conN = conL.item(i);
        if (conN.getNodeType() == Node.ELEMENT_NODE) {
          Element conE = (Element)conN;
          m_connectionIDs[i] = conE.getAttribute("from");
          String weight = conE.getAttribute("weight");
          m_weights[i] = Double.parseDouble(weight);
        }
      }
    }
    
    protected String getID() {
      return m_ID;
    }    
    
    /**
     * Compute the output of this Neuron.
     * 
     * @param incoming a Map of input values. The keys are the IDs
     * of incoming connections (either neural inputs or neurons) and
     * the values are the output values of the neural input/neuron in
     * question.
     * 
     * @return the output of this neuron
     * @throws Exception if any of our incoming connection IDs cannot be
     * located in the Map
     */
    protected double getValue(HashMap<String, Double> incoming) throws Exception {
      
      double z = 0;
      double result = Double.NaN;
      
      double width = (Double.isNaN(m_neuronWidth))
        ? m_layer.getWidth()
        : m_neuronWidth;

      z = m_bias;
      for (int i = 0; i < m_connectionIDs.length; i++) {
        Double inVal = incoming.get(m_connectionIDs[i]);
        if (inVal == null) {
          throw new Exception("[Neuron] unable to find connection " 
              + m_connectionIDs[i] + " in input Map!");
        }

        if (m_layer.getActivationFunction() != ActivationFunction.RADIALBASIS) {
          // multiply with weight
          double inV = inVal.doubleValue() * m_weights[i];
          z += inV;
        } else {
          // Euclidean distance to the center (stored in m_weights)
          double inV = Math.pow((inVal.doubleValue() - m_weights[i]), 2.0);
          z += inV;
        }
      }
      
      // apply the width if necessary
      if (m_layer.getActivationFunction() == ActivationFunction.RADIALBASIS) {
        z /= (2.0 * (width * width));
      }

      double threshold = m_layer.getThreshold();
      double altitude = (Double.isNaN(m_neuronAltitude))
        ? m_layer.getAltitude()
        : m_neuronAltitude;
        
      double fanIn = m_connectionIDs.length;        
      result = m_layer.getActivationFunction().eval(z, threshold, altitude, fanIn);
      
      return result;
    }
    
    public String toString() {
      StringBuffer temp = new StringBuffer();
      temp.append("Nueron (" + m_ID + ") [bias:" + m_bias);
      if (!Double.isNaN(m_neuronWidth)) {
        temp.append(" width:" + m_neuronWidth);
      }
      if (!Double.isNaN(m_neuronAltitude)) {
        temp.append(" altitude:" + m_neuronAltitude);
      }
      temp.append("]\n");
      temp.append("  con. (ID:weight): ");
      for (int i = 0; i < m_connectionIDs.length; i++) {
        temp.append(m_connectionIDs[i] + ":" + Utils.doubleToString(m_weights[i], 2));
        if ((i + 1) % 10 == 0 || i == m_connectionIDs.length - 1) {
          temp.append("\n                    ");
        } else {
          temp.append(", ");
        }
      }
      return temp.toString();
    }
  }
  
  static class NeuralOutputs implements Serializable {
    
    /**
     * For serialization
     */
    private static final long serialVersionUID = -233611113950482952L;

    /** The neurons we are mapping */
    private String[] m_outputNeurons = null;
    
    /**
     *  In the case of a nominal class, the index of the value
     * being predicted by each output neuron
     */
    private int[] m_categoricalIndexes = null;
    
    /** The class attribute we are mapping to */
    private Attribute m_classAttribute = null;
    
    /** Used when the class is numeric */
    private NormContinuous m_regressionMapping = null;
        
    protected NeuralOutputs(Element outputs, MiningSchema miningSchema) throws Exception {
      m_classAttribute = miningSchema.getMiningSchemaAsInstances().classAttribute();
      
      int vals = (m_classAttribute.isNumeric())
        ? 1
        : m_classAttribute.numValues();
      
      m_outputNeurons = new String[vals];
      m_categoricalIndexes = new int[vals];
      
      NodeList outputL = outputs.getElementsByTagName("NeuralOutput");
      if (outputL.getLength() != m_outputNeurons.length) {
        throw new Exception("[NeuralOutputs] the number of neural outputs does not match "
            + "the number expected!");
      }
      
      for (int i = 0; i < outputL.getLength(); i++) {
        Node outputN = outputL.item(i);
        if (outputN.getNodeType() == Node.ELEMENT_NODE) {
          Element outputE = (Element)outputN;
          // get the ID for this output neuron
          m_outputNeurons[i] = outputE.getAttribute("outputNeuron");
          
          if (m_classAttribute.isNumeric()) {
            // get the single norm continuous
            NodeList contL = outputE.getElementsByTagName("NormContinuous");
            if (contL.getLength() != 1) {
              throw new Exception("[NeuralOutputs] Should be exactly one norm continuous element "
                  + "for numeric class!");
            }
            Node normContNode = contL.item(0);
            String attName = ((Element)normContNode).getAttribute("field");
            Attribute dummyTargetDef = new Attribute(attName);
            ArrayList<Attribute> dummyFieldDefs = new ArrayList<Attribute>();
            dummyFieldDefs.add(dummyTargetDef);
            
            m_regressionMapping = new NormContinuous((Element)normContNode, 
                FieldMetaInfo.Optype.CONTINUOUS, dummyFieldDefs);
            break;
          } else {
            // we just need to grab the categorical value (out of the NormDiscrete element)
            // that this output neuron is associated with
            NodeList discL = outputE.getElementsByTagName("NormDiscrete");
            if (discL.getLength() != 1) {
              throw new Exception("[NeuralOutputs] Should be only one norm discrete element "
                  + "per derived field/neural output for a nominal class!");
            }
            Node normDiscNode = discL.item(0);
            String attValue = ((Element)normDiscNode).getAttribute("value");
            int index = m_classAttribute.indexOfValue(attValue);
            if (index < 0) {
              throw new Exception("[NeuralOutputs] Can't find specified target value "
                  + attValue + " in class attribute " + m_classAttribute.name());
            }
            m_categoricalIndexes[i] = index;
          }
        }
      }
    }
    
    /**
     * Compute the output. Either a probability distribution or a single
     * value (regression).
     * 
     * @param incoming the values from the last hidden layer
     * @param preds the array to fill with predicted values
     * @throws Exception if there is a problem computing the output
     */
    protected void getOuput(HashMap<String, Double> incoming, double[] preds) throws Exception {
      
      if (preds.length != m_outputNeurons.length) {
        throw new Exception("[NeuralOutputs] Incorrect number of predictions requested: "
            + preds.length + "requested, " + m_outputNeurons.length + " expected");
      }
      for (int i = 0; i < m_outputNeurons.length; i++) {
        Double neuronOut = incoming.get(m_outputNeurons[i]);
        if (neuronOut == null) {
          throw new Exception("[NeuralOutputs] Unable to find output neuron "
              + m_outputNeurons[i] + " in the incoming HashMap!!");
        }
        if (m_classAttribute.isNumeric()) {
          // will be only one output neuron anyway
          preds[0] = neuronOut.doubleValue();
          
          preds[0] = m_regressionMapping.getResultInverse(preds);
        } else {

          // clip at zero
          // preds[m_categoricalIndexes[i]] = (neuronOut < 0) ? 0.0 : neuronOut;
          preds[m_categoricalIndexes[i]] = neuronOut;
        }
      }
      
      if (m_classAttribute.isNominal() || m_classAttribute.isRanking()) {
        // check for negative values and adjust
        double min = preds[Utils.minIndex(preds)];
        if (min < 0) {
          for (int i = 0; i < preds.length; i++) {
            preds[i] -= min;
          }
        }
        // do a simplemax normalization
        Utils.normalize(preds);
      }
    }
    
    public String toString() {
      StringBuffer temp = new StringBuffer();
      
      for (int i = 0; i < m_outputNeurons.length; i++) {
        temp.append("Output neuron (" + m_outputNeurons[i] + ")\n");
        temp.append("mapping:\n");
        if (m_classAttribute.isNumeric()) {
          temp.append(m_regressionMapping +"\n");
        } else {
          temp.append(m_classAttribute.name() + " = " 
              + m_classAttribute.value(m_categoricalIndexes[i]) + "\n");
        }
      }
      
      return temp.toString();
    }
  }
  
  /**
   * Enumerated type for the mining function
   */
  enum MiningFunction {
    CLASSIFICATION,
    REGRESSION;
  }
  
  /** The mining function */
  protected MiningFunction m_functionType = MiningFunction.CLASSIFICATION;
  
  /**
   * Enumerated type for the activation function.
   */
  enum ActivationFunction {
    THRESHOLD("threshold") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        if (z > threshold) {
          return 1.0;
        }
        return 0.0;
      }
    },
    LOGISTIC("logistic") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return 1.0 / (1.0 + Math.exp(-z));
      }
    },
    TANH("tanh") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        double a = Math.exp( z );
        double b = Math.exp( -z );
        return ((a-b)/(a+b));
        //return (1.0 - Math.exp(-2.0 * z)) / (1.0 + Math.exp(-2.0 * z));
      }
    },
    IDENTITY("identity") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return z;
      }
    },
    EXPONENTIAL("exponential") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return Math.exp(z);
      }
    },
    RECIPROCAL("reciprocal") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return 1.0 / z;
      }
    },
    SQUARE("square") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return  z * z;
      }
    },
    GAUSS("gauss") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return Math.exp(-(z * z));
      }
    },
    SINE("sine") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return Math.sin(z);
      }
    },
    COSINE("cosine") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return Math.cos(z);
      }
    },
    ELLICOT("ellicot") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return z / (1.0 + Math.abs(z));
      }
    },
    ARCTAN("arctan") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return 2.0 * Math.atan(z) / Math.PI;
      }
    },
    RADIALBASIS("radialBasis") {
      double eval(double z, double threshold, double altitude, double fanIn) {
        return Math.exp(fanIn * Math.log(altitude) - z);
      }
    };
    
    abstract double eval(double z, double threshold, double altitude, double fanIn);
    
    private final String m_stringVal;
    
    ActivationFunction(String name) {
      m_stringVal = name;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }
  
  /** The activation function to use */
  protected ActivationFunction m_activationFunction = ActivationFunction.ARCTAN;
  
  /**
   * Enumerated type for the normalization method
   */
  enum Normalization {
    NONE ("none"),
    SIMPLEMAX ("simplemax"),
    SOFTMAX ("softmax");
    
    private final String m_stringVal;
    
    Normalization(String name) {
      m_stringVal = name;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }
    
  /** The normalization method */
  protected Normalization m_normalizationMethod = Normalization.NONE;
  
  /** Threshold activation */
  protected double m_threshold = 0.0; // default = 0
  
  /** Width for radial basis */
  protected double m_width = Double.NaN; // no default
  
  /** Altitude for radial basis */
  protected double m_altitude = 1.0; // default = 1
  
  /** The number of inputs to the network */
  protected int m_numberOfInputs = 0;
  
  /** Number of hidden layers in the network */
  protected int m_numberOfLayers = 0;
  
  /** The inputs to the network */
  protected NeuralInput[] m_inputs = null;
  
  /** A map for storing network input values (computed from an incoming instance) */
  protected HashMap<String, Double> m_inputMap = new HashMap<String, Double>();
    
  /** The hidden layers in the network */
  protected NeuralLayer[] m_layers = null;
  
  /** The outputs of the network */
  protected NeuralOutputs m_outputs = null;
  
  public NeuralNetwork(Element model, Instances dataDictionary,
                       MiningSchema miningSchema) throws Exception {
    
    super(dataDictionary, miningSchema);
    
    String fn = model.getAttribute("functionName");
    if (fn.equals("regression")) {
      m_functionType = MiningFunction.REGRESSION;
    }
    
    String act = model.getAttribute("activationFunction");
    if (act == null || act.length() == 0) {
      throw new Exception("[NeuralNetwork] no activation functon defined");
    }
    
    // get the activation function
    for (ActivationFunction a : ActivationFunction.values()) {
      if (a.toString().equals(act)) {
        m_activationFunction = a;
        break;
      }
    }
    
    // get the normalization method (if specified)
    String norm = model.getAttribute("normalizationMethod");
    if (norm != null && norm.length() > 0) {
      for (Normalization n : Normalization.values()) {
        if (n.toString().equals(norm)) {
          m_normalizationMethod = n;
          break;
        }
      }
    }
    
    String thresh = model.getAttribute("threshold");
    if (thresh != null && thresh.length() > 0) {
      m_threshold = Double.parseDouble(thresh);
    }
    String width = model.getAttribute("width");
    if (width != null && width.length() > 0) {
      m_width = Double.parseDouble(width);
    }
    String alt = model.getAttribute("altitude");
    if (alt != null && alt.length() > 0) {
      m_altitude = Double.parseDouble(alt);
    }
    
    // get all the inputs
    NodeList inputL = model.getElementsByTagName("NeuralInput");
    m_numberOfInputs = inputL.getLength();
    m_inputs = new NeuralInput[m_numberOfInputs];
    for (int i = 0; i < m_numberOfInputs; i++) {
      Node inputN = inputL.item(i);
      if (inputN.getNodeType() == Node.ELEMENT_NODE) {
        NeuralInput nI = new NeuralInput((Element)inputN, m_miningSchema);
        m_inputs[i] = nI;
      }
    }
    
    // get the layers
    NodeList layerL = model.getElementsByTagName("NeuralLayer");
    m_numberOfLayers = layerL.getLength();
    m_layers = new NeuralLayer[m_numberOfLayers];
    for (int i = 0; i < m_numberOfLayers; i++) {
      Node layerN = layerL.item(i);
      if (layerN.getNodeType() == Node.ELEMENT_NODE) {
        NeuralLayer nL = new NeuralLayer((Element)layerN);
        m_layers[i] = nL;
      }
    }
    
    // get the outputs
    NodeList outputL = model.getElementsByTagName("NeuralOutputs");
    if (outputL.getLength() != 1) {
      throw new Exception("[NeuralNetwork] Should be just one NeuralOutputs element defined!");
    }
    
    m_outputs = new NeuralOutputs((Element)outputL.item(0), m_miningSchema);
  }

  /* (non-Javadoc)
   * @see weka.core.RevisionHandler#getRevision()
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
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
        String message = "[NeuralNetwork] WARNING: Instance to predict has missing value(s) but "
          + "there is no missing value handling meta data and no "
          + "prior probabilities/default value to fall back to. No "
          + "prediction will be made (" 
          + ((m_miningSchema.getFieldsAsInstances().classAttribute().isNominal() || m_miningSchema.getFieldsAsInstances().classAttribute().isNominal()
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
      
      // construct the input to the network for this instance
      m_inputMap.clear();
      for (int i = 0; i < m_inputs.length; i++) {
        double networkInVal = m_inputs[i].getValue(incoming);
        String ID = m_inputs[i].getID();
        m_inputMap.put(ID, networkInVal);
      }
      
      // now compute the output of each layer
      HashMap<String, Double> layerOut = m_layers[0].computeOutput(m_inputMap);
      for (int i = 1; i < m_layers.length; i++) {
        layerOut = m_layers[i].computeOutput(layerOut);
      }
      
      // now do the output
      m_outputs.getOuput(layerOut, preds);
    }
    
    return preds;
  }

  public String toString() {
    StringBuffer temp = new StringBuffer();
    
    temp.append("PMML version " + getPMMLVersion());
    if (!getCreatorApplication().equals("?")) {
      temp.append("\nApplication: " + getCreatorApplication());
    }
    temp.append("\nPMML Model: Neural network");
    temp.append("\n\n");
    temp.append(m_miningSchema);
    
    temp.append("Inputs:\n");
    for (int i = 0; i < m_inputs.length; i++) {
      temp.append(m_inputs[i] + "\n");
    }

    for (int i = 0; i < m_layers.length; i++) {
      temp.append("Layer: " + (i+1) + "\n");
      temp.append(m_layers[i] + "\n");
    }
    
    temp.append("Outputs:\n");
    temp.append(m_outputs);
    
    return temp.toString();
  }
}
