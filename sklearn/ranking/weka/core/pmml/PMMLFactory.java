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
 *    PMMLFactory.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.classifiers.AbstractClassifier;
import weka.classifiers.pmml.consumer.GeneralRegression;
import weka.classifiers.pmml.consumer.NeuralNetwork;
import weka.classifiers.pmml.consumer.PMMLClassifier;
import weka.classifiers.pmml.consumer.Regression;
import weka.classifiers.pmml.consumer.RuleSetModel;
import weka.classifiers.pmml.consumer.SupportVectorMachineModel;
import weka.classifiers.pmml.consumer.TreeModel;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;
import weka.gui.Logger;

/**
 * This class is a factory class for reading/writing PMML models
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6467 $
 */
public class PMMLFactory {

  /** for serialization */
  
  protected enum ModelType {
    UNKNOWN_MODEL ("unknown"),
    REGRESSION_MODEL ("Regression"),
    GENERAL_REGRESSION_MODEL ("GeneralRegression"),
    NEURAL_NETWORK_MODEL ("NeuralNetwork"),
    TREE_MODEL ("TreeModel"),
    RULESET_MODEL("RuleSetModel"),
    SVM_MODEL ("SupportVectorMachineModel");
    
    private final String m_stringVal;
    
    ModelType(String name) {
      m_stringVal = name;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }

  /**
   * Read and return a PMML model.
   *
   * @param filename the name of the file to read from
   * @return a PMML model
   * @throws Exception if there is a problem while reading the file
   */
  public static PMMLModel getPMMLModel(String filename) throws Exception {
    return getPMMLModel(filename, null);
  }
  
  /**
   * Read and return a PMML model.
   *
   * @param file a <code>File</code> to read from
   * @return a PMML model
   * @throws Exception if there is a problem while reading the file
   */
  public static PMMLModel getPMMLModel(File file) throws Exception {
    return getPMMLModel(file, null);
  }
  
  /**
   * Read and return a PMML model.
   *
   * @param stream the <code>InputStream</code> to read from
   * @return a PMML model
   * @throws Exception if there is a problem while reading from the stream
   */
  public static PMMLModel getPMMLModel(InputStream stream) throws Exception {
    return getPMMLModel(stream, null);
  }
  
  /**
   * Read and return a PMML model.
   *
   * @param filename the name of the file to read from
   * @param log the logging object to use (or null if none is to be used)
   * @return a PMML model
   * @throws Exception if there is a problem while reading the file
   */
  public static PMMLModel getPMMLModel(String filename, Logger log) throws Exception {
    return getPMMLModel(new File(filename), log);
  }

  /**
   * Read and return a PMML model.
   *
   * @param file a <code>File</code> to read from
   * @param log the logging object to use (or null if none is to be used)
   * @return a PMML model
   * @throws Exception if there is a problem while reading the file
   */
  public static PMMLModel getPMMLModel(File file, Logger log) throws Exception {
    return getPMMLModel(new BufferedInputStream(new FileInputStream(file)), log);
  }
  
  private static boolean isPMML(Document doc) {
    NodeList tempL = doc.getElementsByTagName("PMML");
    if (tempL.getLength() == 0) {
      return false;
    }
    
    return true;
  }

  /**
   * Read and return a PMML model.
   *
   * @param stream the <code>InputStream</code> to read from
   * @param log the logging object to use (or null if none is to be used)
   * @return a PMML model
   * @throws Exception if there is a problem while reading from the stream
   */
  public static PMMLModel getPMMLModel(InputStream stream, Logger log) throws Exception {
    DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
    DocumentBuilder db = dbf.newDocumentBuilder();
    Document doc = db.parse(stream);
    stream.close();
    doc.getDocumentElement().normalize();
    if (!isPMML(doc)) {
      throw new IllegalArgumentException("[PMMLFactory] Source is not a PMML file!!");
    }
    
    //    System.out.println("Root element " + doc.getDocumentElement().getNodeName());

    Instances dataDictionary = getDataDictionaryAsInstances(doc);
    TransformationDictionary transDict = getTransformationDictionary(doc, dataDictionary);
    
    ModelType modelType = getModelType(doc);
    if (modelType == ModelType.UNKNOWN_MODEL) {
      throw new Exception("Unsupported PMML model type");
    }
    Element model = getModelElement(doc, modelType);

    // Construct mining schema and meta data
    MiningSchema ms = new MiningSchema(model, dataDictionary, transDict);
    
    //System.out.println(ms);
    //System.exit(1);
    //    Instances miningSchema = getMiningSchemaAsInstances(model, dataDictionary);
    PMMLModel theModel = getModelInstance(doc, modelType, model, dataDictionary, ms);
    if (log != null) {
      theModel.setLog(log);
    }
    return theModel;
  }
  
  /**
   * Get the transformation dictionary (if there is one).
   * 
   * @param doc the Document containing the PMML model
   * @param dataDictionary the data dictionary as an Instances object
   * @return the transformation dictionary or null if there is none defined in
   * the Document
   * @throws Exception if there is a problem getting the transformation
   * dictionary
   */
  protected static TransformationDictionary getTransformationDictionary(Document doc, 
      Instances dataDictionary) throws Exception {
    TransformationDictionary transDict = null;
    
    NodeList transL = doc.getElementsByTagName("TransformationDictionary");
    // should be of size 0 or 1
    if (transL.getLength() > 0) {
      Node transNode = transL.item(0);
      if (transNode.getNodeType() == Node.ELEMENT_NODE) {
        transDict = new TransformationDictionary((Element)transNode, dataDictionary);
      }
    }
    
    return transDict;
  }
  

  /**
   * Serialize a <code>PMMLModel</code> object that encapsulates a PMML model
   *
   * @param model the <code>PMMLModel</code> to serialize
   * @param filename the name of the file to save to
   * @throws Exception if something goes wrong during serialization
   */
  public static void serializePMMLModel(PMMLModel model, String filename) 
    throws Exception {
    serializePMMLModel(model, new File(filename));
  }

  /**
   * Serialize a <code>PMMLModel</code> object that encapsulates a PMML model
   *
   * @param model the <code>PMMLModel</code> to serialize
   * @param file the <code>File</code> to save to
   * @throws Exception if something goes wrong during serialization
   */
  public static void serializePMMLModel(PMMLModel model, File file)
    throws Exception {
    serializePMMLModel(model, new BufferedOutputStream(new FileOutputStream(file)));
  }

  /**
   * Serialize a <code>PMMLModel</code> object that encapsulates a PMML model
   *
   * @param model the <code>PMMLModel</code> to serialize
   * @param stream the <code>OutputStream</code> to serialize to
   * @throws Exception if something goes wrong during serialization
   */
  public static void serializePMMLModel(PMMLModel model, OutputStream stream)
    throws Exception {
    ObjectOutputStream oo = new ObjectOutputStream(stream);
    Instances header = model.getMiningSchema().getFieldsAsInstances();
    oo.writeObject(header);
    oo.writeObject(model);
    oo.flush();
    oo.close();
  }

  /**
   * Get an instance of a PMMLModel from the supplied Document
   *
   * @param doc the Document holding the pmml
   * @param modelType the type of model
   * @param model the Element encapsulating the model part of the Document
   * @param dataDictionary the data dictionary as an Instances object
   * @param miningSchema the mining schema
   * @return a PMMLModel object
   * @throws Exception if there is a problem constructing the model or
   * if the model type is not supported
   */
  protected static PMMLModel getModelInstance(Document doc, 
                                              ModelType modelType, Element model,
                                              Instances dataDictionary,
                                              MiningSchema miningSchema) throws Exception {
    PMMLModel pmmlM = null;
    switch (modelType) {
    case REGRESSION_MODEL:
      pmmlM = new Regression(model, dataDictionary, miningSchema);
      //System.out.println(pmmlM);
      break;
    case GENERAL_REGRESSION_MODEL:
      pmmlM = new GeneralRegression(model, dataDictionary, miningSchema);
      //System.out.println(pmmlM);
      break;
    case NEURAL_NETWORK_MODEL:
      pmmlM = new NeuralNetwork(model, dataDictionary, miningSchema);
      break;
    case TREE_MODEL:
      pmmlM = new TreeModel(model, dataDictionary, miningSchema);
      break;
    case RULESET_MODEL:
      pmmlM = new RuleSetModel(model, dataDictionary, miningSchema);
      break;
    case SVM_MODEL:
      pmmlM = new SupportVectorMachineModel(model, dataDictionary, miningSchema);
      break;
    default:
      throw new Exception("[PMMLFactory] Unknown model type!!");
    }
    pmmlM.setPMMLVersion(doc);
    pmmlM.setCreatorApplication(doc);
    return pmmlM;
  }

  /**
   * Get the type of model
   *
   * @param doc the Document encapsulating the pmml
   * @return the type of model
   */
  protected static ModelType getModelType(Document doc) {
    NodeList temp = doc.getElementsByTagName("RegressionModel");
    if (temp.getLength() > 0) {
      return ModelType.REGRESSION_MODEL;
    }

    temp = doc.getElementsByTagName("GeneralRegressionModel");
    if (temp.getLength() > 0) {
      return ModelType.GENERAL_REGRESSION_MODEL;
    }
    
    temp = doc.getElementsByTagName("NeuralNetwork");
    if (temp.getLength() > 0) {
      return ModelType.NEURAL_NETWORK_MODEL;
    }
    
    temp = doc.getElementsByTagName("TreeModel");
    if (temp.getLength() > 0) {
      return ModelType.TREE_MODEL;
    }
    
    temp = doc.getElementsByTagName("RuleSetModel");
    if (temp.getLength() > 0) {
      return ModelType.RULESET_MODEL;
    }
    
    temp = doc.getElementsByTagName("SupportVectorMachineModel");
    if (temp.getLength() > 0) {
      return ModelType.SVM_MODEL;
    }

    return ModelType.UNKNOWN_MODEL;
  }

  /**
   * Get the Element that contains the pmml model
   *
   * @param doc the Document encapsulating the pmml
   * @param modelType the type of model
   * @throws Exception if the model type is unsupported/unknown
   */
  protected static Element getModelElement(Document doc, ModelType modelType) 
    throws Exception {
    NodeList temp = null;
    Element model = null;
    switch (modelType) {
    case REGRESSION_MODEL:
      temp = doc.getElementsByTagName("RegressionModel");
      break;
    case GENERAL_REGRESSION_MODEL:
      temp = doc.getElementsByTagName("GeneralRegressionModel");
      break;
    case NEURAL_NETWORK_MODEL:
      temp = doc.getElementsByTagName("NeuralNetwork");
      break;
    case TREE_MODEL:
      temp = doc.getElementsByTagName("TreeModel");
      break;
    case RULESET_MODEL:
      temp = doc.getElementsByTagName("RuleSetModel");
      break;
    case SVM_MODEL:
      temp = doc.getElementsByTagName("SupportVectorMachineModel");
      break;
    default:
      throw new Exception("[PMMLFactory] unknown/unsupported model type.");
    }

    if (temp != null && temp.getLength() > 0) {
      Node modelNode = temp.item(0);
      if (modelNode.getNodeType() == Node.ELEMENT_NODE) {
        model = (Element)modelNode;
      }
    }

    return model;
  }

  /**
   * Get the mining schema as an Instances object
   *
   * @param model the Element containing the pmml model
   * @param dataDictionary the data dictionary as an Instances object
   * @return the mining schema as an Instances object
   * @throws Exception if something goes wrong during reading the mining schema
   * @deprecated Use the MiningSchema class instead
   */
  protected static Instances getMiningSchemaAsInstances(Element model,
                                                        Instances dataDictionary) 
    throws Exception {
    ArrayList<Attribute> attInfo = new ArrayList<Attribute>();
    NodeList fieldList = model.getElementsByTagName("MiningField");
    int classIndex = -1;
    int addedCount = 0;
    for (int i = 0; i < fieldList.getLength(); i++) {
      Node miningField = fieldList.item(i);
      if (miningField.getNodeType() == Node.ELEMENT_NODE) {
        Element miningFieldEl = (Element)miningField;
        String name = miningFieldEl.getAttribute("name");
        String usage = miningFieldEl.getAttribute("usageType");
        // TO-DO: also missing value replacement etc.

        // find this attribute in the dataDictionary
        Attribute miningAtt = dataDictionary.attribute(name);
        if (miningAtt != null) {
          if (usage.length() == 0 || usage.equals("active") || usage.equals("predicted")) {
            attInfo.add(miningAtt);
            addedCount++;
          }
          if (usage.equals("predicted")) {
            classIndex = addedCount - 1;
          }
        } else {
          throw new Exception("Can't find mining field: " + name 
                              + " in the data dictionary.");
        }
      }
    }
    
    Instances insts = new Instances("miningSchema", attInfo, 0);
    //    System.out.println(insts);
    if (classIndex != -1) {
      insts.setClassIndex(classIndex);
    }


    return insts;
  }

  /**
   * Get the data dictionary as an Instances object
   *
   * @param doc the Document encapsulating the pmml
   * @return the data dictionary as an Instances object
   * @throws Exception if there are fields that are not continuous, 
   * ordinal or categorical in the data dictionary
   */
  protected static Instances getDataDictionaryAsInstances(Document doc) 
    throws Exception {

    // TO-DO: definition of missing values (see below)

    ArrayList<Attribute> attInfo = new ArrayList<Attribute>();
    NodeList dataDictionary = doc.getElementsByTagName("DataField");
    for (int i = 0; i < dataDictionary.getLength(); i++) {
      Node dataField = dataDictionary.item(i);
      if (dataField.getNodeType() == Node.ELEMENT_NODE) {
        Element dataFieldEl = (Element)dataField;
        String name = dataFieldEl.getAttribute("name");
        String type = dataFieldEl.getAttribute("optype");
        Attribute tempAtt = null;
        if (name != null && type != null) {
          if (type.equals("continuous")) {
            tempAtt = new Attribute(name);
          } else if (type.equals("categorical") || type.equals("ordinal")) {
            NodeList valueList = dataFieldEl.getElementsByTagName("Value");
            if (valueList == null || valueList.getLength() == 0) {
              // assume that categorical values will be revealed in the actual model.
              // Create a string attribute for now
              ArrayList<String> nullV = null;
              tempAtt = new Attribute(name, (ArrayList<String>)nullV);
            } else {
              // add the values (if defined as "valid")
              ArrayList<String> valueVector = new ArrayList<String>();
              for (int j = 0; j < valueList.getLength(); j++) {
                Node val = valueList.item(j);
                if (val.getNodeType() == Node.ELEMENT_NODE) {
                  // property is optional (default value is "valid")
                  String property = ((Element)val).getAttribute("property");
                  if (property == null || property.length() == 0 || property.equals("valid")) {
                    String value = ((Element)val).getAttribute("value");
                    valueVector.add(value);
                  } else {
                    // Just ignore invalid or missing value definitions for now...
                    // TO-DO: implement Value meta data with missing/invalid value defs.
                  }
                }
              }
              tempAtt = new Attribute(name, valueVector);
            }
          } else {
            throw new Exception("[PMMLFactory] can't handle " + type + "attributes.");
          }
          attInfo.add(tempAtt);
        }
      }
    }

    // TO-DO: check whether certain values are declared to represent
    // missing or invalid values (applies to both categorical and continuous
    // attributes

    // create the Instances structure
    Instances insts = new Instances("dataDictionary", attInfo, 0);
    //    System.out.println(insts);

    return insts;
  }

  public static String applyClassifier(PMMLModel model, Instances test) throws Exception {
    StringBuffer buff = new StringBuffer();
    if (!(model instanceof PMMLClassifier)) {
      throw new Exception("PMML model is not a classifier!");
    }

    double[] preds = null;
    PMMLClassifier classifier = (PMMLClassifier)model;
    for (int i = 0; i < test.numInstances(); i++) {
      buff.append("Actual: ");
      Instance temp = test.instance(i);
      if (temp.classAttribute().isNumeric()) {
        buff.append(temp.value(temp.classIndex()) + " ");
      } else {
        buff.append(temp.classAttribute().value((int)temp.value(temp.classIndex())) + " ");
      }
      preds = classifier.distributionForInstance(temp);
      buff.append(" Predicted: ");
      for (int j = 0; j < preds.length; j++) {
        buff.append("" + preds[j] + " ");
      }
      buff.append("\n");
    }
    return buff.toString();
  }
  
  private static class PMMLClassifierRunner extends AbstractClassifier {
    public double[] distributionForInstance(Instance test) throws Exception {
      throw new Exception("Don't call this method!!");
    }
    
    public void buildClassifier(Instances instances) throws Exception {
      throw new Exception("Don't call this method!!");
    }
    
    public String getRevision() {
      return weka.core.RevisionUtils.extract("$Revision: 6467 $");
    }
    
    public void evaluatePMMLClassifier(String[] options) {
      runClassifier(this, options);
    }
  }

  public static void main(String[] args) {
    try {
      String[] optionsTmp = new String[args.length];
      for (int i = 0; i < args.length; i++) {
        optionsTmp[i] = args[i];
      }
      String pmmlFile = Utils.getOption('l', optionsTmp); 
      if (pmmlFile.length() == 0) {
        throw new Exception("[PMMLFactory] must specify a PMML file using the -l option.");
      }
      // see if it is supported before going any further
      getPMMLModel(pmmlFile, null);
      
      PMMLClassifierRunner pcr = new PMMLClassifierRunner();
      pcr.evaluatePMMLClassifier(args);

      
      /*PMMLModel model = getPMMLModel(args[0], null);
      System.out.println(model);
      if (args.length == 2) {
        // load an arff file
        Instances testData = new Instances(new java.io.BufferedReader(new java.io.FileReader(args[1])));
        Instances miningSchemaI = model.getMiningSchema().getFieldsAsInstances();
        if (miningSchemaI.classIndex() >= 0) {
          String className = miningSchemaI.classAttribute().name();
          for (int i = 0; i < testData.numAttributes(); i++) {
            if (testData.attribute(i).name().equals(className)) {
              testData.setClassIndex(i);
              System.out.println("Found class " + className + " in test data.");
              break;
            }
          }
        }
        System.out.println(applyClassifier(model, testData));
      }*/
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
