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
 *    MiningSchema.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.lang.String;
import java.io.Serializable;
import java.util.ArrayList;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.core.Attribute;
import weka.core.Instances;

/**
 * This class encapsulates the mining schema from
 * a PMML xml file. Specifically, it contains the
 * fields used in the PMML model as an Instances
 * object (just the header). It also contains meta
 * information such as value ranges and how to handle
 * missing values, outliers etc.
 *
 * We also store various other PMML elements here, such as
 * the TransformationDictionary, DerivedFields and Targets 
 * (if defined). They are not part of the mining schema per se, but
 * relate to inputs used by the model and it is convenient to
 * store them here.
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 5987 $
 */
public class MiningSchema implements Serializable {

  /** For serialization */
  private static final long serialVersionUID = 7144380586726330455L;

  /** The structure of all the fields (both mining schema and derived) as Instances */
  protected Instances m_fieldInstancesStructure;
  
  /** Just the mining schema fields as Instances */
  protected Instances m_miningSchemaInstancesStructure;

  /** Meta information about the mining schema fields */
  protected ArrayList<MiningFieldMetaInfo> m_miningMeta = 
    new ArrayList<MiningFieldMetaInfo>();
  
  /** 
   * Meta information about derived fields (those defined in
   * the TransformationDictionary followed by those defined in
   * LocalTransformations)
   */
  protected ArrayList<DerivedFieldMetaInfo> m_derivedMeta = 
    new ArrayList<DerivedFieldMetaInfo>();
  
  /** The transformation dictionary (if defined) */
  protected TransformationDictionary m_transformationDictionary = null;

  /** target meta info (may be null if not defined) */
  protected TargetMetaInfo m_targetMetaInfo = null;
  
  private void getLocalTransformations(Element model) throws Exception {
    NodeList temp = model.getElementsByTagName("LocalTransformations");
    
    if (temp.getLength() > 0) {
      // should be just one LocalTransformations element
      Element localT = (Element)temp.item(0);
      
      // Set up some field defs to pass in
      /*ArrayList<Attribute> fieldDefs = new ArrayList<Attribute>();
      for (int i = 0; i < m_miningSchemaInstancesStructure.numAttributes(); i++) {
        fieldDefs.add(m_miningSchemaInstancesStructure.attribute(i));
      } */
      
      NodeList localDerivedL = localT.getElementsByTagName("DerivedField");
      for (int i = 0; i < localDerivedL.getLength(); i++) {
        Node localDerived = localDerivedL.item(i);
        if (localDerived.getNodeType() == Node.ELEMENT_NODE) {
          DerivedFieldMetaInfo d = 
            new DerivedFieldMetaInfo((Element)localDerived, null /*fieldDefs*/, m_transformationDictionary);
          m_derivedMeta.add(d);
        }
      }
    }
  }

  /**
   * Constructor for MiningSchema.
   *
   * @param model the <code>Element</code> encapsulating the pmml model
   * @param dataDictionary the data dictionary as an Instances object
   * @throws Exception if something goes wrong during construction of the
   * mining schema
   */
  public MiningSchema(Element model, 
                      Instances dataDictionary,
                      TransformationDictionary transDict) throws Exception {
    
    /*// First check for transformation dictionary/local transformations and derived fields.
    // These are not supported yet.
    NodeList temp = model.getElementsByTagName("LocalTransformations");
    if (temp.getLength() > 0) {
      throw new Exception("[MiningSchema] LocalTransformations "
          + "are not supported yet.");
    }*/

    ArrayList<Attribute> attInfo = new ArrayList<Attribute>();
    NodeList fieldList = model.getElementsByTagName("MiningField");
    int classIndex = -1;
    int addedCount = 0;
    for (int i = 0; i < fieldList.getLength(); i++) {
      Node miningField = fieldList.item(i);
      if (miningField.getNodeType() == Node.ELEMENT_NODE) {
        Element miningFieldEl = (Element)miningField;

        MiningFieldMetaInfo mfi = new MiningFieldMetaInfo(miningFieldEl);

        if (mfi.getUsageType() == MiningFieldMetaInfo.Usage.ACTIVE ||
            mfi.getUsageType() == MiningFieldMetaInfo.Usage.PREDICTED) {

          // find this attribute in the dataDictionary
          Attribute miningAtt = dataDictionary.attribute(mfi.getName());
          if (miningAtt != null) {
            mfi.setIndex(addedCount);
            attInfo.add(miningAtt);
            addedCount++;

            if (mfi.getUsageType() == MiningFieldMetaInfo.Usage.PREDICTED) {
              classIndex = addedCount - 1;
            }

            // add to the array list
            m_miningMeta.add(mfi);
          } else {
          throw new Exception("Can't find mining field: " + mfi.getName() 
                              + " in the data dictionary.");
          }
        }
      }
    }

    m_miningSchemaInstancesStructure = new Instances("miningSchema", attInfo, 0);
    
    // set these instances on the MiningFieldMetaInfos so that the
    // toString() method can operate correctly
    for (MiningFieldMetaInfo m : m_miningMeta) {
      m.setMiningSchemaInstances(m_miningSchemaInstancesStructure);
    }
    
    m_transformationDictionary = transDict;
    
    // Handle transformation dictionary and any local transformations
    if (m_transformationDictionary != null) {      
      ArrayList<DerivedFieldMetaInfo> transDerived = transDict.getDerivedFields();
      m_derivedMeta.addAll(transDerived);
    }
    
    // Get any local transformations
    getLocalTransformations(model);
    
    // Set up the full instances structure: combo of mining schema fields and
    // all derived fields
    ArrayList<Attribute> newStructure = new ArrayList<Attribute>();
    for (MiningFieldMetaInfo m : m_miningMeta) {
      newStructure.add(m.getFieldAsAttribute());
    }
    
    for (DerivedFieldMetaInfo d : m_derivedMeta) {
      newStructure.add(d.getFieldAsAttribute());
    }
    m_fieldInstancesStructure = new Instances("FieldStructure", newStructure, 0);
    
    if (m_transformationDictionary != null) {
      // first update the field defs for any derived fields in the transformation dictionary
      // and our complete list of derived fields, now that we have a fixed 
      // ordering for the mining schema + derived attributes (i.e. could
      // be different from the order of attributes in the data dictionary that was
      // used when the transformation dictionary was initially constructed
      m_transformationDictionary.setFieldDefsForDerivedFields(m_fieldInstancesStructure);
    }

    // update the field defs for all our derived fields.
    for (DerivedFieldMetaInfo d : m_derivedMeta) {
      d.setFieldDefs(m_fieldInstancesStructure);
    }
    
    if (classIndex != -1) {
      m_fieldInstancesStructure.setClassIndex(classIndex);
      m_miningSchemaInstancesStructure.setClassIndex(classIndex);
    }

    // do Targets (if any)
    NodeList targetsList = model.getElementsByTagName("Targets");
    if (targetsList.getLength() > 0) {
      if (targetsList.getLength() > 1) {
        throw new Exception("[MiningSchema] Can only handle a single Target");
      } else {
        Node te = targetsList.item(0);
        if (te.getNodeType() == Node.ELEMENT_NODE) {
          m_targetMetaInfo = new TargetMetaInfo((Element)te);

          // fill in any necessary categorical values in the mining schema 
          // class attribute
          if (m_fieldInstancesStructure.classIndex() >= 0 && 
              m_fieldInstancesStructure.classAttribute().isString()) {
            ArrayList<String> targetVals = m_targetMetaInfo.getValues();
            if (targetVals.size() > 0) {
              Attribute classAtt = m_fieldInstancesStructure.classAttribute();
              for (int i = 0; i < targetVals.size(); i++) {
                classAtt.addStringValue(targetVals.get(i));
              }
            }
          }
        }
      }
    }
  }

  /**
   * Apply the missing value treatments (if any) to an incoming instance.
   *
   * @param values an array of doubles in order of the fields in the mining schema
   * that represents the incoming instance (note: use PMMLUtils.instanceToSchema()
   * to generate this).
   * @throws Exception if something goes wrong during missing value handling
   */
  public void applyMissingValuesTreatment(double[] values) throws Exception {
    for (int i = 0; i < m_miningMeta.size(); i++) {
      MiningFieldMetaInfo mfi = m_miningMeta.get(i);
      values[i] = mfi.applyMissingValueTreatment(values[i]);
    }
  }

  /**
   * Apply the outlier treatment methods (if any) to an incoming instance.
   *
   * @param values an array of doubles in order of the fields in the mining schema
   * that represents the incoming instance (note: use PMMLUtils.instanceToSchema()
   * to generate this).
   * @throws Exception if something goes wrong during outlier treatment handling
   */
  public void applyOutlierTreatment(double[] values) throws Exception {
    for (int i = 0; i < m_miningMeta.size(); i++) {
      MiningFieldMetaInfo mfi = m_miningMeta.get(i);
      values[i] = mfi.applyOutlierTreatment(values[i]);
    }
  }

  /**
   * Apply both missing and outlier treatments to an incoming instance.
   * @param values an array of doubles in order of the fields in the mining schema
   * that represents the incoming instance (note: use MappingInfo.instanceToSchema()
   * to generate this).
   * @throws Exception if something goes wrong during this process
   */
  public void applyMissingAndOutlierTreatments(double[] values) throws Exception {
    for (int i = 0; i < m_miningMeta.size(); i++) {
      MiningFieldMetaInfo mfi = m_miningMeta.get(i);
      values[i] = mfi.applyMissingValueTreatment(values[i]);
      values[i] = mfi.applyOutlierTreatment(values[i]);
    }
  }

  /**
   * Get the all the fields (both mining schema and derived) as Instances.
   * Attributes are in order of those in the mining schema, followed by
   * derived attributes from the TransformationDictionary followed by
   * derived attributes from LocalTransformations.
   *
   * @return all the fields as an Instances object
   */
  public Instances getFieldsAsInstances() {
    return m_fieldInstancesStructure;
  }
  
  /**
   * Get the mining schema fields as an Instances object.
   * 
   * @return the mining schema fields as an Instances object.
   */
  public Instances getMiningSchemaAsInstances() {
    return m_miningSchemaInstancesStructure;
  }
  
  /**
   * Get the transformation dictionary .
   * 
   * @return the transformation dictionary or null if none is
   * defined.
   */
  public TransformationDictionary getTransformationDictionary() {
    return m_transformationDictionary;
  }
  
  /**
   * Returns true if there is Target meta data.
   *
   * @return true if there is Target meta data
   */
  public boolean hasTargetMetaData() {
    return (m_targetMetaInfo != null);
  }

  /**
   * Get the Target meta data.
   *
   * @return the Target meta data
   */
  public TargetMetaInfo getTargetMetaData() {
    return m_targetMetaInfo;
  }

  /**
   * Method to convert any string attributes in the mining schema
   * Instances to nominal attributes. This may be necessary if there are
   * no Value elements defined for categorical fields in the data dictionary.
   * In this case, elements in the actual model definition will probably reveal
   * the valid values for categorical fields.
   */
  public void convertStringAttsToNominal() {
    Instances miningSchemaI = getFieldsAsInstances();
    if (miningSchemaI.checkForStringAttributes()) {
      ArrayList<Attribute> attInfo = new ArrayList<Attribute>();
      for (int i = 0; i < miningSchemaI.numAttributes(); i++) {
        Attribute tempA = miningSchemaI.attribute(i);
        if (tempA.isString()) {
          ArrayList<String> valueVector = new ArrayList<String>();
          for (int j = 0; j < tempA.numValues(); j++) {
            valueVector.add(tempA.value(j));
          }
          Attribute newAtt = new Attribute(tempA.name(), valueVector);
          attInfo.add(newAtt);
        } else {
          attInfo.add(tempA);
        }
      }
      Instances newI = new Instances("miningSchema", attInfo, 0);
      if (m_fieldInstancesStructure.classIndex() >= 0) {
        newI.setClassIndex(m_fieldInstancesStructure.classIndex());
      }
      m_fieldInstancesStructure = newI;

      /*      StringToNominal stn = new StringToNominal();
      stn.setInputFormat(miningSchemaI);
      Instances newI = Filter.useFilter(miningSchemaI, stn);
      m_miningSchema = newI; */
    }
  }

  /**
   * Convert a numeric attribute in the mining schema to nominal.
   * 
   * @param index the index of the attribute to convert
   * @param newVals an ArrayList of the values of the nominal attribute
   */
  public void convertNumericAttToNominal(int index, 
                                         ArrayList<String> newVals) {
    Instances miningSchemaI = getFieldsAsInstances();
    if (miningSchemaI.attribute(index).isNominal() || miningSchemaI.attribute(index).isRanking()) {
      throw new IllegalArgumentException("[MiningSchema] convertNumericAttToNominal: attribute is "
                                         + "already nominal!");
    }

    ArrayList<String> newValues = new ArrayList<String>();
    for (int i = 0; i < newVals.size(); i++) {
      newValues.add(newVals.get(i));
    }

    ArrayList<Attribute> attInfo = new ArrayList<Attribute>();
    for (int i = 0; i < miningSchemaI.numAttributes(); i++) {
      Attribute tempA = miningSchemaI.attribute(i);
      if (i == index) {
        Attribute newAtt = new Attribute(tempA.name(), newValues);
        attInfo.add(newAtt);
      } else {
        attInfo.add(tempA);
      }
    }

    Instances newI = new Instances("miningSchema", attInfo, 0);
    if (m_fieldInstancesStructure.classIndex() >= 0) {
      newI.setClassIndex(m_fieldInstancesStructure.classIndex());
    }
    m_fieldInstancesStructure = newI;
  }
  
  public ArrayList<DerivedFieldMetaInfo> getDerivedFields() {
    return m_derivedMeta;
  }
  
  public ArrayList<MiningFieldMetaInfo> getMiningFields() {
    return m_miningMeta;
  }

  /**
   * Get a textual description of the mining schema.
   *
   * @return a textual description of the mining schema
   */
  public String toString() {
    StringBuffer temp = new StringBuffer();
    
    if (m_transformationDictionary != null) {
      temp.append(m_transformationDictionary);
    }
    
    temp.append("Mining schema:\n\n");
    for (MiningFieldMetaInfo m : m_miningMeta) {
      temp.append(m + "\n");
    }
    
    if (m_derivedMeta.size() > 0) {
      temp.append("\nDerived fields:\n\n");
      for (DerivedFieldMetaInfo d : m_derivedMeta) {
        temp.append(d + "\n");
      }
    }
    temp.append("\n");
    return temp.toString();
  }
}
