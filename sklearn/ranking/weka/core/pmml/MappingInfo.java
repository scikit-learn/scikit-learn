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
 *    MappingInfo.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.io.Serializable;
import java.util.ArrayList;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;
import weka.gui.Logger;

/**
 * Class that maintains the mapping between incoming data set structure
 * and that of the mining schema.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com
 * @version $Revision: 6463 $
 */
public class MappingInfo implements Serializable {
  
  /** For serialization */
  
  /** 
   * Index for incoming nominal values that are not defined in the mining
   * schema.
   */
  public static final int UNKNOWN_NOMINAL_VALUE = -1;
  
  /** 
   * Map the incoming attributes to the mining schema attributes.
   * Each entry holds the index of the incoming attribute that
   * corresponds to this mining schema attribute.  
   */
  private int[] m_fieldsMap = null;
  
  /** 
   * Map indexes for nominal values in incoming structure to those
   * in the mining schema. There will be as many entries as there are
   * attributes in this array. Non-nominal attributes will have
   * null entries. Each non-null entry is an array of integer indexes.
   * Each entry in a given array (for a given attribute) holds the index of
   * the mining schema value that corresponds to this incoming value.
   * UNKNOWN_NOMINAL_VALUE is used as the index for those incoming values 
   * that are not defined in the mining schema. 
   */
  private int[][] m_nominalValueMaps = null;
  
  /** Holds a textual description of the fields mapping */
  private String m_fieldsMappingText = null;
  
  /** For logging */
  private Logger m_log = null;
  
  public MappingInfo(Instances dataSet, MiningSchema miningSchema,
                     Logger log) throws Exception {
    m_log = log;
    //miningSchema.convertStringAttsToNominal();
    Instances fieldsI = miningSchema.getMiningSchemaAsInstances();
 
    m_fieldsMap = new int[fieldsI.numAttributes()];
    m_nominalValueMaps = new int[fieldsI.numAttributes()][];
    
    for (int i = 0; i < fieldsI.numAttributes(); i++) {
      String schemaAttName = fieldsI.attribute(i).name();
      boolean found = false;
      for (int j = 0; j < dataSet.numAttributes(); j++) {
        if (dataSet.attribute(j).name().equals(schemaAttName)) {
          Attribute miningSchemaAtt = fieldsI.attribute(i);
          Attribute incomingAtt = dataSet.attribute(j);
          // check type match
          if (miningSchemaAtt.type() != incomingAtt.type()) {
            if (miningSchemaAtt.isString() && (incomingAtt.isNominal() || incomingAtt.isRanking())) {
              // don't worry about String attributes in the mining schema
              // (as long as the corresponding incoming is a String or nominal),
              // since values for the String attributes are more than likely revealed
              // by FieldRef elements in the actual model itself
            } else {
              throw new Exception("[MappingInfo] type mismatch for field " +
                  schemaAttName + ". Mining schema type " + 
                  miningSchemaAtt.toString() + ". Incoming type " +
                  incomingAtt.toString() + ".");
            }
          }

          // check nominal values (number, names...)
          if (miningSchemaAtt.numValues() != incomingAtt.numValues()) {
            String warningString = "[MappingInfo] WARNING: incoming nominal attribute "
              + incomingAtt.name() + " does not have the same "
              + "number of values as the corresponding mining "
              + "schema attribute.";
            if (m_log != null) {
              m_log.logMessage(warningString);
            } else {
              System.err.println(warningString);
            }
          }
          if (miningSchemaAtt.isNominal()|| miningSchemaAtt.isRanking() || miningSchemaAtt.isString()) {
            int[] valuesMap = new int[incomingAtt.numValues()];
            for (int k = 0; k < incomingAtt.numValues(); k++) {
              String incomingNomVal = incomingAtt.value(k);
              int indexInSchema = miningSchemaAtt.indexOfValue(incomingNomVal);
              if (indexInSchema < 0) {
                String warningString = "[MappingInfo] WARNING: incoming nominal attribute "
                  + incomingAtt.name() + " has value "
                  + incomingNomVal + " that doesn't occur in the mining schema.";
                if (m_log != null) {
                  m_log.logMessage(warningString);
                } else {
                  System.err.println(warningString);
                }
                valuesMap[k] = UNKNOWN_NOMINAL_VALUE;
              } else {
                valuesMap[k] = indexInSchema;
              }
            }
            m_nominalValueMaps[i] = valuesMap;
          }

          /*if (miningSchemaAtt.isNominal()) {
            for (int k = 0; k < miningSchemaAtt.numValues(); k++) {
              if (!miningSchemaAtt.value(k).equals(incomingAtt.value(k))) {
                throw new Exception("[PMMLUtils] value " + k + " (" + 
                                    miningSchemaAtt.value(k) + ") does not match " +
                                    "incoming value (" + incomingAtt.value(k) +
                                    ") for attribute " + miningSchemaAtt.name() +
                                    ".");

              }
            }
          }*/
          found = true;
          m_fieldsMap[i] = j;
        }
      }
      if (!found) {
        throw new Exception("[MappingInfo] Unable to find a match for mining schema "
            + "attribute " + schemaAttName + " in the "
            + "incoming instances!");
      }
    }

    // check class attribute (if set)
   if (fieldsI.classIndex() >= 0) {
      if (dataSet.classIndex() < 0) {
	// first see if we can find a matching class
	String className = fieldsI.classAttribute().name();
	Attribute classMatch = dataSet.attribute(className);
	if (classMatch == null) {
	  throw new Exception("[MappingInfo] Can't find match for target field " + className
	      + "in incoming instances!");
	}
	dataSet.setClass(classMatch);
      } else if (!fieldsI.classAttribute().name().equals(dataSet.classAttribute().name())) {
        throw new Exception("[MappingInfo] class attribute in mining schema does not match "
            + "class attribute in incoming instances!");
      }
    }
    
    // Set up the textual description of the mapping
    fieldsMappingString(fieldsI, dataSet);
  }
  
  private void fieldsMappingString(Instances miningSchemaI, Instances incomingI) {
    StringBuffer result = new StringBuffer();
    
    int maxLength = 0;
    for (int i = 0; i < miningSchemaI.numAttributes(); i++) {
      if (miningSchemaI.attribute(i).name().length() > maxLength) {
        maxLength = miningSchemaI.attribute(i).name().length();
      }
    }
    maxLength += 12; // length of " (nominal)"/" (numeric)"
    
    int minLength = 13; // "Mining schema".length()                                                        
    String headerS = "Mining schema";
    String sep = "-------------";

    if (maxLength < minLength) {
      maxLength = minLength;
    }
    
    headerS = PMMLUtils.pad(headerS, " ", maxLength, false);
    sep = PMMLUtils.pad(sep, "-", maxLength, false);
    
    sep += "\t    ----------------\n";
    headerS += "\t    Incoming fields\n";
    result.append(headerS);
    result.append(sep);
    
    for (int i = 0; i < miningSchemaI.numAttributes(); i++) {
      Attribute temp = miningSchemaI.attribute(i);
      String attName = "("
        + ((temp.isNumeric())
           ? "numeric)"
           : "nominal)")
        + " " + temp.name();
      attName = PMMLUtils.pad(attName, " ", maxLength, false);
      attName +=  "\t--> ";
      result.append(attName);
      
      Attribute incoming = incomingI.attribute(m_fieldsMap[i]);
      String fieldName = "" + (m_fieldsMap[i] + 1) + " ("
        + ((incoming.isNumeric())
            ? "numeric)"
            : "nominal)");
      fieldName += " " + incoming.name();
      result.append(fieldName + "\n");
    }
    
    m_fieldsMappingText = result.toString();
  }
  
  /**
   * Convert an <code>Instance</code> to an array of values that matches the
   * format of the mining schema. First maps raw attribute values and then
   * applies rules for missing values, outliers etc.
   *
   * @param inst the <code>Instance</code> to convert
   * @param miningSchema the mining schema
   * incoming instance attributes
   * @return an array of doubles that are values from the incoming Instances,
   * correspond to the format of the mining schema and have had missing values,
   * outliers etc. dealt with.
   * @throws Exception if something goes wrong
   */
  public double[] instanceToSchema(Instance inst, 
                                   MiningSchema miningSchema) throws Exception {
    Instances miningSchemaI = miningSchema.getMiningSchemaAsInstances();
    
    // allocate enough space for both mining schema fields and any derived fields
    double[] result = new double[miningSchema.getFieldsAsInstances().numAttributes()];

    // Copy over the values
    for (int i = 0; i < miningSchemaI.numAttributes(); i++) {
      //if (miningSchemaI.attribute(i).isNumeric()) {
      result[i] = inst.value(m_fieldsMap[i]);
      if (miningSchemaI.attribute(i).isNominal() || 
    	  miningSchemaI.attribute(i).isRanking() ||
          miningSchemaI.attribute(i).isString()) {
        // If not missing, look up the index of this incoming categorical value in
        // the mining schema
        if (!Utils.isMissingValue(inst.value(m_fieldsMap[i]))) {
          int[] valueMap = m_nominalValueMaps[i];
          int index = valueMap[(int)inst.value(m_fieldsMap[i])];
          String incomingAttValue = 
            inst.attribute(m_fieldsMap[i]).value((int)inst.value(m_fieldsMap[i]));
          /*int index = miningSchemaI.attribute(i).indexOfValue(incomingAttValue); */
          if (index >= 0) {
            result[i] = index;
          } else {
            // set this to "unknown" (-1) for nominal valued attributes
            result[i] = UNKNOWN_NOMINAL_VALUE;
            String warningString = "[MappingInfo] WARNING: Can't match nominal value "
              + incomingAttValue;
            if (m_log != null) {
              m_log.logMessage(warningString);
            } else {
              System.err.println(warningString);
            }
          }
        }
      }
    }

    // Now deal with missing values and outliers...
    miningSchema.applyMissingAndOutlierTreatments(result);
    //    printInst(result);
    
    // now fill in any derived values
    ArrayList<DerivedFieldMetaInfo> derivedFields = miningSchema.getDerivedFields();
    for (int i = 0; i < derivedFields.size(); i++) {
      DerivedFieldMetaInfo temp = derivedFields.get(i);
//      System.err.println("Applying : " + temp);
      double r = temp.getDerivedValue(result);
      result[i + miningSchemaI.numAttributes()] = r;
    }
    
    /*System.err.print("==> ");
    for (int i = 0; i < result.length; i++) {
      System.err.print(" " + result[i]);
    }
    System.err.println();*/
    
    return result;
  }
  
  /**
   * Get a textual description of them mapping between mining schema
   * fields and incoming data fields.
   * 
   * @return a description of the fields mapping as a String
   */
  public String getFieldsMappingString() {
    if (m_fieldsMappingText == null) {
      return "No fields mapping constructed!";
    }
    return m_fieldsMappingText;
  }
}

