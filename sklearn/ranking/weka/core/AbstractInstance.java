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
 *    DenseInstance.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.ArrayList;

import weka.core.labelranking.PreferenceAttribute;
import weka.core.labelranking.PreferenceDenseInstance;

/**
 * Abstract class providing common functionality for the original
 * instance implementations.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5987 $ 
 */
public abstract class AbstractInstance
  implements Instance, Serializable, RevisionHandler {
  
  /** for serialization */
  static final long serialVersionUID = 1482635194499365155L;

  /** 
   * The dataset the instance has access to.  Null if the instance
   * doesn't have access to any dataset.  Only if an instance has
   * access to a dataset, it knows about the actual attribute types.  
   */
  protected /*@spec_public@*/ Instances m_Dataset;

  /** The instance's attribute values. */
  protected /*@spec_public non_null@*/ double[] m_AttValues;

  /** The instance's weight. */
  protected double m_Weight;

  /**
   * Returns the attribute with the given index.
   *
   * @param index the attribute's index
   * @return the attribute at the given position
   * @throws UnassignedDatasetException if instance doesn't have access to a
   * dataset
   */ 
  //@ requires m_Dataset != null;
  public /*@pure@*/ Attribute attribute(int index) {
   
    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    }
    return m_Dataset.attribute(index);
  }

  /**
   * Returns the attribute with the given index in the sparse representation.
   *
   * @param indexOfIndex the index of the attribute's index 
   * @return the attribute at the given position
   * @throws UnassignedDatasetException if instance doesn't have access to a
   * dataset
   */ 
  //@ requires m_Dataset != null;
  public /*@pure@*/ Attribute attributeSparse(int indexOfIndex) {
   
    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    }
    return m_Dataset.attribute(index(indexOfIndex));
  }

  /**
   * Returns class attribute.
   *
   * @return the class attribute
   * @throws UnassignedDatasetException if the class is not set or the
   * instance doesn't have access to a dataset
   */
  //@ requires m_Dataset != null;
  public /*@pure@*/ Attribute classAttribute() {

    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    }
    return m_Dataset.classAttribute();
  }

  /**
   * Returns the class attribute's index.
   *
   * @return the class index as an integer 
   * @throws UnassignedDatasetException if instance doesn't have access to a dataset 
   */
  //@ requires m_Dataset != null;
  //@ ensures  \result == m_Dataset.classIndex();
  public /*@pure@*/ int classIndex() {
    
    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    }
    return m_Dataset.classIndex();
  }

  /**
   * Tests if an instance's class is missing.
   *
   * @return true if the instance's class is missing
   * @throws UnassignedClassException if the class is not set or the instance doesn't
   * have access to a dataset
   */
  //@ requires classIndex() >= 0;
  public /*@pure@*/ boolean classIsMissing() {

    if (classIndex() < 0) {
      throw new UnassignedClassException("Class is not set!");
    }
    return isMissing(classIndex());
  }

  /**
   * Returns an instance's class value in internal format. (ie. as a
   * floating-point number)
   *
   * @return the corresponding value as a double (If the 
   * corresponding attribute is nominal (or a string) then it returns the 
   * value's index as a double).
   * @throws UnassignedClassException if the class is not set or the instance doesn't
   * have access to a dataset 
   */
  //@ requires classIndex() >= 0;
  public /*@pure@*/ double classValue() {
    
    if (classIndex() < 0) {
      throw new UnassignedClassException("Class is not set!");
    }
    return value(classIndex());
  }

  /**
   * Returns the dataset this instance has access to. (ie. obtains
   * information about attribute types from) Null if the instance
   * doesn't have access to a dataset.
   *
   * @return the dataset the instance has accesss to
   */
  //@ ensures \result == m_Dataset;
  public /*@pure@*/ Instances dataset() {

    return m_Dataset;
  }

  /**
   * Deletes an attribute at the given position (0 to 
   * numAttributes() - 1). Only succeeds if the instance does not
   * have access to any dataset because otherwise inconsistencies
   * could be introduced.
   *
   * @param position the attribute's position
   * @throws RuntimeException if the instance has access to a
   * dataset 
   */
  //@ requires m_Dataset != null;
  public void deleteAttributeAt(int position) {

    if (m_Dataset != null) {
      throw new RuntimeException("DenseInstance has access to a dataset!");
    }
    forceDeleteAttributeAt(position);
  }

  /**
   * Returns an enumeration of all the attributes.
   *
   * @return enumeration of all the attributes
   * @throws UnassignedDatasetException if the instance doesn't
   * have access to a dataset 
   */
  //@ requires m_Dataset != null;
  public /*@pure@*/ Enumeration enumerateAttributes() {

    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    }
    return m_Dataset.enumerateAttributes();
  }

  /**
   * Tests if the headers of two instances are equivalent.
   *
   * @param inst another instance
   * @return true if the header of the given instance is 
   * equivalent to this instance's header
   * @throws UnassignedDatasetException if instance doesn't have access to any
   * dataset
   */
  //@ requires m_Dataset != null;
  public /*@pure@*/ boolean equalHeaders(Instance inst) {

    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    }
    return m_Dataset.equalHeaders(inst.dataset());
  }

  /**
   * Checks if the headers of two instances are equivalent. 
   * If not, then returns a message why they differ.
   *
   * @param dataset 	another instance
   * @return 		null if the header of the given instance is equivalent 
   * 			to this instance's header, otherwise a message with details on
   * 			why they differ
   */
  public String equalHeadersMsg(Instance inst) {
    if (m_Dataset == null)
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");

    return m_Dataset.equalHeadersMsg(inst.dataset());
  }

  /**
   * Tests whether an instance has a missing value. Skips the class attribute if set.
   * @return true if instance has a missing value.
   * @throws UnassignedDatasetException if instance doesn't have access to any
   * dataset
   */
  //@ requires m_Dataset != null;
  public /*@pure@*/ boolean hasMissingValue() {
    
    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    }
    for (int i = 0; i < numValues(); i++) {
      if (index(i) != classIndex()) {
	if (isMissingSparse(i)) {
	  return true;
	}
      }
    }
    return false;
  }

  /**
   * Inserts an attribute at the given position (0 to 
   * numAttributes()). Only succeeds if the instance does not
   * have access to any dataset because otherwise inconsistencies
   * could be introduced.
   *
   * @param position the attribute's position
   * @throws RuntimeException if the instance has accesss to a
   * dataset
   * @throws IllegalArgumentException if the position is out of range
   */
  //@ requires m_Dataset == null;
  //@ requires 0 <= position && position <= numAttributes();
  public void insertAttributeAt(int position) {

    if (m_Dataset != null) {
      throw new RuntimeException("DenseInstance has accesss to a dataset!");
    }
    if ((position < 0) ||
	(position > numAttributes())) {
      throw new IllegalArgumentException("Can't insert attribute: index out "+
                                         "of range");
    }
    forceInsertAttributeAt(position);
  }

  /**
   * Tests if a specific value is "missing".
   *
   * @param attIndex the attribute's index
   * @return true if the value is "missing"
   */
  public /*@pure@*/ boolean isMissing(int attIndex) {

    if (Utils.isMissingValue(value(attIndex))) {
      return true;
    }
    return false;
  }

  /**
   * Tests if a specific value is "missing", given
   * an index in the sparse representation.
   *
   * @param indexOfIndex the index of the attribute's index 
   * @return true if the value is "missing"
   */
  public /*@pure@*/ boolean isMissingSparse(int indexOfIndex) {

    if (Utils.isMissingValue(valueSparse(indexOfIndex))) {
      return true;
    }
    return false;
  }

  /**
   * Tests if a specific value is "missing".
   * The given attribute has to belong to a dataset.
   *
   * @param att the attribute
   * @return true if the value is "missing"
   */
  public /*@pure@*/ boolean isMissing(Attribute att) {

    return isMissing(att.index());
  }

  /**
   * Returns the number of class labels.
   *
   * @return the number of class labels as an integer if the 
   * class attribute is nominal, 1 otherwise.
   * @throws UnassignedDatasetException if instance doesn't have access to any
   * dataset
   */
  //@ requires m_Dataset != null;
  public /*@pure@*/ int numClasses() {
    
    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    }
    return m_Dataset.numClasses();
  }

  /**
   * Sets the class value of an instance to be "missing". A deep copy of
   * the vector of attribute values is performed before the
   * value is set to be missing.
   *
   * @throws UnassignedClassException if the class is not set
   * @throws UnassignedDatasetException if the instance doesn't
   * have access to a dataset
   */
  //@ requires classIndex() >= 0;
  public void setClassMissing() {

    if (classIndex() < 0) {
      throw new UnassignedClassException("Class is not set!");
    }
    setMissing(classIndex());
  }

  /**
   * Sets the class value of an instance to the given value (internal
   * floating-point format).  A deep copy of the vector of attribute
   * values is performed before the value is set.
   *
   * @param value the new attribute value (If the corresponding
   * attribute is nominal (or a string) then this is the new value's
   * index as a double).  
   * @throws UnassignedClassException if the class is not set
   * @throws UnaddignedDatasetException if the instance doesn't
   * have access to a dataset 
   */
  //@ requires classIndex() >= 0;
  public void setClassValue(double value) {

    if (classIndex() < 0) {
      throw new UnassignedClassException("Class is not set!");
    }
    setValue(classIndex(), value);
  }

  /**
   * Sets the class value of an instance to the given value. A deep
   * copy of the vector of attribute values is performed before the
   * value is set.
   *
   * @param value the new class value (If the class
   * is a string attribute and the value can't be found,
   * the value is added to the attribute).
   * @throws UnassignedClassException if the class is not set
   * @throws UnassignedDatasetException if the dataset is not set
   * @throws IllegalArgumentException if the attribute is not
   * nominal or a string, or the value couldn't be found for a nominal
   * attribute 
   */
  //@ requires classIndex() >= 0;
  public final void setClassValue(String value) {

    if (classIndex() < 0) {
      throw new UnassignedClassException("Class is not set!");
    }
    setValue(classIndex(), value);
  }

  /**
   * Sets the reference to the dataset. Does not check if the instance
   * is compatible with the dataset. Note: the dataset does not know
   * about this instance. If the structure of the dataset's header
   * gets changed, this instance will not be adjusted automatically.
   *
   * @param instances the reference to the dataset 
   */
  public final void setDataset(Instances instances) {
    
    m_Dataset = instances;
  }

  /**
   * Sets a specific value to be "missing". Performs a deep copy
   * of the vector of attribute values before the value is set to
   * be missing.
   *
   * @param attIndex the attribute's index
   */
  public final void setMissing(int attIndex) {

    setValue(attIndex, Utils.missingValue());
  }

  /**
   * Sets a specific value to be "missing". Performs a deep copy
   * of the vector of attribute values before the value is set to
   * be missing. The given attribute has to belong to a dataset.
   *
   * @param att the attribute
   */
  public final void setMissing(Attribute att) {

    setMissing(att.index());
  }

  /**
   * Sets a value of a nominal or string attribute to the given
   * value. Performs a deep copy of the vector of attribute values
   * before the value is set.
   *
   * @param attIndex the attribute's index
   * @param value the new attribute value (If the attribute
   * is a string attribute and the value can't be found,
   * the value is added to the attribute).
   * @throws UnassignedDatasetException if the dataset is not set
   * @throws IllegalArgumentException if the selected
   * attribute is not nominal or a string, or the supplied value couldn't 
   * be found for a nominal attribute 
   */
  //@ requires m_Dataset != null;
  public final void setValue(int attIndex, String value) {
    
    int valIndex;

    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    }
    if (!attribute(attIndex).isNominal()&& !attribute(attIndex).isRanking() &&
	!attribute(attIndex).isString()) {
      throw new IllegalArgumentException("Attribute neither nominal nor string!");
    }
    valIndex = attribute(attIndex).indexOfValue(value);
    if (valIndex == -1) {
    //RANKING BEGIN
      if (attribute(attIndex).isNominal() || attribute(attIndex).isRanking()) {
    	  if(!(attribute(attIndex)instanceof PreferenceAttribute || m_Dataset.instance(0) instanceof PreferenceDenseInstance))
	throw new IllegalArgumentException("Value not defined for given nominal attribute!");
	//RANKING END
      } 
    else {
	attribute(attIndex).forceAddValue(value);
	valIndex = attribute(attIndex).indexOfValue(value);
      }
    }
    setValue(attIndex, (double)valIndex); 
  }

  /**
   * Sets a specific value in the instance to the given value
   * (internal floating-point format). Performs a deep copy of the
   * vector of attribute values before the value is set, so if you are
   * planning on calling setValue many times it may be faster to
   * create a new instance using toDoubleArray.  The given attribute
   * has to belong to a dataset.
   *
   * @param att the attribute 
   * @param value the new attribute value (If the corresponding
   * attribute is nominal (or a string) then this is the new value's
   * index as a double).  
   */
  public final void setValue(Attribute att, double value) {

    setValue(att.index(), value);
  }

  /**
   * Sets a value of an nominal or string attribute to the given
   * value. Performs a deep copy of the vector of attribute values
   * before the value is set, so if you are planning on calling setValue many
   * times it may be faster to create a new instance using toDoubleArray.
   * The given attribute has to belong to a dataset.
   *
   * @param att the attribute
   * @param value the new attribute value (If the attribute
   * is a string attribute and the value can't be found,
   * the value is added to the attribute).
   * @throws IllegalArgumentException if the the attribute is not
   * nominal or a string, or the value couldn't be found for a nominal
   * attribute 
   */
  public final void setValue(Attribute att, String value) {

    if (!att.isNominal() && !att.isRanking() &&
	!att.isString()) {
      throw new IllegalArgumentException("Attribute neither nominal nor string!");
    }
    int valIndex = att.indexOfValue(value);
    if (valIndex == -1) {
      if (att.isNominal() || att.isRanking()) {
	throw new IllegalArgumentException("Value not defined for given nominal attribute!");
      } else {
	att.forceAddValue(value);
	valIndex = att.indexOfValue(value);
      }
    }
    setValue(att.index(), (double)valIndex);
  }
  
  /**
   * Sets the weight of an instance.
   *
   * @param weight the weight
   */
  public final void setWeight(double weight) {

    m_Weight = weight;
  }

  /** 
   * Returns the relational value of a relational attribute.
   *
   * @param attIndex the attribute's index
   * @return the corresponding relation as an Instances object
   * @throws IllegalArgumentException if the attribute is not a
   * relation-valued attribute
   * @throws UnassignedDatasetException if the instance doesn't belong
   * to a dataset.
   */
  //@ requires m_Dataset != null;
  public final /*@pure@*/ Instances relationalValue(int attIndex) {

    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    } 
    return relationalValue(m_Dataset.attribute(attIndex));
  }


  /** 
   * Returns the relational value of a relational attribute.
   *
   * @param att the attribute
   * @return the corresponding relation as an Instances object
   * @throws IllegalArgumentException if the attribute is not a
   * relation-valued attribute
   * @throws UnassignedDatasetException if the instance doesn't belong
   * to a dataset.
   */
  public final /*@pure@*/ Instances relationalValue(Attribute att) {

    int attIndex = att.index();
    if (att.isRelationValued()) {
      return att.relation((int) value(attIndex));
    } else {
      throw new IllegalArgumentException("Attribute isn't relation-valued!");
    }
  }

  /** 
   * Returns the value of a nominal, string, date, or relational attribute
   * for the instance as a string.
   *
   * @param attIndex the attribute's index
   * @return the value as a string
   * @throws IllegalArgumentException if the attribute is not a nominal,
   * string, date, or relation-valued attribute.
   * @throws UnassignedDatasetException if the instance doesn't belong
   * to a dataset.
   */
  //@ requires m_Dataset != null;
  public final /*@pure@*/ String stringValue(int attIndex) {

    if (m_Dataset == null) {
      throw new UnassignedDatasetException("DenseInstance doesn't have access to a dataset!");
    } 
    return stringValue(m_Dataset.attribute(attIndex));
  }


  /** 
   * Returns the value of a nominal, string, date, or relational attribute
   * for the instance as a string.
   *
   * @param att the attribute
   * @return the value as a string
   * @throws IllegalArgumentException if the attribute is not a nominal,
   * string, date, or relation-valued attribute.
   * @throws UnassignedDatasetException if the instance doesn't belong
   * to a dataset.
   */
  public final /*@pure@*/ String stringValue(Attribute att) {

    int attIndex = att.index();
    switch (att.type()) {
    case Attribute.NOMINAL:
    case PreferenceAttribute.RANKING:
    case Attribute.STRING:
      return att.value((int) value(attIndex));
    case Attribute.DATE:
      return att.formatDate(value(attIndex));
    case Attribute.RELATIONAL:
      return att.relation((int) value(attIndex)).stringWithoutHeader();
    default:
      throw new IllegalArgumentException("Attribute isn't nominal, string or date!");
    }
  }

  /**
   * Returns the description of one instance. If the instance
   * doesn't have access to a dataset, it returns the internal
   * floating-point values. Quotes string
   * values that contain whitespace characters.
   *
   * @return the instance's description as a string
   */
  public String toString() {

    StringBuffer text = new StringBuffer(toStringNoWeight());

    if (m_Weight != 1.0) {
      text.append(",{" + Utils.doubleToString(m_Weight, 6) + "}");
    }

    return text.toString();
  }

  /**
   * Returns the description of one value of the instance as a 
   * string. If the instance doesn't have access to a dataset, it 
   * returns the internal floating-point value. Quotes string
   * values that contain whitespace characters, or if they
   * are a question mark.
   *
   * @param attIndex the attribute's index
   * @return the value's description as a string
   */
  public final /*@pure@*/ String toString(int attIndex) {

   StringBuffer text = new StringBuffer();
   
   if (isMissing(attIndex)) {
     text.append("?");
   } else {
     if (m_Dataset == null) {
       text.append(Utils.doubleToString(value(attIndex),6));
     } else {
       switch (m_Dataset.attribute(attIndex).type()) {
       case Attribute.NOMINAL:
       case PreferenceAttribute.RANKING:
       case Attribute.STRING:
       case Attribute.DATE:
       case Attribute.RELATIONAL:
         text.append(Utils.quote(stringValue(attIndex)));
         break;
       case Attribute.NUMERIC:
	 text.append(Utils.doubleToString(value(attIndex),6));
         break;
       default:
         throw new IllegalStateException("Unknown attribute type");
       }
     }
   }
   return text.toString();
  }

  /**
   * Returns the description of one value of the instance as a 
   * string. If the instance doesn't have access to a dataset it 
   * returns the internal floating-point value. Quotes string
   * values that contain whitespace characters, or if they
   * are a question mark.
   * The given attribute has to belong to a dataset.
   *
   * @param att the attribute
   * @return the value's description as a string
   */
  public final String toString(Attribute att) {
   
   return toString(att.index());
  }

  /**
   * Returns an instance's attribute value in internal format.
   * The given attribute has to belong to a dataset.
   *
   * @param att the attribute
   * @return the specified value as a double (If the corresponding
   * attribute is nominal (or a string) then it returns the value's index as a
   * double).
   */
  public /*@pure@*/ double value(Attribute att) {

    return value(att.index());
  }

  /**
   * Returns an instance's attribute value in internal format, given
   * an index in the sparse representation.
   *
   * @param indexOfIndex the index of the attribute's index
   * @return the specified value as a double (If the corresponding
   * attribute is nominal (or a string) then it returns the value's index as a 
   * double).
   */
  public /*@pure@*/ double valueSparse(int indexOfIndex) {

    return m_AttValues[indexOfIndex];
  }  

  /**
   * Returns the instance's weight.
   *
   * @return the instance's weight as a double
   */
  public final /*@pure@*/ double weight() {

    return m_Weight;
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
   * Deletes an attribute at the given position (0 to 
   * numAttributes() - 1).
   *
   * @param position the attribute's position
   */
  protected abstract void forceDeleteAttributeAt(int position);

  /**
   * Inserts an attribute at the given position
   * (0 to numAttributes()) and sets its value to be missing. 
   *
   * @param position the attribute's position
   */
  protected abstract void forceInsertAttributeAt(int position);
}
