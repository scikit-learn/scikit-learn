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
 *    Instance.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.util.Enumeration;

/**
 * Interface representing an instance. All values (numeric, date,
 * nominal, string or relational) are internally stored as
 * floating-point numbers in the original concrete class
 * implementations (now called DenseInstance.java and
 * SparseInstance.java), and the methods in this interface reflect
 * this. If an attribute is nominal (or a string or relational), the
 * stored value is the index of the corresponding nominal (or string
 * or relational) value in the attribute's definition. We have chosen
 * this approach in favor of a more elegant object-oriented approach
 * because it is much faster. <p>
 *
 * Typical usage (code from the main() method of this class): <p>
 *
 * <code>
 * ... <br>
 *      
 * // Create empty instance with three attribute values <br>
 * Instance inst = new DenseInstance(3); <br><br>
 *     
 * // Set instance's values for the attributes "length", "weight", and "position"<br>
 * inst.setValue(length, 5.3); <br>
 * inst.setValue(weight, 300); <br>
 * inst.setValue(position, "first"); <br><br>
 *   
 * // Set instance's dataset to be the dataset "race" <br>
 * inst.setDataset(race); <br><br>
 *   
 * // Print the instance <br>
 * System.out.println("The instance: " + inst); <br>
 *
 * ... <br>
 * </code><p>
 *
 * All methods that change an instance's attribute values must be
 * safe, ie. a change of an instance's attribute values must not
 * affect any other instances.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5987 $ 
 */
public interface Instance extends Copyable {

  /**
   * Returns the attribute with the given index.
   *
   * @param index the attribute's index
   * @return the attribute at the given position
   * @throws UnassignedDatasetException if instance doesn't have access to a
   * dataset
   */ 
  public Attribute attribute(int index);

  /**
   * Returns the attribute with the given index in the sparse representation.
   * Same as attribute(int) for a DenseInstance.
   *
   * @param indexOfIndex the index of the attribute's index 
   * @return the attribute at the given position
   * @throws UnassignedDatasetException if instance doesn't have access to a
   * dataset
   */ 
  public Attribute attributeSparse(int indexOfIndex);

  /**
   * Returns class attribute.
   *
   * @return the class attribute
   * @throws UnassignedDatasetException if the class is not set or the
   * instance doesn't have access to a dataset
   */
  public Attribute classAttribute();

  /**
   * Returns the class attribute's index.
   *
   * @return the class index as an integer 
   * @throws UnassignedDatasetException if instance doesn't have access to a dataset 
   */
  public int classIndex();

  /**
   * Tests if an instance's class is missing.
   *
   * @return true if the instance's class is missing
   * @throws UnassignedClassException if the class is not set or the instance doesn't
   * have access to a dataset
   */
  public boolean classIsMissing();

  /**
   * Returns an instance's class value as a floating-point number.
   *
   * @return the corresponding value as a double (If the 
   * corresponding attribute is nominal (or a string) then it returns the 
   * value's index as a double).
   * @throws UnassignedClassException if the class is not set or the instance doesn't
   * have access to a dataset 
   */
  public double classValue();

  /**
   * Returns the dataset this instance has access to. (ie. obtains
   * information about attribute types from) Null if the instance
   * doesn't have access to a dataset.
   *
   * @return the dataset the instance has accesss to
   */
  public Instances dataset();

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
  public void deleteAttributeAt(int position);

  /**
   * Returns an enumeration of all the attributes.
   *
   * @return enumeration of all the attributes
   * @throws UnassignedDatasetException if the instance doesn't
   * have access to a dataset 
   */
  public Enumeration enumerateAttributes();

  /**
   * Tests if the headers of two instances are equivalent.
   *
   * @param inst another instance
   * @return true if the header of the given instance is 
   * equivalent to this instance's header
   * @throws UnassignedDatasetException if instance doesn't have access to any
   * dataset
   */
  public boolean equalHeaders(Instance inst);

  /**
   * Checks if the headers of two instances are equivalent. 
   * If not, then returns a message why they differ.
   *
   * @param dataset 	another instance
   * @return 		null if the header of the given instance is equivalent 
   * 			to this instance's header, otherwise a message with details on
   * 			why they differ
   */
  public String equalHeadersMsg(Instance inst);

  /**
   * Tests whether an instance has a missing value. Skips the class attribute if set.
   * @return true if instance has a missing value.
   * @throws UnassignedDatasetException if instance doesn't have access to any
   * dataset
   */
  public boolean hasMissingValue();

  /**
   * Returns the index of the attribute stored at the given position in the sparse
   * representation. Identify function for an instance of type DenseInstance.
   *
   * @param position the position 
   * @return the index of the attribute stored at the given position
   */
  public int index(int position);

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
  public void insertAttributeAt(int position);

  /**
   * Tests if a specific value is "missing".
   *
   * @param attIndex the attribute's index
   * @return true if the value is "missing"
   */
  public boolean isMissing(int attIndex);

  /**
   * Tests if a specific value is "missing" in the sparse
   * representation. Samse as isMissing(int) for a DenseInstance.
   *
   * @param indexOfIndex the index of the attribute's index 
   * @return true if the value is "missing"
   */
  public boolean isMissingSparse(int indexOfIndex);

  /**
   * Tests if a specific value is "missing".
   * The given attribute has to belong to a dataset.
   *
   * @param att the attribute
   * @return true if the value is "missing"
   */
  public boolean isMissing(Attribute att);

  /**
   * Merges this instance with the given instance and returns
   * the result. Dataset is set to null. The returned instance
   * is of the same type as this instance.
   *
   * @param inst the instance to be merged with this one
   * @return the merged instances
   */
  public Instance mergeInstance(Instance inst);

  /**
   * Returns the number of attributes.
   *
   * @return the number of attributes as an integer
   */
  public int numAttributes();

  /**
   * Returns the number of class labels.
   *
   * @return the number of class labels as an integer if the 
   * class attribute is nominal, 1 otherwise.
   * @throws UnassignedDatasetException if instance doesn't have access to any
   * dataset
   */
  public int numClasses();

  /**
   * Returns the number of values present in a sparse representation.
   *
   * @return the number of values
   */
  public int numValues();

  /** 
   * Replaces all missing values in the instance with the
   * values contained in the given array. A deep copy of
   * the vector of attribute values is performed before the
   * values are replaced.
   *
   * @param array containing the means and modes
   * @throws IllegalArgumentException if numbers of attributes are unequal
   */
  public void replaceMissingValues(double[] array);

  /**
   * Sets the class value of an instance to be "missing". A deep copy of
   * the vector of attribute values is performed before the
   * value is set to be missing.
   *
   * @throws UnassignedClassException if the class is not set
   * @throws UnassignedDatasetException if the instance doesn't
   * have access to a dataset
   */
  public void setClassMissing();

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
  public void setClassValue(double value);

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
  public void setClassValue(String value);

  /**
   * Sets the reference to the dataset. Does not check if the instance
   * is compatible with the dataset. Note: the dataset does not know
   * about this instance. If the structure of the dataset's header
   * gets changed, this instance will not be adjusted automatically.
   *
   * @param instances the reference to the dataset 
   */
  public void setDataset(Instances instances);

  /**
   * Sets a specific value to be "missing". Performs a deep copy
   * of the vector of attribute values before the value is set to
   * be missing.
   *
   * @param attIndex the attribute's index
   */
  public void setMissing(int attIndex);

  /**
   * Sets a specific value to be "missing". Performs a deep copy
   * of the vector of attribute values before the value is set to
   * be missing. The given attribute has to belong to a dataset.
   *
   * @param att the attribute
   */
  public void setMissing(Attribute att);

  /**
   * Sets a specific value in the instance to the given value 
   * (internal floating-point format). Performs a deep copy
   * of the vector of attribute values before the value is set.
   *
   * @param attIndex the attribute's index 
   * @param value the new attribute value (If the corresponding
   * attribute is nominal (or a string) then this is the new value's
   * index as a double).  
   */
  public void setValue(int attIndex, double value);

  /**
   * Sets a specific value in the instance to the given value
   * (internal floating-point format), given an index into the sparse
   * representation. Performs a deep copy of the vector of attribute
   * values before the value is set. Same as setValue(int, double)
   * for a DenseInstance.
   *
   * @param indexOfIndex the index of the attribute's index 
   * @param value the new attribute value (If the corresponding
   * attribute is nominal (or a string) then this is the new value's
   * index as a double).  
   */
  public void setValueSparse(int indexOfIndex, double value);

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
  public void setValue(int attIndex, String value);

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
  public void setValue(Attribute att, double value);

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
  public void setValue(Attribute att, String value);
  
  /**
   * Sets the weight of an instance.
   *
   * @param weight the weight
   */
  public void setWeight(double weight);

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
  public Instances relationalValue(int attIndex);


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
  public Instances relationalValue(Attribute att);

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
  public String stringValue(int attIndex);

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
  public String stringValue(Attribute att);

  /**
   * Returns the values of each attribute as an array of doubles.
   *
   * @return an array containing all the instance attribute values
   */
  public double[] toDoubleArray();

  /**
   * Returns the description of one instance (without weight
   * appended). If the instance
   * doesn't have access to a dataset, it returns the internal
   * floating-point values. Quotes string
   * values that contain whitespace characters.
   *
   * This method is used by getRandomNumberGenerator() in
   * Instances.java in order to maintain backwards compatibility
   * with weka 3.4.
   *
   * @return the instance's description as a string
   */
  public String toStringNoWeight();

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
  public String toString(int attIndex);

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
  public String toString(Attribute att);

  /**
   * Returns an instance's attribute value in internal format.
   *
   * @param attIndex the attribute's index
   * @return the specified value as a double (If the corresponding
   * attribute is nominal (or a string) then it returns the value's index as a 
   * double).
   */
  public double value(int attIndex);

  /**
   * Returns an instance's attribute value in internal format, given
   * an index in the sparse representation. Same as value(int) for
   * a DenseInstance.
   *
   * @param indexOfIndex the index of the attribute's index
   * @return the specified value as a double (If the corresponding
   * attribute is nominal (or a string) then it returns the value's index as a 
   * double).
   */
  public double valueSparse(int indexOfIndex);

  /**
   * Returns an instance's attribute value in internal format.
   * The given attribute has to belong to a dataset.
   *
   * @param att the attribute
   * @return the specified value as a double (If the corresponding
   * attribute is nominal (or a string) then it returns the value's index as a
   * double).
   */
  public double value(Attribute att);

  /**
   * Returns the instance's weight.
   *
   * @return the instance's weight as a double
   */
  public double weight();
}
