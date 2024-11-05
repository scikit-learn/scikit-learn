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
 *    SparseInstance.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.ArrayList;

/**
 * Class for storing an instance as a sparse vector. A sparse instance
 * only requires storage for those attribute values that are non-zero.
 * Since the objective is to reduce storage requirements for datasets
 * with large numbers of default values, this also includes nominal
 * attributes -- the first nominal value (i.e. that which has index 0)
 * will not require explicit storage, so rearrange your nominal attribute
 * value orderings if necessary. Missing values will be stored
 * explicitly.
 *
 * @author Eibe Frank
 * @version $Revision: 5987 $
 */
public class SparseInstance extends AbstractInstance {

  /** for serialization */
  private static final long serialVersionUID = -3579051291332630149L;

  /** The index of the attribute associated with each stored value. */
  protected int[] m_Indices;

  /** The maximum number of values that can be stored. */
  protected int m_NumAttributes;

  /**
   * Constructor that generates a sparse instance from the given
   * instance. Reference to the dataset is set to null.
   * (ie. the instance doesn't have access to information about the
   * attribute types)
   *
   * @param instance the instance from which the attribute values
   * and the weight are to be copied
   */
  public SparseInstance(Instance instance) {
    
    m_Weight = instance.weight();
    m_Dataset = null;
    m_NumAttributes = instance.numAttributes();
    if (instance instanceof SparseInstance) {
      m_AttValues = ((SparseInstance)instance).m_AttValues;
      m_Indices = ((SparseInstance)instance).m_Indices;
    } else {
      double[] tempValues = new double[instance.numAttributes()];
      int[] tempIndices = new int[instance.numAttributes()];
      int vals = 0;
      for (int i = 0; i < instance.numAttributes(); i++) {
	if (instance.value(i) != 0) {
	  tempValues[vals] = instance.value(i);
	  tempIndices[vals] = i;
	  vals++;
	}
      }
      m_AttValues = new double[vals];
      m_Indices = new int[vals];
      System.arraycopy(tempValues, 0, m_AttValues, 0, vals);
      System.arraycopy(tempIndices, 0, m_Indices, 0, vals);
    }
  }

  /**
   * Constructor that copies the info from the given instance. 
   * Reference to the dataset is set to null.
   * (ie. the instance doesn't have access to information about the
   * attribute types)
   *
   * @param instance the instance from which the attribute
   * info is to be copied 
   */
  public SparseInstance(SparseInstance instance) {
    
    m_AttValues = instance.m_AttValues;
    m_Indices = instance.m_Indices;
    m_Weight = instance.m_Weight;
    m_NumAttributes = instance.m_NumAttributes;
    m_Dataset = null;
  }

  /**
   * Constructor that generates a sparse instance from the given
   * parameters. Reference to the dataset is set to null.
   * (ie. the instance doesn't have access to information about the
   * attribute types)
   *
   * @param weight the instance's weight
   * @param attValues a vector of attribute values 
   */
  public SparseInstance(double weight, double[] attValues) {
    
    m_Weight = weight;
    m_Dataset = null;
    m_NumAttributes = attValues.length;
    double[] tempValues = new double[m_NumAttributes];
    int[] tempIndices = new int[m_NumAttributes];
    int vals = 0;
    for (int i = 0; i < m_NumAttributes; i++) {
      if (attValues[i] != 0) {
	tempValues[vals] = attValues[i];
	tempIndices[vals] = i;
	vals++;
      }
    }
    m_AttValues = new double[vals];
    m_Indices = new int[vals];
    System.arraycopy(tempValues, 0, m_AttValues, 0, vals);
    System.arraycopy(tempIndices, 0, m_Indices, 0, vals);
  }
  
  /**
   * Constructor that inititalizes instance variable with given
   * values. Reference to the dataset is set to null. (ie. the instance
   * doesn't have access to information about the attribute types)
   * Note that the indices need to be sorted in ascending order. Otherwise
   * things won't work properly.
   *
   * @param weight the instance's weight
   * @param attValues a vector of attribute values (just the ones to be stored)
   * @param indices the indices of the given values in the full vector (need to
   * be sorted in ascending order)
   * @param maxNumValues the maximium number of values that can be stored
   */
  public SparseInstance(double weight, double[] attValues,
			int[] indices, int maxNumValues){
    
    int vals = 0; 
    m_AttValues = new double [attValues.length];
    m_Indices = new int [indices.length];
    for (int i = 0; i < attValues.length; i++) {
      if (attValues[i] != 0) {
        m_AttValues[vals] = attValues[i];
        m_Indices[vals] = indices[i];
        vals++;
      }
    }
    if (vals != attValues.length) {
      // Need to truncate.
      double [] newVals = new double[vals];
      System.arraycopy(m_AttValues, 0, newVals, 0, vals);
      m_AttValues = newVals;
      int [] newIndices = new int[vals];
      System.arraycopy(m_Indices, 0, newIndices, 0, vals);
      m_Indices = newIndices;
    }
    m_Weight = weight;
    m_NumAttributes = maxNumValues;
    m_Dataset = null;
  }

  /**
   * Constructor of an instance that sets weight to one, all values to
   * be missing, and the reference to the dataset to null. (ie. the instance
   * doesn't have access to information about the attribute types)
   *
   * @param numAttributes the size of the instance 
   */
  public SparseInstance(int numAttributes) {
    
    m_AttValues = new double[numAttributes];
    m_NumAttributes = numAttributes;
    m_Indices = new int[numAttributes];
    for (int i = 0; i < m_AttValues.length; i++) {
      m_AttValues[i] = Utils.missingValue();
      m_Indices[i] = i;
    }
    m_Weight = 1;
    m_Dataset = null;
  }

  /**
   * Produces a shallow copy of this instance. The copy has
   * access to the same dataset. (if you want to make a copy
   * that doesn't have access to the dataset, use 
   * <code>new SparseInstance(instance)</code>
   *
   * @return the shallow copy
   */
  public Object copy() {

    SparseInstance result = new SparseInstance(this);
    result.m_Dataset = m_Dataset;
    return result;
  }

  /**
   * Returns the index of the attribute stored at the given position.
   *
   * @param position the position 
   * @return the index of the attribute stored at the given position
   */
  public int index(int position) {

    return m_Indices[position];
  }

  /**
   * Locates the greatest index that is not greater than the
   * given index.
   *
   * @return the internal index of the attribute index. Returns
   * -1 if no index with this property could be found
   */
  public int locateIndex(int index) {

    int min = 0, max = m_Indices.length - 1;

    if (max == -1) {
      return -1;
    }

    // Binary search
    while ((m_Indices[min] <= index) && (m_Indices[max] >= index)) {
      int current = (max + min) / 2;
      if (m_Indices[current] > index) {
	max = current - 1;
      } else if (m_Indices[current] < index) {
	min = current + 1;
      } else {
	return current;
      }
    }
    if (m_Indices[max] < index) {
      return max;
    } else {
      return min - 1;
    }
  }

  /**
   * Merges this instance with the given instance and returns
   * the result. Dataset is set to null.
   *
   * @param inst the instance to be merged with this one
   * @return the merged instances
   */
  public Instance mergeInstance(Instance inst) {

    double[] values = new double[numValues() + inst.numValues()];
    int[] indices = new int[numValues() + inst.numValues()];

    int m = 0;
    for (int j = 0; j < numValues(); j++, m++) {
      values[m] = valueSparse(j);
      indices[m] = index(j);
    }
    for (int j = 0; j < inst.numValues(); j++, m++) {
      values[m] = inst.valueSparse(j);
      indices[m] = numAttributes() + inst.index(j);
    }
    
    return new SparseInstance(1.0, values, indices, numAttributes() +
			      inst.numAttributes());
  }

  /**
   * Returns the number of attributes.
   *
   * @return the number of attributes as an integer
   */
  public int numAttributes() {

    return m_NumAttributes;
  }

  /**
   * Returns the number of values in the sparse vector.
   *
   * @return the number of values
   */
  public int numValues() {

    return m_Indices.length;
  }

  /** 
   * Replaces all missing values in the instance with the 
   * values contained in the given array. A deep copy of
   * the vector of attribute values is performed before the
   * values are replaced.
   *
   * @param array containing the means and modes
   * @exception IllegalArgumentException if numbers of attributes are unequal
   */
  public void replaceMissingValues(double[] array) {
	 
    if ((array == null) || (array.length != m_NumAttributes)) {
      throw new IllegalArgumentException("Unequal number of attributes!");
    }
    double[] tempValues = new double[m_AttValues.length];
    int[] tempIndices = new int[m_AttValues.length];
    int vals = 0;
    for (int i = 0; i < m_AttValues.length; i++) {
      if (isMissingSparse(i)) {
	if (array[m_Indices[i]] != 0) {
	  tempValues[vals] = array[m_Indices[i]];
	  tempIndices[vals] = m_Indices[i];
	  vals++;
	} 
      } else {
	tempValues[vals] = m_AttValues[i];
	tempIndices[vals] = m_Indices[i];
	vals++;
      }
    }
    m_AttValues = new double[vals];
    m_Indices = new int[vals];
    System.arraycopy(tempValues, 0, m_AttValues, 0, vals);
    System.arraycopy(tempIndices, 0, m_Indices, 0, vals);
  }

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
  public void setValue(int attIndex, double value) {

    int index = locateIndex(attIndex);
    
    if ((index >= 0) && (m_Indices[index] == attIndex)) {
      if (value != 0) {
	double[] tempValues = new double[m_AttValues.length];
	System.arraycopy(m_AttValues, 0, tempValues, 0, m_AttValues.length);
	tempValues[index] = value;
	m_AttValues = tempValues;
      } else {
	double[] tempValues = new double[m_AttValues.length - 1];
	int[] tempIndices = new int[m_Indices.length - 1];
	System.arraycopy(m_AttValues, 0, tempValues, 0, index);
	System.arraycopy(m_Indices, 0, tempIndices, 0, index);
	System.arraycopy(m_AttValues, index + 1, tempValues, index, 
			 m_AttValues.length - index - 1);
	System.arraycopy(m_Indices, index + 1, tempIndices, index, 
			 m_Indices.length - index - 1);
	m_AttValues = tempValues;
	m_Indices = tempIndices;
      }
    } else {
      if (value != 0) {
	double[] tempValues = new double[m_AttValues.length + 1];
	int[] tempIndices = new int[m_Indices.length + 1];
	System.arraycopy(m_AttValues, 0, tempValues, 0, index + 1);
	System.arraycopy(m_Indices, 0, tempIndices, 0, index + 1);
	tempIndices[index + 1] = attIndex;
	tempValues[index + 1] = value;
	System.arraycopy(m_AttValues, index + 1, tempValues, index + 2, 
			 m_AttValues.length - index - 1);
	System.arraycopy(m_Indices, index + 1, tempIndices, index + 2, 
			 m_Indices.length - index - 1);
	m_AttValues = tempValues;
	m_Indices = tempIndices;
      }
    }
  }

  /**
   * Sets a specific value in the instance to the given value 
   * (internal floating-point format). Performs a deep copy
   * of the vector of attribute values before the value is set.
   *
   * @param indexOfIndex the index of the attribute's index 
   * @param value the new attribute value (If the corresponding
   * attribute is nominal (or a string) then this is the new value's
   * index as a double).  
   */
  public void setValueSparse(int indexOfIndex, double value) {

    if (value != 0) {
      double[] tempValues = new double[m_AttValues.length];
      System.arraycopy(m_AttValues, 0, tempValues, 0, m_AttValues.length);
      m_AttValues = tempValues;
      m_AttValues[indexOfIndex] = value;
    } else {
      double[] tempValues = new double[m_AttValues.length - 1];
      int[] tempIndices = new int[m_Indices.length - 1];
      System.arraycopy(m_AttValues, 0, tempValues, 0, indexOfIndex);
      System.arraycopy(m_Indices, 0, tempIndices, 0, indexOfIndex);
      System.arraycopy(m_AttValues, indexOfIndex + 1, tempValues, indexOfIndex, 
		       m_AttValues.length - indexOfIndex - 1);
      System.arraycopy(m_Indices, indexOfIndex + 1, tempIndices, indexOfIndex, 
		       m_Indices.length - indexOfIndex - 1);
      m_AttValues = tempValues;
      m_Indices = tempIndices;
    }
  }
  
  /**
   * Returns the values of each attribute as an array of doubles.
   *
   * @return an array containing all the instance attribute values
   */
  public double[] toDoubleArray() {

    double[] newValues = new double[m_NumAttributes];
    for (int i = 0; i < m_AttValues.length; i++) {
      newValues[m_Indices[i]] = m_AttValues[i];
    }
    return newValues;
  }

  /**
   * Returns the description of one instance in sparse format. 
   * If the instance doesn't have access to a dataset, it returns the 
   * internal floating-point values. Quotes string values that contain 
   * whitespace characters.
   *
   * @return the instance's description as a string
   */
  public String toStringNoWeight() {

    StringBuffer text = new StringBuffer();
    
    text.append('{');
    for (int i = 0; i < m_Indices.length; i++) {
      if (i > 0) text.append(",");
      if (isMissingSparse(i)) {
	text.append(m_Indices[i] + " ?");
      } else {
	if (m_Dataset == null) {
	  text.append(m_Indices[i] + " " + 
		      Utils.doubleToString(m_AttValues[i],6));
	} else {
	  if (m_Dataset.attribute(m_Indices[i]).isNominal() || 
		  m_Dataset.attribute(m_Indices[i]).isRanking() ||
	      m_Dataset.attribute(m_Indices[i]).isString() ||
	      m_Dataset.attribute(m_Indices[i]).isDate()) {
	    try {
	      text.append(m_Indices[i] + " " +
		  Utils.quote(stringValue(m_Indices[i])));
	    } catch (Exception e) {
              e.printStackTrace();
              System.err.println(new Instances(m_Dataset, 0));
              System.err.println("Att:" + m_Indices[i] + " Val:" + valueSparse(i));
	      throw new Error("This should never happen!");
	    }
	  } else if (m_Dataset.attribute(m_Indices[i]).isRelationValued()) {
	    try {
	      text.append(m_Indices[i] + " " +
			  Utils.quote(m_Dataset.attribute(m_Indices[i]).
				      relation((int)valueSparse(i)).
                                      stringWithoutHeader()));
            } catch (Exception e) {
              e.printStackTrace();
              System.err.println(new Instances(m_Dataset, 0));
              System.err.println("Att:" + m_Indices[i] + " Val:" + valueSparse(i));
	      throw new Error("This should never happen!");
	    }
	  } else {
	    text.append(m_Indices[i] + " " +
			Utils.doubleToString(m_AttValues[i],6));
	  }
	}
      }
    }
    text.append('}');

    return text.toString();
  }

  /**
   * Returns an instance's attribute value in internal format.
   *
   * @param attIndex the attribute's index
   * @return the specified value as a double (If the corresponding
   * attribute is nominal (or a string) then it returns the value's index as a 
   * double).
   */
  public double value(int attIndex) {

    int index = locateIndex(attIndex);
    if ((index >= 0) && (m_Indices[index] == attIndex)) {
      return m_AttValues[index];
    } else {
      return 0.0;
    }
  }  

  /**
   * Deletes an attribute at the given position (0 to 
   * numAttributes() - 1).
   *
   * @param pos the attribute's position
   */
  protected void forceDeleteAttributeAt(int position) {

    int index = locateIndex(position);

    m_NumAttributes--;
    if ((index >= 0) && (m_Indices[index] == position)) {
      int[] tempIndices = new int[m_Indices.length - 1];
      double[] tempValues = new double[m_AttValues.length - 1];
      System.arraycopy(m_Indices, 0, tempIndices, 0, index);
      System.arraycopy(m_AttValues, 0, tempValues, 0, index);
      for (int i = index; i < m_Indices.length - 1; i++) {
	tempIndices[i] = m_Indices[i + 1] - 1;
	tempValues[i] = m_AttValues[i + 1];
      }
      m_Indices = tempIndices;
      m_AttValues = tempValues;
    } else {
      int[] tempIndices = new int[m_Indices.length];
      double[] tempValues = new double[m_AttValues.length];
      System.arraycopy(m_Indices, 0, tempIndices, 0, index + 1);
      System.arraycopy(m_AttValues, 0, tempValues, 0, index + 1);
      for (int i = index + 1; i < m_Indices.length; i++) {
	tempIndices[i] = m_Indices[i] - 1;
	tempValues[i] = m_AttValues[i];
      }
      m_Indices = tempIndices;
      m_AttValues = tempValues;
    }
  }

  /**
   * Inserts an attribute at the given position
   * (0 to numAttributes()) and sets its value to be missing. 
   *
   * @param pos the attribute's position
   */
  protected void forceInsertAttributeAt(int position)  {

    int index = locateIndex(position);

    m_NumAttributes++;
    if ((index >= 0) && (m_Indices[index] == position)) {
      int[] tempIndices = new int[m_Indices.length + 1];
      double[] tempValues = new double[m_AttValues.length + 1];
      System.arraycopy(m_Indices, 0, tempIndices, 0, index);
      System.arraycopy(m_AttValues, 0, tempValues, 0, index);
      tempIndices[index] = position;
      tempValues[index] = Utils.missingValue();
      for (int i = index; i < m_Indices.length; i++) {
	tempIndices[i + 1] = m_Indices[i] + 1;
	tempValues[i + 1] = m_AttValues[i];
      }
      m_Indices = tempIndices;
      m_AttValues = tempValues;
    } else {
      int[] tempIndices = new int[m_Indices.length + 1];
      double[] tempValues = new double[m_AttValues.length + 1];
      System.arraycopy(m_Indices, 0, tempIndices, 0, index + 1);
      System.arraycopy(m_AttValues, 0, tempValues, 0, index + 1);
      tempIndices[index + 1] = position;
      tempValues[index + 1] = Utils.missingValue();
      for (int i = index + 1; i < m_Indices.length; i++) {
	tempIndices[i + 1] = m_Indices[i] + 1;
	tempValues[i + 1] = m_AttValues[i];
      }
      m_Indices = tempIndices;
      m_AttValues = tempValues;
    }
  }

  /**
   * Constructor for sub classes.
   */
  protected SparseInstance() {};

  /**
   * Main method for testing this class.
   */
  public static void main(String[] options) {

    try {

      // Create numeric attributes "length" and "weight"
      Attribute length = new Attribute("length");
      Attribute weight = new Attribute("weight");
      
      // Create vector to hold nominal values "first", "second", "third" 
      ArrayList<String> my_nominal_values = new ArrayList<String>(3); 
      my_nominal_values.add("first"); 
      my_nominal_values.add("second"); 
      my_nominal_values.add("third"); 
      
      // Create nominal attribute "position" 
      Attribute position = new Attribute("position", my_nominal_values);
      
      // Create vector of the above attributes 
      ArrayList<Attribute> attributes = new ArrayList<Attribute>(3);
      attributes.add(length);
      attributes.add(weight);
      attributes.add(position);
      
      // Create the empty dataset "race" with above attributes
      Instances race = new Instances("race", attributes, 0);
      
      // Make position the class attribute
      race.setClassIndex(position.index());
      
      // Create empty instance with three attribute values
      SparseInstance inst = new SparseInstance(3);
      
      // Set instance's values for the attributes "length", "weight", and "position"
      inst.setValue(length, 5.3);
      inst.setValue(weight, 300);
      inst.setValue(position, "first");
      
      // Set instance's dataset to be the dataset "race"
      inst.setDataset(race);
      
      // Print the instance
      System.out.println("The instance: " + inst);
      
      // Print the first attribute
      System.out.println("First attribute: " + inst.attribute(0));
      
      // Print the class attribute
      System.out.println("Class attribute: " + inst.classAttribute());
      
      // Print the class index
      System.out.println("Class index: " + inst.classIndex());
      
      // Say if class is missing
      System.out.println("Class is missing: " + inst.classIsMissing());
      
      // Print the instance's class value in internal format
      System.out.println("Class value (internal format): " + inst.classValue());
      
      // Print a shallow copy of this instance
      SparseInstance copy = (SparseInstance) inst.copy();
      System.out.println("Shallow copy: " + copy);
      
      // Set dataset for shallow copy
      copy.setDataset(inst.dataset());
      System.out.println("Shallow copy with dataset set: " + copy);

      // Print out all values in internal format
      System.out.print("All stored values in internal format: ");
      for (int i = 0; i < inst.numValues(); i++) {
	if (i > 0) {
	  System.out.print(",");
	}
	System.out.print(inst.valueSparse(i));
      }
      System.out.println();

      // Set all values to zero
      System.out.print("All values set to zero: ");
      while (inst.numValues() > 0) {
	inst.setValueSparse(0, 0);
      }
      for (int i = 0; i < inst.numValues(); i++) {
	if (i > 0) {
	  System.out.print(",");
	}
	System.out.print(inst.valueSparse(i));
      }
      System.out.println();

      // Set all values to one
      System.out.print("All values set to one: ");
      for (int i = 0; i < inst.numAttributes(); i++) {
	inst.setValue(i, 1);
      }
      for (int i = 0; i < inst.numValues(); i++) {
	if (i > 0) {
	  System.out.print(",");
	}
	System.out.print(inst.valueSparse(i));
      }
      System.out.println();

      // Unset dataset for copy, delete first attribute, and insert it again
      copy.setDataset(null);
      copy.deleteAttributeAt(0);
      copy.insertAttributeAt(0);
      copy.setDataset(inst.dataset());
      System.out.println("Copy with first attribute deleted and inserted: " + copy); 

      // Same for second attribute
      copy.setDataset(null);
      copy.deleteAttributeAt(1);
      copy.insertAttributeAt(1);
      copy.setDataset(inst.dataset());
      System.out.println("Copy with second attribute deleted and inserted: " + copy); 

      // Same for last attribute
      copy.setDataset(null);
      copy.deleteAttributeAt(2);
      copy.insertAttributeAt(2);
      copy.setDataset(inst.dataset());
      System.out.println("Copy with third attribute deleted and inserted: " + copy); 
      
      // Enumerate attributes (leaving out the class attribute)
      System.out.println("Enumerating attributes (leaving out class):");
      Enumeration enu = inst.enumerateAttributes();
      while (enu.hasMoreElements()) {
	Attribute att = (Attribute) enu.nextElement();
	System.out.println(att);
      }
      
      // Headers are equivalent?
      System.out.println("Header of original and copy equivalent: " +
			 inst.equalHeaders(copy));

      // Test for missing values
      System.out.println("Length of copy missing: " + copy.isMissing(length));
      System.out.println("Weight of copy missing: " + copy.isMissing(weight.index()));
      System.out.println("Length of copy missing: " + 
			 Utils.isMissingValue(copy.value(length)));

      // Prints number of attributes and classes
      System.out.println("Number of attributes: " + copy.numAttributes());
      System.out.println("Number of classes: " + copy.numClasses());

      // Replace missing values
      double[] meansAndModes = {2, 3, 0};
      copy.replaceMissingValues(meansAndModes);
      System.out.println("Copy with missing value replaced: " + copy);

      // Setting and getting values and weights
      copy.setClassMissing();
      System.out.println("Copy with missing class: " + copy);
      copy.setClassValue(0);
      System.out.println("Copy with class value set to first value: " + copy);
      copy.setClassValue("third");
      System.out.println("Copy with class value set to \"third\": " + copy);
      copy.setMissing(1);
      System.out.println("Copy with second attribute set to be missing: " + copy);
      copy.setMissing(length);
      System.out.println("Copy with length set to be missing: " + copy);
      copy.setValue(0, 0);
      System.out.println("Copy with first attribute set to 0: " + copy);
      copy.setValue(weight, 1);
      System.out.println("Copy with weight attribute set to 1: " + copy);
      copy.setValue(position, "second");
      System.out.println("Copy with position set to \"second\": " + copy);
      copy.setValue(2, "first");
      System.out.println("Copy with last attribute set to \"first\": " + copy);
      System.out.println("Current weight of instance copy: " + copy.weight());
      copy.setWeight(2);
      System.out.println("Current weight of instance copy (set to 2): " + copy.weight());
      System.out.println("Last value of copy: " + copy.toString(2));
      System.out.println("Value of position for copy: " + copy.toString(position));
      System.out.println("Last value of copy (internal format): " + copy.value(2));
      System.out.println("Value of position for copy (internal format): " + 
			 copy.value(position));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }
}
