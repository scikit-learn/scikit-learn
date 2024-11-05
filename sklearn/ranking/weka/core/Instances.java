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
 *    Instances.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import weka.core.converters.ArffLoader.ArffReader;
import weka.core.converters.ConverterUtils.DataSource;
import weka.core.labelranking.PreferenceAttribute;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.core.labelranking.RankUtilities;

/**
 * Class for handling an ordered set of weighted instances. <p>
 *
 * Typical usage: <p>
 * <pre>
 * import weka.core.converters.ConverterUtils.DataSource;
 * ...
 * 
 * // Read all the instances in the file (ARFF, CSV, XRFF, ...)
 * DataSource source = new DataSource(filename);
 * Instances instances = source.getDataSet();
 *
 * // Make the last attribute be the class
 * instances.setClassIndex(instances.numAttributes() - 1);
 * 
 * // Print header and instances.
 * System.out.println("\nDataset:\n");
 * System.out.println(instances);
 * 
 * ...
 * </pre><p>
 *
 * All methods that change a set of instances are safe, ie. a change
 * of a set of instances does not affect any other sets of
 * instances. All methods that change a datasets's attribute
 * information clone the dataset before it is changed.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6354 $ 
 */
public class Instances extends AbstractList<Instance>
  implements Serializable, RevisionHandler {
	//RANKING BEGIN
	//Two dimensional string arrays are included. First position: the label name,
	//second position: the label index.
	public List<String[]> labels = new ArrayList<String[]>();
	//sum of all pairwise ranking weights, used for calculating a dotplot.
	public int[][] weightMatrix;
	//sum of all pairwise rankings, used for calculating a dotplot.
	public int[][] rankMatrix;
    //RANKING END
	
  /** for serialization */
  static final long serialVersionUID = -19412345060742748L;
  
  /** The filename extension that should be used for arff files */
  public final static String FILE_EXTENSION = ".arff";

  /** The filename extension that should be used for bin. serialized instances files */
  public final static String SERIALIZED_OBJ_FILE_EXTENSION = ".bsi";

  /** The keyword used to denote the start of an arff header */
  public final static String ARFF_RELATION = "@relation";

  /** The keyword used to denote the start of the arff data section */
  public final static String ARFF_DATA = "@data";

  /** The dataset's name. */
  protected /*@spec_public non_null@*/ String m_RelationName;         

  /** The attribute information. */
  protected /*@spec_public non_null@*/ ArrayList<Attribute> m_Attributes;
  /*  public invariant (\forall int i; 0 <= i && i < m_Attributes.size(); 
                    m_Attributes.get(i) != null);
  */

  /** The instances. */
  protected /*@spec_public non_null@*/ ArrayList<Instance> m_Instances;

  /** The class attribute's index */
  protected int m_ClassIndex;
  //@ protected invariant classIndex() == m_ClassIndex;

  /** The lines read so far in case of incremental loading. Since the 
   * StreamTokenizer will be re-initialized with every instance that is read,
   * we have to keep track of the number of lines read so far. 
   * @see #readInstance(Reader) */
  protected int m_Lines = 0;
  
  /**
   * Reads an ARFF file from a reader, and assigns a weight of
   * one to each instance. Lets the index of the class 
   * attribute be undefined (negative).
   *
   * @param reader the reader
   * @throws IOException if the ARFF file is not read 
   * successfully
   */
  public Instances(/*@non_null@*/Reader reader) throws IOException {
    ArffReader arff = new ArffReader(reader);
    Instances dataset = arff.getData();
    initialize(dataset, dataset.numInstances());
    dataset.copyInstances(0, this, dataset.numInstances());
    compactify();
  }
 
  /**
   * Reads the header of an ARFF file from a reader and 
   * reserves space for the given number of instances. Lets
   * the class index be undefined (negative).
   *
   * @param reader the reader
   * @param capacity the capacity
   * @throws IllegalArgumentException if the header is not read successfully
   * or the capacity is negative.
   * @throws IOException if there is a problem with the reader.
   * @deprecated instead of using this method in conjunction with the
   * <code>readInstance(Reader)</code> method, one should use the 
   * <code>ArffLoader</code> or <code>DataSource</code> class instead.
   * @see weka.core.converters.ArffLoader
   * @see weka.core.converters.ConverterUtils.DataSource
   */
  //@ requires capacity >= 0;
  //@ ensures classIndex() == -1;
  @Deprecated public Instances(/*@non_null@*/Reader reader, int capacity)
    throws IOException {
//	labels = new ArrayList<String[]>();
    ArffReader arff = new ArffReader(reader, 0);
    Instances header = arff.getStructure();
    initialize(header, capacity);
    m_Lines = arff.getLineNo();
  }

  /**
   * Constructor copying all instances and references to
   * the header information from the given set of instances.
   *
   * @param dataset the set to be copied
   */
  public Instances(/*@non_null@*/Instances dataset) {

    this(dataset, dataset.numInstances());
//    labels = new ArrayList<String[]>();
    this.setLabels(dataset.getLabels());
    dataset.copyInstances(0, this, dataset.numInstances());
  }

  /**
   * Constructor creating an empty set of instances. Copies references
   * to the header information from the given set of instances. Sets
   * the capacity of the set of instances to 0 if its negative.
   *
   * @param dataset the instances from which the header 
   * information is to be taken
   * @param capacity the capacity of the new dataset 
   */
  public Instances(/*@non_null@*/Instances dataset, int capacity) {
    initialize(dataset, capacity);
    this.setLabels(dataset.getLabels());
//    labels = new ArrayList<String[]>();
  }

  /**
   * initializes with the header information of the given dataset and sets
   * the capacity of the set of instances.
   * 
   * @param dataset the dataset to use as template
   * @param capacity the number of rows to reserve
   */
  protected void initialize(Instances dataset, int capacity) {
    if (capacity < 0)
      capacity = 0;
    
    // Strings only have to be "shallow" copied because
    // they can't be modified.
    m_ClassIndex   = dataset.m_ClassIndex;
    m_RelationName = dataset.m_RelationName;
    m_Attributes   = dataset.m_Attributes;
    this.setLabels(dataset.getLabels());
    m_Instances    = new ArrayList<Instance>(capacity);
    
  }
  
  /**
   * Creates a new set of instances by copying a 
   * subset of another set.
   *
   * @param source the set of instances from which a subset 
   * is to be created
   * @param first the index of the first instance to be copied
   * @param toCopy the number of instances to be copied
   * @throws IllegalArgumentException if first and toCopy are out of range
   */
  //@ requires 0 <= first;
  //@ requires 0 <= toCopy;
  //@ requires first + toCopy <= source.numInstances();
  public Instances(/*@non_null@*/Instances source, int first, int toCopy) {
    
    this(source, toCopy);
//    labels = new ArrayList<String[]>();
//    this.setLabels(source.getLabels(););
    if ((first < 0) || ((first + toCopy) > source.numInstances())) {
      throw new IllegalArgumentException("Parameters first and/or toCopy out "+
                                         "of range");
    }
    source.copyInstances(first, this, toCopy);
    this.setLabels(source.getLabels());
  }

  /**
   * Creates an empty set of instances. Uses the given
   * attribute information. Sets the capacity of the set of 
   * instances to 0 if its negative. Given attribute information
   * must not be changed after this constructor has been used.
   *
   * @param name the name of the relation
   * @param attInfo the attribute information
   * @param capacity the capacity of the set
   * @throws IllegalArgumentException if attribute names are not unique
   */
  public Instances(/*@non_null@*/String name, 
		   /*@non_null@*/ArrayList<Attribute> attInfo, int capacity) {

    // check whether the attribute names are unique
    HashSet<String> names = new HashSet<String>();
//    labels = new ArrayList<String[]>();
    StringBuffer nonUniqueNames = new StringBuffer();
    for (Attribute att: attInfo) {
      if (names.contains(att.name())) {
        nonUniqueNames.append("'" + att.name() +"' ");
      }
      names.add(att.name());
    }
    if (names.size() != attInfo.size())
      throw new IllegalArgumentException("Attribute names are not unique!" +
      		" Causes: " + nonUniqueNames.toString());
    names.clear();
    
    m_RelationName = name;
    m_ClassIndex = -1;
    m_Attributes = attInfo;
    for (int i = 0; i < numAttributes(); i++) {
      attribute(i).setIndex(i);
    }
    m_Instances = new ArrayList<Instance>(capacity);
    
  }

  /**
   * Create a copy of the structure if the data has string or
   * relational attributes, "cleanses" string types (i.e. doesn't
   * contain references to the strings seen in the past) and all
   * relational attributes.
   *
   * @return a copy of the instance structure.
   */
  public Instances stringFreeStructure() {

    ArrayList<Attribute> newAtts = new ArrayList<Attribute>();
    for (int i = 0 ; i < m_Attributes.size(); i++) {
      Attribute att = (Attribute)m_Attributes.get(i);
      if (att.type() == Attribute.STRING) {
        newAtts.add(new Attribute(att.name(), (List<String>)null, i));
      } else if (att.type() == Attribute.RELATIONAL) {
        newAtts.add(new Attribute(att.name(), new Instances(att.relation(), 0), i));
      }
    }
    if (newAtts.size() == 0) {
      return new Instances(this, 0);
    }
    ArrayList<Attribute> atts = Utils.cast(m_Attributes.clone());
    for (int i = 0; i < newAtts.size(); i++) {
      atts.set(((Attribute)newAtts.get(i)).index(), newAtts.get(i));
    }
    Instances result = new Instances(this, 0);
    result.m_Attributes = atts;
    return result;
  }

  /**
   * Adds one instance to the end of the set. 
   * Shallow copies instance before it is added. Increases the
   * size of the dataset if it is not large enough. Does not
   * check if the instance is compatible with the dataset.
   * Note: String or relational values are not transferred.
   *
   * @param instance the instance to be added
   */
  public boolean add(/*@non_null@*/ Instance instance) {

    Instance newInstance = (Instance)instance.copy();

    newInstance.setDataset(this);
    m_Instances.add(newInstance);

    return true;
  }

  /**
   * Adds one instance to the end of the set. 
   * Shallow copies instance before it is added. Increases the
   * size of the dataset if it is not large enough. Does not
   * check if the instance is compatible with the dataset.
   * Note: String or relational values are not transferred.
   *
   * @param index position where instance is to be inserted
   * @param instance the instance to be added
   */
  //@ requires 0 <= index;
  //@ requires index < m_Instances.size();
  public void add(int index, /*@non_null@*/ Instance instance) {

    Instance newInstance = (Instance)instance.copy();

    newInstance.setDataset(this);
    m_Instances.add(index, newInstance);
  }

  /**
   * Returns an attribute.
   *
   * @param index the attribute's index (index starts with 0)
   * @return the attribute at the given position
   */
  //@ requires 0 <= index;
  //@ requires index < m_Attributes.size();
  //@ ensures \result != null;
  public /*@pure@*/ Attribute attribute(int index) {
    
    return (Attribute) m_Attributes.get(index);
  }

  /**
   * Returns an attribute given its name. If there is more than
   * one attribute with the same name, it returns the first one.
   * Returns null if the attribute can't be found.
   *
   * @param name the attribute's name
   * @return the attribute with the given name, null if the
   * attribute can't be found
   */ 
  public /*@pure@*/ Attribute attribute(String name) {
    
    for (int i = 0; i < numAttributes(); i++) {
      if (attribute(i).name().equals(name)) {
	return attribute(i);
      }
    }
    return null;
  }

  /**
   * Checks for attributes of the given type in the dataset
   *
   * @param attType  the attribute type to look for
   * @return         true if attributes of the given type are present
   */
  public boolean checkForAttributeType(int attType) {
    
    int i = 0;
    
    while (i < m_Attributes.size()) {
      if (attribute(i++).type() == attType) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks for string attributes in the dataset
   *
   * @return true if string attributes are present, false otherwise
   */
  public /*@pure@*/ boolean checkForStringAttributes() {
    return checkForAttributeType(Attribute.STRING);
  }

  /**
   * Checks if the given instance is compatible
   * with this dataset. Only looks at the size of
   * the instance and the ranges of the values for 
   * nominal and string attributes.
   *
   * @param instance the instance to check
   * @return true if the instance is compatible with the dataset 
   */
  public /*@pure@*/ boolean checkInstance(Instance instance) {

    if (instance.numAttributes() != numAttributes()) {
      return false;
    }
    for (int i = 0; i < numAttributes(); i++) {
      if (instance.isMissing(i)) {
	continue;
      } else if (attribute(i).isNominal() || attribute(i).isRanking()||
		 attribute(i).isString()) {
	if (!(Utils.eq(instance.value(i),
		       (double)(int)instance.value(i)))) {
	  return false;
	} else if (Utils.sm(instance.value(i), 0) ||
		   Utils.gr(instance.value(i),
			    attribute(i).numValues())) {
	  return false;
	}
      }
    }
    return true;
  }
	
  /**
   * Returns the class attribute.
   *
   * @return the class attribute
   * @throws UnassignedClassException if the class is not set
   */
  //@ requires classIndex() >= 0;
  public /*@pure@*/ Attribute classAttribute() {

    if (m_ClassIndex < 0) {
      throw new UnassignedClassException("Class index is negative (not set)!");
    }
    return attribute(m_ClassIndex);
  }

  /**
   * Returns the class attribute's index. Returns negative number
   * if it's undefined.
   *
   * @return the class index as an integer
   */
  // ensures \result == m_ClassIndex;
  public /*@pure@*/ int classIndex() {
    
    return m_ClassIndex;
  }
 
  /**
   * Compactifies the set of instances. Decreases the capacity of
   * the set so that it matches the number of instances in the set.
   */
  public void compactify() {

    m_Instances.trimToSize();
  }

  /**
   * Removes all instances from the set.
   */
  public void delete() {
    
    m_Instances = new ArrayList<Instance>();
  }

  /**
   * Removes an instance at the given position from the set.
   *
   * @param index the instance's position (index starts with 0)
   */
  //@ requires 0 <= index && index < numInstances();
  public void delete(int index) {
    
    m_Instances.remove(index);
  }

  /**
   * Deletes an attribute at the given position 
   * (0 to numAttributes() - 1). A deep copy of the attribute
   * information is performed before the attribute is deleted.
   *
   * @param position the attribute's position (position starts with 0)
   * @throws IllegalArgumentException if the given index is out of range 
   *            or the class attribute is being deleted
   */
  //@ requires 0 <= position && position < numAttributes();
  //@ requires position != classIndex();
  public void deleteAttributeAt(int position) {
	 
    if ((position < 0) || (position >= m_Attributes.size())) {
      throw new IllegalArgumentException("Index out of range");
    }
    if (position == m_ClassIndex) {
      throw new IllegalArgumentException("Can't delete class attribute");
    }
    freshAttributeInfo();
    if (m_ClassIndex > position) {
      m_ClassIndex--;
    }
    m_Attributes.remove(position);
    for (int i = position; i < m_Attributes.size(); i++) {
      Attribute current = (Attribute)m_Attributes.get(i);
      current.setIndex(current.index() - 1);
    }
    for (int i = 0; i < numInstances(); i++) {
      instance(i).setDataset(null);
      instance(i).deleteAttributeAt(position); 
      instance(i).setDataset(this);
    }
  }

  /**
   * Deletes all attributes of the given type in the dataset. A deep copy of 
   * the attribute information is performed before an attribute is deleted.
   *
   * @param attType the attribute type to delete
   * @throws IllegalArgumentException if attribute couldn't be 
   * successfully deleted (probably because it is the class attribute).
   */
  public void deleteAttributeType(int attType) {
    int i = 0;
    while (i < m_Attributes.size()) {
      if (attribute(i).type() == attType) {
        deleteAttributeAt(i);
      } else {
        i++;
      }
    }
  }

  /**
   * Deletes all string attributes in the dataset. A deep copy of the attribute
   * information is performed before an attribute is deleted.
   *
   * @throws IllegalArgumentException if string attribute couldn't be 
   * successfully deleted (probably because it is the class attribute).
   * @see #deleteAttributeType(int)
   */
  public void deleteStringAttributes() {
    deleteAttributeType(Attribute.STRING);
  }

  /**
   * Removes all instances with missing values for a particular
   * attribute from the dataset.
   *
   * @param attIndex the attribute's index (index starts with 0)
   */
  //@ requires 0 <= attIndex && attIndex < numAttributes();
  public void deleteWithMissing(int attIndex) {

    ArrayList<Instance> newInstances = new ArrayList<Instance>(numInstances());

    for (int i = 0; i < numInstances(); i++) {
      if (!instance(i).isMissing(attIndex)) {
	newInstances.add(instance(i));
      }
    }
    m_Instances = newInstances;
  }

  /**
   * Removes all instances with missing values for a particular
   * attribute from the dataset.
   *
   * @param att the attribute
   */
  public void deleteWithMissing(/*@non_null@*/ Attribute att) {

    deleteWithMissing(att.index());
  }

  /**
   * Removes all instances with a missing class value
   * from the dataset.
   *
   * @throws UnassignedClassException if class is not set
   */
  public void deleteWithMissingClass() {

    if (m_ClassIndex < 0) {
      throw new UnassignedClassException("Class index is negative (not set)!");
    }
    deleteWithMissing(m_ClassIndex);
  }

  /**
   * Returns an enumeration of all the attributes.
   *
   * @return enumeration of all the attributes.
   */
  public /*@non_null pure@*/ Enumeration enumerateAttributes() {

    return new WekaEnumeration(m_Attributes, m_ClassIndex);
  }

  /**
   * Returns an enumeration of all instances in the dataset.
   *
   * @return enumeration of all instances in the dataset
   */
  public /*@non_null pure@*/ Enumeration enumerateInstances() {

    return new WekaEnumeration(m_Instances);
  }

  /**
   * Checks if two headers are equivalent. If not, then returns a message why
   * they differ.
   *
   * @param dataset 	another dataset
   * @return 		null if the header of the given dataset is equivalent 
   * 			to this header, otherwise a message with details on
   * 			why they differ
   */
  public String equalHeadersMsg(Instances dataset) {
    // Check class and all attributes
    if (m_ClassIndex != dataset.m_ClassIndex)
      return "Class index differ: " + (m_ClassIndex+1) + " != " + (dataset.m_ClassIndex+1);

    if (m_Attributes.size() != dataset.m_Attributes.size())
      return "Different number of attributes: " + m_Attributes.size() + " != " + dataset.m_Attributes.size();
    
    for (int i = 0; i < m_Attributes.size(); i++) {
      String msg = attribute(i).equalsMsg(dataset.attribute(i));
      if (msg != null)
	return "Attributes differ at position " + (i+1) + ":\n" + msg;
    }
    
    return null;
  }

  /**
   * Checks if two headers are equivalent.
   *
   * @param dataset another dataset
   * @return true if the header of the given dataset is equivalent 
   * to this header
   */
  public /*@pure@*/ boolean equalHeaders(Instances dataset){
    return (equalHeadersMsg(dataset) == null);
  }
 
  /**
   * Returns the first instance in the set.
   *
   * @return the first instance in the set
   */
  //@ requires numInstances() > 0;
  public /*@non_null pure@*/ Instance firstInstance() {
    
    return (Instance)m_Instances.get(0);
  }

  /**
   * Returns a random number generator. The initial seed of the random
   * number generator depends on the given seed and the hash code of
   * a string representation of a instances chosen based on the given
   * seed. 
   *
   * @param seed the given seed
   * @return the random number generator
   */
  public Random getRandomNumberGenerator(long seed) {

    Random r = new Random(seed);
    r.setSeed(instance(r.nextInt(numInstances())).toStringNoWeight().hashCode() + seed);
    return r;
  }
 
  /**
   * Inserts an attribute at the given position (0 to 
   * numAttributes()) and sets all values to be missing.
   * Shallow copies the attribute before it is inserted, and performs
   * a deep copy of the existing attribute information.
   *
   * @param att the attribute to be inserted
   * @param position the attribute's position (position starts with 0)
   * @throws IllegalArgumentException if the given index is out of range
   */
  //@ requires 0 <= position;
  //@ requires position <= numAttributes();
  public void insertAttributeAt(/*@non_null@*/ Attribute att, int position) {
	 
    if ((position < 0) ||
	(position > m_Attributes.size())) {
      throw new IllegalArgumentException("Index out of range");
    }
    att = (Attribute)att.copy();
    freshAttributeInfo();
    att.setIndex(position);
    m_Attributes.add(position, att);
    for (int i = position + 1; i < m_Attributes.size(); i++) {
      Attribute current = (Attribute)m_Attributes.get(i);
      current.setIndex(current.index() + 1);
    }
    for (int i = 0; i < numInstances(); i++) {
      instance(i).setDataset(null);
      instance(i).insertAttributeAt(position);
      instance(i).setDataset(this);
    }
    if (m_ClassIndex >= position) {
      m_ClassIndex++;
    }
  }

  /**
   * Returns the instance at the given position.
   *
   * @param index the instance's index (index starts with 0)
   * @return the instance at the given position
   */
  //@ requires 0 <= index;
  //@ requires index < numInstances();
  public /*@non_null pure@*/ Instance instance(int index) {

    return m_Instances.get(index);
  }

  /**
   * Returns the instance at the given position.
   *
   * @param index the instance's index (index starts with 0)
   * @return the instance at the given position
   */
  //@ requires 0 <= index;
  //@ requires index < numInstances();
  public /*@non_null pure@*/ Instance get(int index) {

    return m_Instances.get(index);
  }

  /**
   * Returns the kth-smallest attribute value of a numeric attribute.
   * Note that calling this method will change the order of the data!
   *
   * @param att the Attribute object
   * @param k the value of k
   * @return the kth-smallest value
   */
  public double kthSmallestValue(Attribute att, int k) {

    return kthSmallestValue(att.index(), k);
  }

  /**
   * Returns the kth-smallest attribute value of a numeric attribute.
   * Note that calling this method will change the order of the data!
   * The number of non-missing values in the data must be as least
   * as last as k for this to work.
   *
   * @param attIndex the attribute's index
   * @param k the value of k
   * @return the kth-smallest value
   */
  public double kthSmallestValue(int attIndex, int k) {
    
    if (!attribute(attIndex).isNumeric()) {
      throw new IllegalArgumentException("Instances: attribute must be numeric to compute kth-smallest value.");
    }

    int i,j;

    // move all instances with missing values to end
    j = numInstances() - 1;
    i = 0;
    while (i <= j) {
      if (instance(j).isMissing(attIndex)) {
	j--;
      } else {
	if (instance(i).isMissing(attIndex)) {
	  swap(i,j);
	  j--;
	}
	i++;
      }
    }

    if ((k < 1) || (k > j+1)) {
      throw new IllegalArgumentException("Instances: value for k for computing kth-smallest value too large.");
    }

    return instance(select(attIndex, 0, j, k)).value(attIndex);
  }

  /**
   * Returns the last instance in the set.
   *
   * @return the last instance in the set
   */
  //@ requires numInstances() > 0;
  public /*@non_null pure@*/ Instance lastInstance() {
    
    return (Instance)m_Instances.get(m_Instances.size() - 1);
  }

  /**
   * Returns the mean (mode) for a numeric (nominal) attribute as
   * a floating-point value. Returns 0 if the attribute is neither nominal nor 
   * numeric. If all values are missing it returns zero.
   *
   * @param attIndex the attribute's index (index starts with 0)
   * @return the mean or the mode
   */
  public /*@pure@*/ double meanOrMode(int attIndex) {

    double result, found;
    int [] counts;

    if (attribute(attIndex).isNumeric()) {
      result = found = 0;
      for (int j = 0; j < numInstances(); j++) {
	if (!instance(j).isMissing(attIndex)) {
	  found += instance(j).weight();
	  result += instance(j).weight()*instance(j).value(attIndex);
	}
      }
      if (found <= 0) {
	return 0;
      } else {
	return result / found;
      }
    } else if (attribute(attIndex).isNominal() || attribute(attIndex).isRanking()) {
      counts = new int[attribute(attIndex).numValues()];
      for (int j = 0; j < numInstances(); j++) {
	if (!instance(j).isMissing(attIndex)) {
	  counts[(int) instance(j).value(attIndex)] += instance(j).weight();
	}
      }
      return (double)Utils.maxIndex(counts);
    } else {
      return 0;
    }
  }

  /**
   * Returns the mean (mode) for a numeric (nominal) attribute as a
   * floating-point value.  Returns 0 if the attribute is neither
   * nominal nor numeric.  If all values are missing it returns zero.
   *
   * @param att the attribute
   * @return the mean or the mode 
   */
  public /*@pure@*/ double meanOrMode(Attribute att) {

    return meanOrMode(att.index());
  }

  /**
   * Returns the number of attributes.
   *
   * @return the number of attributes as an integer
   */
  //@ ensures \result == m_Attributes.size();
  public /*@pure@*/ int numAttributes() {

    return m_Attributes.size();
  }

  /**
   * Returns the number of class labels.
   *
   * @return the number of class labels as an integer if the class 
   * attribute is nominal, 1 otherwise.
   * @throws UnassignedClassException if the class is not set
   */
  //@ requires classIndex() >= 0;
  public /*@pure@*/ int numClasses() {
    
    if (m_ClassIndex < 0) {
      throw new UnassignedClassException("Class index is negative (not set)!");
    }
    if (!classAttribute().isNominal() && !classAttribute().isRanking()) {
      return 1;
    } else {
      return classAttribute().numValues();
    }
  }

  /**
   * Returns the number of distinct values of a given attribute.
   * Returns the number of instances if the attribute is a
   * string attribute. The value 'missing' is not counted.
   *
   * @param attIndex the attribute (index starts with 0)
   * @return the number of distinct values of a given attribute
   */
  //@ requires 0 <= attIndex;
  //@ requires attIndex < numAttributes();
  public /*@pure@*/ int numDistinctValues(int attIndex) {

    if (attribute(attIndex).isNumeric()) {
      double [] attVals = attributeToDoubleArray(attIndex);
      int [] sorted = Utils.sort(attVals);
      double prev = 0;
      int counter = 0;
      for (int i = 0; i < sorted.length; i++) {
	Instance current = instance(sorted[i]);
	if (current.isMissing(attIndex)) {
	  break;
	}
	if ((i == 0) || 
	    (current.value(attIndex) > prev)) {
	  prev = current.value(attIndex);
	  counter++;
	}
      }
      return counter;
    } else {
      return attribute(attIndex).numValues();
    }
  }

  /**
   * Returns the number of distinct values of a given attribute.
   * Returns the number of instances if the attribute is a
   * string attribute. The value 'missing' is not counted.
   *
   * @param att the attribute
   * @return the number of distinct values of a given attribute
   */
  public /*@pure@*/ int numDistinctValues(/*@non_null@*/Attribute att) {

    return numDistinctValues(att.index());
  }
  
  /**
   * Returns the number of instances in the dataset.
   *
   * @return the number of instances in the dataset as an integer
   */
  //@ ensures \result == m_Instances.size();
  public /*@pure@*/ int numInstances() {

    return m_Instances.size();
  }
  
  /**
   * Returns the number of instances in the dataset.
   *
   * @return the number of instances in the dataset as an integer
   */
  //@ ensures \result == m_Instances.size();
  public /*@pure@*/ int size() {

    return m_Instances.size();
  }

  /**
   * Shuffles the instances in the set so that they are ordered 
   * randomly.
   *
   * @param random a random number generator
   */
  public void randomize(Random random) {

    for (int j = numInstances() - 1; j > 0; j--)
      swap(j, random.nextInt(j+1));
  }
  
  /**
   * Reads a single instance from the reader and appends it
   * to the dataset.  Automatically expands the dataset if it
   * is not large enough to hold the instance. This method does
   * not check for carriage return at the end of the line.
   *
   * @param reader the reader 
   * @return false if end of file has been reached
   * @throws IOException if the information is not read 
   * successfully
   * @deprecated instead of using this method in conjunction with the
   * <code>readInstance(Reader)</code> method, one should use the 
   * <code>ArffLoader</code> or <code>DataSource</code> class instead.
   * @see weka.core.converters.ArffLoader
   * @see weka.core.converters.ConverterUtils.DataSource
   */ 
  @Deprecated public boolean readInstance(Reader reader) throws IOException {

    ArffReader arff = new ArffReader(reader, this, m_Lines, 1);
    Instance inst = arff.readInstance(arff.getData(), false);
    m_Lines = arff.getLineNo();
    if (inst != null) {
      add(inst);
      return true;
    }
    else {
      return false;
    }
  }    

  /**
   * Returns the relation's name.
   *
   * @return the relation's name as a string
   */
  //@ ensures \result == m_RelationName;
  public /*@pure@*/ String relationName() {

    return m_RelationName;
  }

  /**
   * Removes the instance at the given position.
   *
   * @param index the instance's index (index starts with 0)
   * @return the instance at the given position
   */
  //@ requires 0 <= index;
  //@ requires index < numInstances();
  public Instance remove(int index) {

    return m_Instances.remove(index);
  }

  /**
   * Renames an attribute. This change only affects this
   * dataset.
   *
   * @param att the attribute's index (index starts with 0)
   * @param name the new name
   */
  public void renameAttribute(int att, String name) {

    Attribute newAtt = attribute(att).copy(name);
    ArrayList<Attribute> newVec = new ArrayList<Attribute>(numAttributes());

    for (int i = 0; i < numAttributes(); i++) {
      if (i == att) {
	newVec.add(newAtt);
      } else {
	newVec.add(attribute(i));
      }
    }
    m_Attributes = newVec;
  }

  /**
   * Renames an attribute. This change only affects this
   * dataset.
   *
   * @param att the attribute
   * @param name the new name
   */
  public void renameAttribute(Attribute att, String name) {

    renameAttribute(att.index(), name);
  }

  /**
   * Renames the value of a nominal (or string) attribute value. This
   * change only affects this dataset.
   *
   * @param att the attribute's index (index starts with 0)
   * @param val the value's index (index starts with 0)
   * @param name the new name 
   */
  public void renameAttributeValue(int att, int val, String name) {

    Attribute newAtt = (Attribute)attribute(att).copy();
    ArrayList<Attribute> newVec = new ArrayList<Attribute>(numAttributes());

    newAtt.setValue(val, name);
    for (int i = 0; i < numAttributes(); i++) {
      if (i == att) {
	newVec.add(newAtt);
      } else {
	newVec.add(attribute(i));
      }
    }
    m_Attributes = newVec;
  }

  /**
   * Renames the value of a nominal (or string) attribute value. This
   * change only affects this dataset.
   *
   * @param att the attribute
   * @param val the value
   * @param name the new name
   */
  public void renameAttributeValue(Attribute att, String val, 
                                         String name) {

    int v = att.indexOfValue(val);
    if (v == -1) throw new IllegalArgumentException(val + " not found");
    renameAttributeValue(att.index(), v, name);
  }

  /**
   * Creates a new dataset of the same size using random sampling
   * with replacement.
   *
   * @param random a random number generator
   * @return the new dataset
   */
  public Instances resample(Random random) {

    Instances newData = new Instances(this, numInstances());
    while (newData.numInstances() < numInstances()) {
      newData.add(instance(random.nextInt(numInstances())));
    }
    return newData;
  }

  /**
   * Creates a new dataset of the same size using random sampling
   * with replacement according to the current instance weights. The
   * weights of the instances in the new dataset are set to one.
   *
   * @param random a random number generator
   * @return the new dataset
   */
  public Instances resampleWithWeights(Random random) {

    double [] weights = new double[numInstances()];
    for (int i = 0; i < weights.length; i++) {
      weights[i] = instance(i).weight();
    }
    return resampleWithWeights(random, weights);
  }


  /**
   * Creates a new dataset of the same size using random sampling
   * with replacement according to the given weight vector. The
   * weights of the instances in the new dataset are set to one.
   * The length of the weight vector has to be the same as the
   * number of instances in the dataset, and all weights have to
   * be positive.
   *
   * @param random a random number generator
   * @param weights the weight vector
   * @return the new dataset
   * @throws IllegalArgumentException if the weights array is of the wrong
   * length or contains negative weights.
   */
  public Instances resampleWithWeights(Random random, 
					     double[] weights) {

    if (weights.length != numInstances()) {
      throw new IllegalArgumentException("weights.length != numInstances.");
    }
    Instances newData = new Instances(this, numInstances());
    if (numInstances() == 0) {
      return newData;
    }
    double[] probabilities = new double[numInstances()];
    double sumProbs = 0, sumOfWeights = Utils.sum(weights);
    for (int i = 0; i < numInstances(); i++) {
      sumProbs += random.nextDouble();
      probabilities[i] = sumProbs;
    }
    Utils.normalize(probabilities, sumProbs / sumOfWeights);

    // Make sure that rounding errors don't mess things up
    probabilities[numInstances() - 1] = sumOfWeights;
    int k = 0; int l = 0;
    sumProbs = 0;
    while ((k < numInstances() && (l < numInstances()))) {
      if (weights[l] < 0) {
	throw new IllegalArgumentException("Weights have to be positive.");
      }
      sumProbs += weights[l];
      while ((k < numInstances()) &&
	     (probabilities[k] <= sumProbs)) { 
	newData.add(instance(l));
	newData.instance(k).setWeight(1);
	k++;
      }
      l++;
    }
    return newData;
  }

  /**
   * Replaces the instance at the given position.
   * Shallow copies instance before it is added. Does not
   * check if the instance is compatible with the dataset.
   * Note: String or relational values are not transferred.
   *
   * @param index position where instance is to be inserted
   * @param instance the instance to be inserted
   * @return the instance previously at that position
   */
  //@ requires 0 <= index;
  //@ requires index < m_Instances.size();
  public Instance set(int index, /*@non_null@*/ Instance instance) {

    Instance newInstance = (Instance)instance.copy();
    Instance oldInstance = m_Instances.get(index);

    newInstance.setDataset(this);
    m_Instances.set(index, newInstance);

    return oldInstance;
  }

  /** 
   * Sets the class attribute.
   *
   * @param att attribute to be the class
   */
  public void setClass(Attribute att) {
    m_ClassIndex = att.index();
  }

  /** 
   * Sets the class index of the set.
   * If the class index is negative there is assumed to be no class.
   * (ie. it is undefined)
   *
   * @param classIndex the new class index (index starts with 0)
   * @throws IllegalArgumentException if the class index is too big or < 0
   */
  public void setClassIndex(int classIndex) {

    if (classIndex >= numAttributes()) {
      throw new IllegalArgumentException("Invalid class index: " + classIndex);
    }
    m_ClassIndex = classIndex;
  }

  /**
   * Sets the relation's name.
   *
   * @param newName the new relation name.
   */
  public void setRelationName(/*@non_null@*/String newName) {
    
    m_RelationName = newName;
  }

  /**
   * Sorts the instances based on an attribute. For numeric attributes, 
   * instances are sorted in ascending order. For nominal attributes, 
   * instances are sorted based on the attribute label ordering 
   * specified in the header. Instances with missing values for the 
   * attribute are placed at the end of the dataset.
   *
   * @param attIndex the attribute's index (index starts with 0)
   */
  public void sort(int attIndex) {

    int i,j;

    // move all instances with missing values to end
    j = numInstances() - 1;
    i = 0;
    while (i <= j) {
      if (instance(j).isMissing(attIndex)) {
	j--;
      } else {
	if (instance(i).isMissing(attIndex)) {
	  swap(i,j);
	  j--;
	}
	i++;
      }
    }
    quickSort(attIndex, 0, j);
  }

  /**
   * Sorts the instances based on an attribute. For numeric attributes, 
   * instances are sorted into ascending order. For nominal attributes, 
   * instances are sorted based on the attribute label ordering 
   * specified in the header. Instances with missing values for the 
   * attribute are placed at the end of the dataset.
   *
   * @param att the attribute
   */
  public void sort(Attribute att) {

    sort(att.index());
  }

  /**
   * Stratifies a set of instances according to its class values 
   * if the class attribute is nominal (so that afterwards a 
   * stratified cross-validation can be performed).
   *
   * @param numFolds the number of folds in the cross-validation
   * @throws UnassignedClassException if the class is not set
   */
  public void stratify(int numFolds) {
    
    if (numFolds <= 1) {
      throw new IllegalArgumentException("Number of folds must be greater than 1");
    }
    if (m_ClassIndex < 0) {
      throw new UnassignedClassException("Class index is negative (not set)!");
    }
    if (classAttribute().isNominal() || classAttribute().isRanking()) {

      // sort by class
      int index = 1;
      while (index < numInstances()) {
	Instance instance1 = instance(index - 1);
	for (int j = index; j < numInstances(); j++) {
	  Instance instance2 = instance(j);
	  if ((instance1.classValue() == instance2.classValue()) ||
	      (instance1.classIsMissing() && 
	       instance2.classIsMissing())) {
	    swap(index,j);
	    index++;
	  }
	}
	index++;
      }
      stratStep(numFolds);
    }
  }
 
  /**
   * Computes the sum of all the instances' weights.
   *
   * @return the sum of all the instances' weights as a double
   */
  public /*@pure@*/ double sumOfWeights() {
    
    double sum = 0;

    for (int i = 0; i < numInstances(); i++) {
      sum += instance(i).weight();
    }
    return sum;
  }

  /**
   * Creates the test set for one fold of a cross-validation on 
   * the dataset.
   *
   * @param numFolds the number of folds in the cross-validation. Must
   * be greater than 1.
   * @param numFold 0 for the first fold, 1 for the second, ...
   * @return the test set as a set of weighted instances
   * @throws IllegalArgumentException if the number of folds is less than 2
   * or greater than the number of instances.
   */
  //@ requires 2 <= numFolds && numFolds < numInstances();
  //@ requires 0 <= numFold && numFold < numFolds;
  public Instances testCV(int numFolds, int numFold) {

    int numInstForFold, first, offset;
    Instances test;
    
    if (numFolds < 2) {
      throw new IllegalArgumentException("Number of folds must be at least 2!");
    }
    if (numFolds > numInstances()) {
      throw new IllegalArgumentException("Can't have more folds than instances!");
    }
    numInstForFold = numInstances() / numFolds;
    if (numFold < numInstances() % numFolds){
      numInstForFold++;
      offset = numFold;
    }else
      offset = numInstances() % numFolds;
    test = new Instances(this, numInstForFold);
    first = numFold * (numInstances() / numFolds) + offset;
    copyInstances(first, test, numInstForFold);
    return test;
  }
 
  /**
   * Returns the dataset as a string in ARFF format. Strings
   * are quoted if they contain whitespace characters, or if they
   * are a question mark.
   *
   * @return the dataset in ARFF format as a string
   */
  public String toString() {
    
    StringBuffer text = new StringBuffer();
    
    text.append(ARFF_RELATION).append(" ").
      append(Utils.quote(m_RelationName)).append("\n\n");
    for (int i = 0; i < numAttributes(); i++) {
    	if(attribute(i).toString().contains("RANKING"))
    		text.append(attribute(i).toString().replace('\'', ' ')).append("\n");
    	else
    		text.append(attribute(i)).append("\n");
      
    }
    text.append("\n").append(ARFF_DATA).append("\n");

    text.append(stringWithoutHeader());
    return text.toString();
  }

  /**
   * Returns the instances in the dataset as a string in ARFF format. Strings
   * are quoted if they contain whitespace characters, or if they
   * are a question mark.
   *
   * @return the dataset in ARFF format as a string
   */
  protected String stringWithoutHeader() {
    
    StringBuffer text = new StringBuffer();

    for (int i = 0; i < numInstances(); i++) {
      //RANKING BEGIN
    	String rnk="";
    	
      if(instance(i) instanceof PreferenceDenseInstance){
    	  double[] ranking;
    	  PreferenceDenseInstance temp = (PreferenceDenseInstance) instance(i);
    	  for(int j=0; ;j++){
    		  if(temp.getHashMap().get(j)!=null){
    			  RankUtilities.noPartialRanking(temp.getHashMap().get(j));
    			  ranking = RankUtilities.triangleLabels;
    			  break;
    		  }
    	  }
    	  if(getLabels().size()==0){
    		  setLabels(RankUtilities.labels);
    	  }
    	  
    	  //modified by erez
    	  for(int o=0; o<ranking.length; o++) {
	    	  for(int k=0; k<ranking.length; k++){
	    		  if (ranking[k]==o) {
		    		  for(int l=0; l<getLabels().size(); l++){
		    			  if(k==Integer.parseInt(getLabels().get(l)[1])){
		    				  rnk+=getLabels().get(l)[0]+">";
		    			  }
		    		  }
	    		  }
	    	  }
    	  }
    	  rnk = rnk.substring(0, rnk.length()-1); 
    	 for(int x=0; x<instance(i).numAttributes()-1; x++){
    		 text.append(instance(i).toString(x));
    		 text.append(", ");
    	 }
    	 text.append(rnk);
      }
      //RANKING END
      else
    	  text.append(instance(i));
      
      
      if (i < numInstances() - 1) {
	text.append('\n');
      }
    }
    return text.toString();
  }

  /**
   * Creates the training set for one fold of a cross-validation 
   * on the dataset. 
   *
   * @param numFolds the number of folds in the cross-validation. Must
   * be greater than 1.
   * @param numFold 0 for the first fold, 1 for the second, ...
   * @return the training set 
   * @throws IllegalArgumentException if the number of folds is less than 2
   * or greater than the number of instances.
   */
  //@ requires 2 <= numFolds && numFolds < numInstances();
  //@ requires 0 <= numFold && numFold < numFolds;
  public Instances trainCV(int numFolds, int numFold) {

    int numInstForFold, first, offset;
    Instances train;
 
    if (numFolds < 2) {
      throw new IllegalArgumentException("Number of folds must be at least 2!");
    }
    if (numFolds > numInstances()) {
      throw new IllegalArgumentException("Can't have more folds than instances!");
    }
    numInstForFold = numInstances() / numFolds;
    if (numFold < numInstances() % numFolds) {
      numInstForFold++;
      offset = numFold;
    }else
      offset = numInstances() % numFolds;
    train = new Instances(this, numInstances() - numInstForFold);
    first = numFold * (numInstances() / numFolds) + offset;
    copyInstances(0, train, first);
    copyInstances(first + numInstForFold, train,
		  numInstances() - first - numInstForFold);

    return train;
  }

  /**
   * Creates the training set for one fold of a cross-validation 
   * on the dataset. The data is subsequently randomized based
   * on the given random number generator.
   *
   * @param numFolds the number of folds in the cross-validation. Must
   * be greater than 1.
   * @param numFold 0 for the first fold, 1 for the second, ...
   * @param random the random number generator
   * @return the training set 
   * @throws IllegalArgumentException if the number of folds is less than 2
   * or greater than the number of instances.
   */
  //@ requires 2 <= numFolds && numFolds < numInstances();
  //@ requires 0 <= numFold && numFold < numFolds;
  public Instances trainCV(int numFolds, int numFold, Random random) {

    Instances train = trainCV(numFolds, numFold);
    train.randomize(random);
    return train;
  }

  /**
   * Computes the variance for a numeric attribute.
   *
   * @param attIndex the numeric attribute (index starts with 0)
   * @return the variance if the attribute is numeric
   * @throws IllegalArgumentException if the attribute is not numeric
   */
  public /*@pure@*/ double variance(int attIndex) {
  
    double sum = 0, sumSquared = 0, sumOfWeights = 0;

    if (!attribute(attIndex).isNumeric()) {
      throw new IllegalArgumentException("Can't compute variance because attribute is " +
			  "not numeric!");
    }
    for (int i = 0; i < numInstances(); i++) {
      if (!instance(i).isMissing(attIndex)) {
	sum += instance(i).weight() * 
	  instance(i).value(attIndex);
	sumSquared += instance(i).weight() * 
	  instance(i).value(attIndex) *
	  instance(i).value(attIndex);
	sumOfWeights += instance(i).weight();
      }
    }
    if (sumOfWeights <= 1) {
      return 0;
    }
    double result = (sumSquared - (sum * sum / sumOfWeights)) / 
      (sumOfWeights - 1);

    // We don't like negative variance
    if (result < 0) {
      return 0;
    } else {
      return result;
    }
  }

  /**
   * Computes the variance for a numeric attribute.
   *
   * @param att the numeric attribute
   * @return the variance if the attribute is numeric
   * @throws IllegalArgumentException if the attribute is not numeric
   */
  public /*@pure@*/ double variance(Attribute att) {
    
    return variance(att.index());
  }
  
  /**
   * Calculates summary statistics on the values that appear in this
   * set of instances for a specified attribute.
   *
   * @param index the index of the attribute to summarize (index starts with 0)
   * @return an AttributeStats object with it's fields calculated.
   */
  //@ requires 0 <= index && index < numAttributes();
  public AttributeStats attributeStats(int index) {

    AttributeStats result = new AttributeStats();
    if (attribute(index).isNominal() || attribute(index).isRanking()) {
      result.nominalCounts = new int [attribute(index).numValues()];
      result.nominalWeights = new double[attribute(index).numValues()];
    }
    if (attribute(index).isNumeric()) {
      result.numericStats = new weka.experiment.Stats();
    }
    result.totalCount = numInstances();

    double [] attVals = attributeToDoubleArray(index);
    int [] sorted = Utils.sort(attVals);
    int currentCount = 0;
    double currentWeight = 0;
    double prev = Double.NaN;
    for (int j = 0; j < numInstances(); j++) {
      Instance current = instance(sorted[j]);
      if (current.isMissing(index)) {
	result.missingCount = numInstances() - j;
	break;
      }
      
      if(!(attribute(index) instanceof PreferenceAttribute)){
    	  if (current.value(index) == prev) {
    		  currentCount++;
    		  currentWeight += current.weight();
    	  } else {
    		  result.addDistinct(prev, currentCount, currentWeight);
    		  currentCount = 1;
    		  currentWeight = current.weight();
    		  prev = current.value(index);
    	  }
    
    	  result.addDistinct(prev, currentCount, currentWeight);
    	  result.distinctCount--; // So we don't count "missing" as a value
      }
      //For a PreferenceAttribute, we obtain the weights of pairwise rankings inside
      //of all preferences per instance.
      else{
    	  if(weightMatrix==null){
    		  weightMatrix = new int[labels.size()][labels.size()];
    	  }
    	  
    	  PreferenceDenseInstance pdi = (PreferenceDenseInstance) current;
    	  int[][] tmp;
    	  for(int c=0; ;c++){    		  
    		  tmp = pdi.getHashMap().get(c);
    		  if(tmp!=null)break;
    	  }
    	  
    	  for(int a=0; a<weightMatrix.length; a++){
    		  for(int b=0; b<weightMatrix.length; b++){
    			  if(tmp[a][b]==1){
    				  weightMatrix[a][b] += pdi.weight();
    			  }
    		  }
    	  }
      }
    }
    result.pairRankWeights = weightMatrix;
    weightMatrix=null;
    return result;
  }
  
  /**
   * Gets the value of all instances in this dataset for a particular
   * attribute. Useful in conjunction with Utils.sort to allow iterating
   * through the dataset in sorted order for some attribute.
   *
   * @param index the index of the attribute.
   * @return an array containing the value of the desired attribute for
   * each instance in the dataset. 
   */
  //@ requires 0 <= index && index < numAttributes();
  public /*@pure@*/ double [] attributeToDoubleArray(int index) {

    double [] result = new double[numInstances()];
    for (int i = 0; i < result.length; i++) {
      result[i] = instance(i).value(index);
    }
    return result;
  }

  /**
   * Generates a string summarizing the set of instances. Gives a breakdown
   * for each attribute indicating the number of missing/discrete/unique
   * values and other information.
   *
   * @return a string summarizing the dataset
   */
  public String toSummaryString() {

    StringBuffer result = new StringBuffer();
    result.append("Relation Name:  ").append(relationName()).append('\n');
    result.append("Num Instances:  ").append(numInstances()).append('\n');
    result.append("Num Attributes: ").append(numAttributes()).append('\n');
    result.append('\n');

    result.append(Utils.padLeft("", 5)).append(Utils.padRight("Name", 25));
    result.append(Utils.padLeft("Type", 5)).append(Utils.padLeft("Nom", 5));
    result.append(Utils.padLeft("Int", 5)).append(Utils.padLeft("Real", 5));
    result.append(Utils.padLeft("Missing", 12));
    result.append(Utils.padLeft("Unique", 12));
    result.append(Utils.padLeft("Dist", 6)).append('\n');
    for (int i = 0; i < numAttributes(); i++) {
      Attribute a = attribute(i);
      AttributeStats as = attributeStats(i);
      result.append(Utils.padLeft("" + (i + 1), 4)).append(' ');
      result.append(Utils.padRight(a.name(), 25)).append(' ');
      long percent;
      switch (a.type()) {
      
      case PreferenceAttribute.RANKING:
    result.append(Utils.padLeft("Rnk", 4)).append(' ');
    percent = Math.round(100.0 * as.intCount / as.totalCount);
    result.append(Utils.padLeft("" + percent, 3)).append("% ");
    result.append(Utils.padLeft("" + 0, 3)).append("% ");
    percent = Math.round(100.0 * as.realCount / as.totalCount);
    result.append(Utils.padLeft("" + percent, 3)).append("% ");
    break;
      case Attribute.NOMINAL:
	result.append(Utils.padLeft("Nom", 4)).append(' ');
	percent = Math.round(100.0 * as.intCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	result.append(Utils.padLeft("" + 0, 3)).append("% ");
	percent = Math.round(100.0 * as.realCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	break;
      case Attribute.NUMERIC:
	result.append(Utils.padLeft("Num", 4)).append(' ');
	result.append(Utils.padLeft("" + 0, 3)).append("% ");
	percent = Math.round(100.0 * as.intCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	percent = Math.round(100.0 * as.realCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	break;
      case Attribute.DATE:
	result.append(Utils.padLeft("Dat", 4)).append(' ');
	result.append(Utils.padLeft("" + 0, 3)).append("% ");
	percent = Math.round(100.0 * as.intCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	percent = Math.round(100.0 * as.realCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	break;
      case Attribute.STRING:
	result.append(Utils.padLeft("Str", 4)).append(' ');
	percent = Math.round(100.0 * as.intCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	result.append(Utils.padLeft("" + 0, 3)).append("% ");
	percent = Math.round(100.0 * as.realCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	break;
      case Attribute.RELATIONAL:
	result.append(Utils.padLeft("Rel", 4)).append(' ');
	percent = Math.round(100.0 * as.intCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	result.append(Utils.padLeft("" + 0, 3)).append("% ");
	percent = Math.round(100.0 * as.realCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	break;
      default:
	result.append(Utils.padLeft("???", 4)).append(' ');
	result.append(Utils.padLeft("" + 0, 3)).append("% ");
	percent = Math.round(100.0 * as.intCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	percent = Math.round(100.0 * as.realCount / as.totalCount);
	result.append(Utils.padLeft("" + percent, 3)).append("% ");
	break;
      }
      result.append(Utils.padLeft("" + as.missingCount, 5)).append(" /");
      percent = Math.round(100.0 * as.missingCount / as.totalCount);
      result.append(Utils.padLeft("" + percent, 3)).append("% ");
      result.append(Utils.padLeft("" + as.uniqueCount, 5)).append(" /");
      percent = Math.round(100.0 * as.uniqueCount / as.totalCount);
      result.append(Utils.padLeft("" + percent, 3)).append("% ");
      result.append(Utils.padLeft("" + as.distinctCount, 5)).append(' ');
      result.append('\n');
    }
    return result.toString();
  }

  /**
   * Copies instances from one set to the end of another 
   * one.
   *
   * @param from the position of the first instance to be copied
   * @param dest the destination for the instances
   * @param num the number of instances to be copied
   */
  //@ requires 0 <= from && from <= numInstances() - num;
  //@ requires 0 <= num;
  protected void copyInstances(int from, /*@non_null@*/ Instances dest, int num) {
    
    for (int i = 0; i < num; i++) {
      dest.add(instance(from + i));
    }
  }
  
  /**
   * Replaces the attribute information by a clone of
   * itself.
   */
  protected void freshAttributeInfo() {

    ArrayList<Attribute> newList = new ArrayList<Attribute>(m_Attributes.size());
    for (Attribute att : m_Attributes) {
      newList.add((Attribute)att.copy());
    }
    m_Attributes = newList;
  }
 
  /**
   * Returns string including all instances, their weights and
   * their indices in the original dataset.
   *
   * @return description of instance and its weight as a string
   */
  protected /*@pure@*/ String instancesAndWeights(){

    StringBuffer text = new StringBuffer();

    for (int i = 0; i < numInstances(); i++) {
      text.append(instance(i) + " " + instance(i).weight());
      if (i < numInstances() - 1) {
	text.append("\n");
      }
    }
    return text.toString();
  }
  
  /**
   * Partitions the instances around a pivot. Used by quicksort and
   * kthSmallestValue.
   *
   * @param attIndex the attribute's index (index starts with 0)
   * @param l the first index of the subset (index starts with 0)
   * @param r the last index of the subset (index starts with 0)
   *
   * @return the index of the middle element
   */
  //@ requires 0 <= attIndex && attIndex < numAttributes();
  //@ requires 0 <= left && left <= right && right < numInstances();
  protected int partition(int attIndex, int l, int r) {
    
    double pivot = instance((l + r) / 2).value(attIndex);

    while (l < r) {
      while ((instance(l).value(attIndex) < pivot) && (l < r)) {
        l++;
      }
      while ((instance(r).value(attIndex) > pivot) && (l < r)) {
        r--;
      }
      if (l < r) {
        swap(l, r);
        l++;
        r--;
      }
    }
    if ((l == r) && (instance(r).value(attIndex) > pivot)) {
      r--;
    } 

    return r;
  }
  
  /**
   * Implements quicksort according to Manber's "Introduction to
   * Algorithms".
   *
   * @param attIndex the attribute's index (index starts with 0)
   * @param left the first index of the subset to be sorted (index starts with 0)
   * @param right the last index of the subset to be sorted (index starts with 0)
   */
  //@ requires 0 <= attIndex && attIndex < numAttributes();
  //@ requires 0 <= first && first <= right && right < numInstances();
  protected void quickSort(int attIndex, int left, int right) {

    if (left < right) {
      int middle = partition(attIndex, left, right);
      quickSort(attIndex, left, middle);
      quickSort(attIndex, middle + 1, right);
    }
  }
  
  /**
   * Implements computation of the kth-smallest element according
   * to Manber's "Introduction to Algorithms".
   *
   * @param attIndex the attribute's index (index starts with 0)
   * @param left the first index of the subset (index starts with 0)
   * @param right the last index of the subset (index starts with 0)
   * @param k the value of k
   *
   * @return the index of the kth-smallest element
   */
  //@ requires 0 <= attIndex && attIndex < numAttributes();
  //@ requires 0 <= first && first <= right && right < numInstances();
  protected int select(int attIndex, int left, int right, int k) {
    
    if (left == right) {
      return left;
    } else {
      int middle = partition(attIndex, left, right);
      if ((middle - left + 1) >= k) {
        return select(attIndex, left, middle, k);
      } else {
        return select(attIndex, middle + 1, right, k - (middle - left + 1));
      }
    }
  }

  /**
   * Help function needed for stratification of set.
   *
   * @param numFolds the number of folds for the stratification
   */
  protected void stratStep (int numFolds){
    
    ArrayList<Instance> newVec = new ArrayList<Instance>(m_Instances.size());
    int start = 0, j;

    // create stratified batch
    while (newVec.size() < numInstances()) {
      j = start;
      while (j < numInstances()) {
	newVec.add(instance(j));
	j = j + numFolds;
      }
      start++;
    }
    m_Instances = newVec;
  }
  
  /**
   * Swaps two instances in the set.
   *
   * @param i the first instance's index (index starts with 0)
   * @param j the second instance's index (index starts with 0)
   */
  //@ requires 0 <= i && i < numInstances();
  //@ requires 0 <= j && j < numInstances();
  public void swap(int i, int j){
    
    Instance in = m_Instances.get(i);
    m_Instances.set(i, m_Instances.get(j));
    m_Instances.set(j, in);
  }

  /**
   * Merges two sets of Instances together. The resulting set will have
   * all the attributes of the first set plus all the attributes of the 
   * second set. The number of instances in both sets must be the same.
   *
   * @param first the first set of Instances
   * @param second the second set of Instances
   * @return the merged set of Instances
   * @throws IllegalArgumentException if the datasets are not the same size
   */
  public static Instances mergeInstances(Instances first, Instances second) {

    if (first.numInstances() != second.numInstances()) {
      throw new IllegalArgumentException("Instance sets must be of the same size");
    }

    // Create the vector of merged attributes
    ArrayList<Attribute> newAttributes = new ArrayList<Attribute>();
    for (int i = 0; i < first.numAttributes(); i++) {
      newAttributes.add(first.attribute(i));
    }
    for (int i = 0; i < second.numAttributes(); i++) {
      newAttributes.add(second.attribute(i));
    }
    
    // Create the set of Instances
    Instances merged = new Instances(first.relationName() + '_'
				     + second.relationName(), 
				     newAttributes, 
				     first.numInstances());
    // Merge each instance
    for (int i = 0; i < first.numInstances(); i++) {
      merged.add(first.instance(i).mergeInstance(second.instance(i)));
    }
    return merged;
  }

  /**
   * Method for testing this class.
   *
   * @param argv should contain one element: the name of an ARFF file
   */
  //@ requires argv != null;
  //@ requires argv.length == 1;
  //@ requires argv[0] != null;
  public static void test(String [] argv) {

    Instances instances, secondInstances, train, test, empty;
    Random random = new Random(2);
    Reader reader;
    int start, num;
    ArrayList<Attribute> testAtts;
    ArrayList<String> testVals;
    int i,j;
    
    try{
      if (argv.length > 1) {
	throw (new Exception("Usage: Instances [<filename>]"));
      }
      
      // Creating set of instances from scratch
      testVals = new ArrayList<String>(2);
      testVals.add("first_value");
      testVals.add("second_value");
      testAtts = new ArrayList<Attribute>(2);
      testAtts.add(new Attribute("nominal_attribute", testVals));
      testAtts.add(new Attribute("numeric_attribute"));
      instances = new Instances("test_set", testAtts, 10);
      instances.add(new DenseInstance(instances.numAttributes()));
      instances.add(new DenseInstance(instances.numAttributes()));
      instances.add(new DenseInstance(instances.numAttributes()));
      instances.setClassIndex(0);
      System.out.println("\nSet of instances created from scratch:\n");
      System.out.println(instances);
      
      if (argv.length == 1) {
	String filename = argv[0];
	reader = new FileReader(filename);
	
	// Read first five instances and print them
	System.out.println("\nFirst five instances from file:\n");
	instances = new Instances(reader, 1);
	instances.setClassIndex(instances.numAttributes() - 1);
	i = 0;
	while ((i < 5) && (instances.readInstance(reader))) {
	  i++;
	}
	System.out.println(instances);

	// Read all the instances in the file
	reader = new FileReader(filename);
	instances = new Instances(reader);

	// Make the last attribute be the class 
	instances.setClassIndex(instances.numAttributes() - 1);
	
	// Print header and instances.
	System.out.println("\nDataset:\n");
	System.out.println(instances);
	System.out.println("\nClass index: "+instances.classIndex());
      }
      
      // Test basic methods based on class index.
      System.out.println("\nClass name: "+instances.classAttribute().name());
      System.out.println("\nClass index: "+instances.classIndex());
      System.out.println("\nClass is nominal: " +
			 instances.classAttribute().isNominal());
      System.out.println("\nClass is a ranking: " +
 			 instances.classAttribute().isRanking());
      System.out.println("\nClass is numeric: " +
			 instances.classAttribute().isNumeric());
      System.out.println("\nClasses:\n");
      for (i = 0; i < instances.numClasses(); i++) {
	System.out.println(instances.classAttribute().value(i));
      }
      System.out.println("\nClass values and labels of instances:\n");
      for (i = 0; i < instances.numInstances(); i++) {
	Instance inst = instances.instance(i);
	System.out.print(inst.classValue() + "\t");
	System.out.print(inst.toString(inst.classIndex()));
	if (instances.instance(i).classIsMissing()) {
	  System.out.println("\tis missing");
	} else {
	  System.out.println();
	}
      }
      
      // Create random weights.
      System.out.println("\nCreating random weights for instances.");
      for (i = 0; i < instances.numInstances(); i++) {
	instances.instance(i).setWeight(random.nextDouble()); 
      }
      
      // Print all instances and their weights (and the sum of weights).
      System.out.println("\nInstances and their weights:\n");
      System.out.println(instances.instancesAndWeights());
      System.out.print("\nSum of weights: ");
      System.out.println(instances.sumOfWeights());
      
      // Insert an attribute
      secondInstances = new Instances(instances);
      Attribute testAtt = new Attribute("Inserted");
      secondInstances.insertAttributeAt(testAtt, 0);
      System.out.println("\nSet with inserted attribute:\n");
      System.out.println(secondInstances);
      System.out.println("\nClass name: "
			 + secondInstances.classAttribute().name());
      
      // Delete the attribute
      secondInstances.deleteAttributeAt(0);
      System.out.println("\nSet with attribute deleted:\n");
      System.out.println(secondInstances);
      System.out.println("\nClass name: "
			 + secondInstances.classAttribute().name());
      
      // Test if headers are equal
      System.out.println("\nHeaders equal: "+
			 instances.equalHeaders(secondInstances) + "\n");
      
      // Print data in internal format.
      System.out.println("\nData (internal values):\n");
      for (i = 0; i < instances.numInstances(); i++) {
	for (j = 0; j < instances.numAttributes(); j++) {
	  if (instances.instance(i).isMissing(j)) {
	    System.out.print("? ");
	  } else {
	    System.out.print(instances.instance(i).value(j) + " ");
	  }
	}
	System.out.println();
      }
      
      // Just print header
      System.out.println("\nEmpty dataset:\n");
      empty = new Instances(instances, 0);
      System.out.println(empty);
      System.out.println("\nClass name: "+empty.classAttribute().name());

      // Create copy and rename an attribute and a value (if possible)
      if (empty.classAttribute().isNominal() || empty.classAttribute().isRanking()) {
	Instances copy = new Instances(empty, 0);
	copy.renameAttribute(copy.classAttribute(), "new_name");
	copy.renameAttributeValue(copy.classAttribute(), 
				  copy.classAttribute().value(0), 
				  "new_val_name");
	System.out.println("\nDataset with names changed:\n" + copy);
	System.out.println("\nOriginal dataset:\n" + empty);
      }

      // Create and prints subset of instances.
      start = instances.numInstances() / 4;
      num = instances.numInstances() / 2;
      System.out.print("\nSubset of dataset: ");
      System.out.println(num + " instances from " + (start + 1) 
			 + ". instance");
      secondInstances = new Instances(instances, start, num);
      System.out.println("\nClass name: "
			 + secondInstances.classAttribute().name());

      // Print all instances and their weights (and the sum of weights).
      System.out.println("\nInstances and their weights:\n");
      System.out.println(secondInstances.instancesAndWeights());
      System.out.print("\nSum of weights: ");
      System.out.println(secondInstances.sumOfWeights());
      
      // Create and print training and test sets for 3-fold
      // cross-validation.
      System.out.println("\nTrain and test folds for 3-fold CV:");
      if (instances.classAttribute().isNominal() || instances.classAttribute().isRanking()) {
	instances.stratify(3);
      }
      for (j = 0; j < 3; j++) {
        train = instances.trainCV(3,j, new Random(1));
	test = instances.testCV(3,j);
                      
	// Print all instances and their weights (and the sum of weights).
	System.out.println("\nTrain: ");
	System.out.println("\nInstances and their weights:\n");
	System.out.println(train.instancesAndWeights());
	System.out.print("\nSum of weights: ");
	System.out.println(train.sumOfWeights());
	System.out.println("\nClass name: "+train.classAttribute().name());
	System.out.println("\nTest: ");
	System.out.println("\nInstances and their weights:\n");
	System.out.println(test.instancesAndWeights());
	System.out.print("\nSum of weights: ");
	System.out.println(test.sumOfWeights());
	System.out.println("\nClass name: "+test.classAttribute().name());
      }

      // Randomize instances and print them.
      System.out.println("\nRandomized dataset:");
      instances.randomize(random);
      
      // Print all instances and their weights (and the sum of weights).
      System.out.println("\nInstances and their weights:\n");
      System.out.println(instances.instancesAndWeights());
      System.out.print("\nSum of weights: ");
      System.out.println(instances.sumOfWeights());

      // Sort instances according to first attribute and
      // print them.
      System.out.print("\nInstances sorted according to first attribute:\n ");
      instances.sort(0);
        
      // Print all instances and their weights (and the sum of weights).
      System.out.println("\nInstances and their weights:\n");
      System.out.println(instances.instancesAndWeights());
      System.out.print("\nSum of weights: ");
      System.out.println(instances.sumOfWeights());
    } catch (Exception e) {
      e.printStackTrace(); 
    }
  }

  /**
   * Main method for this class. The following calls are possible:
   * <ul>
   *   <li>
   *     <code>weka.core.Instances</code> help<br/>
   *     prints a short list of possible commands.
   *   </li>
   *   <li>
   *     <code>weka.core.Instances</code> &lt;filename&gt;<br/>
   *     prints a summary of a set of instances.
   *   </li>
   *   <li>
   *     <code>weka.core.Instances</code> merge &lt;filename1&gt; &lt;filename2&gt;<br/>
   *     merges the two datasets (must have same number of instances) and
   *     outputs the results on stdout.
   *   </li>
   *   <li>
   *     <code>weka.core.Instances</code> append &lt;filename1&gt; &lt;filename2&gt;<br/>
   *     appends the second dataset to the first one (must have same headers) and
   *     outputs the results on stdout.
   *   </li>
   *   <li>
   *     <code>weka.core.Instances</code> headers &lt;filename1&gt; &lt;filename2&gt;<br/>
   *     Compares the headers of the two datasets and prints whether they match
   *     or not.
   *   </li>
   *   <li>
   *     <code>weka.core.Instances</code> randomize &lt;seed&gt; &lt;filename&gt;<br/>
   *     randomizes the dataset with the given seed and outputs the result on stdout.
   *   </li>
   * </ul>
   *
   * @param args 	the commandline parameters
   */
  public static void main(String[] args) {

    try {
      Instances i;
      // read from stdin and print statistics
      if (args.length == 0) {
	DataSource source = new DataSource(System.in);
	i = source.getDataSet();
	System.out.println(i.toSummaryString());
      }
      // read file and print statistics
      else if ((args.length == 1) && (!args[0].equals("-h")) && (!args[0].equals("help"))) {
	DataSource source = new DataSource(args[0]);
	i = source.getDataSet();
	System.out.println(i.toSummaryString());
      }
      // read two files, merge them and print result to stdout
      else if ((args.length == 3) && (args[0].toLowerCase().equals("merge"))) {
	DataSource source1 = new DataSource(args[1]);
	DataSource source2 = new DataSource(args[2]);
	i = Instances.mergeInstances(source1.getDataSet(), source2.getDataSet());
	System.out.println(i);
      }
      // read two files, append them and print result to stdout
      else if ((args.length == 3) && (args[0].toLowerCase().equals("append"))) {
	DataSource source1 = new DataSource(args[1]);
	DataSource source2 = new DataSource(args[2]);
	String msg = source1.getStructure().equalHeadersMsg(source2.getStructure());
	if (msg != null)
	  throw new Exception("The two datasets have different headers:\n" + msg);
	Instances structure = source1.getStructure();
	System.out.println(source1.getStructure());
	while (source1.hasMoreElements(structure))
	  System.out.println(source1.nextElement(structure));
	structure = source2.getStructure();
	while (source2.hasMoreElements(structure))
	  System.out.println(source2.nextElement(structure));
      }
      // read two files and compare their headers
      else if ((args.length == 3) && (args[0].toLowerCase().equals("headers"))) {
	DataSource source1 = new DataSource(args[1]);
	DataSource source2 = new DataSource(args[2]);
	String msg = source1.getStructure().equalHeadersMsg(source2.getStructure());
	if (msg == null)
	  System.out.println("Headers match");
	else
	  System.out.println("Headers don't match:\n" + msg);
      }
      // read file and seed value, randomize data and print result to stdout
      else if ((args.length == 3) && (args[0].toLowerCase().equals("randomize"))) {
	DataSource source = new DataSource(args[2]);
	i = source.getDataSet();
	i.randomize(new Random(Integer.parseInt(args[1])));
	System.out.println(i);
      }
      // wrong parameters or help
      else {
	System.err.println(
	    "\nUsage:\n"
	    // help
	    + "\tweka.core.Instances help\n"
	    + "\t\tPrints this help\n"
	    // stats
	    + "\tweka.core.Instances <filename>\n"
	    + "\t\tOutputs dataset statistics\n"
	    // merge
	    + "\tweka.core.Instances merge <filename1> <filename2>\n"
	    + "\t\tMerges the datasets (must have same number of rows).\n"
	    + "\t\tGenerated dataset gets output on stdout.\n"
	    // append
	    + "\tweka.core.Instances append <filename1> <filename2>\n"
	    + "\t\tAppends the second dataset to the first (must have same number of attributes).\n"
	    + "\t\tGenerated dataset gets output on stdout.\n"
	    // headers
	    + "\tweka.core.Instances headers <filename1> <filename2>\n"
	    + "\t\tCompares the structure of the two datasets and outputs whether they\n"
	    + "\t\tdiffer or not.\n"
	    // randomize
	    + "\tweka.core.Instances randomize <seed> <filename>\n"
	    + "\t\tRandomizes the dataset and outputs it on stdout.\n"
	);
      }
    }
    catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6354 $");
  }
  
  //RANKING BEGIN
  /**
   * Returns a list of labels with all the labels read from the input file.
   * The first item of the String array contains the label's name, the second item is an index.
   * 
   * @return an indexed list of the single labels appearing.
   */
  public List<String[]> getLabels(){
  	return labels;
  }
  
  /**
   * Setting the list of labels.
   * 
   * @param input
   */
  public void setLabels(List<String[]> input){
  	labels = input;
  }
  
  /**
   * Returns a label out of the label arrayList for a specified index.
   * @param index
   * @return
   */
  public String getLabelFromIndex(int index){
	  String indx = String.valueOf(index);
	  for(int i=0; i<getLabels().size(); i++){
		  if(getLabels().get(i)[1].equals(indx)){
			  return getLabels().get(i)[0];
		  }
	  }
	  return null;
  }
  //RANKING END
}
