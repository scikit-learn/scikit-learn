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
 * JSONInstances.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.core.json;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.DenseInstance;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSource;
import weka.core.labelranking.PreferenceAttribute;

import java.util.ArrayList;

/**
 * Class for transforming Instances objects into <a href="http://www.json.org/" target="_blank">JSON</a>
 * objects and vice versa.
 * 
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class JSONInstances {

  /** the header section. */
  public final static String HEADER = "header";

  /** the data section. */
  public final static String DATA = "data";

  /** the relation name. */
  public final static String RELATION = "relation";

  /** the attributes object. */
  public final static String ATTRIBUTES = "attributes";

  /** the name attribute. */
  public final static String NAME = "name";

  /** the type attribute. */
  public final static String TYPE = "type";

  /** the class attribute indicator. */
  public final static String CLASS = "class";

  /** the labels attribute. */
  public final static String LABELS = "labels";

  /** the weight attribute. */
  public final static String WEIGHT = "weight";

  /** the dateformat attribute. */
  public final static String DATEFORMAT = "dateformat";

  /** the sparse attribute. */
  public final static String SPARSE = "sparse";

  /** the values attribute. */
  public final static String VALUES = "values";

  /** the separator for index/value in case of sparse instances. */
  public final static String SPARSE_SEPARATOR = ":";
  
  /**
   * Turns the JSON object into an Attribute, if possible.
   * 
   * @param att		the JSON object to turn into an Attribute
   * @param classAtt	for storing whether the attribute is the class attribute
   * @return		the Attribute, null in case of an error
   */
  protected static Attribute toAttribute(JSONNode att, boolean[] classAtt) {
    Attribute	result;
    String	name;
    String	type;
    String	dateformat;
    JSONNode	labels;
    ArrayList<String>	values;
    int		i;
    double	weight;
    
    name   = (String) att.getChild(NAME).getValue("noname");
    type   = (String) att.getChild(TYPE).getValue("");
    weight = (Double) att.getChild(WEIGHT).getValue(new Double(1.0));
    if (type.equals(Attribute.typeToString(Attribute.NUMERIC))) {
      result = new Attribute(name);
    }
    else if (type.equals(Attribute.typeToString(Attribute.NOMINAL)) || type.equals(Attribute.typeToString(PreferenceAttribute.RANKING))) {
      labels = att.getChild(LABELS);
      values = new ArrayList<String>();
      for (i = 0; i < labels.getChildCount(); i++)
	values.add((String)((JSONNode) labels.getChildAt(i)).getValue());
      result = new Attribute(name, values);
    }
    else if (type.equals(Attribute.typeToString(Attribute.DATE))) {
      dateformat = (String) att.getChild(TYPE).getValue("yyyy-MM-dd'T'HH:mm:ss");
      result     = new Attribute(name, dateformat);
    }
    else if (type.equals(Attribute.typeToString(Attribute.STRING))) {
      result = new Attribute(name, (ArrayList<String>) null);
    }
    else {
      System.err.println("Unhandled attribute type '" + type + "'!");
      return null;
    }
    result.setWeight(weight);
    
    return result;
  }

  /**
   * Turns the JSON Object into an Instance, if possible.
   * 
   * @param inst	the JSON object to turn into an Instance
   * @param data	the data so far (only used for header information)
   * @return		the Instance, null in case of an error
   */
  protected static Instance toInstance(JSONNode inst, Instances data) {
    Instance	result;
    boolean	sparse;
    double	weight;
    JSONNode	values;
    int		i;
    int		index;
    int		pos;
    String	value;
    double[]	vals;

    sparse = (Boolean) inst.getChild(SPARSE).getValue(new Boolean(false));
    weight = (Double) inst.getChild(WEIGHT).getValue(new Double(1.0));
    values = inst.getChild(VALUES);
    vals   = new double[data.numAttributes()];
    for (i = 0; i < values.getChildCount(); i++) {
      if (sparse) {
	value = "" + ((JSONNode) values.getChildAt(i)).getValue();
	pos   = value.indexOf(SPARSE_SEPARATOR);
	index = Integer.parseInt(value.substring(0, pos));
	value = value.substring(pos + 1);
      }
      else {
	index = i;
	value = "" + ((JSONNode) values.getChildAt(i)).getValue();
      }
      
      try {
	if (data.attribute(index).isNumeric()) {
	  vals[index] = Double.parseDouble(value);
	}
	else if (data.attribute(index).isNominal() || data.attribute(index).isRanking()) {
	  vals[index] = data.attribute(index).indexOfValue(value);
	  if ((vals[index] == -1) && value.startsWith("'") && value.endsWith("'"))
	    vals[index] = data.attribute(index).indexOfValue(Utils.unquote(value));
	  if (vals[index] == -1) {
	    System.err.println("Unknown label '" + value + "' for attribute #" + (index+1) + "!");
	    return null;
	  }
	}
	else if (data.attribute(index).isDate()) {
	  vals[index] = data.attribute(index).parseDate(value);
	}
	else if (data.attribute(index).isString()) {
	  vals[index] = data.attribute(index).addStringValue(value);
	}
	else {
	  System.err.println("Unhandled attribute type '" + Attribute.typeToString(data.attribute(index).type()) + "'!");
	  return null;
	}
      }
      catch (Exception e) {
	System.err.println("Error parsing value #" + (index+1) + ": " + e.toString());
	return null;
      }
    }

    result = new DenseInstance(weight, vals);
    result.setDataset(data);
      
    return result;
  }
  
  /**
   * Turns a JSON object, if possible, into an Instances object.
   * 
   * @param json	the JSON object to convert
   * @param onlyHeader	whether to retrieve only the header
   * @return		the generated Instances object, null if not possible
   */
  protected static Instances toInstances(JSONNode json, boolean onlyHeader) {
    Instances	result;
    JSONNode	header;
    JSONNode	attributes;
    JSONNode	data;
    ArrayList<Attribute>	atts;
    Attribute	att;
    Instance	inst;
    int		i;
    int		classIndex;
    boolean[]	classAtt;
    
    header = json.getChild(HEADER);
    if (header == null) {
      System.err.println("No '" + HEADER + "' section!");
      return null;
    }
    data = json.getChild(DATA);
    if (data == null) {
      System.err.println("No '" + DATA + "' section!");
      return null;
    }
    
    // attributes
    attributes = header.getChild(ATTRIBUTES);
    if (attributes == null) {
      System.err.println("No '" + ATTRIBUTES + "' array!");
      return null;
    }
    atts       = new ArrayList<Attribute>();
    classAtt   = new boolean[1];
    classIndex = -1;
    for (i = 0; i < attributes.getChildCount(); i++) {
      att = toAttribute((JSONNode) attributes.getChildAt(i), classAtt);
      if (att == null) {
	System.err.println("Could not convert attribute #" + (i+1) + "!");
	return null;
      }
      if (classAtt[0])
	classIndex = i;
      atts.add(att);
    }
    result = new Instances(
	header.getChild(RELATION).getValue("unknown").toString(), 
	atts, 
	(onlyHeader ? 0 : data.getChildCount()));
    result.setClassIndex(classIndex);
    
    // data
    if (!onlyHeader) {
      for (i = 0; i < data.getChildCount(); i++) {
	inst = toInstance((JSONNode) data.getChildAt(i), result);
	if (inst == null) {
	  System.err.println("Could not convert instance #" + (i+1) + "!");
	  return null;
	}
	result.add(inst);
      }
    }
    
    return result;
  }
  
  /**
   * Turns a JSON object, if possible, into an Instances object.
   * 
   * @param json	the JSON object to convert
   * @return		the generated Instances object, null if not possible
   */
  public static Instances toInstances(JSONNode json) {
    return toInstances(json, false);
  }
  
  /**
   * Turns a JSON object, if possible, into an Instances object (only header).
   * 
   * @param json	the JSON object to convert
   * @return		the generated Instances header object, null if not possible
   */
  public static Instances toHeader(JSONNode json) {
    return toInstances(json, true);
  }
  
  /**
   * Turns the Attribute into a JSON object.
   * 
   * @param inst	the corresponding dataset
   * @param att		the attribute to convert
   * @return		the JSON object
   */
  protected static JSONNode toJSON(Instances inst, Attribute att) {
    JSONNode	result;
    JSONNode	labels;
    int		i;
    
    result = new JSONNode();

    result.addPrimitive(NAME, att.name());
    result.addPrimitive(TYPE, Attribute.typeToString(att));
    result.addPrimitive(CLASS, (att.index() == inst.classIndex()));
    result.addPrimitive(WEIGHT, att.weight());
    if (att.isNominal()|| att.isRanking()) {
      labels = result.addArray(LABELS);
      for (i = 0; i < att.numValues(); i++)
	labels.addArrayElement(att.value(i));
    }
    if (att.isDate())
      result.addPrimitive(DATEFORMAT, att.getDateFormat());
    
    return result;
  }
  
  /**
   * Turns the Instance into a JSON object.
   * 
   * @param inst	the Instance to convert
   * @return		the JSON object
   */
  protected static JSONNode toJSON(Instance inst) {
    JSONNode	result;
    JSONNode	values;
    int		i;
    boolean	sparse;
    
    result = new JSONNode();
    
    sparse = (inst instanceof SparseInstance);
    result.addPrimitive(SPARSE, sparse);
    result.addPrimitive(WEIGHT, inst.weight());
    values = result.addArray(VALUES);
    if (sparse) {
      for (i = 0; i < inst.numValues(); i++)
	values.addArrayElement(inst.index(i) + SPARSE_SEPARATOR + inst.toString(inst.index(i)));
    }
    else {
      for (i = 0; i < inst.numAttributes(); i++)
	values.addArrayElement(inst.toString(i));
    }
    
    return result;
  }
  
  /**
   * Turns the Instances object into a JSON object.
   * 
   * @param inst	the Instances to turn into a JSON object
   * @return		the JSON object
   */
  public static JSONNode toJSON(Instances inst) {
    JSONNode	result;
    JSONNode	header;
    JSONNode	atts;
    JSONNode	data;
    int		i;
    
    result = new JSONNode();
    
    // header
    header = result.addObject(HEADER);
    header.addPrimitive(RELATION, inst.relationName());
    atts = header.addArray(ATTRIBUTES);
    for (i = 0; i < inst.numAttributes(); i++)
      atts.add(toJSON(inst, inst.attribute(i)));
    
    // data
    data = result.addArray(DATA);
    for (i = 0; i < inst.numInstances(); i++)
      data.add(toJSON(inst.instance(i)));
    
    return result;
  }
  
  /**
   * For testing only.
   * 
   * @param args	expects a dataset as first parameter
   * @throws Exception	if something goes wrong
   */
  public static void main(String[] args) throws Exception {
    if (args.length != 1) {
      System.err.println("No dataset supplied!");
      System.exit(1);
    }

    // load dataset
    Instances data = DataSource.read(args[0]);
    
    // turn Instances into JSON object and output it
    JSONNode json = toJSON(data);
    StringBuffer buffer = new StringBuffer();
    json.toString(buffer);
    System.out.println(buffer.toString());
    
    // turn JSON object back into Instances and output it
    Instances inst = toInstances(json);
    System.out.println(inst);
  }
}
