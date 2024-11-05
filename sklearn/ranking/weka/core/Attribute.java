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
 *    Attribute.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.io.IOException;
import java.io.Serializable;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Properties;
import java.util.List;
import java.util.ArrayList;
import javax.swing.JOptionPane;

import weka.core.labelranking.PreferenceAttribute;
import weka.core.labelranking.RankUtilities;

/** 
 * Class for handling an attribute. Once an attribute has been created,
 * it can't be changed. <p>
 *
 * The following attribute types are supported:
 * <ul>
 *    <li> numeric: <br/>
 *         This type of attribute represents a floating-point number.
 *    </li>
 *    <li> nominal: <br/>
 *         This type of attribute represents a fixed set of nominal values.
 *    </li>
 *    <li> string: <br/>
 *         This type of attribute represents a dynamically expanding set of
 *         nominal values. Usually used in text classification.
 *    </li>
 *    <li> date: <br/>
 *         This type of attribute represents a date, internally represented as 
 *         floating-point number storing the milliseconds since January 1, 
 *         1970, 00:00:00 GMT. The string representation of the date must be
 *         <a href="http://www.iso.org/iso/en/prods-services/popstds/datesandtime.html" target="_blank">
 *         ISO-8601</a> compliant, the default is <code>yyyy-MM-dd'T'HH:mm:ss</code>.
 *    </li>
 *    <li> relational: <br/>
 *         This type of attribute can contain other attributes and is, e.g., 
 *         used for representing Multi-Instance data. (Multi-Instance data
 *         consists of a nominal attribute containing the bag-id, then a 
 *         relational attribute with all the attributes of the bag, and 
 *         finally the class attribute.)
 *    </li>
 * </ul>
 * 
 * Typical usage (code from the main() method of this class): <p>
 *
 * <code>
 * ... <br>
 *
 * // Create numeric attributes "length" and "weight" <br>
 * Attribute length = new Attribute("length"); <br>
 * Attribute weight = new Attribute("weight"); <br><br>
 * 
 * // Create list to hold nominal values "first", "second", "third" <br>
 * List<String> my_nominal_values = new ArrayList<String>(3); <br>
 * my_nominal_values.add("first"); <br>
 * my_nominal_values.add("second"); <br>
 * my_nominal_values.add("third"); <br><br>
 *
 * // Create nominal attribute "position" <br>
 * Attribute position = new Attribute("position", my_nominal_values);<br>
 *
 * ... <br>
 * </code><p>
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5987 $
 */
public class Attribute
  implements Copyable, Serializable, RevisionHandler {


  /** for serialization */
  static final long serialVersionUID = -742180568732916383L;
  
  /** Constant set for numeric attributes. */
  public static final int NUMERIC = 0;

  /** Constant set for nominal attributes. */
  public static final int NOMINAL = 1;

  /** Constant set for attributes with string values. */
  public static final int STRING = 2;

  /** Constant set for attributes with date values. */
  public static final int DATE = 3;

  /** Constant set for relation-valued attributes. */
  public static final int RELATIONAL = 4;

  /** Constant set for symbolic attributes. */
  public static final int ORDERING_SYMBOLIC = 0;

  /** Constant set for ordered attributes. */
  public static final int ORDERING_ORDERED  = 1;

  /** Constant set for modulo-ordered attributes. */
  public static final int ORDERING_MODULO   = 2;

  /** The keyword used to denote the start of an arff attribute declaration */
  public final static String ARFF_ATTRIBUTE = "@attribute";

  /** A keyword used to denote a numeric attribute */
  public final static String ARFF_ATTRIBUTE_INTEGER = "integer";

  /** A keyword used to denote a numeric attribute */
  public final static String ARFF_ATTRIBUTE_REAL = "real";

  /** A keyword used to denote a numeric attribute */
  public final static String ARFF_ATTRIBUTE_NUMERIC = "numeric";

  /** The keyword used to denote a string attribute */
  public final static String ARFF_ATTRIBUTE_STRING = "string";

  /** The keyword used to denote a date attribute */
  public final static String ARFF_ATTRIBUTE_DATE = "date";

  /** The keyword used to denote a relation-valued attribute */
  public final static String ARFF_ATTRIBUTE_RELATIONAL = "relational";

  /** The keyword used to denote the end of the declaration of a subrelation */
  public final static String ARFF_END_SUBRELATION = "@end";

  /** Strings longer than this will be stored compressed. */
  private static final int STRING_COMPRESS_THRESHOLD = 200;

  /** The attribute's name. */
  private /*@ spec_public non_null @*/ String m_Name;

  /** The attribute's type. */
  protected /*@ spec_public @*/ int m_Type;
  /*@ invariant m_Type == NUMERIC || 
                m_Type == DATE || 
                m_Type == STRING || 
                m_Type == NOMINAL ||
                m_Type == RELATIONAL ||
                m_Type == RANKING
  */

  /** The attribute's values (if nominal or string). */
  private /*@ spec_public @*/ ArrayList<Object> m_Values;

  /** Mapping of values to indices (if nominal or string). */
  private Hashtable<Object,Integer> m_Hashtable;

  /** The header information for a relation-valued attribute. */
  private Instances m_Header;

  /** Date format specification for date attributes */
  private SimpleDateFormat m_DateFormat;

  /** The attribute's index. */
  private /*@ spec_public @*/ int m_Index;

  /** The attribute's metadata. */
  private ProtectedProperties m_Metadata;

  /** The attribute's ordering. */
  private int m_Ordering;

  /** Whether the attribute is regular. */
  private boolean m_IsRegular;

  /** Whether the attribute is averagable. */
  private boolean m_IsAveragable;

  /** Whether the attribute has a zeropoint. */
  private boolean m_HasZeropoint;

  /** The attribute's weight. */
  private double m_Weight;

  /** The attribute's lower numeric bound. */
  private double m_LowerBound;

  /** Whether the lower bound is open. */
  private boolean m_LowerBoundIsOpen;

  /** The attribute's upper numeric bound. */
  private double m_UpperBound;

  /** Whether the upper bound is open */
  private boolean m_UpperBoundIsOpen;

  /**
   * Manipulates the attribute's type.
   * @param type
   */ 
  public void setM_type (int type){
	  m_Type = type;
  }
  
  /**
   * Constructor for a numeric attribute.
   *
   * @param attributeName the name for the attribute
   */
  //@ requires attributeName != null;
  //@ ensures  m_Name == attributeName;
  public Attribute(String attributeName) {

    this(attributeName, new ProtectedProperties(new Properties()));
//    instances = new ArrayList<Double>();
  }

  /**
   * Constructor for a numeric attribute, where metadata is supplied.
   *
   * @param attributeName the name for the attribute
   * @param metadata the attribute's properties
   */
  //@ requires attributeName != null;
  //@ requires metadata != null;
  //@ ensures  m_Name == attributeName;
  public Attribute(String attributeName, ProtectedProperties metadata) {
//	instances = new ArrayList<Double>();
    m_Name = attributeName;
    m_Index = -1;
    m_Values = null;
    m_Hashtable = null;
    m_Header = null;
    m_Type = NUMERIC;
    setMetadata(metadata);
  }

  /**
   * Constructor for a date attribute.
   *
   * @param attributeName the name for the attribute
   * @param dateFormat a string suitable for use with
   * SimpleDateFormatter for parsing dates.
   */
  //@ requires attributeName != null;
  //@ requires dateFormat != null;
  //@ ensures  m_Name == attributeName;
  public Attribute(String attributeName, String dateFormat) {

    this(attributeName, dateFormat,
	 new ProtectedProperties(new Properties()));
//    instances = new ArrayList<Double>();
  }

  /**
   * Constructor for a date attribute, where metadata is supplied.
   *
   * @param attributeName the name for the attribute
   * @param dateFormat a string suitable for use with
   * SimpleDateFormatter for parsing dates.
   * @param metadata the attribute's properties
   */
  //@ requires attributeName != null;
  //@ requires dateFormat != null;
  //@ requires metadata != null;
  //@ ensures  m_Name == attributeName;
  public Attribute(String attributeName, String dateFormat,
		   ProtectedProperties metadata) {

    m_Name = attributeName;
    m_Index = -1;
    m_Values = null;
    m_Hashtable = null;
    m_Header = null;
    m_Type = DATE;
    if (dateFormat != null) {
      m_DateFormat = new SimpleDateFormat(dateFormat);
    } else {
      m_DateFormat = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");
    }
    m_DateFormat.setLenient(false);
    setMetadata(metadata);
//    instances = new ArrayList<Double>();
  }

  /**
   * Constructor for nominal attributes and string attributes.
   * If a null vector of attribute values is passed to the method,
   * the attribute is assumed to be a string.
   *
   * @param attributeName the name for the attribute
   * @param attributeValues a vector of strings denoting the 
   * attribute values. Null if the attribute is a string attribute.
   */
  //@ requires attributeName != null;
  //@ ensures  m_Name == attributeName;
  public Attribute(String attributeName, 
		   List<String> attributeValues) {

    this(attributeName, attributeValues,
	 new ProtectedProperties(new Properties()));
//    instances = new ArrayList<Double>();
  }

  /**
   * Constructor for nominal attributes and string attributes, where
   * metadata is supplied. If a null vector of attribute values is passed
   * to the method, the attribute is assumed to be a string.
   *
   * @param attributeName the name for the attribute
   * @param attributeValues a vector of strings denoting the 
   * attribute values. Null if the attribute is a string attribute.
   * @param metadata the attribute's properties
   */
  //@ requires attributeName != null;
  //@ requires metadata != null;
  /*@ ensures  m_Name == attributeName;
      ensures  m_Index == -1;
      ensures  attributeValues == null && m_Type == STRING
            || attributeValues != null && m_Type == NOMINAL 
                  && m_Values.size() == attributeValues.size();
      signals (IllegalArgumentException ex) 
                 (* if duplicate strings in attributeValues *);
  */
  public Attribute(String attributeName, 
		   List<String> attributeValues,
		   ProtectedProperties metadata) {

    m_Name = attributeName;
    m_Index = -1;
    if (attributeValues == null) {
      m_Values = new ArrayList<Object>();
      m_Hashtable = new Hashtable<Object,Integer>();
      m_Header = null;
      m_Type = STRING;
    } else {

      //RANKING BEGIN 
        boolean preference = false;
        m_Values = new ArrayList<Object>(attributeValues.size());
        m_Hashtable = new Hashtable<Object,Integer>(attributeValues.size());
        m_Header = null;
    	  for (int i = 0; i < attributeValues.size(); i++) {
    		  Object store = attributeValues.get(i);
    		  if (((String)store).length() > STRING_COMPRESS_THRESHOLD) {
    			  try {
    				  store = new SerializedObject(attributeValues.get(i), true);
    			  } catch (Exception ex) {
    				  System.err.println("Couldn't compress nominal attribute value -"
    						  + " storing uncompressed.");
    			  }
    		  }

    		  if (m_Hashtable.containsKey(store) && !store.toString().contains(">")) {
    			  throw new IllegalArgumentException("A nominal attribute (" +
					     attributeName + ") cannot"
					     + " have duplicate labels (" + store + ").");
    		  }
    		  
    		  m_Values.add(store);
    		  m_Hashtable.put(store, new Integer(i));
    	  }
 //     }
    	  
    	  if(attributeName.toUpperCase().contains("RANKING")) 
    		  m_Type = PreferenceAttribute.RANKING;
    	  else 
    		  m_Type = NOMINAL;
    	//RANKING END 
    }
    setMetadata(metadata);
//    instances = new ArrayList<Double>();
  }


  /**
   * Constructor for relation-valued attributes.
   *
   * @param attributeName the name for the attribute
   * @param header an Instances object specifying the header of the relation.
   */
  public Attribute(String attributeName, Instances header) {

    this(attributeName, header,
	 new ProtectedProperties(new Properties()));
//    instances = new ArrayList<Double>();
  }

  /**
   * Constructor for relation-valued attributes.
   *
   * @param attributeName the name for the attribute
   * @param header an Instances object specifying the header of the relation.
   * @param metadata the attribute's properties
   */
  public Attribute(String attributeName, 
		   Instances header,
		   ProtectedProperties metadata) {

    if (header.numInstances() > 0) {
      throw new IllegalArgumentException("Header for relation-valued " +
                                         "attribute should not contain " +
                                         "any instances");
    }
    m_Name = attributeName;
    m_Index = -1;
    m_Values = new ArrayList<Object>();
    m_Hashtable = new Hashtable<Object,Integer>();
    m_Header = header;
    m_Type = RELATIONAL;
    setMetadata(metadata);
//    instances = new ArrayList<Double>();
  }

  /**
   * Produces a shallow copy of this attribute.
   *
   * @return a copy of this attribute with the same index
   */
  //@ also ensures \result instanceof Attribute;
  public /*@ pure non_null @*/ Object copy() {

    Attribute copy = new Attribute(m_Name);

    copy.m_Index = m_Index;
    copy.m_Type = m_Type;
    copy.m_Values = m_Values;
    copy.m_Hashtable = m_Hashtable;
    copy.m_DateFormat = m_DateFormat;
    copy.m_Header = m_Header;
    copy.setMetadata(m_Metadata);
 
    return copy;
  }

  /**
   * Returns an enumeration of all the attribute's values if the
   * attribute is a ranking, nominal, string, or relation-valued, null otherwise.
   *
   * @return enumeration of all the attribute's values
   */
  public final /*@ pure @*/ Enumeration enumerateValues() {
	  //RANKING BEGIN
    if (isRanking() || isNominal() || isString()) {
    	//RANKING END
      final Enumeration ee = new WekaEnumeration(m_Values);
      return new Enumeration () {
          public boolean hasMoreElements() {
            return ee.hasMoreElements();
          }
          public Object nextElement() {
            Object oo = ee.nextElement();
            if (oo instanceof SerializedObject) {
              return ((SerializedObject)oo).getObject();
            } else {
              return oo;
            }
          }
        };
    }
    return null;
  }

  /**
   * Tests if given attribute is equal to this attribute.
   *
   * @param other the Object to be compared to this attribute
   * @return true if the given attribute is equal to this attribute
   */
  public final /*@ pure @*/ boolean equals(Object other) {
    return (equalsMsg(other) == null);
  }

  /**
   * Tests if given attribute is equal to this attribute. If they're not
   * the same a message detailing why they differ will be returned, otherwise
   * null.
   *
   * @param other 	the Object to be compared to this attribute
   * @return 		null if the given attribute is equal to this attribute
   */
  public final String equalsMsg(Object other) {
    if (other == null)
      return "Comparing with null object";
    
    if (!(other.getClass().equals(this.getClass())))
      return "Object has wrong class";
    
    Attribute att = (Attribute) other;
    if (!m_Name.equals(att.m_Name))
      return "Names differ: " + m_Name + " != " + att.m_Name;

    if ((isRanking() && att.isRanking()) || (isNominal() && att.isNominal())) {
      if (m_Values.size() != att.m_Values.size())
        return "Different number of labels: " + m_Values.size() + " != " + att.m_Values.size();
      
      for (int i = 0; i < m_Values.size(); i++) {
        if (!m_Values.get(i).equals(att.m_Values.get(i)))
          return "Labels differ at position " + (i+1) + ": " + m_Values.get(i) + " != " + att.m_Values.get(i);
      }
      
      return null;
    } 
    
    if (isRelationValued() && att.isRelationValued())
      return m_Header.equalHeadersMsg(att.m_Header);
    
    if ((type() != att.type()))
      return "Types differ: " + typeToString(this) + " != " + typeToString(att);
    
    return null;
  }
  
  /**
   * Returns a string representation of the attribute type.
   * 
   * @param att		the attribute to return the type string for
   * @return		the string representation of the attribute type
   */
  public static String typeToString(Attribute att) {
    return typeToString(att.type());
  }
  
  /**
   * Returns a string representation of the attribute type.
   * 
   * @param type	the type of the attribute
   * @return		the string representation of the attribute type
   */
  public static String typeToString(int type) {
    String	result;
    
    switch(type) {
      case NUMERIC:
	result = "numeric";
	break;
	
      case NOMINAL:    	  
	result = "nominal";
	break;
	
      case STRING:
	result = "string";
	break;
	
      case DATE:
	result = "date";
	break;
	
      case RELATIONAL:
	result = "relational";
	break;
	
      case PreferenceAttribute.RANKING:
    result = "ranking";
    break;
    
      default:
	result = "unknown(" + type + ")";
    }
    
    return result;
  }

  /**
   * Returns the index of this attribute.
   *
   * @return the index of this attribute
   */
  //@ ensures \result == m_Index;
  public final /*@ pure @*/ int index() {

    return m_Index;
  }

  /**
   * Returns the index of a given attribute value. (The index of
   * the first occurence of this value.)
   *
   * @param value the value for which the index is to be returned
   * @return the index of the given attribute value if attribute
   * is nominal or a string, -1 if it is not or the value 
   * can't be found
   */
  public final int indexOfValue(String value) {
  //Ranking attribute added.
    if (!isRanking() && !isNominal() && !isString())
      return -1;
    Object store = value;
    if (value.length() > STRING_COMPRESS_THRESHOLD) {
      try {
        store = new SerializedObject(value, true);
      } catch (Exception ex) {
        System.err.println("Couldn't compress string attribute value -"
                           + " searching uncompressed.");
      }
    }
    Integer val = (Integer)m_Hashtable.get(store);
    if (val == null) return -1;
    else return val.intValue();
  }

  /**
   * Test if the attribute is nominal.
   *
   * @return true if the attribute is nominal
   */
  //@ ensures \result <==> (m_Type == NOMINAL);
  public final /*@ pure @*/ boolean isNominal() {

    return (m_Type == NOMINAL);
  }
  
  //RANKING BEGIN
  /**
   * Test if the attribute is a ranking.
   *
   * @return true if the attribute is a ranking.
   */
  //@ ensures \result <==> (m_Type == PreferenceAttribute.RANKING);
  public final boolean isRanking(){
	  return (m_Type == PreferenceAttribute.RANKING);
  }
  //RANKING END
  

  /**
   * Tests if the attribute is numeric.
   *
   * @return true if the attribute is numeric
   */
  //@ ensures \result <==> ((m_Type == NUMERIC) || (m_Type == DATE));
  public final /*@ pure @*/ boolean isNumeric() {

    return ((m_Type == NUMERIC) || (m_Type == DATE));
  }

  /**
   * Tests if the attribute is relation valued.
   *
   * @return true if the attribute is relation valued
   */
  //@ ensures \result <==> (m_Type == RELATIONAL);
  public final /*@ pure @*/ boolean isRelationValued() {

    return (m_Type == RELATIONAL);
  }

  /**
   * Tests if the attribute is a string.
   *
   * @return true if the attribute is a string
   */
  //@ ensures \result <==> (m_Type == STRING);
  public final /*@ pure @*/ boolean isString() {

    return (m_Type == STRING);
  }

  /**
   * Tests if the attribute is a date type.
   *
   * @return true if the attribute is a date type
   */
  //@ ensures \result <==> (m_Type == DATE);
  public final /*@ pure @*/ boolean isDate() {

    return (m_Type == DATE);
  }

  /**
   * Returns the attribute's name.
   *
   * @return the attribute's name as a string
   */
  //@ ensures \result == m_Name;
  public final /*@ pure @*/ String name() {

    return m_Name;
  }
  
  /**
   * Returns the number of attribute values. Returns 0 for 
   * attributes that are not either nominal, string, or
   * relation-valued.
   *
   * @return the number of attribute values
   */
  public final /*@ pure @*/ int numValues() {

    if (!isRanking() && !isNominal() && !isString() && !isRelationValued()) {
      return 0;
    } else {
      return m_Values.size();
    }
  }

  /**
   * Returns a description of this attribute in ARFF format. Quotes
   * strings if they contain whitespace characters, or if they
   * are a question mark.
   *
   * @return a description of this attribute as a string
   */
  public final String toString() {
    
    StringBuffer text = new StringBuffer();
    
    text.append(ARFF_ATTRIBUTE).append(" ").append(Utils.quote(m_Name)).append(" ");
    switch (m_Type) {
    case NOMINAL:
      text.append('{');
      Enumeration enu = enumerateValues();
      while (enu.hasMoreElements()) {
	text.append(Utils.quote((String) enu.nextElement()));
	if (enu.hasMoreElements())
	  text.append(',');
      }
      text.append('}');
      break;
    case NUMERIC:
      text.append(ARFF_ATTRIBUTE_NUMERIC);
      break;
    case STRING:
      text.append(ARFF_ATTRIBUTE_STRING);
      break;
    case DATE:
      text.append(ARFF_ATTRIBUTE_DATE).append(" ").append(Utils.quote(m_DateFormat.toPattern()));
      break;
    case RELATIONAL:
      text.append(ARFF_ATTRIBUTE_RELATIONAL).append("\n");
      Enumeration enm = m_Header.enumerateAttributes();
      while (enm.hasMoreElements()) {
        text.append(enm.nextElement()).append("\n");
      }
      text.append(ARFF_END_SUBRELATION).append(" ").append(Utils.quote(m_Name));
      break;
    default:
      text.append("UNKNOWN");
      break;
    }
    return text.toString();
  }

  /**
   * Returns the attribute's type as an integer.
   *
   * @return the attribute's type.
   */
  //@ ensures \result == m_Type;
  public final /*@ pure @*/ int type() {

    return m_Type;
  }
  
  /**
   * Returns the Date format pattern in case this attribute is of type DATE,
   * otherwise an empty string.
   * 
   * @return the date format pattern
   * @see SimpleDateFormat
   */
  public final String getDateFormat() {
    if (isDate())
      return m_DateFormat.toPattern();
    else
      return "";
  }
  
  /**
   * Returns a value of a nominal or string attribute.  Returns an
   * empty string if the attribute is neither a string nor a nominal
   * attribute.
   *
   * @param valIndex the value's index
   * @return the attribute's value as a string
   */
  public /*@ non_null pure @*/ String value(int valIndex) {
    
    if ( !isNominal() && !isString() && !isRanking()) {
      return "";
    } else {
      
      Object val;

        val = m_Values.get(valIndex);

      // If we're storing strings compressed, uncompress it.
      if (val instanceof SerializedObject) {
        val = ((SerializedObject)val).getObject();
      }
      return (String) val;
    }
  }
  
  /**
   * Returns the header info for a relation-valued attribute,
   * null if the attribute is not relation-valued.
   *
   * @return the attribute's value as an Instances object
   */
  public final /*@ non_null pure @*/ Instances relation() {
    
    if (!isRelationValued()) {
      return null;
    } else {
      return m_Header;
    }
  }

  /**
   * Returns a value of a relation-valued attribute. Returns
   * null if the attribute is not relation-valued.
   *
   * @param valIndex the value's index
   * @return the attribute's value as an Instances object
   */
  public final /*@ non_null pure @*/ Instances relation(int valIndex) {
    
    if (!isRelationValued()) {
      return null;
    } else {
      return (Instances) m_Values.get(valIndex);
    }
  }

  /**
   * Constructor for a numeric attribute with a particular index.
   *
   * @param attributeName the name for the attribute
   * @param index the attribute's index
   */
  //@ requires attributeName != null;
  //@ requires index >= 0;
  //@ ensures  m_Name == attributeName;
  //@ ensures  m_Index == index;
  public Attribute(String attributeName, int index) {

    this(attributeName);
    m_Index = index;
  }

  /**
   * Constructor for date attributes with a particular index.
   *
   * @param attributeName the name for the attribute
   * @param dateFormat a string suitable for use with
   * SimpleDateFormatter for parsing dates.  Null for a default format
   * string.
   * @param index the attribute's index
   */
  //@ requires attributeName != null;
  //@ requires index >= 0;
  //@ ensures  m_Name == attributeName;
  //@ ensures  m_Index == index;
  public Attribute(String attributeName, String dateFormat, 
	    int index) {

    this(attributeName, dateFormat);
    m_Index = index;
  }

  /**
   * Constructor for nominal attributes and string attributes with
   * a particular index.
   * If a null vector of attribute values is passed to the method,
   * the attribute is assumed to be a string.
   *
   * @param attributeName the name for the attribute
   * @param attributeValues a vector of strings denoting the attribute values.
   * Null if the attribute is a string attribute.
   * @param index the attribute's index
   */
  //@ requires attributeName != null;
  //@ requires index >= 0;
  //@ ensures  m_Name == attributeName;
  //@ ensures  m_Index == index;
  public Attribute(String attributeName, List<String> attributeValues, 
	    int index) {

    this(attributeName, attributeValues);
    m_Index = index;
  }

  /**
   * Constructor for a relation-valued attribute with a particular index.
   *
   * @param attributeName the name for the attribute
   * @param header the header information for this attribute
   * @param index the attribute's index
   */
  //@ requires attributeName != null;
  //@ requires index >= 0;
  //@ ensures  m_Name == attributeName;
  //@ ensures  m_Index == index;
  public Attribute(String attributeName, Instances header,
	    int index) {

    this(attributeName, header);
    m_Index = index;
  }

  /**
   * Adds a string value to the list of valid strings for attributes
   * of type STRING and returns the index of the string.
   *
   * @param value The string value to add
   * @return the index assigned to the string, or -1 if the attribute is not
   * of type Attribute.STRING 
   */
  /*@ requires value != null;
      ensures  isString() && 0 <= \result && \result < m_Values.size() ||
             ! isString() && \result == -1;
  */
  public int addStringValue(String value) {

    if (!isString()) {
      return -1;
    }
    Object store = value;

    if (value.length() > STRING_COMPRESS_THRESHOLD) {
      try {
        store = new SerializedObject(value, true);
      } catch (Exception ex) {
        System.err.println("Couldn't compress string attribute value -"
                           + " storing uncompressed.");
      }
    }
    Integer index = (Integer)m_Hashtable.get(store);
    if (index != null) {
      return index.intValue();
    } else {
      int intIndex = m_Values.size();
      m_Values.add(store);
      m_Hashtable.put(store, new Integer(intIndex));
      return intIndex;
    }
  }

  /**
   * Adds a string value to the list of valid strings for attributes
   * of type STRING and returns the index of the string. This method is
   * more efficient than addStringValue(String) for long strings.
   *
   * @param src The Attribute containing the string value to add.
   * @param index the index of the string value in the source attribute.
   * @return the index assigned to the string, or -1 if the attribute is not
   * of type Attribute.STRING 
   */
  /*@ requires src != null;
      requires 0 <= index && index < src.m_Values.size();
      ensures  isString() && 0 <= \result && \result < m_Values.size() ||
             ! isString() && \result == -1;
  */
  public int addStringValue(Attribute src, int index) {

    if (!isString()) {
      return -1;
    }
    Object store = src.m_Values.get(index);
    Integer oldIndex = (Integer)m_Hashtable.get(store);
    if (oldIndex != null) {
      return oldIndex.intValue();
    } else {
      int intIndex = m_Values.size();
      m_Values.add(store);
      m_Hashtable.put(store, new Integer(intIndex));
      return intIndex;
    }
  }

  /**
   * Adds a relation to a relation-valued attribute.
   *
   * @param value The value to add
   * @return the index assigned to the value, or -1 if the attribute is not
   * of type Attribute.RELATIONAL 
   */
  public int addRelation(Instances value) {

    if (!isRelationValued()) {
      return -1;
    }
    if (!m_Header.equalHeaders(value)) {
      throw new IllegalArgumentException("Incompatible value for " +
                                         "relation-valued attribute.\n" + 
                                         m_Header.equalHeadersMsg(value));
    }
    Integer index = (Integer)m_Hashtable.get(value);
    if (index != null) {
      return index.intValue();
    } else {
      int intIndex = m_Values.size();
      m_Values.add(value);
      m_Hashtable.put(value, new Integer(intIndex));
      return intIndex;
    }
  }

  /**
   * Adds an attribute value. Creates a fresh list of attribute
   * values before adding it.
   *
   * @param value the attribute value
   */
  final void addValue(String value) {

    m_Values = Utils.cast(m_Values.clone());
    m_Hashtable = Utils.cast(m_Hashtable.clone());
    forceAddValue(value);
  }

  /**
   * Produces a shallow copy of this attribute with a new name.
   *
   * @param newName the name of the new attribute
   * @return a copy of this attribute with the same index
   */
  //@ requires newName != null;
  //@ ensures \result.m_Name  == newName;
  //@ ensures \result.m_Index == m_Index;
  //@ ensures \result.m_Type  == m_Type;
  public final /*@ pure non_null @*/ Attribute copy(String newName) {

    Attribute copy = new Attribute(newName);

    copy.m_Index = m_Index;
    copy.m_DateFormat = m_DateFormat;
    copy.m_Type = m_Type;
    copy.m_Values = m_Values;
    copy.m_Hashtable = m_Hashtable;
    copy.m_Header = m_Header;
    copy.setMetadata(m_Metadata);
 
    return copy;
  }

  /**
   * Removes a value of a nominal, string, or relation-valued
   * attribute. Creates a fresh list of attribute values before
   * removing it.
   *
   * @param index the value's index
   * @throws IllegalArgumentException if the attribute is not 
   * of the correct type
   */
  //@ requires isRanking() || isNominal() || isString() || isRelationValued();
  //@ requires 0 <= index && index < m_Values.size();
  final void delete(int index) {
    if (!isRanking() && !isNominal() && !isString() && !isRelationValued()) 
      throw new IllegalArgumentException("Can only remove value of " +
                                         "nominal, string or relation-" +
                                         " valued attribute!");
    else {
      m_Values = Utils.cast(m_Values.clone());
      m_Values.remove(index);
      if (!isRelationValued()) {
        Hashtable<Object,Integer> hash = new Hashtable<Object,Integer>(m_Hashtable.size());
        Enumeration enu = m_Hashtable.keys();
        while (enu.hasMoreElements()) {
          Object string = enu.nextElement();
          Integer valIndexObject = (Integer)m_Hashtable.get(string);
          int valIndex = valIndexObject.intValue();
          if (valIndex > index) {
            hash.put(string, new Integer(valIndex - 1));
          } else if (valIndex < index) {
            hash.put(string, valIndexObject);
          }
        }
        m_Hashtable = hash;
      }
    }
  }

  /**
   * Adds an attribute value.
   *
   * @param value the attribute value
   */
  //@ requires value != null;
  //@ ensures  m_Values.size() == \old(m_Values.size()) + 1;
  final void forceAddValue(String value) {

    Object store = value;
    if (value.length() > STRING_COMPRESS_THRESHOLD) {
      try {
        store = new SerializedObject(value, true);
      } catch (Exception ex) {
        System.err.println("Couldn't compress string attribute value -"
                           + " storing uncompressed.");
      }
    }
    m_Values.add(store);
    m_Hashtable.put(store, new Integer(m_Values.size() - 1));
  }

  /**
   * Sets the index of this attribute.
   *
   * @param index the index of this attribute
   */
  //@ requires 0 <= index;
  //@ assignable m_Index;
  //@ ensures m_Index == index;
  final void setIndex(int index) {

    m_Index = index;
  }

  /**
   * Sets a value of a nominal attribute or string attribute.
   * Creates a fresh list of attribute values before it is set.
   *
   * @param index the value's index
   * @param string the value
   * @throws IllegalArgumentException if the attribute is not nominal or 
   * string.
   */
  //@ requires string != null;
  //@ requires isNominal() || isString() || isRanking();
  //@ requires 0 <= index && index < m_Values.size();
  final void setValue(int index, String string) {
    
    switch (m_Type) {
    case PreferenceAttribute.RANKING:
    case NOMINAL:
    case STRING:
      m_Values = Utils.cast(m_Values.clone());
      m_Hashtable = Utils.cast(m_Hashtable.clone());
      Object store = string;
      if (string.length() > STRING_COMPRESS_THRESHOLD) {
        try {
          store = new SerializedObject(string, true);
        } catch (Exception ex) {
          System.err.println("Couldn't compress string attribute value -"
                             + " storing uncompressed.");
        }
      }
      m_Hashtable.remove(m_Values.get(index));
      m_Values.set(index, store);
      m_Hashtable.put(store, new Integer(index));
      break;
    default:
      throw new IllegalArgumentException("Can only set values for nominal"
                                         + " or string attributes!");
    }
  }

  /**
   * Sets a value of a relation-valued attribute.
   * Creates a fresh list of attribute values before it is set.
   *
   * @param index the value's index
   * @param data the value
   * @throws IllegalArgumentException if the attribute is not 
   * relation-valued.
   */
  final void setValue(int index, Instances data) {
    
    if (isRelationValued()) { 
      if (!data.equalHeaders(m_Header)) {
        throw new IllegalArgumentException("Can't set relational value. " +
                                           "Headers not compatible.\n" +
                                           data.equalHeadersMsg(m_Header));
      }
      m_Values = Utils.cast(m_Values.clone());
      m_Values.set(index, data);
    } else {
      throw new IllegalArgumentException("Can only set value for"
                                         + " relation-valued attributes!");
    }
  }

  /**
   * Returns the given amount of milliseconds formatted according to the
   * current Date format.
   * 
   * @param date 	the date, represented in milliseconds since 
   * 			January 1, 1970, 00:00:00 GMT, to return as string
   * @return 		the formatted date
   */
  //@ requires isDate();
  public /*@pure@*/ String formatDate(double date) {
    switch (m_Type) {
    case DATE:
      return m_DateFormat.format(new Date((long)date));
    default:
      throw new IllegalArgumentException("Can only format date values for date"
                                         + " attributes!");
    }
  }

  /**
   * Parses the given String as Date, according to the current format and
   * returns the corresponding amount of milliseconds.
   * 
   * @param string the date to parse
   * @return the date in milliseconds since January 1, 1970, 00:00:00 GMT
   * @throws ParseException if parsing fails
   */
  //@ requires isDate();
  //@ requires string != null;
  public double parseDate(String string) throws ParseException {
    switch (m_Type) {
    case DATE:
      long time = m_DateFormat.parse(string).getTime();
      // TODO put in a safety check here if we can't store the value in a double.
      return (double)time;
    default:
      throw new IllegalArgumentException("Can only parse date values for date"
                                         + " attributes!");
    }
  }

  /**
   * Returns the properties supplied for this attribute.
   *
   * @return metadata for this attribute
   */  
  public final /*@ pure @*/ ProtectedProperties getMetadata() {

    return m_Metadata;
  }

  /**
   * Returns the ordering of the attribute. One of the following:
   * 
   * ORDERING_SYMBOLIC - attribute values should be treated as symbols.
   * ORDERING_ORDERED  - attribute values have a global ordering.
   * ORDERING_MODULO   - attribute values have an ordering which wraps.
   *
   * @return the ordering type of the attribute
   */
  public final /*@ pure @*/ int ordering() {

    return m_Ordering;
  }

  /**
   * Returns whether the attribute values are equally spaced.
   *
   * @return whether the attribute is regular or not
   */
  public final /*@ pure @*/ boolean isRegular() {

    return m_IsRegular;
  }

  /**
   * Returns whether the attribute can be averaged meaningfully.
   *
   * @return whether the attribute can be averaged or not
   */
  public final /*@ pure @*/ boolean isAveragable() {

    return m_IsAveragable;
  }

  /**
   * Returns whether the attribute has a zeropoint and may be
   * added meaningfully.
   *
   * @return whether the attribute has a zeropoint or not
   */
  public final /*@ pure @*/ boolean hasZeropoint() {

    return m_HasZeropoint;
  }

  /**
   * Returns the attribute's weight.
   *
   * @return the attribute's weight as a double
   */
  public final /*@ pure @*/ double weight() {

    return m_Weight;
  }

  /**
   * Sets the new attribute's weight
   * 
   * @param value	the new weight
   */
  public void setWeight(double value) {
    Properties	props;
    Enumeration names;
    String	name;
    
    m_Weight = value;

    // generate new metadata object
    props = new Properties();
    names = m_Metadata.propertyNames();
    while (names.hasMoreElements()) {
      name = (String) names.nextElement();
      if (!name.equals("weight"))
	props.setProperty(name, m_Metadata.getProperty(name));
    }
    props.setProperty("weight", "" + m_Weight);
    m_Metadata = new ProtectedProperties(props);
  }
  
  /**
   * Returns the lower bound of a numeric attribute.
   *
   * @return the lower bound of the specified numeric range
   */
  public final /*@ pure @*/ double getLowerNumericBound() {

    return m_LowerBound;
  }

  /**
   * Returns whether the lower numeric bound of the attribute is open.
   *
   * @return whether the lower numeric bound is open or not (closed)
   */
  public final /*@ pure @*/ boolean lowerNumericBoundIsOpen() {

    return m_LowerBoundIsOpen;
  }

  /**
   * Returns the upper bound of a numeric attribute.
   *
   * @return the upper bound of the specified numeric range
   */
  public final /*@ pure @*/ double getUpperNumericBound() {

    return m_UpperBound;
  }

  /**
   * Returns whether the upper numeric bound of the attribute is open.
   *
   * @return whether the upper numeric bound is open or not (closed)
   */
  public final /*@ pure @*/ boolean upperNumericBoundIsOpen() {

    return m_UpperBoundIsOpen;
  }

  /**
   * Determines whether a value lies within the bounds of the attribute.
   *
   * @param value the value to check
   * @return whether the value is in range
   */
  public final /*@ pure @*/ boolean isInRange(double value) {

    // dates and missing values are a special case 
    if (m_Type == DATE || Utils.isMissingValue(value)) return true;
    if (m_Type != NUMERIC) {
      // do label range check
      int intVal = (int) value;
      if (intVal < 0 || intVal >= m_Hashtable.size()) return false;
    } else {
      // do numeric bounds check
      if (m_LowerBoundIsOpen) {
	if (value <= m_LowerBound) return false;
      } else {
	if (value < m_LowerBound) return false;
      }
      if (m_UpperBoundIsOpen) {
	if (value >= m_UpperBound) return false;
      } else {
	if (value > m_UpperBound) return false;
      }
    }
    return true;
  }

  /**
   * Sets the metadata for the attribute. Processes the strings stored in the
   * metadata of the attribute so that the properties can be set up for the
   * easy-access metadata methods. Any strings sought that are omitted will
   * cause default values to be set.
   * 
   * The following properties are recognised:
   * ordering, averageable, zeropoint, regular, weight, and range.
   *
   * All other properties can be queried and handled appropriately by classes
   * calling the getMetadata() method.
   *
   * @param metadata the metadata
   * @throws IllegalArgumentException if the properties are not consistent
   */
  //@ requires metadata != null;
  protected void setMetadata(ProtectedProperties metadata) {
    
    m_Metadata = metadata;

    if (m_Type == DATE) {
      m_Ordering = ORDERING_ORDERED;
      m_IsRegular = true;
      m_IsAveragable = false;
      m_HasZeropoint = false;
    } else {

      // get ordering
      String orderString = m_Metadata.getProperty("ordering","");
      
      // numeric ordered attributes are averagable and zeropoint by default
      String def;
      if (m_Type == NUMERIC
	  && orderString.compareTo("modulo") != 0
	  && orderString.compareTo("symbolic") != 0)
	def = "true";
      else def = "false";
      
      // determine boolean states
      m_IsAveragable =
	(m_Metadata.getProperty("averageable",def).compareTo("true") == 0);
      m_HasZeropoint =
	(m_Metadata.getProperty("zeropoint",def).compareTo("true") == 0);
      // averagable or zeropoint implies regular
      if (m_IsAveragable || m_HasZeropoint) def = "true";
      m_IsRegular =
	(m_Metadata.getProperty("regular",def).compareTo("true") == 0);
      
      // determine ordering
      if (orderString.compareTo("symbolic") == 0)
	m_Ordering = ORDERING_SYMBOLIC;
      else if (orderString.compareTo("ordered") == 0)
	m_Ordering = ORDERING_ORDERED;
      else if (orderString.compareTo("modulo") == 0)
	m_Ordering = ORDERING_MODULO;
      else {
	if (m_Type == NUMERIC || m_IsAveragable || m_HasZeropoint)
	  m_Ordering = ORDERING_ORDERED;
	else m_Ordering = ORDERING_SYMBOLIC;
      }
    }

    // consistency checks
    if (m_IsAveragable && !m_IsRegular)
      throw new IllegalArgumentException("An averagable attribute must be"
					 + " regular");
    if (m_HasZeropoint && !m_IsRegular)
      throw new IllegalArgumentException("A zeropoint attribute must be"
					 + " regular");
    if (m_IsRegular && m_Ordering == ORDERING_SYMBOLIC)
      throw new IllegalArgumentException("A symbolic attribute cannot be"
					 + " regular");
    if (m_IsAveragable && m_Ordering != ORDERING_ORDERED)
      throw new IllegalArgumentException("An averagable attribute must be"
					 + " ordered");
    if (m_HasZeropoint && m_Ordering != ORDERING_ORDERED)
      throw new IllegalArgumentException("A zeropoint attribute must be"
					 + " ordered");

    // determine weight
    m_Weight = 1.0;
    String weightString = m_Metadata.getProperty("weight");
    if (weightString != null) {
      try{
	m_Weight = Double.valueOf(weightString).doubleValue();
      } catch (NumberFormatException e) {
	// Check if value is really a number
	throw new IllegalArgumentException("Not a valid attribute weight: '" 
					   + weightString + "'");
      }
    }

    // determine numeric range
    if (m_Type == NUMERIC) setNumericRange(m_Metadata.getProperty("range"));
  }

  /**
   * Sets the numeric range based on a string. If the string is null the range
   * will default to [-inf,+inf]. A square brace represents a closed interval, a
   * curved brace represents an open interval, and 'inf' represents infinity.
   * Examples of valid range strings: "[-inf,20)","(-13.5,-5.2)","(5,inf]"
   *
   * @param rangeString the string to parse as the attribute's numeric range
   * @throws IllegalArgumentException if the range is not valid
   */
  //@ requires rangeString != null;
  private void setNumericRange(String rangeString)
  {
    // set defaults
    m_LowerBound = Double.NEGATIVE_INFINITY;
    m_LowerBoundIsOpen = false;
    m_UpperBound = Double.POSITIVE_INFINITY;
    m_UpperBoundIsOpen = false;

    if (rangeString == null) return;

    // set up a tokenzier to parse the string
    StreamTokenizer tokenizer =
      new StreamTokenizer(new StringReader(rangeString));
    tokenizer.resetSyntax();         
    tokenizer.whitespaceChars(0, ' ');    
    tokenizer.wordChars(' '+1,'\u00FF');
    tokenizer.ordinaryChar('[');
    tokenizer.ordinaryChar('(');
    tokenizer.ordinaryChar(',');
    tokenizer.ordinaryChar(']');
    tokenizer.ordinaryChar(')');

    try {

      // get opening brace
      tokenizer.nextToken();
    
      if (tokenizer.ttype == '[') m_LowerBoundIsOpen = false;
      else if (tokenizer.ttype == '(') m_LowerBoundIsOpen = true;
      else throw new IllegalArgumentException("Expected opening brace on range,"
					      + " found: "
					      + tokenizer.toString());

      // get lower bound
      tokenizer.nextToken();
      if (tokenizer.ttype != tokenizer.TT_WORD)
	throw new IllegalArgumentException("Expected lower bound in range,"
					   + " found: "
					   + tokenizer.toString());
      if (tokenizer.sval.compareToIgnoreCase("-inf") == 0)
	m_LowerBound = Double.NEGATIVE_INFINITY;
      else if (tokenizer.sval.compareToIgnoreCase("+inf") == 0)
	m_LowerBound = Double.POSITIVE_INFINITY;
      else if (tokenizer.sval.compareToIgnoreCase("inf") == 0)
	m_LowerBound = Double.NEGATIVE_INFINITY;
      else try {
	m_LowerBound = Double.valueOf(tokenizer.sval).doubleValue();
      } catch (NumberFormatException e) {
	throw new IllegalArgumentException("Expected lower bound in range,"
					   + " found: '" + tokenizer.sval + "'");
      }

      // get separating comma
      if (tokenizer.nextToken() != ',')
	throw new IllegalArgumentException("Expected comma in range,"
					   + " found: "
					   + tokenizer.toString());

      // get upper bound
      tokenizer.nextToken();
      if (tokenizer.ttype != tokenizer.TT_WORD)
	throw new IllegalArgumentException("Expected upper bound in range,"
					   + " found: "
					   + tokenizer.toString());
      if (tokenizer.sval.compareToIgnoreCase("-inf") == 0)
	m_UpperBound = Double.NEGATIVE_INFINITY;
      else if (tokenizer.sval.compareToIgnoreCase("+inf") == 0)
	m_UpperBound = Double.POSITIVE_INFINITY;
      else if (tokenizer.sval.compareToIgnoreCase("inf") == 0)
	m_UpperBound = Double.POSITIVE_INFINITY;
      else try {
	m_UpperBound = Double.valueOf(tokenizer.sval).doubleValue();
      } catch (NumberFormatException e) {
	throw new IllegalArgumentException("Expected upper bound in range,"
					   + " found: '" + tokenizer.sval + "'");
      }

      // get closing brace
      tokenizer.nextToken();
    
      if (tokenizer.ttype == ']') m_UpperBoundIsOpen = false;
      else if (tokenizer.ttype == ')') m_UpperBoundIsOpen = true;
      else throw new IllegalArgumentException("Expected closing brace on range,"
					      + " found: "
					      + tokenizer.toString());

      // check for rubbish on end
      if (tokenizer.nextToken() != tokenizer.TT_EOF)
	throw new IllegalArgumentException("Expected end of range string,"
					   + " found: "
					   + tokenizer.toString());

    } catch (IOException e) {
      throw new IllegalArgumentException("IOException reading attribute range"
					 + " string: " + e.getMessage());
    }

    if (m_UpperBound < m_LowerBound)
      throw new IllegalArgumentException("Upper bound (" + m_UpperBound
					 + ") on numeric range is"
					 + " less than lower bound ("
					 + m_LowerBound + ")!");
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
   * Simple main method for testing this class.
   * 
   * @param ops the commandline options
   */
  //@ requires ops != null;
  //@ requires \nonnullelements(ops);
  public static void main(String[] ops) {

    try {
      
      // Create numeric attributes "length" and "weight"
      Attribute length = new Attribute("length");
      Attribute weight = new Attribute("weight");

      // Create date attribute "date"
      Attribute date = new Attribute("date", "yyyy-MM-dd HH:mm:ss");

      System.out.println(date);
      double dd = date.parseDate("2001-04-04 14:13:55");
      System.out.println("Test date = " + dd);
      System.out.println(date.formatDate(dd));

      dd = new Date().getTime();
      System.out.println("Date now = " + dd);
      System.out.println(date.formatDate(dd));
      
      // Create vector to hold nominal values "first", "second", "third" 
      List<String> my_nominal_values = new ArrayList<String>(3); 
      my_nominal_values.add("first"); 
      my_nominal_values.add("second"); 
      my_nominal_values.add("third"); 
      
      // Create nominal attribute "position" 
      Attribute position = new Attribute("position", my_nominal_values);

      // Print the name of "position"
      System.out.println("Name of \"position\": " + position.name());

      // Print the values of "position"
      Enumeration attValues = position.enumerateValues();
      while (attValues.hasMoreElements()) {
	String string = (String)attValues.nextElement();
	System.out.println("Value of \"position\": " + string);
      }

      // Shallow copy attribute "position"
      Attribute copy = (Attribute) position.copy();

      // Test if attributes are the same
      System.out.println("Copy is the same as original: " + copy.equals(position));

      // Print index of attribute "weight" (should be unset: -1)
      System.out.println("Index of attribute \"weight\" (should be -1): " + 
			 weight.index());

      // Print index of value "first" of attribute "position"
      System.out.println("Index of value \"first\" of \"position\" (should be 0): " +
			 position.indexOfValue("first"));

      // Tests type of attribute "position"
      System.out.println("\"position\" is numeric: " + position.isNumeric());
      System.out.println("\"position\" is ranking: " + position.isRanking());
      System.out.println("\"position\" is nominal: " + position.isNominal());      
      System.out.println("\"position\" is string: " + position.isString());

      // Prints name of attribute "position"
      System.out.println("Name of \"position\": " + position.name());
    
      // Prints number of values of attribute "position"
      System.out.println("Number of values for \"position\": " + position.numValues());

      // Prints the values (againg)
      for (int i = 0; i < position.numValues(); i++) {
	System.out.println("Value " + i + ": " + position.value(i));
      }

      // Prints the attribute "position" in ARFF format
      System.out.println(position);

      // Checks type of attribute "position" using constants
      switch (position.type()) {
      case Attribute.NUMERIC:
	System.out.println("\"position\" is numeric");
	break;
      case Attribute.NOMINAL:
	System.out.println("\"position\" is nominal");
	break;
      case PreferenceAttribute.RANKING:
    System.out.println("\"position\" is a ranking");
    break;
      case Attribute.STRING:
	System.out.println("\"position\" is string");
	break;
      case Attribute.DATE:
	System.out.println("\"position\" is date");
	break;
      case Attribute.RELATIONAL:
	System.out.println("\"position\" is relation-valued");
	break;
      default:
	System.out.println("\"position\" has unknown type");
      }

      ArrayList<Attribute> atts = new ArrayList<Attribute>(1);
      atts.add(position);
      Instances relation = new Instances("Test", atts, 0);
      Attribute relationValuedAtt = new Attribute("test", relation);
      System.out.println(relationValuedAtt);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
  
