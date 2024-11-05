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
 *    Constant.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.util.ArrayList;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.core.Attribute;

/**
 * Class encapsulating a Constant Expression.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com
 * @version $Revision 1.0 $
 */
public class Constant extends Expression {
  
  /**
   * For serialization 
   */
  private static final long serialVersionUID = -304829687822452424L;
  
  protected String m_categoricalConst = null;
  protected double m_continuousConst = Double.NaN;

  /**
   * Construct an new Constant Expression.
   * 
   * @param constant the xml Element containing the Constant
   * @param opType the optype for the Constant
   * @param fieldDefs an ArrayList of Attributes for the fields that this
   * Expression might need to access (not needed for a constant!)
   * @throws Exception if the optype is specified as continuous
   * and there is a problem parsing the value of the Constant
   */
  public Constant(Element constant, FieldMetaInfo.Optype opType, ArrayList<Attribute> fieldDefs) 
    throws Exception {
    super(opType, fieldDefs);
    
    NodeList constL = constant.getChildNodes();
    String c = constL.item(0).getNodeValue();
    
    if (m_opType == FieldMetaInfo.Optype.CATEGORICAL ||
        m_opType == FieldMetaInfo.Optype.ORDINAL) {
      m_categoricalConst = c;
    } else {
      try {
        m_continuousConst = Double.parseDouble(c);
      } catch (IllegalArgumentException ex) {
        throw new Exception("[Constant] Unable to parse continuous constant: "
            + c);
      }
    }
  }
  
  /**
   * Return the structure of the result of applying this Expression
   * as an Attribute.
   * 
   * @return the structure of the result of applying this Expression as an
   * Attribute.
   */
  protected Attribute getOutputDef() {
    
    if (m_opType == FieldMetaInfo.Optype.CONTINUOUS) {
      return new Attribute("Constant: " + m_continuousConst);
    }
    
    ArrayList<String> nom = new ArrayList<String>();
    nom.add(m_categoricalConst);
    return new Attribute("Constant: " + m_categoricalConst, nom);
  }
  
  /**
   * Get the result of evaluating the expression. In the case
   * of a continuous optype, a real number is returned; in
   * the case of a categorical/ordinal optype, the index of the nominal
   * value is returned as a double.
   * 
   * @param incoming the incoming parameter values
   * @return the result of evaluating the expression
   */
  public double getResult(double[] incoming) {
    if (m_opType == FieldMetaInfo.Optype.CONTINUOUS) {
      return m_continuousConst;
    }
    return 0; // constant (first and only value of a nominal attribute)
  }
  
  /**
   * Gets the result of evaluating the expression when the
   * optype is categorical or ordinal as the actual String
   * value.
   * 
   * @param incoming the incoming parameter values 
   * @return the result of evaluating the expression
   * @throws Exception if the optype is continuous
   */
  public String getResultCategorical(double[] incoming) 
    throws Exception {
    if (m_opType == FieldMetaInfo.Optype.CONTINUOUS) {
      throw new IllegalArgumentException("[Constant] Cant't return result as "
          +"categorical/ordinal as optype is continuous!");
    }
    return m_categoricalConst;
  }
  
  public static void main(String[] args) {
    try {
      java.io.File f = new java.io.File(args[0]);
      javax.xml.parsers.DocumentBuilderFactory dbf = javax.xml.parsers.DocumentBuilderFactory.newInstance();
      javax.xml.parsers.DocumentBuilder db = dbf.newDocumentBuilder();
      org.w3c.dom.Document doc = db.parse(f);
      doc.getDocumentElement().normalize();
      NodeList constL = doc.getElementsByTagName("Constant");
      Node c = constL.item(0);
      
      if (c.getNodeType() == Node.ELEMENT_NODE) {
        Constant constC = new Constant((Element)c, FieldMetaInfo.Optype.CONTINUOUS, null);
        System.err.println("Value of first constant: " + constC.getResult(null));
      }
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
  
  public String toString(String pad) {
    return pad + "Constant: " + ((m_categoricalConst != null)
        ? m_categoricalConst
        : "" + m_continuousConst); 
  }
}
