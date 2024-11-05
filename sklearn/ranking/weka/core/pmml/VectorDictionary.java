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
 *    VectorDictionary.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.core.Attribute;
import weka.core.Instances;

/**
 * Class encapsulating the PMML VectorDictionary construct.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6466 $
 */
public class VectorDictionary implements Serializable {
  
  /**
   * For serialization
   */
  private static final long serialVersionUID = -5538024467333813123L;

  /** The number of support vectors in the dictionary */
  protected int m_numberOfVectors;
  
  /** The fields accessed by the support vectors **/
  protected List<FieldRef> m_vectorFields = new ArrayList<FieldRef>();
  
  /** The vectors in the dictionary */
  protected Map<String, VectorInstance> m_vectorInstances = 
    new HashMap<String, VectorInstance>();
  
  /**
   * Returns a new VectorDictionary constructed from the supplied XML container
   * 
   * @param container the containing XML
   * @param ms the mining schema
   * @return a VectorDictionary
   * @throws Exception if the VectorDictionary can't be read from the XML container
   */
  public static VectorDictionary getVectorDictionary(Element container, MiningSchema ms) 
    throws Exception {

    VectorDictionary vectDict = null;

    NodeList vecL = container.getElementsByTagName("VectorDictionary");
    if (vecL.getLength() > 0) {
      Node vecNode = vecL.item(0);
      if (vecNode.getNodeType() == Node.ELEMENT_NODE) {
        vectDict = new VectorDictionary((Element)vecNode, ms);
      }
    }
    
    return vectDict;
  }
  
  /**
   * Convert an incoming instance to an array of values that corresponds
   * to the fields referenced by the support vectors in the vector dictionary
   * 
   * @param incoming an incoming instance
   * @return an array of values from the incoming instance that corresponds
   * to just the fields referenced by the support vectors
   * @throws Exception if this array cant be constructed for some reason
   */
  public double[] incomingInstanceToVectorFieldVals(double[] incoming) 
    throws Exception {
    double[] newInst = new double[m_vectorFields.size()];
    
    for (int i = 0; i < m_vectorFields.size(); i++) {
      FieldRef fr = m_vectorFields.get(i);
      newInst[i] = fr.getResult(incoming);
    }
    
    return newInst;
  }
  
  
  /**
   * Constructor.
   * 
   * @param vectNode the XML containing the VectorDictionary
   * @param ms the mining schema
   * @throws Exception if something goes wrong
   */
  public VectorDictionary(Element vectNode, MiningSchema ms) throws Exception {
    NodeList vecFieldsL = vectNode.getElementsByTagName("VectorFields");
    if (vecFieldsL.getLength() == 0) {
      throw new Exception("[VectorDictionary] there are no VectorFields defined!!");
    }
    
    Instances fullStructure = ms.getFieldsAsInstances();
    ArrayList<Attribute> fieldDefs = new ArrayList<Attribute>();
    for (int i = 0; i < fullStructure.numAttributes(); i++) {
      fieldDefs.add(fullStructure.attribute(i));
    }
    
    // should be just one VectorFields element
    Node fieldsNode = vecFieldsL.item(0);
    
    NodeList fieldRefsL = ((Element)fieldsNode).getElementsByTagName("FieldRef");
    
    // build the list of field refs
    for (int i = 0; i < fieldRefsL.getLength(); i++) {
      Element fieldR = (Element)fieldRefsL.item(i);
      String fieldName = fieldR.getAttribute("field");
      Attribute a = fullStructure.attribute(fieldName);
      
      if (a == null) {
        throw new Exception("[VectorDictionary] can't find field '" + fieldName 
            + "' in the mining schema/derived fields!");
      }
      
      FieldMetaInfo.Optype fieldOpt = (a.isNumeric())
        ? FieldMetaInfo.Optype.CONTINUOUS
        : FieldMetaInfo.Optype.CATEGORICAL;
      
      FieldRef fr = new FieldRef(fieldR, fieldOpt, fieldDefs);
      m_vectorFields.add(fr);
    }
    
    // now get the support vectors
    NodeList vecInstL = vectNode.getElementsByTagName("VectorInstance");
    
    if (vecInstL.getLength() == 0) {
      throw new Exception("[VectorDictionary] no VectorInstances defined!");
    }
    
    for (int i = 0; i < vecInstL.getLength(); i++) {
      Element vecInstEl = (Element)vecInstL.item(i);
      VectorInstance temp = new VectorInstance(vecInstEl, m_vectorFields);
      String id = temp.getID();
      if (m_vectorInstances.get(id) != null) {
        throw new Exception("[VectorDictionary] : There is already a vector with ID " 
            + id + " in the dictionary!");
      }
      m_vectorInstances.put(id, temp);
    }
  }
  
  /**
   * Gets a vector from the dictionary corresponding to the supplied ID
   * 
   * @param ID the ID of the vector to retrieve
   * @return the vector with the given ID or null if no vector with
   * that ID exists in the dictionary
   */
  public VectorInstance getVector(String ID) {
    return m_vectorInstances.get(ID);
  }
}
