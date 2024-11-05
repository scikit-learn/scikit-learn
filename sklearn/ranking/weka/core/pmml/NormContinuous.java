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
 *    NormContinuous.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.util.ArrayList;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Utils;


/**
 * Class encapsulating a NormContinuous Expression.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision 1.0 $
 */
public class NormContinuous extends Expression {
  
  /**
   * For serialization
   */
  private static final long serialVersionUID = 4714332374909851542L;

  /** The name of the field to use */
  protected String m_fieldName;
  
  /** The index of the field */
  protected int m_fieldIndex;
  
  /** True if a replacement for missing values has been specified */
  protected boolean m_mapMissingDefined = false;
  
  /** The value of the missing value replacement (if defined) */
  protected double m_mapMissingTo;
  
  /** Outlier treatment method (default = asIs) */
  protected MiningFieldMetaInfo.Outlier m_outlierTreatmentMethod =
    MiningFieldMetaInfo.Outlier.ASIS;
  
  /** original values for the LinearNorm entries */
  protected double[] m_linearNormOrig;
  
  /** norm values for the LinearNorm entries */
  protected double[] m_linearNormNorm;
  
  public NormContinuous(Element normCont, FieldMetaInfo.Optype opType, ArrayList<Attribute> fieldDefs) 
    throws Exception {
    super(opType, fieldDefs);
    
    if (opType != FieldMetaInfo.Optype.CONTINUOUS) {
      throw new Exception("[NormContinuous] can only have a continuous optype");
    }
    
    m_fieldName = normCont.getAttribute("field");
    
    String mapMissing = normCont.getAttribute("mapMissingTo");
    if (mapMissing != null && mapMissing.length() > 0) {
      m_mapMissingTo = Double.parseDouble(mapMissing);
      m_mapMissingDefined = true;
    }
    
    String outliers = normCont.getAttribute("outliers");
    if (outliers != null && outliers.length() > 0) {
      for (MiningFieldMetaInfo.Outlier o : MiningFieldMetaInfo.Outlier.values()) {
        if (o.toString().equals(outliers)) {
          m_outlierTreatmentMethod = o;
          break;
        }
      }
    }
    
    // get the LinearNorm elements
    NodeList lnL = normCont.getElementsByTagName("LinearNorm");
    if (lnL.getLength() < 2) {
      throw new Exception("[NormContinuous] Must be at least 2 LinearNorm elements!");
    }
    m_linearNormOrig = new double[lnL.getLength()];
    m_linearNormNorm = new double[lnL.getLength()];
    
    for (int i = 0; i < lnL.getLength(); i++) {
      Node lnN = lnL.item(i);
      if (lnN.getNodeType() == Node.ELEMENT_NODE) {
        Element lnE = (Element)lnN;
        
        String orig = lnE.getAttribute("orig");
        m_linearNormOrig[i] = Double.parseDouble(orig);
        
        String norm = lnE.getAttribute("norm");
        m_linearNormNorm[i] = Double.parseDouble(norm);
      }
    }
    
    if (fieldDefs != null) {
      setUpField();
    }
  }
  
  /**
   * Set the field definitions for this Expression to use
   * 
   * @param fieldDefs the field definitions to use
   * @throws Exception if there is a problem setting the field definitions
   */
  public void setFieldDefs(ArrayList<Attribute> fieldDefs) throws Exception {
    super.setFieldDefs(fieldDefs);
    setUpField();
  }
  
  private void setUpField() throws Exception {
    m_fieldIndex = -1;
    
    if (m_fieldDefs != null) {
      m_fieldIndex = getFieldDefIndex(m_fieldName);
//      System.err.println("NormCont... index of " + m_fieldName + " " + m_fieldIndex);
      if (m_fieldIndex < 0) {
        throw new Exception("[NormContinuous] Can't find field " + m_fieldName
            + " in the supplied field definitions.");
      }
      
      Attribute field = m_fieldDefs.get(m_fieldIndex);
      if (!field.isNumeric()) {
        throw new Exception("[NormContinuous] reference field " + m_fieldName
            +" must be continuous.");
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
    return new Attribute(m_fieldName + "_normContinuous");
  }

  /**
   * Get the result of evaluating the expression. In the case
   * of a continuous optype, a real number is returned; in
   * the case of a categorical/ordinal optype, the index of the nominal
   * value is returned as a double.
   * 
   * @param incoming the incoming parameter values
   * @return the result of normalizing the input field
   * @throws Exception if there is a problem computing the result
   */
  public double getResult(double[] incoming) throws Exception {
    
    double[] a = m_linearNormOrig;
    double[] b = m_linearNormNorm;
    
    return computeNorm(a, b, incoming);
  }
  
  /**
   * Compute the inverse of the normalization (i.e. map back to a unormalized value).
   * 
   * @param incoming the incoming parameter values
   * @return the unormalized value
   */
  public double getResultInverse(double[] incoming) {
    double[] a = m_linearNormNorm;
    double[] b = m_linearNormOrig;
    
    return computeNorm(a, b, incoming);
  }
  
  private double computeNorm(double[] a, double[] b, double[] incoming) {
    double result = 0.0;
    
    if (Utils.isMissingValue(incoming[m_fieldIndex])) {
      if (m_mapMissingDefined) {
        result = m_mapMissingTo;
      } else {
        result = incoming[m_fieldIndex]; // just return the missing value
      }
    } else {
      double x = incoming[m_fieldIndex];
      /*System.err.println("NormCont (index): " + m_fieldIndex);
      System.err.println("NormCont (input val): " + x); */
      
      // boundary cases first
      if (x < a[0]) {
        if (m_outlierTreatmentMethod == MiningFieldMetaInfo.Outlier.ASIS) {
          double slope = (b[1] - b[0]) /
            (a[1] - a[0]);
          double offset = b[0] - (slope * a[0]);
          result = slope * x + offset;
        } else if (m_outlierTreatmentMethod == MiningFieldMetaInfo.Outlier.ASEXTREMEVALUES) {
          result = b[0];
        } else {
          // map to missing replacement value
          result = m_mapMissingTo;
        }
      } else if (x > a[a.length - 1]) {
        int length = a.length;
        if (m_outlierTreatmentMethod == MiningFieldMetaInfo.Outlier.ASIS) {
          double slope = (b[length - 1] - b[length - 2]) /
            (a[length - 1] - a[length - 2]);
          double offset = b[length - 1] - (slope * a[length - 1]);
          result = slope * x + offset;
        } else if (m_outlierTreatmentMethod == MiningFieldMetaInfo.Outlier.ASEXTREMEVALUES) {
          result = b[length - 1];
        } else {
          // map to missing replacement value
          result = m_mapMissingTo;
        }
      } else {
        // find the segment that this value falls in to
        for (int i = 1; i < a.length; i++) {
          if (x <= a[i]) {
            result = b[i - 1];
            result += ((x - a[i - 1])/(a[i] - a[i - 1]) * 
                        (b[i] - b[i - 1]));
            break;
          }
        }
      }
    }
    return result;
  }

  /**
   * Always throws an Exception since the result of NormContinuous must
   * be continuous.
   * 
   * @param incoming the incoming parameter values
   * @throws Exception always
   */
  public String getResultCategorical(double[] incoming) throws Exception {
    throw new Exception("[NormContinuous] Can't return the result as a categorical value!");
  }
  
  public String toString(String pad) {
    StringBuffer buff = new StringBuffer();
    
    buff.append(pad + "NormContinuous (" + m_fieldName + "):\n" + pad + "linearNorm: ");
    for (int i = 0; i < m_linearNormOrig.length; i++) {
      buff.append("" + m_linearNormOrig[i] + ":" + m_linearNormNorm[i] + " ");
    }
    buff.append("\n" + pad);
    buff.append("outlier treatment: " + m_outlierTreatmentMethod.toString());
    if (m_mapMissingDefined) {
      buff.append("\n" + pad);
      buff.append("map missing values to: " + m_mapMissingTo);
    }
    
    return buff.toString();
  }

}
