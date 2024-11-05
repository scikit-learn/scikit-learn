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
 *    TargetMetaInfo.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.io.Serializable;
import java.util.ArrayList;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Utils;

/**
 * Class to encapsulate information about a Target.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision 1.0 $
 */
public class TargetMetaInfo extends FieldMetaInfo implements Serializable {

  /** For serialization */
  private static final long serialVersionUID = 863500462237904927L;

  /** min and max */
  protected double m_min = Double.NaN;
  protected double m_max = Double.NaN;

  /** re-scaling of target value (if defined) */
  protected double m_rescaleConstant = 0;
  protected double m_rescaleFactor = 1.0;

  /** cast integers (default no casting) */
  protected String m_castInteger = "";

  // -------------------------------------------------------

  /** default value (numeric) or prior distribution (categorical) */
  protected double[] m_defaultValueOrPriorProbs;

  /**  for categorical values. Actual values */
  protected ArrayList<String> m_values = new ArrayList<String>();
  
  /** corresponding display values */
  protected ArrayList<String> m_displayValues = new ArrayList<String>();


  // TODO: toString method.

  /**
   * Constructor.
   *
   * @param target the <code>Element</code> encapsulating a Target
   * @throws Exception if there is a problem reading the Target
   */
  protected TargetMetaInfo(Element target) throws Exception {
    super(target);

    // check for an OPTYPE
    /*String op = target.getAttribute("optype");
    if (op != null && op.length() > 0) {
      for (int i = 0; i < Optype.values().length; i++) {
        if (op.equals(Optype.values()[i].toString())) {
          m_optype = Optype.values()[i];
          break;
        }
      }
    }*/

    // min and max (if defined)
    String min = target.getAttribute("min");
    if (min != null && min.length() > 0) {
      try {
        m_min = Double.parseDouble(min);
      } catch (IllegalArgumentException ex) {
        throw new Exception("[TargetMetaInfo] can't parse min value for target field "
                            + m_fieldName);
      }
    }

    String max = target.getAttribute("max");
    if (max != null && max.length() > 0) {
      try {
        m_max = Double.parseDouble(max);
      } catch (IllegalArgumentException ex) {
        throw new Exception("[TargetMetaInfo] can't parse max value for target field "
                            + m_fieldName);
      }
    }

    // Re-scaling (if any)
    String rsc = target.getAttribute("rescaleConstant");
    if (rsc != null && rsc.length() > 0) {
      try {
        m_rescaleConstant = Double.parseDouble(rsc);
      } catch (IllegalArgumentException ex) {
        throw new Exception("[TargetMetaInfo] can't parse rescale constant value for "
                            + "target field " + m_fieldName);
      }
    }
    String rsf = target.getAttribute("rescaleFactor");
    if (rsf != null && rsf.length() > 0) {
      try {
        m_rescaleFactor = Double.parseDouble(rsf);
      } catch (IllegalArgumentException ex) {
        throw new Exception("[TargetMetaInfo] can't parse rescale factor value for "
                            + "target field " + m_fieldName);
      }
    }

    // Cast integers
    String cstI = target.getAttribute("castInteger");
    if (cstI != null && cstI.length() > 0) {
      m_castInteger = cstI;
    }
    
    // Get the target value(s). Apparently, there doesn't have to
    // be any target values defined.
    NodeList vals = target.getElementsByTagName("TargetValue");
    if (vals.getLength() > 0) {
      m_defaultValueOrPriorProbs = new double[vals.getLength()];
      
      for (int i = 0; i < vals.getLength(); i++) {
        Node value = vals.item(i);
        if (value.getNodeType() == Node.ELEMENT_NODE) {
          Element valueE = (Element)value;
          String valueName = valueE.getAttribute("value");
          if (valueName != null && valueName.length() > 0) {
            // we have a categorical value - set optype if it's not
            // already set
            if (m_optype != Optype.CATEGORICAL &&
                m_optype != Optype.NONE) {
              throw new Exception("[TargetMetaInfo] TargetValue element has categorical value but "
                                  + "optype is not categorical!");
            }

            if (m_optype == Optype.NONE) {
              m_optype = Optype.CATEGORICAL;
            }

            m_values.add(valueName);
            // get display value (if any)
            String displayValue = valueE.getAttribute("displayValue");
            if (displayValue != null && displayValue.length() > 0) {
              m_displayValues.add(displayValue);
            } else {
              // use the value as the display value
              m_displayValues.add(valueName);
            }

            // get prior probability (should be defined!!)
            String prior = valueE.getAttribute("priorProbability");
            if (prior != null && prior.length() > 0) {
              try {
                m_defaultValueOrPriorProbs[i] = Double.parseDouble(prior);
              } catch (IllegalArgumentException ex) {
                throw new Exception("[TargetMetaInfo] Can't parse probability from "
                                    + "TargetValue element.");
              }
            } else {
              throw new Exception("[TargetMetaInfo] No prior probability defined for value "
                                  + valueName);
            }
          } else {
            // we have a numeric field
            // check the optype
            if (m_optype != Optype.CONTINUOUS &&
                m_optype != Optype.NONE) {
              throw new Exception("[TargetMetaInfo] TargetValue element has continuous value but "
                                  + "optype is not continuous!");
            }

            if (m_optype == Optype.NONE) {
              m_optype = Optype.CONTINUOUS;
            }

            // get the default value
            String defaultV = valueE.getAttribute("defaultValue");
            if (defaultV != null && defaultV.length() > 0) {
              try {
                m_defaultValueOrPriorProbs[i] = Double.parseDouble(defaultV);
              } catch (IllegalArgumentException ex) {
                throw new Exception("[TargetMetaInfo] Can't parse default value from "
                                    + "TargetValue element.");
              }
            } else {
              throw new Exception("[TargetMetaInfo] No default value defined for target "
                                  + m_fieldName);
            }
          }
        }
      }
    }
  }

  /**
   * Get the prior probability for the supplied value.
   * 
   * @param value the value to get the probability for
   * @return the probability
   * @throws Exception if there are no TargetValues defined or
   * if the supplied value is not in the list of TargetValues
   */
  public double getPriorProbability(String value) throws Exception {
    if (m_defaultValueOrPriorProbs == null) {
      throw new Exception("[TargetMetaInfo] no TargetValues defined (getPriorProbability)");
    }
    double result = Double.NaN;
    boolean found = false;
    for (int i = 0; i < m_values.size(); i++) {
      if (value.equals(m_values.get(i))) {
        found = true;
        result = m_defaultValueOrPriorProbs[i];
        break;
      }
    }
    if (!found) {
      throw new Exception("[TargetMetaInfo] couldn't find value " + value 
                          + "(getPriorProbability)");
    }
    return result;
  }

  /**
   * Get the default value (numeric target)
   *
   * @return the default value
   * @throws Exception if there is no TargetValue defined
   */
  public double getDefaultValue() throws Exception {
    if (m_defaultValueOrPriorProbs == null) {
      throw new Exception("[TargetMetaInfo] no TargetValues defined (getPriorProbability)");
    }
    return m_defaultValueOrPriorProbs[0];
  }

  /**
   * Get the values (discrete case only) for this Target. Note: the
   * list may be empty if the pmml doesn't specify any values.
   *
   * @return the values of this Target
   */
  public ArrayList<String> getValues() {
    return new ArrayList<String>(m_values);
  }

  /**
   * Apply min and max, rescaleFactor, rescaleConstant and castInteger - in
   * that order (where defined).
   *
   * @param prediction the prediction to apply these modification to
   * @return the modified prediction
   * @throws Exception if this target is not a continuous one
   */
  public double applyMinMaxRescaleCast(double prediction) throws Exception {
    if (m_optype != Optype.CONTINUOUS) {
      throw new Exception("[TargetMetaInfo] target must be continuous!");
    }

    if (!Utils.isMissingValue(m_min) && prediction < m_min) {
      prediction = m_min;
    }
    if (!Utils.isMissingValue(m_max) && prediction > m_max) {
      prediction = m_max;
    }

    prediction *= m_rescaleFactor;
    prediction += m_rescaleConstant;

    if (m_castInteger.length() > 0) {
      if (m_castInteger.equals("round")) {
        prediction = Math.round(prediction);
      } else if (m_castInteger.equals("ceiling")) {
        prediction = Math.ceil(prediction);
      } else if (m_castInteger.equals("floor")) {
        prediction = Math.floor(prediction);
      } else {
        throw new Exception("[TargetMetaInfo] unknown castInteger value "
                            + m_castInteger);
      }
    }
    
    return prediction;
  }
  
  /**
   * Return this field as an Attribute.
   * 
   * @return an Attribute for this field.
   */
  public Attribute getFieldAsAttribute() {
    if (m_optype == Optype.CONTINUOUS) {
      return new Attribute(m_fieldName);
    }
    if (m_values.size() == 0) {
      // return a String attribute
      return new Attribute(m_fieldName, (ArrayList<String>)null);
    }
    
    ArrayList<String> values = new ArrayList<String>();
    for (String val : m_values) {
      values.add(val);
    }
    return new Attribute(m_fieldName, values);
  }
}
