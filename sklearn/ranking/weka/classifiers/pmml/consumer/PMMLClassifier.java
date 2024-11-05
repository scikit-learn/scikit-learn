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
 *    PMMLClassifier.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.pmml.consumer;

import java.io.Serializable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Instances;
import weka.core.pmml.*;
import weka.gui.Logger;

/**
 * Abstract base class for all PMML classifiers.
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 5928 $
 */
public abstract class PMMLClassifier extends AbstractClassifier
  implements Serializable, PMMLModel {
  
  /** For serialization */
  private static final long serialVersionUID = -5371600590320702971L;

  /** PMML version */
  protected String m_pmmlVersion = "?";
  
  /** Creator application */
  protected String m_creatorApplication = "?";
  
  /** Logger */
  protected Logger m_log = null;

  /** The data dictionary */
  protected Instances m_dataDictionary;

  /** The fields and meta data used by the model */
  protected MiningSchema m_miningSchema;

  /** The mapping between mining schema fields and incoming instance
      attributes */
  protected transient MappingInfo m_fieldsMap;

  /** Has the classifier been initialized (i.e. have we established
      a mapping between the mining schema and the incoming instances)? */
  protected transient boolean m_initialized = false;

  /**
   * Constructor.
   *
   * @param dataDictionary the data dictionary
   * @param miningSchema the mining schema
   */
   PMMLClassifier(Instances dataDictionary,
                        MiningSchema miningSchema) {
    m_dataDictionary = dataDictionary;
    m_miningSchema = miningSchema;
  }

  /**
   * Set the version of PMML used for this model.
   *
   * @param doc the Document encapsulating the pmml
   */
  public void setPMMLVersion(Document doc) {
    NodeList tempL = doc.getElementsByTagName("PMML");
    Node pmml = tempL.item(0);
    if (pmml.getNodeType() == Node.ELEMENT_NODE) {
      String version = ((Element)pmml).getAttribute("version");
      if (version.length() > 0) {
        m_pmmlVersion = version;
      }
    }
  }
  
  /**
   * Set the name of the application (if specified) that created this
   * model
   * 
   * @param doc the Document encapsulating the pmml
   */
  public void setCreatorApplication(Document doc) {
    NodeList tempL = doc.getElementsByTagName("Header");
    Node header = tempL.item(0);
    if (header.getNodeType() == Node.ELEMENT_NODE) {
      NodeList appL = ((Element)header).getElementsByTagName("Application");
      if (appL.getLength() > 0) {
        Node app = appL.item(0);
        if (app.getNodeType() == Node.ELEMENT_NODE) {
          String appName = ((Element)app).getAttribute("name");
          if (appName != null && appName.length() > 0) {
            String version = ((Element)app).getAttribute("version");
            if (version != null && version.length() > 0) {
              appName += " v. " + version;
            }
            m_creatorApplication = appName;
          }
        }
      }
    }
  }

  /**
   * Get the data dictionary.
   *
   * @return the data dictionary
   */
  public Instances getDataDictionary() {
    return m_dataDictionary;
  }

  /**
   * Get the mining schema for this model.
   *
   * @return the mining schema
   */
  public MiningSchema getMiningSchema() {
    return m_miningSchema;
  }

  /**
   * Get the PMML version used for this model.
   *
   * @return the PMML version
   */
  public String getPMMLVersion() {
    return m_pmmlVersion;
  }
  
  /**
   * Get the name of the application that created this model
   * 
   * @return the name of the creating application or null
   * if not specified in the pmml.
   */
  public String getCreatorApplication() {
    return m_creatorApplication;
  }
  
  /**
   * Set a logger to use.
   * 
   * @param log the logger to use
   */
  public void setLog(Logger log) {
    m_log = log;
  }
  
  /**
   * Get the logger.
   * 
   * @return the logger (or null if none is being used)
   */
  public Logger getLog() {
    return m_log;
  }

  /**
   * Throw an exception - PMML models are pre-built.
   *
   * @param data the Instances to learn from
   * @throws Exception if something goes wrong
   */
  public void buildClassifier(Instances data) throws Exception {
    throw new Exception("[PMMLClassifier] PMML models are pre-built "
                        + "and static!");
  }
  
  /**
   * Signal that a scoring run has been completed. Resets
   * the initialized state to false so that a subsequent
   * scoring run will trigger the mapping of the mining
   * schema to incoming instances. If not called after a
   * scoring run, then the classifier will assume that
   * the current mapping is still valid.
   */
  public void done() {
    m_initialized = false;
    m_fieldsMap = null;
  }

  /**
   * Map mining schema to incoming instances.
   *
   * @param dataSet the structure of the incoming Instances
   * @throws Exception if something goes wrong
   */
  public void mapToMiningSchema(Instances dataSet) throws Exception {
    if (m_fieldsMap == null) {
      // PMMLUtils.mapToMiningSchema(dataSet, m_miningSchema);
      m_fieldsMap = new MappingInfo(dataSet, m_miningSchema, m_log);
      m_initialized = true;
    }
  }
  
  /**
   * Get a textual description of the mapping between mining schema
   * fields and incoming data fields.
   * 
   * @return a description of the fields mapping as a String or null if
   * no mapping has been constructed yet.
   */
  public String getFieldsMappingString() {
    if (!m_initialized) {
      return null;
    }
    return m_fieldsMap.getFieldsMappingString();
  }
}
