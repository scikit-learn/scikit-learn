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
 * XMLBeans.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.beans.xml;

import weka.core.converters.ConverterUtils;
import weka.core.xml.XMLBasicSerialization;
import weka.core.xml.XMLDocument;
import weka.core.Environment;
import weka.core.EnvironmentHandler;
import weka.gui.beans.BeanConnection;
import weka.gui.beans.BeanInstance;
import weka.gui.beans.BeanVisual;
import weka.gui.beans.MetaBean;
import weka.gui.beans.Visible;
import weka.gui.beans.BeanCommon;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Point;
import java.beans.BeanInfo;
import java.beans.EventSetDescriptor;
import java.beans.Introspector;
import java.beans.beancontext.BeanContextSupport;
import java.io.File;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.plaf.ColorUIResource;
import javax.swing.plaf.FontUIResource;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * This class serializes and deserializes a KnowledgeFlow setup to and fro XML.
 * <br>
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6392 $
 */
public class XMLBeans 
  extends XMLBasicSerialization {

  /** the value of the id property */
  public final static String VAL_ID = "id";

  /** the value of the x property */
  public final static String VAL_X = "x";

  /** the value of the y property */
  public final static String VAL_Y = "y";

  /** the value of the bean property */
  public final static String VAL_BEAN = "bean";

  /** the value of the customName property */
  public final static String VAL_CUSTOM_NAME = "custom_name";
 
  /** the value of the source property */
  public final static String VAL_SOURCEID = "source_id";
 
  /** the value of the target property */
  public final static String VAL_TARGETID = "target_id";
  
  /** the value of the eventname property */
  public final static String VAL_EVENTNAME = "eventname";
  
  /** the value of the hidden property */
  public final static String VAL_HIDDEN = "hidden";
  
  /** the value of the file property */
  public final static String VAL_FILE = "file";
  
  /** the value of the dir property */
  public final static String VAL_DIR = "dir";
  
  /** the value of the prefix property */
  public final static String VAL_PREFIX = "prefix";

  public final static String VAL_RELATIVE_PATH = "useRelativePath";
  
  /** the value of the options property */
  public final static String VAL_OPTIONS = "options";
  
  /** the value of the saver property */
  public final static String VAL_SAVER = "wrappedAlgorithm";
  
  /** the value of the loader property */
  public final static String VAL_LOADER = "wrappedAlgorithm";
  
  /** the value of the text property */
  public final static String VAL_TEXT = "text";
  
  /** the value of the beanContext property */
  public final static String VAL_BEANCONTEXT = "beanContext";

  /** the value of the width property */
  public final static String VAL_WIDTH = "width";

  /** the value of the height property */
  public final static String VAL_HEIGHT = "height";

  /** the value of the red property */
  public final static String VAL_RED = "red";

  /** the value of the green property */
  public final static String VAL_GREEN = "green";
  
  /** the value of the blue property */
  public final static String VAL_BLUE = "blue";

  /** the value of the value property */
  public final static String VAL_NAME = "name";

  /** the value of the style property */
  public final static String VAL_STYLE = "style";
 
  /** the value of the location property */
  public final static String VAL_LOCATION = "location";
 
  /** the value of the size property */
  public final static String VAL_SIZE = "size";
 
  /** the value of the color property */
  public final static String VAL_COLOR = "color";
  
  /** the value of the font property */
  public final static String VAL_FONT = "font";
  
  /** the value of the iconpath property */
  public final static String VAL_ICONPATH = "iconPath";
  
  /** the value of the animatedIconPath property */
  public final static String VAL_ANIMATEDICONPATH = "animatedIconPath";
  
  /** the value of the associatedConnections property */
  public final static String VAL_ASSOCIATEDCONNECTIONS = "associatedConnections";
  
  /** the value of the input property */
  public final static String VAL_INPUTS = "inputs";
  
  /** the value of the input id property */
  public final static String VAL_INPUTSID = "inputs_id";
  
  /** the value of the outputs id property */
  public final static String VAL_OUTPUTS = "outputs";
  
  /** the value of the outputs property */
  public final static String VAL_OUTPUTSID = "outputs_id";
  
  /** the value of the subFlow property */
  public final static String VAL_SUBFLOW = "subFlow";
  
  /** the value of the originalCoords property */
  public final static String VAL_ORIGINALCOORDS = "originalCoords";
  
  /** the value of the relationNameForFilename property (Saver) */
  public final static String VAL_RELATIONNAMEFORFILENAME = "relationNameForFilename";

  /** the index in the Vector, where the BeanInstances are stored 
   * (Instances and Connections are stored in a Vector and then serialized) */
  public final static int INDEX_BEANINSTANCES = 0;

  /** the index in the Vector, where the BeanConnections are stored
  * (Instances and Connections are stored in a Vector and then serialized) */
  public final static int INDEX_BEANCONNECTIONS = 1;
  
  /** the component that manages the layout of the beans */
  protected JComponent m_BeanLayout;
  
  /** keeps track of the BeanInstances read so far, used for the BeanConnections */
  protected Vector m_BeanInstances;
  
  /** keeps track of the BeanInstances read so far, used for the BeanConnections */
  protected Vector m_BeanInstancesID;

  /** whether to ignore the BeanConnection */
  protected boolean m_IgnoreBeanConnections;
  
  /** the current MetaBean (for the BeanConnections) */
  protected MetaBean m_CurrentMetaBean;

  /** the identifier for regular BeanConnections */
  protected final static String REGULAR_CONNECTION = "regular_connection";
  
  /** the relation between Bean and connection, MetaBean BeanConnections
   * are stored under the reference of the MetaBean, regular connections 
   * are stored under REGULAR_CONNECTION. The relation has the following 
   * format (is a string): sourcePos,targetPos,event,hidden
   * @see #REGULAR_CONNECTION */
  protected Hashtable m_BeanConnectionRelation;
  
  /** the data that is about to be read/written contains a complete layout 
   * @see #m_DataType */
  public final static int DATATYPE_LAYOUT = 0;
  
  /** the data that is about to be read/written contains user-components, i.e., 
   * Metabeans 
   * @see #m_DataType */
  public final static int DATATYPE_USERCOMPONENTS = 1;
  
  /** the type of data that is be read/written
   * @see #DATATYPE_LAYOUT
   * @see #DATATYPE_USERCOMPONENTS */
  protected int m_DataType = DATATYPE_LAYOUT;
  
  /** the beancontext to use for loading from XML and the beancontext is
   * null in the bean */
  protected BeanContextSupport m_BeanContextSupport = null;
  
  /**
   * initializes the serialization for layouts
   * 
   * @param layout      the component that manages the layout
   * @param context     the bean context support to use
   * @throws Exception  if initialization fails
   */
  public XMLBeans(JComponent layout, BeanContextSupport context) throws Exception {
    this(layout, context, DATATYPE_LAYOUT);
  }
  
  /**
   * initializes the serialization for different types of data
   * 
   * @param layout      the component that manages the layout
   * @param context     the bean context support to use
   * @param datatype    the type of data to read/write
   * @throws Exception  if initialization fails
   */
  public XMLBeans(JComponent layout, BeanContextSupport context, int datatype) throws Exception {
    super();
    
    m_BeanLayout = layout;
    m_BeanContextSupport = context;
    setDataType(datatype);
  }
  
  /**
   * sets what kind of data is to be read/written
   * @param value       the type of data
   * @see #m_DataType
   */
  public void setDataType(int value) {
    if (value == DATATYPE_LAYOUT)
      m_DataType = value;
    else if (value == DATATYPE_USERCOMPONENTS)
      m_DataType = value;
    else
      System.out.println("DataType '" + value + "' is unknown!");
  }
  
  /**
   * returns the type of data that is to be read/written
   * @return the type of data
   * @see #m_DataType
   */
  public int getDataType() {
    return m_DataType;
  }

  /**
   * generates internally a new XML document and clears also the IgnoreList and
   * the mappings for the Read/Write-Methods
   * 
   * @throws Exception if something goes wrong
   */
  public void clear() throws Exception {
    Vector<String>	classnames;
    int			i;
    
    super.clear();
    
    // ignore: suppress unnecessary GUI stuff 
    // needs to be checked for new Java versions (might introduce new properties) 
    // - works with Java 1.5
    m_Properties.addIgnored("UI");
    m_Properties.addIgnored("actionMap");
    m_Properties.addIgnored("alignmentX");
    m_Properties.addIgnored("alignmentY");
    m_Properties.addIgnored("autoscrolls");
    m_Properties.addIgnored("background");
    m_Properties.addIgnored("border");
    m_Properties.addIgnored("componentPopupMenu");
    m_Properties.addIgnored("debugGraphicsOptions");
    m_Properties.addIgnored("doubleBuffered");
    m_Properties.addIgnored("enabled");
    m_Properties.addIgnored("focusCycleRoot");
    m_Properties.addIgnored("focusTraversalPolicy");
    m_Properties.addIgnored("focusTraversalPolicyProvider");
    m_Properties.addIgnored("focusable");
    m_Properties.addIgnored("font");
    m_Properties.addIgnored("foreground");
    m_Properties.addIgnored("inheritsPopupMenu");
    m_Properties.addIgnored("inputVerifier");
    m_Properties.addIgnored("layout");
    m_Properties.addIgnored("locale");
    m_Properties.addIgnored("maximumSize");
    m_Properties.addIgnored("minimumSize");
    m_Properties.addIgnored("nextFocusableComponent");
    m_Properties.addIgnored("opaque");
    m_Properties.addIgnored("preferredSize");
    m_Properties.addIgnored("requestFocusEnabled");
    m_Properties.addIgnored("toolTipText");
    m_Properties.addIgnored("transferHandler");
    m_Properties.addIgnored("verifyInputWhenFocusTarget");
    m_Properties.addIgnored("visible");

    // special ignore
    m_Properties.addIgnored("size");  // otherwise you get an endless loop with Dimension!
    m_Properties.addIgnored("location");  // otherwise you get an endless loop with Point!

    // allow
    m_Properties.addAllowed(weka.gui.beans.BeanInstance.class, "x");
    m_Properties.addAllowed(weka.gui.beans.BeanInstance.class, "y");
    m_Properties.addAllowed(weka.gui.beans.BeanInstance.class, "bean");
    m_Properties.addAllowed(weka.gui.beans.Saver.class, "wrappedAlgorithm");
    m_Properties.addAllowed(weka.gui.beans.Loader.class, "wrappedAlgorithm");
    m_Properties.addAllowed(weka.gui.beans.Saver.class, "relationNameForFilename");
    if (getDataType() == DATATYPE_LAYOUT)
      m_Properties.addAllowed(weka.gui.beans.Loader.class, "beanContext");
    else
      m_Properties.addIgnored(weka.gui.beans.Loader.class, "beanContext");   // TODO: more classes???
    m_Properties.addAllowed(weka.gui.beans.Filter.class, "filter");
    m_Properties.addAllowed(weka.gui.beans.Classifier.class, "wrappedAlgorithm");
    m_Properties.addAllowed(weka.gui.beans.Clusterer.class, "wrappedAlgorithm");
    m_Properties.addAllowed(weka.gui.beans.Classifier.class, "executionSlots");
    m_Properties.addAllowed(weka.gui.beans.Classifier.class, "blockOnLastFold");

    m_Properties.addAllowed(weka.classifiers.Classifier.class, "debug");
    m_Properties.addAllowed(weka.classifiers.Classifier.class, "options");
    m_Properties.addAllowed(weka.filters.Filter.class, "options");
    
    m_Properties.addAllowed(weka.core.converters.DatabaseSaver.class, "options");
    m_Properties.addAllowed(weka.core.converters.DatabaseLoader.class, "options");
    m_Properties.addAllowed(weka.core.converters.TextDirectoryLoader.class, "options");

    // we assume that classes implementing SplitEvaluator also implement OptionHandler
    m_Properties.addAllowed(weka.experiment.SplitEvaluator.class, "options");
    // we assume that classes implementing ResultProducer also implement OptionHandler
    m_Properties.addAllowed(weka.experiment.ResultProducer.class, "options");

    // read/write methods
    m_CustomMethods.register(this, Color.class, "Color");
    m_CustomMethods.register(this, Dimension.class, "Dimension");
    m_CustomMethods.register(this, Font.class, "Font");
    m_CustomMethods.register(this, Point.class, "Point");
    m_CustomMethods.register(this, ColorUIResource.class, "ColorUIResource");
    m_CustomMethods.register(this, FontUIResource.class, "FontUIResource");

    m_CustomMethods.register(this, weka.gui.beans.BeanInstance.class, "BeanInstance");
    m_CustomMethods.register(this, weka.gui.beans.BeanConnection.class, "BeanConnection");
    m_CustomMethods.register(this, weka.gui.beans.BeanVisual.class, "BeanVisual");
    m_CustomMethods.register(this, weka.gui.beans.Saver.class, "BeanSaver");
    m_CustomMethods.register(this, weka.gui.beans.MetaBean.class, "MetaBean");

    classnames = ConverterUtils.getFileLoaders();
    for (i = 0; i < classnames.size(); i++)
      m_CustomMethods.register(this, Class.forName(classnames.get(i)), "Loader");
    classnames = ConverterUtils.getFileSavers();
    for (i = 0; i < classnames.size(); i++)
      m_CustomMethods.register(this, Class.forName(classnames.get(i)), "Saver");
    
    // other variables
    m_BeanInstances          = null;
    m_BeanInstancesID        = null;
    m_CurrentMetaBean        = null;
    m_IgnoreBeanConnections  = true;
    m_BeanConnectionRelation = null;
  }
  
  /**
   * traverses over all BeanInstances (or MetaBeans) and stores them in a vector 
   * (recurses into MetaBeans, since the sub-BeanInstances are not visible)
   * @param list       the BeanInstances/MetaBeans to traverse
   */
  protected void addBeanInstances(Vector list) {
    int             i;
    BeanInstance    beaninst;
    
    for (i = 0; i < list.size(); i++) {
      if (list.get(i) instanceof BeanInstance) {
        beaninst = (BeanInstance) list.get(i);
        
        m_BeanInstancesID.add(new Integer(m_BeanInstances.size()));
        m_BeanInstances.add(beaninst);
        
        if (beaninst.getBean() instanceof MetaBean)
          addBeanInstances(((MetaBean) beaninst.getBean()).getBeansInSubFlow());
      }
      else if (list.get(i) instanceof MetaBean) {
        addBeanInstances(((MetaBean) list.get(i)).getBeansInSubFlow());
      }
      else {
        System.out.println("addBeanInstances does not support Vectors of class '" + list.get(i) + "'!");
      }
    }
  }
  
  /**
   * enables derived classes to due some pre-processing on the objects, that's
   * about to be serialized. Right now it only returns the object.
   * 
   * @param o the object that is serialized into XML
   * @return the possibly altered object
   * @throws Exception if post-processing fails
   */
  protected Object writePreProcess(Object o) throws Exception {
    o = super.writePreProcess(o);
    
    // gather all BeanInstances, also the ones in MetaBeans
    m_BeanInstances   = new Vector();
    m_BeanInstancesID = new Vector();
    
    switch (getDataType()) {
      case DATATYPE_LAYOUT:
        addBeanInstances(BeanInstance.getBeanInstances());
        break;

      case DATATYPE_USERCOMPONENTS:
        addBeanInstances((Vector) o);
        break;
        
      default:
        System.out.println("writePreProcess: data type '" + getDataType() + "' is not recognized!");
        break;
    }
    
    return o;
  }
  
  /**
   * enables derived classes to add other properties to the DOM tree, e.g.
   * ones that do not apply to the get/set convention of beans. only implemented
   * with empty method body.
   * 
   * @param o the object that is serialized into XML
   * @throws Exception if post-processing fails
   */
  protected void writePostProcess(Object o) throws Exception {
    Element         root;
    NodeList        list;
    Element         conns;
    Element         child;
    int             i;

    // since not all BeanConnections get saved to XML (e.g., MetaBeans in the 
    // UserToolBar) if one saves a layout, the numbering in the Vector of the 
    // BeanConnections is not correct. The "name" attribute of the nodes has
    // to be modified
    if (getDataType() == DATATYPE_LAYOUT) {
      root  = m_Document.getDocument().getDocumentElement();
      conns = (Element) root.getChildNodes().item(INDEX_BEANCONNECTIONS);
      list  = conns.getChildNodes();
      for (i = 0; i < list.getLength(); i++) {
        child = (Element) list.item(i);
        child.setAttribute(ATT_NAME, "" + i);
      }
    }
  }
  
  /**
   * additional pre-processing can happen in derived classes before the 
   * actual reading from XML (working on the raw XML). right now it does 
   * nothing with the document, only empties the help-vector for the
   * BeanInstances and reads the IDs for the BeanInstances, s.t. the correct
   * references can be set again
   * 
   * @param document the document to pre-process
   * @return the processed object
   * @throws Exception if post-processing fails
   * @see #m_BeanInstances
   */
  protected Document readPreProcess(Document document) throws Exception {
    NodeList        list;
    int             i;
    Element         node;
    String          clsName;
    Vector          children;
    int             id;
    int             n;
    Element         child;
    
    m_BeanInstances   = new Vector();
    m_BeanInstancesID = new Vector();
    
    // get all BeanInstance nodes
    list    = document.getElementsByTagName("*");
    clsName = BeanInstance.class.getName();
    for (i = 0; i < list.getLength(); i++) {
      node = (Element) list.item(i);
      
      // is it a BeanInstance?
      if (node.getAttribute(ATT_CLASS).equals(clsName)) {
        children = XMLDocument.getChildTags(node);
        id       = m_BeanInstancesID.size();
        
        // get id-tag (if available)
        for (n = 0; n < children.size(); n++) {
          child = (Element) children.get(n);
          if (child.getAttribute(ATT_NAME).equals(VAL_ID))
            id = readIntFromXML((Element) child);
        }
        
        m_BeanInstancesID.add(new Integer(id));
      }
    }
    
    m_BeanInstances.setSize(m_BeanInstancesID.size());
    
    // set MetaBean to null
    m_CurrentMetaBean = null;

    // no BeanConnections -> see readPostProcess(Object)
    m_IgnoreBeanConnections = true;
    
    // reset BeanConnection-Relations
    m_BeanConnectionRelation = new Hashtable();
    
    return document;
  }
  
  /**
   * puts the given BeanConnection onto the next null in the given Vector,
   * or at the end of the list, if no null is found.
   * (during the de-serializing, no BeanConnections are set, only nulls)
   * @param conn      the connection to add to the list
   * @param list      the list to add the BeanConnection to
   */
  protected void setBeanConnection(BeanConnection conn, Vector list) {
    int         i;
    boolean     added;
    
    added = false;
    for (i = 0; i < list.size(); i++) {
      if (list.get(i) == null) {
        list.set(i, conn);
        added = true;
        break;
      }
    }
    
    if (!added)
      list.add(conn);
  }
  
  /**
   * generates a connection based on the given parameters
   * @param sourcePos the source position in the m_BeanInstances vector
   * @param targetPos the target position in the m_BeanInstances vector
   * @param event the name of the event, i.e., the connection
   * @param hidden true if the connection is hidden
   * @return the generated BeanConnection
   * @throws Exception if something goes wrong
   */
  protected BeanConnection createBeanConnection(int sourcePos, int targetPos, String event, boolean hidden) throws Exception {
    BeanConnection          result;
    BeanInfo                compInfo;
    EventSetDescriptor[]    esds;
    int                     i;
    BeanInstance            instSource;
    BeanInstance            instTarget;

    result = null;
    
    // was there a connection?
    if ( (sourcePos == -1) || (targetPos == -1) )
      return result;
    
    instSource = (BeanInstance) m_BeanInstances.get(sourcePos);
    instTarget = (BeanInstance) m_BeanInstances.get(targetPos);
    
    compInfo = Introspector.getBeanInfo(((BeanInstance) m_BeanInstances.get(sourcePos)).getBean().getClass());
    esds     = compInfo.getEventSetDescriptors();

    for (i = 0; i < esds.length; i++) {
      if (esds[i].getName().equals(event)) {
        result = new BeanConnection(instSource, instTarget, esds[i]);
        ((BeanConnection) result).setHidden(hidden);
        break;
      }
    }
    
    return result;
  }
  
  /**
   * rebuilds all the connections for a certain key in the hashtable.
   * for the ones being part of a MetaBean, no new instance is built, but only
   * the reference to the actual BeanConnection set.
   * @param deserialized    the deserialized knowledgeflow
   * @param key             the key of the hashtable to rebuild all connections for
   * @throws Exception if something goes wrong
   */
  protected void rebuildBeanConnections(Vector deserialized, Object key) throws Exception {
    int                     i;
    int                     n;
    int                     sourcePos;
    int                     targetPos;
    String                  event;
    boolean                 hidden;
    Vector                  conns;
    BeanConnection          conn;
    StringTokenizer         tok;
    Vector                  beanconns;

    conns = (Vector) m_BeanConnectionRelation.get(key);
    
    // no connections?
    if (conns == null)
      return;
    
    for (n = 0; n < conns.size(); n++) {
      tok       = new StringTokenizer(conns.get(n).toString(), ",");
      conn      = null;
      sourcePos = Integer.parseInt(tok.nextToken());
      targetPos = Integer.parseInt(tok.nextToken());
      event     = tok.nextToken();
      hidden    = stringToBoolean(tok.nextToken());

      // regular connection? -> new instance
      // or MetaBean from user toolbar
      if ( (!(key instanceof MetaBean)) || (getDataType() == DATATYPE_USERCOMPONENTS)) {
        conn = createBeanConnection(sourcePos, targetPos, event, hidden);
      }
      // MetaBean? -> find BeanConnection 
      else {
        beanconns = BeanConnection.getConnections();
        
        for (i = 0; i < beanconns.size(); i++) {
          conn = (BeanConnection) beanconns.get(i);
          if (    (conn.getSource() == (BeanInstance) m_BeanInstances.get(sourcePos))
               && (conn.getTarget() == (BeanInstance) m_BeanInstances.get(targetPos))
               && (conn.getEventName().equals(event)) ) {
            break;
          }
          conn = null;
        }
      }
      
      // add the connection to the corresponding list/MetaBean
      if (key instanceof MetaBean)
        setBeanConnection(conn, ((MetaBean) key).getAssociatedConnections());
      else
        setBeanConnection(conn, (Vector) deserialized.get(INDEX_BEANCONNECTIONS));
    }
  }
  
  /**
   * removes the given meta beans from the layout, since they're only listed
   * in the user toolbar
   * 
   * @param metabeans         the list of MetaBeans in the user toolbar
   */
  protected void removeUserToolBarBeans(Vector metabeans) {
    int           i;
    int           n;
    MetaBean      meta;
    Vector        subflow;
    BeanInstance  beaninst;
    
    for (i = 0; i < metabeans.size(); i++) {
      meta    = (MetaBean) metabeans.get(i);
      subflow = meta.getSubFlow();
      
      for (n = 0; n < subflow.size(); n++) {
        beaninst = (BeanInstance) subflow.get(n);
        beaninst.removeBean(m_BeanLayout);
      }
    }
  }
  
  /**
   * additional post-processing can happen in derived classes after reading 
   * from XML. re-builds the BeanConnections.
   * 
   * @param o the object to perform some additional processing on
   * @return the processed object
   * @throws Exception if post-processing fails
   */
  protected Object readPostProcess(Object o) throws Exception {
    Enumeration             enm;
    Vector                  deserialized;
    Object                  key;

    deserialized = (Vector) super.readPostProcess(o);
    
    // rebuild the actual connections
    rebuildBeanConnections(deserialized, REGULAR_CONNECTION);

    // rebuild the references in the MetaBeans
    enm = m_BeanConnectionRelation.keys();
    while (enm.hasMoreElements()) {
      key = enm.nextElement();
      
      // skip the regular connections
      if (!(key instanceof MetaBean))
        continue;
      
      rebuildBeanConnections(deserialized, key);
    }

    // remove MetaBean and subflow from BeanInstance (not part of the flow!)
    if (getDataType() == DATATYPE_USERCOMPONENTS)
      removeUserToolBarBeans(deserialized);
    
    return deserialized;
  }

  /**
   * returns the relation for the given MetaBean, for the regular connections,
   * null has to be used
   * @param meta      the MetaBean (or null for regular connections) to retrieve
   *                  the connections for
   * @return          the associated connections
   * @see #REGULAR_CONNECTION
   */
  protected Vector getBeanConnectionRelation(MetaBean meta) {
    Vector      result;
    Object      key;
    
    if (meta == null)
      key = REGULAR_CONNECTION;
    else
      key = meta;
    
    // not yet in there?
    if (!m_BeanConnectionRelation.containsKey(key)) {
      m_BeanConnectionRelation.put(key, new Vector());
    }

    result = (Vector) m_BeanConnectionRelation.get(key);
    
    return result;
  }
  
  /**
   * adds the given connection-relation for the specified MetaBean (or null in 
   * case of regular connections)
   * @param meta        the MetaBean (or null for regular connections) to add
   *                    the relationship for
   * @param connection  the connection relation to add
   */
  protected void addBeanConnectionRelation(MetaBean meta, String connection) {
    Vector      relations;
    Object      key;
    
    relations = getBeanConnectionRelation(meta);
    
    // add relation
    relations.add(connection);
    
    // update
    if (meta == null)
      key = REGULAR_CONNECTION;
    else
      key = meta;
    m_BeanConnectionRelation.put(key, relations);
  }
  
  /**
   * adds the given Color to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeColor(Element parent, Object o, String name)
    throws Exception {
    
    Element     node;
    Color       color;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    color = (Color) o;
    node  = addElement(parent, name, color.getClass().getName(), false);

    writeIntToXML(node, color.getRed(), VAL_RED);
    writeIntToXML(node, color.getGreen(), VAL_GREEN);
    writeIntToXML(node, color.getBlue(), VAL_BLUE);
    
    return node;
  }

  /**
   * builds the Color from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readColor(Element node) throws Exception {
    Object      result;
    Vector      children;
    Element     child;
    int         i;
    int         red;
    int         green;
    int         blue;
    String      name;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = null;
    children = XMLDocument.getChildTags(node);
    red      = 0;
    green    = 0;
    blue     = 0;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_RED))
        red = readIntFromXML(child);
      else if (name.equals(VAL_GREEN))
        green = readIntFromXML(child);
      else if (name.equals(VAL_BLUE))
        blue = readIntFromXML(child);
      else
        System.out.println("WARNING: '" + name
            + "' is not a recognized name for " + node.getAttribute(ATT_NAME) + "!");
    }
    
    result = new Color(red, green, blue);

    return result;
  }

  /**
   * adds the given Dimension to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeDimension(Element parent, Object o, String name)
      throws Exception {
    
    Element     node;
    Dimension   dim;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    dim = (Dimension) o;
    node = addElement(parent, name, dim.getClass().getName(), false);

    writeDoubleToXML(node, dim.getWidth(), VAL_WIDTH);
    writeDoubleToXML(node, dim.getHeight(), VAL_HEIGHT);
    
    return node;
  }

  /**
   * builds the Dimension from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readDimension(Element node) throws Exception {
    Object      result;
    Vector      children;
    Element     child;
    int         i;
    double      width;
    double      height;
    String      name;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = null;
    children = XMLDocument.getChildTags(node);
    width    = 0;
    height   = 0;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_WIDTH))
        width = readDoubleFromXML(child);
      else if (name.equals(VAL_HEIGHT))
        height = readDoubleFromXML(child);
      else
        System.out.println("WARNING: '" + name
            + "' is not a recognized name for " + node.getAttribute(ATT_NAME) + "!");
    }
    
    result = new Dimension();
    ((Dimension) result).setSize(width, height);

    return result;
  }

  /**
   * adds the given Font to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeFont(Element parent, Object o, String name)
      throws Exception {
    
    Element     node;
    Font        font;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    font = (Font) o;
    node = addElement(parent, name, font.getClass().getName(), false);

    invokeWriteToXML(node, font.getName(), VAL_NAME);
    writeIntToXML(node, font.getStyle(), VAL_STYLE);
    writeIntToXML(node, font.getSize(), VAL_SIZE);
    
    return node;
  }

  /**
   * builds the Font from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readFont(Element node) throws Exception {
    Object      result;
    Vector      children;
    Element     child;
    int         i;
    int         style;
    int         size;
    String      name;
    String      fontname;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = null;
    children = XMLDocument.getChildTags(node);
    fontname = "";
    style    = 0;
    size     = 0;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_NAME))
        name = (String) invokeReadFromXML(child);
      else if (name.equals(VAL_STYLE))
        style = readIntFromXML(child);
      else if (name.equals(VAL_SIZE))
        size = readIntFromXML(child);
      else
        System.out.println("WARNING: '" + name
            + "' is not a recognized name for " + node.getAttribute(ATT_NAME) + "!");
    }
    
    result = new Font(fontname, style, size);

    return result;
  }

  /**
   * adds the given Point to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writePoint(Element parent, Object o, String name)
      throws Exception {
    
    Element     node;
    Point       p;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    p    = (Point) o;
    node = addElement(parent, name, p.getClass().getName(), false);

    writeDoubleToXML(node, p.getX(), VAL_X);
    writeDoubleToXML(node, p.getY(), VAL_Y);
    
    return node;
  }

  /**
   * builds the Point from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readPoint(Element node) throws Exception {
    Object      result;
    Vector      children;
    Element     child;
    int         i;
    double      x;
    double      y;
    String      name;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = null;
    children = XMLDocument.getChildTags(node);
    x        = 0;
    y        = 0;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_X))
        x = readDoubleFromXML(child);
      else if (name.equals(VAL_Y))
        y = readDoubleFromXML(child);
      else
        System.out.println("WARNING: '" + name
            + "' is not a recognized name for " + node.getAttribute(ATT_NAME) + "!");
    }
    
    result = new Point();
    ((Point) result).setLocation(x, y);

    return result;
  }

  /**
   * adds the given ColorUIResource to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeColorUIResource(Element parent, Object o, String name)
      throws Exception {
    
    Element           node;
    ColorUIResource   resource;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    resource = (ColorUIResource) o;
    node     = addElement(parent, name, resource.getClass().getName(), false);
    invokeWriteToXML(node, new Color(resource.getRGB()), VAL_COLOR);
    
    return node;
  }

  /**
   * builds the ColorUIResource from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readColorUIResource(Element node) throws Exception {
    Object      result;
    Vector      children;
    Element     child;
    int         i;
    String      name;
    Color       color;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = null;
    children = XMLDocument.getChildTags(node);
    color    = null;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_COLOR))
        color = (Color) invokeReadFromXML(child);
      else
        System.out.println("WARNING: '" + name
            + "' is not a recognized name for " + node.getAttribute(ATT_NAME) + "!");
    }
    
    result = new ColorUIResource(color);

    return result;
  }

  /**
   * adds the given FontUIResource to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeFontUIResource(Element parent, Object o, String name)
      throws Exception {
    
    Element           node;
    FontUIResource    resource;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    resource = (FontUIResource) o;
    node     = addElement(parent, name, resource.getClass().getName(), false);
    invokeWriteToXML(node, new Font(resource.getName(), resource.getStyle(), resource.getSize()), VAL_COLOR);
    
    return node;
  }

  /**
   * builds the FontUIResource from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readFontUIResource(Element node) throws Exception {
    Object      result;
    Vector      children;
    Element     child;
    int         i;
    String      name;
    Font        font;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = null;
    children = XMLDocument.getChildTags(node);
    font     = null;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_FONT))
        font = (Font) invokeReadFromXML(child);
      else
        System.out.println("WARNING: '" + name
            + "' is not a recognized name for " + node.getAttribute(ATT_NAME) + "!");
    }
    
    result = new FontUIResource(font);

    return result;
  }

  /**
   * adds the given BeanInstance to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeBeanInstance(Element parent, Object o, String name)
      throws Exception {
    
    Element         node;
    BeanInstance    beaninst;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    beaninst = (BeanInstance) o;
    node     = addElement(parent, name, beaninst.getClass().getName(), false);

    writeIntToXML(node, m_BeanInstances.indexOf(beaninst), VAL_ID);
    writeIntToXML(node, beaninst.getX() + beaninst.getWidth()  / 2, VAL_X);   // x is thought to be in the center?
    writeIntToXML(node, beaninst.getY() + beaninst.getHeight() / 2, VAL_Y);   // y is thought to be in the center?
    if (beaninst.getBean() instanceof BeanCommon) {
      // write the custom name of this bean
      String custName = ((BeanCommon)beaninst.getBean()).getCustomName();
      invokeWriteToXML(node, custName, VAL_CUSTOM_NAME);
    }
    invokeWriteToXML(node, beaninst.getBean(), VAL_BEAN);
    
    return node;
  }

  /**
   * builds the BeanInstance from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readBeanInstance(Element node) throws Exception {
    Object          result;
    Vector          children;
    Element         child;
    String          name;
    int             i;
    int             x;
    int             y;
    int             id;
    Object          bean;
    BeanVisual      visual;
    BeanInstance    beaninst;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = null;
    children = XMLDocument.getChildTags(node);
    id       = -1;
    x        = 0;
    y        = 0;
    bean     = null;
    String customName = null;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_ID)) {
        id = readIntFromXML(child);
      } else if (name.equals(VAL_X)) {
        x = readIntFromXML(child);
      } else if (name.equals(VAL_Y)) {
        y = readIntFromXML(child);
      } else if (name.equals(VAL_CUSTOM_NAME)) {
        customName = (String)invokeReadFromXML(child);
      } else if (name.equals(VAL_BEAN)) {
        bean = invokeReadFromXML(child);
      } else {
        System.out.println("WARNING: '" + name
            + "' is not a recognized name for " + node.getAttribute(ATT_NAME) + "!");
      }
    }
    
    result   = new BeanInstance(m_BeanLayout, bean, x, y);
    beaninst = (BeanInstance) result;
    
    // set parent of BeanVisual
    if (beaninst.getBean() instanceof weka.gui.beans.Visible) {
      visual = ((Visible) beaninst.getBean()).getVisual();
      visual.setSize(visual.getPreferredSize());
      if (visual.getParent() == null) {
        ((JPanel) beaninst.getBean()).add(visual);
      }
    }

    if (beaninst.getBean() instanceof BeanCommon &&
        customName != null) {
      ((BeanCommon)beaninst.getBean()).setCustomName(customName);
    }
    
    // no IDs -> get next null position
    if (id == -1) {
      for (i = 0; i < m_BeanInstances.size(); i++) {
        if (m_BeanInstances.get(i) == null) {
          id = ((Integer) m_BeanInstancesID.get(i)).intValue();
          break;
        }
      }
    }
    // get position for id
    i = m_BeanInstancesID.indexOf(new Integer(id));

    // keep track of the BeanInstances for reading the connections later on
    m_BeanInstances.set(i, result);

    // no current MetaBean
    m_CurrentMetaBean = null;
    
    return result;
  }

  /**
   * adds the given BeanConncetion to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeBeanConnection(Element parent, Object o, String name)
    throws Exception {
    
    Element           node;
    BeanConnection    beanconn;
    int               source;
    int               target;
    int               sourcePos;
    int               targetPos;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    beanconn = (BeanConnection) o;
    node     = null;

    // get position
    sourcePos = m_BeanInstances.indexOf(beanconn.getSource());
    targetPos = m_BeanInstances.indexOf(beanconn.getTarget());
    
    // get id (if Connection is from a Bean in the UserToolBar, it's not listed! -> ignore it)
    if ( (sourcePos > -1) && (targetPos > -1) ) {
      source = ((Integer) m_BeanInstancesID.get(sourcePos)).intValue();
      target = ((Integer) m_BeanInstancesID.get(targetPos)).intValue();
    }
    else {
       source = -1;
       target = -1;
    }
    
    // connection exists in the layout?
    if ( (source > -1) && (target > -1) ) {
      node = addElement(parent, name, beanconn.getClass().getName(), false);
  
      writeIntToXML(node, source, VAL_SOURCEID);
      writeIntToXML(node, target, VAL_TARGETID);
      invokeWriteToXML(node, beanconn.getEventName(), VAL_EVENTNAME);
      writeBooleanToXML(node, beanconn.isHidden(), VAL_HIDDEN);
    }
    
    return node;
  }

  /**
   * builds the BeanConnection from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readBeanConnection(Element node) throws Exception {
    Object                  result;
    Vector                  children;
    Element                 child;
    String                  name;
    int                     i;
    int                     source;
    int                     target;
    int                     sourcePos;
    int                     targetPos;
    String                  event;
    boolean                 hidden;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = null;
    children = XMLDocument.getChildTags(node);
    source   = 0;
    target   = 0;
    event    = "";
    hidden   = false;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_SOURCEID))
        source = readIntFromXML(child);
      else if (name.equals(VAL_TARGETID))
        target = readIntFromXML(child);
      else if (name.equals(VAL_EVENTNAME))
        event = (String) invokeReadFromXML(child);
      else if (name.equals(VAL_HIDDEN))
        hidden = readBooleanFromXML(child);
      else
        System.out.println("WARNING: '" + name
            + "' is not a recognized name for " + node.getAttribute(ATT_NAME) + "!");
    }

    // get position of id
    sourcePos = m_BeanInstancesID.indexOf(new Integer(source));
    targetPos = m_BeanInstancesID.indexOf(new Integer(target));
    
    // do we currently ignore the connections?
    // Note: necessary because of the MetaBeans
    if (m_IgnoreBeanConnections) {
      addBeanConnectionRelation(m_CurrentMetaBean, sourcePos + "," + targetPos + "," + event + "," + hidden);
      return result;
    }

    // generate it normally
    result = createBeanConnection(sourcePos, targetPos, event, hidden);

    return result;
  }
  
  /**
   * adds the given Loader (a bean) to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeBeanLoader(Element parent, Object o, String name)
      throws Exception {
    
    Element                 node;
    weka.gui.beans.Loader   loader;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    loader = (weka.gui.beans.Loader) o;
    node   = addElement(parent, name, loader.getClass().getName(), false);

    invokeWriteToXML(node, loader.getLoader(), VAL_LOADER);
    invokeWriteToXML(node, loader.getBeanContext(), VAL_BEANCONTEXT);
    
    return node;
  }

  /**
   * adds the given Saver (a bean) to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeBeanSaver(Element parent, Object o, String name)
      throws Exception {
    
    Element                 node;
    weka.gui.beans.Saver    saver;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    saver = (weka.gui.beans.Saver) o;
    node   = addElement(parent, name, saver.getClass().getName(), false);
    invokeWriteToXML(node, saver.getRelationNameForFilename(), VAL_RELATIONNAMEFORFILENAME);

    invokeWriteToXML(node, saver.getSaverTemplate(), VAL_SAVER);
    
    return node;
  }
  
  /**
   * adds the given Loader to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeLoader(Element parent, Object o, String name)
      throws Exception {
    
    Element                       node;
    weka.core.converters.Loader   loader;
    File                          file;
    boolean                       known;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    loader = (weka.core.converters.Loader) o;
    node   = addElement(parent, name, loader.getClass().getName(), false);
    known  = true;
    file   = null;

    // file
    if (loader instanceof weka.core.converters.AbstractFileLoader)
      file = ((weka.core.converters.AbstractFileLoader) loader).retrieveFile();
    else
      known = false;

    if (!known)
      System.out.println("WARNING: unknown loader class '" + loader.getClass().getName() + "' - cannot retrieve file!");

    Boolean relativeB = null;
    if (loader instanceof weka.core.converters.FileSourcedConverter) {
      boolean relative = ((weka.core.converters.FileSourcedConverter)loader).getUseRelativePath();
      relativeB = new Boolean(relative);
    }
    
    // only save it, if it's a real file!
    if ( (file == null) || (file.isDirectory()) ) {
      invokeWriteToXML(node, "", VAL_FILE);
    } else {
      boolean notAbsolute = 
        (((weka.core.converters.AbstractFileLoader) loader).getUseRelativePath() ||
        (loader instanceof EnvironmentHandler 
            && Environment.containsEnvVariables(file.getPath())));
      
      String path = (notAbsolute)
        ? file.getPath()
        : file.getAbsolutePath();
      // Replace any windows file separators with forward slashes (Java under windows can
      // read paths with forward slashes (apparantly)
      path = path.replace('\\', '/');
      invokeWriteToXML(node, path, VAL_FILE);
    }
    if (relativeB != null) {
      invokeWriteToXML(node, relativeB.toString(), VAL_RELATIVE_PATH);
    }
    
    return node;
  }

  /**
   * builds the Loader from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readLoader(Element node) throws Exception {
    Object      result;
    Vector      children;
    Element     child;
    int         i;
    String      name;
    String      file;
    File        fl;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = Class.forName(node.getAttribute(ATT_CLASS)).newInstance();
    children = XMLDocument.getChildTags(node);
    file     = "";
    Object relativeB = null;
    boolean relative = false;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_FILE)) {
        file = (String) invokeReadFromXML(child);
      } else if (name.equals(VAL_RELATIVE_PATH)) {
        relativeB = readFromXML(child);
        if (relativeB instanceof Boolean) {
          relative = ((Boolean)relativeB).booleanValue();
        }
      } else {
        readFromXML(result, name, child);
      }
    }

    if (result instanceof weka.core.converters.FileSourcedConverter) {
      ((weka.core.converters.FileSourcedConverter)result).setUseRelativePath(relative);
    }

    if (file.equals(""))
      file = null;

    // set file only, if it exists
    if (file != null) {
      String tempFile = file;

      boolean containsEnv = false;
      containsEnv = Environment.containsEnvVariables(file);
      
      fl = new File(file);      
      // only test for existence if the path does not contain environment vars
      // (trust that after they are resolved that everything is hunky dory)
      if (containsEnv || fl.exists()) {
        ((weka.core.converters.AbstractFileLoader) result).setSource(new File(file));
      } else {
        System.out.println("WARNING: The file '" + tempFile + "' does not exist!");
      }
    }
    
    return result;
  }

  /**
   * adds the given Saver to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeSaver(Element parent, Object o, String name)
      throws Exception {
    
    Element                     node;
    weka.core.converters.Saver  saver;
    File                        file;
    String			prefix;
    String			dir;
    boolean                     known;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    saver  = (weka.core.converters.Saver) o;
    node   = addElement(parent, name, saver.getClass().getName(), false);
    known  = true;
    file   = null;
    prefix = "";
    dir    = "";

    // file
    if (saver instanceof weka.core.converters.AbstractFileSaver) {
      file   = ((weka.core.converters.AbstractFileSaver) saver).retrieveFile();
      prefix = ((weka.core.converters.AbstractFileSaver) saver).filePrefix();
      dir    = ((weka.core.converters.AbstractFileSaver) saver).retrieveDir();
      // Replace any windows file separators with forward slashes (Java under windows can
      // read paths with forward slashes (apparantly)
      dir = dir.replace('\\', '/');
    }
    else {
      known = false;
    }
    
    if (!known)
      System.out.println("WARNING: unknown saver class '" + saver.getClass().getName() + "' - cannot retrieve file!");

    Boolean relativeB = null;
    if (saver instanceof weka.core.converters.FileSourcedConverter) {
      boolean relative = ((weka.core.converters.FileSourcedConverter)saver).getUseRelativePath();
      relativeB = new Boolean(relative);
    }
    

//    if ( (file == null) || (file.isDirectory()) ) {
      invokeWriteToXML(node, "",     VAL_FILE);
      invokeWriteToXML(node, dir,    VAL_DIR);
      invokeWriteToXML(node, prefix, VAL_PREFIX);
/*    }
    else {
      String path = (((weka.core.converters.FileSourcedConverter) saver).getUseRelativePath())
        ? file.getPath()
        : file.getAbsolutePath();
      // Replace any windows file separators with forward slashes (Java under windows can
      // read paths with forward slashes (apparantly)
      path = path.replace('\\', '/');
      invokeWriteToXML(node, path, VAL_FILE);
      invokeWriteToXML(node, dir, VAL_DIR);
      invokeWriteToXML(node, prefix, VAL_PREFIX);
    }*/

    if (relativeB != null) {
      invokeWriteToXML(node, relativeB.toString(), VAL_RELATIVE_PATH);
    }
    
    return node;
  }

  /**
   * builds the Saver from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readSaver(Element node) throws Exception {
    Object      result;
    Vector      children;
    Element     child;
    int         i;
    String      name;
    String      file;
    String	dir;
    String	prefix;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = Class.forName(node.getAttribute(ATT_CLASS)).newInstance();
    children = XMLDocument.getChildTags(node);
    file     = null;
    dir      = null;
    prefix   = null;

    Object relativeB = null;
    boolean relative = false;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_FILE)) {
        file = (String) invokeReadFromXML(child);
      } else if (name.equals(VAL_DIR)) {
        dir = (String) invokeReadFromXML(child);
      } else if (name.equals(VAL_PREFIX)) {
        prefix = (String) invokeReadFromXML(child);
      } else if (name.equals(VAL_RELATIVE_PATH)) {
        relativeB = readFromXML(child);
        if (relativeB instanceof Boolean) {
          relative = ((Boolean)relativeB).booleanValue();
        }
      } else {
        readFromXML(result, name, child);
      }
    }

    if ( (file != null) && (file.length() == 0) )
      file = null;

    // savers only get directory and prefix, not file (KnowledgeFlow sets the
    // file/destination based on the relation, dir and prefix)
    if ( (dir != null) && (prefix != null) ) {
      ((weka.core.converters.AbstractFileSaver) result).setDir(dir);
      ((weka.core.converters.AbstractFileSaver) result).setFilePrefix(prefix);
    }

    if (result instanceof weka.core.converters.FileSourcedConverter) {
      ((weka.core.converters.FileSourcedConverter)result).setUseRelativePath(relative);
    }
    
    return result;
  }

  /**
   * adds the given BeanVisual to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeBeanVisual(Element parent, Object o, String name)
      throws Exception {
    
    Element         node;
    BeanVisual      visual;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    visual = (BeanVisual) o;
    node   = writeToXML(parent, o, name);

    // add icon paths
    invokeWriteToXML(node, visual.getIconPath(),         VAL_ICONPATH);
    invokeWriteToXML(node, visual.getAnimatedIconPath(), VAL_ANIMATEDICONPATH);
    
    return node;
  }

  /**
   * builds the BeanVisual from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readBeanVisual(Element node) throws Exception {
    Object      result;
    Vector      children;
    Element     child;
    int         i;
    String      name;
    String      text;
    String      iconPath;
    String      animIconPath;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result       = null;
    children     = XMLDocument.getChildTags(node);
    text         = "";
    iconPath     = "";
    animIconPath = "";

    // find text
    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_TEXT))
        text = (String) invokeReadFromXML(child);
      else if (name.equals(VAL_ICONPATH))
        iconPath = (String) invokeReadFromXML(child);
      else if (name.equals(VAL_ANIMATEDICONPATH))
        animIconPath = (String) invokeReadFromXML(child);
    }

    result = new BeanVisual(text, iconPath, animIconPath);
    
    // set rest of properties
    for (i = 0; i < children.size(); i++)
      readFromXML(result, node.getAttribute(ATT_NAME), (Element) children.get(i));
    
    return result;
  }
  
  /**
   * returns the IDs for the given BeanInstances, i.e., the stored IDs
   * in m_BeanInstancesID, based on m_BeanInstances
   * 
   * @param beans       the beans to retrieve the IDs for
   * @return            the IDs for the given BeanInstances
   * @see #m_BeanInstances
   * @see #m_BeanInstancesID
   */
  protected Vector getIDsForBeanInstances(Vector beans) {
    Vector        result;
    int           i;
    int           pos;
    
    result = new Vector();

    for (i = 0; i < beans.size(); i++) {
      pos = m_BeanInstances.indexOf(beans.get(i));
      result.add(m_BeanInstancesID.get(pos));
    }
    
    return result;
  }
  
  /**
   * adds the given MetaBean to a DOM structure.
   * 
   * @param parent the parent of this object, e.g. the class this object is a member of
   * @param o the Object to describe in XML
   * @param name the name of the object
   * @return the node that was created
   * @throws Exception if the DOM creation fails
   */
  public Element writeMetaBean(Element parent, Object o, String name)
      throws Exception {
    
    Element         node;
    MetaBean        meta;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), name);
    
    m_CurrentNode = parent;
    
    meta = (MetaBean) o;
    node = writeToXML(parent, o, name);

    invokeWriteToXML(node, getIDsForBeanInstances(meta.getBeansInInputs()), VAL_INPUTSID);
    invokeWriteToXML(node, getIDsForBeanInstances(meta.getBeansInOutputs()), VAL_OUTPUTSID);
    
    return node;
  }
  
  /**
   * returns a vector with references to BeanInstances according to the IDs
   * in the given Vector.
   * @param ids       contains the IDs of the BeanInstances
   * @return          the corresponding BeanInstances
   * @see #m_BeanInstances
   * @see #m_BeanInstancesID
   */
  protected Vector getBeanInstancesForIDs(Vector ids) {
    Vector        result;
    int           i;
    int           pos;
    
    result = new Vector();
    
    for (i = 0; i < ids.size(); i++) {
      pos = m_BeanInstancesID.indexOf(ids.get(i));
      result.add(m_BeanInstances.get(pos));
    }
    
    return result;
  }

  /**
   * builds the MetaBean from the given DOM node.
   * 
   * @param node the associated XML node
   * @return the instance created from the XML description
   * @throws Exception if instantiation fails
   */
  public Object readMetaBean(Element node) throws Exception {
    Object          result;
    Vector          children;
    Element         child;
    int             i;
    String          name;
    Vector          inputs;
    Vector          outputs;
    Vector          coords;
    MetaBean        bean;

    // for debugging only
    if (DEBUG)
       trace(new Throwable(), node.getAttribute(ATT_NAME));

    m_CurrentNode = node;
    
    result   = new MetaBean();
    children = XMLDocument.getChildTags(node);
    inputs   = new Vector();
    outputs  = new Vector();
    coords   = new Vector();
    
    // the current MetaBean
    m_CurrentMetaBean = (MetaBean) result;

    for (i = 0; i < children.size(); i++) {
      child = (Element) children.get(i);
      name  = child.getAttribute(ATT_NAME);

      if (name.equals(VAL_ASSOCIATEDCONNECTIONS))
        ((MetaBean) result).setAssociatedConnections((Vector) invokeReadFromXML(child));
      else if (name.equals(VAL_INPUTSID))
        inputs = (Vector) invokeReadFromXML(child);
      else if (name.equals(VAL_OUTPUTSID))
        outputs = (Vector) invokeReadFromXML(child);
      else if (name.equals(VAL_SUBFLOW))
        ((MetaBean) result).setSubFlow((Vector) invokeReadFromXML(child));
      else if (name.equals(VAL_ORIGINALCOORDS))
        coords = (Vector) invokeReadFromXML(child);
      else if (name.equals(VAL_INPUTS))
        System.out.println("INFO: '" + name + "' will be restored later.");
      else if (name.equals(VAL_OUTPUTS))
        System.out.println("INFO: '" + name + "' will be restored later.");
      else
        readFromXML(result, name, child);
    }

    bean = (MetaBean) result;
    
    // set inputs and outputs, after the beans have been instantiated
    bean.setInputs(getBeanInstancesForIDs(inputs));
    bean.setOutputs(getBeanInstancesForIDs(outputs));
    bean.setOriginalCoords(coords);
    
    return result;
  }
}
