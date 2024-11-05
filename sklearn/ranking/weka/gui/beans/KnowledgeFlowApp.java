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
 *    KnowledgeFlowApp.java
 *    Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.ClassloaderUtil;
import weka.core.Copyright;
import weka.core.Environment;
import weka.core.EnvironmentHandler;
import weka.core.Memory;
import weka.core.SerializedObject;
import weka.core.Utils;
import weka.core.xml.KOML;
import weka.core.xml.XStream;
import weka.gui.ExtensionFileFilter;
import weka.gui.GenericObjectEditor;
import weka.gui.GenericPropertiesCreator;
import weka.gui.HierarchyPropertyParser;
import weka.gui.LookAndFeel;
import weka.gui.beans.xml.XMLBeans;
import weka.gui.visualize.PrintablePanel;

import java.awt.BorderLayout;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.MenuItem;
import java.awt.Point;
import java.awt.PopupMenu;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.image.BufferedImage;
import java.beans.BeanInfo;
import java.beans.Beans;
import java.beans.Customizer;
import java.beans.EventSetDescriptor;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.MethodDescriptor;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.beancontext.BeanContextChild;
import java.beans.beancontext.BeanContextSupport;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.lang.reflect.Method;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Properties;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.JWindow;
import javax.swing.filechooser.FileFilter;

/**
 * Main GUI class for the KnowledgeFlow. Modifications to allow interoperability
 * with swt provided by Davide Zerbetto (davide dot zerbetto at eng dot it).
 *
 * @author Mark Hall
 * @version  $Revision: 6624 $
 * @since 1.0
 * @see JPanel
 * @see PropertyChangeListener
 */
public class KnowledgeFlowApp
  extends JPanel
  implements PropertyChangeListener {

  /** for serialization */
  private static final long serialVersionUID = -7064906770289728431L;

  /**
   * Location of the property file for the KnowledgeFlowApp
   */
  protected static String PROPERTY_FILE = "weka/gui/beans/Beans.props";

  /** Contains the editor properties */
  protected static Properties BEAN_PROPERTIES;

  private static ArrayList<Properties> BEAN_PLUGINS_PROPERTIES = 
    new ArrayList<Properties>();

  /**
   * Holds the details needed to construct button bars for various supported
   * classes of weka algorithms/tools 
   */
  private static Vector TOOLBARS = new Vector();
  
  public static void addToPluginBeanProps(File beanPropsFile) throws Exception {
    Properties tempP = new Properties();
    
    tempP.load(new FileInputStream(beanPropsFile));
    if (!BEAN_PLUGINS_PROPERTIES.contains(tempP)) {
      BEAN_PLUGINS_PROPERTIES.add(tempP);
    }
  }
  
  public static void removeFromPluginBeanProps(File beanPropsFile) throws Exception {
    Properties tempP = new Properties();
    
    tempP.load(new FileInputStream(beanPropsFile));
    if (BEAN_PLUGINS_PROPERTIES.contains(tempP)) {
      BEAN_PLUGINS_PROPERTIES.remove(tempP);
    }
  }

  /**
   * Loads KnowledgeFlow properties and any plugins (adds jars to
   * the classpath)
   */
  public static void loadProperties() {
    if (BEAN_PROPERTIES == null) {
      weka.core.WekaPackageManager.loadPackages(false);
      System.out.println("[KnowledgeFlow] Loading properties and plugins...");
      /** Loads the configuration property file */
      //  static {
      // Allow a properties file in the current directory to override
      try {
        BEAN_PROPERTIES = Utils.readProperties(PROPERTY_FILE);
        java.util.Enumeration keys =
          (java.util.Enumeration)BEAN_PROPERTIES.propertyNames();
        if (!keys.hasMoreElements()) {
          throw new Exception( "Could not read a configuration file for the bean\n"
                               +"panel. An example file is included with the Weka distribution.\n"
                               +"This file should be named \"" + PROPERTY_FILE + "\" and\n"
                               +"should be placed either in your user home (which is set\n"
                               + "to \"" + System.getProperties().getProperty("user.home") + "\")\n"
                               + "or the directory that java was started from\n");
        }
      } catch (Exception ex) {
        JOptionPane.showMessageDialog(null,
                                      ex.getMessage(),
                                      "KnowledgeFlow",
                                      JOptionPane.ERROR_MESSAGE);
      }


      // try and load any plugin beans properties
      File pluginDir = new File(System.getProperty("user.home")
                                +File.separator+".knowledgeFlow"
                                +File.separator+"plugins");
      /*if (pluginDir.exists() && pluginDir.isDirectory()) {
        BEAN_PLUGINS_PROPERTIES = new ArrayList<Properties>();
        // How many sub-dirs are there?
        File[] contents = pluginDir.listFiles();
        for (int i = 0; i < contents.length; i++) {
          if (contents[i].isDirectory() && 
              contents[i].listFiles().length > 0) {
            try {      
              Properties tempP = new Properties();
              File propFile = new File(contents[i].getPath()
                                       + File.separator
                                       + "Beans.props");
              tempP.load(new FileInputStream(propFile));
              BEAN_PLUGINS_PROPERTIES.add(tempP);

              // Now try and add all jar files in this directory to the classpath
              File anyJars[] = contents[i].listFiles();
              for (int j = 0; j < anyJars.length; j++) {
                if (anyJars[j].getPath().endsWith(".jar")) {
                  System.out.println("[KnowledgeFlow] Plugins: adding "+anyJars[j].getPath()
                                     +" to classpath...");
                  ClassloaderUtil.addFile(anyJars[j].getPath());
                }
              }
            } catch (Exception ex) {
              // Don't make a fuss
              System.err.println("[KnowledgeFlow] Warning: Unable to load bean properties for plugin "
                                 +"directory: " + contents[i].getPath());
            }
          }
          //        BEAN_PLUGINS_PROPERTIES = new Properties();
          //        BEAN_PLUGINS_PROPERTIES.load(new FileInputStream(pluginDir));
        }
      } else {
        // make the plugin directory for the user
        pluginDir.mkdir();
      } */
    }
  }
  
  public static void reInitialize() {
    if (BEAN_PROPERTIES == null) {
      loadProperties();
    }
    init();
  }

  /**
   * Initializes the temporary files necessary to construct the toolbars
   * from.
   */
  private static void init() {
    System.out.println("[KnowledgeFlow] Initializing KF...");

    try {
      TOOLBARS = new Vector();
      
      TreeMap wrapList = new TreeMap();
      GenericPropertiesCreator creator = new GenericPropertiesCreator();
      Properties GEOProps = null;

      if (creator.useDynamic()) {
        creator.execute(false);
        /* now process the keys in the GenericObjectEditor.props. For each
           key that has an entry in the Beans.props associating it with a
           bean component a button tool bar will be created */
        GEOProps = creator.getOutputProperties();
      } else {
        // Read the static information from the GenericObjectEditor.props
        GEOProps = Utils.readProperties("weka/gui/GenericObjectEditor.props");
      }
      Enumeration en = GEOProps.propertyNames();
      while (en.hasMoreElements()) {
	String geoKey = (String)en.nextElement();

	// try to match this key with one in the Beans.props file
	String beanCompName = BEAN_PROPERTIES.getProperty(geoKey);
	if (beanCompName != null) {
	  // add details necessary to construct a button bar for this class
	  // of algorithms
	  Vector newV = new Vector();
	  // check for a naming alias for this toolbar
	  String toolBarNameAlias = 
	    BEAN_PROPERTIES.getProperty(geoKey+".alias");
	  String toolBarName = (toolBarNameAlias != null) ?
	    toolBarNameAlias :
	    geoKey.substring(geoKey.lastIndexOf('.')+1, geoKey.length());

          // look for toolbar ordering information for this wrapper type
          String order = 
            BEAN_PROPERTIES.getProperty(geoKey+".order");
          Integer intOrder = (order != null) ?
            new Integer(order) :
            new Integer(0);
            
	  // Name for the toolbar (name of weka algorithm class)
	  newV.addElement(toolBarName);
	  // Name of bean capable of handling this class of algorithm
	  newV.addElement(beanCompName);

	  // add the root package for this key
	  String rootPackage = geoKey.substring(0, geoKey.lastIndexOf('.'));

	  newV.addElement(rootPackage);

	  // All the weka algorithms of this class of algorithm
	  String wekaAlgs = GEOProps.getProperty(geoKey);

          Hashtable roots = GenericObjectEditor.sortClassesByRoot(wekaAlgs);
          Hashtable hpps = new Hashtable();
          Enumeration enm = roots.keys();
          while (enm.hasMoreElements()) {
            String root = (String) enm.nextElement();
            String classes = (String) roots.get(root);
            weka.gui.HierarchyPropertyParser hpp = 
              new weka.gui.HierarchyPropertyParser();
            hpp.build(classes, ", ");
            //System.err.println(hpp.showTree());
            hpps.put(root, hpp);
          }

	  //------ test the HierarchyPropertyParser
          /*  weka.gui.HierarchyPropertyParser hpp = 
	    new weka.gui.HierarchyPropertyParser();
	  hpp.build(wekaAlgs, ", ");

	  System.err.println(hpp.showTree()); */
	  // ----- end test the HierarchyPropertyParser
          //	  newV.addElement(hpp); // add the hierarchical property parser
	  newV.addElement(hpps); // add the hierarchical property parser

	  StringTokenizer st = new StringTokenizer(wekaAlgs, ", ");
	  while (st.hasMoreTokens()) {
	    String current = st.nextToken().trim();
	    newV.addElement(current);
	  }
          wrapList.put(intOrder, newV);
          //	  TOOLBARS.addElement(newV);
	}
      }
      Iterator keysetIt = wrapList.keySet().iterator();
      while (keysetIt.hasNext()) {
        Integer key = (Integer)keysetIt.next();
        Vector newV = (Vector)wrapList.get(key);
        if (newV != null) {
          TOOLBARS.addElement(newV);
        }
      }
    } catch (Exception ex) {
      JOptionPane.showMessageDialog(null,
          "Could not read a configuration file for the generic objecte editor"
         +". An example file is included with the Weka distribution.\n"
         +"This file should be named \"GenericObjectEditor.props\" and\n"
         +"should be placed either in your user home (which is set\n"
         + "to \"" + System.getProperties().getProperty("user.home") + "\")\n"
         + "or the directory that java was started from\n",
         "KnowledgeFlow",
         JOptionPane.ERROR_MESSAGE);
    }

    try {
      String standardToolBarNames = 
	BEAN_PROPERTIES.
	getProperty("weka.gui.beans.KnowledgeFlow.standardToolBars");
      StringTokenizer st = new StringTokenizer(standardToolBarNames, ", ");
       while (st.hasMoreTokens()) {
	 String tempBarName = st.nextToken().trim();
	 // construct details for this toolbar
	 Vector newV = new Vector();
	 // add the name of the toolbar
	 newV.addElement(tempBarName);

	 // indicate that this is a standard toolbar (no wrapper bean)
	 newV.addElement("null");
	 String toolBarContents = 
	   BEAN_PROPERTIES.
	   getProperty("weka.gui.beans.KnowledgeFlow."+tempBarName);
	 StringTokenizer st2 = new StringTokenizer(toolBarContents, ", ");
	 while (st2.hasMoreTokens()) {
	   String tempBeanName = st2.nextToken().trim();
	   newV.addElement(tempBeanName);
	 }
	 TOOLBARS.addElement(newV);
       }       
    } catch (Exception ex) {
      JOptionPane.showMessageDialog(null,
				    ex.getMessage(),
				    "KnowledgeFlow",
				    JOptionPane.ERROR_MESSAGE);
    }
  } 
  
  /**
   * Used for displaying the bean components and their visible
   * connections
   *
   * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
   * @version $Revision: 6624 $
   * @since 1.0
   * @see PrintablePanel
   */
  protected class BeanLayout
    extends PrintablePanel {

    /** for serialization */
    private static final long serialVersionUID = -146377012429662757L;

    public void paintComponent(Graphics gx) {
      super.paintComponent(gx);
      BeanInstance.paintLabels(gx);
      BeanConnection.paintConnections(gx);
      //      BeanInstance.paintConnections(gx);
      if (m_mode == CONNECTING) {
	gx.drawLine(m_startX, m_startY, m_oldX, m_oldY);
      } else if (m_mode == SELECTING) {
        gx.drawRect((m_startX < m_oldX) ? m_startX : m_oldX, 
                    (m_startY < m_oldY) ? m_startY : m_oldY, 
                    Math.abs(m_oldX-m_startX), Math.abs(m_oldY-m_startY));
      }
    }

    public void doLayout() {
      super.doLayout();
      Vector comps = BeanInstance.getBeanInstances();
      for (int i = 0; i < comps.size(); i++) {
	BeanInstance bi = (BeanInstance)comps.elementAt(i);
	JComponent c = (JComponent)bi.getBean();
	Dimension d = c.getPreferredSize();
	c.setBounds(bi.getX(), bi.getY(), d.width, d.height);
	c.revalidate();
      }
    }
  }

  // Used for measuring and splitting icon labels
  // over multiple lines
  FontMetrics m_fontM;

  // constants for operations in progress
  protected static final int NONE = 0;
  protected static final int MOVING = 1;
  protected static final int CONNECTING = 2;
  protected static final int ADDING = 3;
  protected static final int SELECTING = 4;

  // which operation is in progress
  private int m_mode = NONE;

  /** the extension for the user components, when serialized to XML */
  protected final static String USERCOMPONENTS_XML_EXTENSION = ".xml";
  
  /**
   * Button group to manage all toolbar buttons
   */
  private ButtonGroup m_toolBarGroup = new ButtonGroup();

  /**
   * Holds the selected toolbar bean
   */
  private Object m_toolBarBean;

  /**
   * The layout area
   */
  private BeanLayout m_beanLayout = new BeanLayout();

  /**
   * Tabbed pane to hold tool bars
   */
  private JTabbedPane m_toolBars = new JTabbedPane();

  /**
   * Stuff relating to plugin beans
   */
  private JToolBar m_pluginsToolBar = null;
  private Box m_pluginsBoxPanel = null;
  
  /**
   * Stuff relating to user created meta beans
   */
  private JToolBar m_userToolBar = null;
  private Box m_userBoxPanel = null;
  private Vector m_userComponents = new Vector();
  private boolean m_firstUserComponentOpp = true;

  private JToggleButton m_pointerB;
  private JButton m_saveB;
  private JButton m_loadB;
  private JButton m_stopB;
  private JButton m_helpB;
  private JButton m_newB;

  /**
   * Reference to bean being manipulated
   */
  private BeanInstance m_editElement;

  /**
   * Event set descriptor for the bean being manipulated
   */
  private EventSetDescriptor m_sourceEventSetDescriptor;

  /**
   * Used to record screen coordinates during move, select and connect
   * operations
   */
  private int m_oldX, m_oldY;
  private int m_startX, m_startY;
  
  /** The file chooser for selecting layout files */
  protected JFileChooser m_FileChooser 
    = new JFileChooser(new File(System.getProperty("user.dir")));

  protected LogPanel m_logPanel = new LogPanel();//new LogPanel(null, true);

  protected BeanContextSupport m_bcSupport = new BeanContextSupport();

  /** the extension for the serialized setups (Java serialization) */
  public final static String FILE_EXTENSION = ".kf";

  /** the extension for the serialized setups (Java serialization) */
  public final static String FILE_EXTENSION_XML = ".kfml";
  
  /** A filter to ensure only KnowledgeFlow files in binary format get shown in
      the chooser */
  protected FileFilter m_KfFilter = 
    new ExtensionFileFilter(FILE_EXTENSION, 
                            "Binary KnowledgeFlow configuration files (*" 
                            + FILE_EXTENSION + ")");

  /** A filter to ensure only KnowledgeFlow files in KOML format 
      get shown in the chooser */
  protected FileFilter m_KOMLFilter = 
    new ExtensionFileFilter(KOML.FILE_EXTENSION + "kf", 
                            "XML KnowledgeFlow configuration files (*" 
                            + KOML.FILE_EXTENSION + "kf)");

  /** A filter to ensure only KnowledgeFlow files in XStream format 
      get shown in the chooser */
  protected FileFilter m_XStreamFilter = 
    new ExtensionFileFilter(XStream.FILE_EXTENSION + "kf", 
                            "XML KnowledgeFlow configuration files (*" 
                            + XStream.FILE_EXTENSION + "kf)");

  /** A filter to ensure only KnowledgeFlow layout files in XML format get 
      shown in the chooser */
  protected FileFilter m_XMLFilter = 
    new ExtensionFileFilter(FILE_EXTENSION_XML, 
                            "XML KnowledgeFlow layout files (*" 
                            + FILE_EXTENSION_XML + ")");

  /** the scrollbar increment of the layout scrollpane */
  protected int m_ScrollBarIncrementLayout = 20;

  /** the scrollbar increment of the components scrollpane */
  protected int m_ScrollBarIncrementComponents = 50;

  /** the flow layout width */
  protected int m_FlowWidth = 1024;

  /** the flow layout height */
  protected int m_FlowHeight = 768;

  /** the preferred file extension */
  protected String m_PreferredExtension = FILE_EXTENSION;
  
  /** whether to store the user components in XML or in binary format */
  protected boolean m_UserComponentsInXML = false;
  
  /** Environment variables for the current flow */
  protected Environment m_flowEnvironment = new Environment();
  
  /**
   * Set the environment variables to use. NOTE: loading a new layout
   * resets back to the default set of variables
   * 
   * @param env
   */
  public void setEnvironment(Environment env) {
    m_flowEnvironment = env;
    setEnvironment();
  }
  
  private void setEnvironment() {
    // pass m_flowEnvironment to all components
    // that implement EnvironmentHandler
    Vector beans = BeanInstance.getBeanInstances();
    for (int i = 0; i < beans.size(); i++) {
      Object temp = ((BeanInstance) beans.elementAt(i)).getBean();

      if (temp instanceof EnvironmentHandler) {
        ((EnvironmentHandler) temp).setEnvironment(m_flowEnvironment);
      }
    }
  }
  
  /**
   * Creates a new <code>KnowledgeFlowApp</code> instance.
   */
  // modifications by Zerbetto
  //public KnowledgeFlowApp() {
  public KnowledgeFlowApp(boolean showFileMenu) {
    if (BEAN_PROPERTIES == null) {
      loadProperties();
      init();
    }

    m_showFileMenu = showFileMenu;

    // end modifications by Zerbetto
    // Grab a fontmetrics object
    JWindow temp = new JWindow();
    temp.setVisible(true);
    temp.getGraphics().setFont(new Font(null, Font.PLAIN, 9));
    m_fontM = temp.getGraphics().getFontMetrics();
    temp.setVisible(false);

    // some GUI defaults
    try {
      m_ScrollBarIncrementLayout = Integer.parseInt(
          BEAN_PROPERTIES.getProperty(
            "ScrollBarIncrementLayout", "" + m_ScrollBarIncrementLayout));
      m_ScrollBarIncrementComponents = Integer.parseInt(
          BEAN_PROPERTIES.getProperty(
            "ScrollBarIncrementComponents", "" + m_ScrollBarIncrementComponents));
      m_FlowWidth = Integer.parseInt(
          BEAN_PROPERTIES.getProperty(
            "FlowWidth", "" + m_FlowWidth));
      m_FlowHeight = Integer.parseInt(
          BEAN_PROPERTIES.getProperty(
            "FlowHeight", "" + m_FlowHeight));
      m_PreferredExtension = BEAN_PROPERTIES.getProperty(
          "PreferredExtension", m_PreferredExtension);
      m_UserComponentsInXML = Boolean.valueOf(
          BEAN_PROPERTIES.getProperty(
            "UserComponentsInXML", "" + m_UserComponentsInXML)).booleanValue();
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }

    // FileChooser
    m_FileChooser.addChoosableFileFilter(m_KfFilter);
    if (KOML.isPresent()) {
      m_FileChooser.addChoosableFileFilter(m_KOMLFilter);
    }
    if (XStream.isPresent()) {
      m_FileChooser.addChoosableFileFilter(m_XStreamFilter);
    }

    m_FileChooser.addChoosableFileFilter(m_XMLFilter);

    if (m_PreferredExtension.equals(FILE_EXTENSION_XML)) {
      m_FileChooser.setFileFilter(m_XMLFilter);
    } else if (KOML.isPresent() && m_PreferredExtension.equals(KOML.FILE_EXTENSION + "kf")) {
      m_FileChooser.setFileFilter(m_KOMLFilter);
    } else if (XStream.isPresent() && m_PreferredExtension.equals(XStream.FILE_EXTENSION + "kf")) {
      m_FileChooser.setFileFilter(m_XStreamFilter);
    } else {
      m_FileChooser.setFileFilter(m_KfFilter);
    }
    m_FileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

    m_bcSupport.setDesignTime(true);
    m_beanLayout.setLayout(null);
    
    // handle mouse events
    m_beanLayout.addMouseListener(new MouseAdapter() {

	public void mousePressed(MouseEvent me) {
	  if (m_toolBarBean == null) {
	    if (((me.getModifiers() & InputEvent.BUTTON1_MASK)
		 == InputEvent.BUTTON1_MASK) && m_mode == NONE) {
	      BeanInstance bi = BeanInstance.findInstance(me.getPoint());
	      JComponent bc = null;
	      if (bi != null) {
		bc = (JComponent)(bi.getBean());
	      }
	      if (bc != null && (bc instanceof Visible)) {
		m_editElement = bi;
		m_oldX = me.getX();
		m_oldY = me.getY();
		m_mode = MOVING;
	      }
              if (m_mode != MOVING) {
                m_mode = SELECTING;
                m_oldX = me.getX();
                m_oldY = me.getY();
                m_startX = m_oldX;
                m_startY = m_oldY;
                Graphics2D gx = (Graphics2D)m_beanLayout.getGraphics();
                gx.setXORMode(java.awt.Color.white);
                //                gx.drawRect(m_oldX, m_oldY, m_oldX, m_oldY);
                //                gx.drawLine(m_startX, m_startY, m_startX, m_startY);
                gx.dispose();
                m_mode = SELECTING;
              }
	    }
	  }
	}

	public void mouseReleased(MouseEvent me) {
	  if (m_editElement != null && m_mode == MOVING) {
	    m_editElement = null;
	    revalidate();
	    m_beanLayout.repaint();
	    m_mode = NONE;
	  }
          if (m_mode == SELECTING) {
            revalidate();
            m_beanLayout.repaint();
            m_mode = NONE;
                        
            checkSubFlow(m_startX, m_startY, me.getX(), me.getY());
          }
	}

	public void mouseClicked(MouseEvent me) {
	  BeanInstance bi = BeanInstance.findInstance(me.getPoint());
	  if (m_mode == ADDING || m_mode == NONE) {
	    // try and popup a context sensitive menu if we have
	    // been clicked over a bean.
	    if (bi != null) {
	      JComponent bc = (JComponent)bi.getBean();
              // if we've been double clicked, then popup customizer
              // as long as we're not a meta bean
              if (me.getClickCount() == 2 && !(bc instanceof MetaBean)) {
                try {
                  Class custClass = 
                    Introspector.getBeanInfo(bc.getClass()).getBeanDescriptor().getCustomizerClass();
                  if (custClass != null) {
                    if (bc instanceof BeanCommon) {
                      if (!((BeanCommon)bc).
                          isBusy()) {
                        popupCustomizer(custClass, bc);
                      }
                    } else {
                      popupCustomizer(custClass, bc);
                    }
                  }
                } catch (IntrospectionException ex) {
                  ex.printStackTrace();
                }
              } else if (((me.getModifiers() & InputEvent.BUTTON1_MASK)
                          != InputEvent.BUTTON1_MASK) || me.isAltDown()) {
		doPopup(me.getPoint(), bi, me.getX(), me.getY());
	      }
	    } else {
	      if (((me.getModifiers() & InputEvent.BUTTON1_MASK)
		   != InputEvent.BUTTON1_MASK) || me.isAltDown()) {
		// find connections if any close to this point
		int delta = 10;
		deleteConnectionPopup(BeanConnection.
		      getClosestConnections(new Point(me.getX(), me.getY()), 
					    delta), me.getX(), me.getY());
	      } else if (m_toolBarBean != null) {
		// otherwise, if a toolbar button is active then 
		// add the component
		addComponent(me.getX(), me.getY());
	      }
	    }
	  }
	
	  if (m_mode == CONNECTING) {
	    // turn off connecting points and remove connecting line
	    m_beanLayout.repaint();
	    Vector beanInstances = BeanInstance.getBeanInstances();
	    for (int i = 0; i < beanInstances.size(); i++) {
	      JComponent bean = 
		(JComponent)((BeanInstance)beanInstances.elementAt(i)).
		getBean();
	      if (bean instanceof Visible) {
		((Visible)bean).getVisual().setDisplayConnectors(false);
	      }
	    }

	    if (bi != null) {
	      boolean doConnection = false;
	      if (!(bi.getBean() instanceof BeanCommon)) {
		doConnection = true;
	      } else {
		// Give the target bean a chance to veto the proposed
		// connection
		if (((BeanCommon)bi.getBean()).
		    //connectionAllowed(m_sourceEventSetDescriptor.getName())) {
		    connectionAllowed(m_sourceEventSetDescriptor)) {
		  doConnection = true;
		}
	      }
	      if (doConnection) {
		// attempt to connect source and target beans

                if (bi.getBean() instanceof MetaBean) {
                  BeanConnection.doMetaConnection(m_editElement, bi,
                                                  m_sourceEventSetDescriptor,
                                                  m_beanLayout);
                } else {
                  BeanConnection bc = 
                    new BeanConnection(m_editElement, bi, 
                                       m_sourceEventSetDescriptor);
                }
	      }
	      m_beanLayout.repaint();
	    }
	    m_mode = NONE;
	    m_editElement = null;
	    m_sourceEventSetDescriptor = null;
	  }
	}
      });
    
     m_beanLayout.addMouseMotionListener(new MouseMotionAdapter() {

	public void mouseDragged(MouseEvent me) {
	  if (m_editElement != null && m_mode == MOVING) {
	    ImageIcon ic = ((Visible)m_editElement.getBean()).
	      getVisual().getStaticIcon();
	    int width = ic.getIconWidth() / 2;
	    int height = ic.getIconHeight() / 2;

            /*	    m_editElement.setX(m_oldX-width);
                    m_editElement.setY(m_oldY-height); */

            m_editElement.setXY(m_oldX-width,
                                m_oldY-height);
	    m_beanLayout.repaint();
	    
	    // note the new points
	    m_oldX = me.getX(); m_oldY = me.getY();
	  }
          if (m_mode == SELECTING) {
            m_beanLayout.repaint();
            m_oldX = me.getX(); m_oldY = me.getY();
          }
	}

	 public void mouseMoved(MouseEvent e) {
	   if (m_mode == CONNECTING) {
	     m_beanLayout.repaint();
	     // note the new coordinates
	     m_oldX = e.getX(); m_oldY = e.getY();
	   }
	 }
       });
     
     String date = (new SimpleDateFormat("EEEE, d MMMM yyyy"))
       .format(new Date());
     m_logPanel.logMessage("Weka Knowledge Flow was written by Mark Hall");
     m_logPanel.logMessage("Weka Knowledge Flow");
     m_logPanel.logMessage("(c) 2002-" + Copyright.getToYear() + " " 
	 + Copyright.getOwner() + ", " + Copyright.getAddress());
     m_logPanel.logMessage("web: " + Copyright.getURL());
     m_logPanel.logMessage( date);
     m_logPanel.statusMessage("[KnowledgeFlow]|Welcome to the Weka Knowledge Flow");
     m_logPanel.getStatusTable().addMouseListener(new MouseAdapter() {
       public void mouseClicked(MouseEvent e) {
         if (m_logPanel.getStatusTable().rowAtPoint(e.getPoint()) == 0) {
           if (((e.getModifiers() & InputEvent.BUTTON1_MASK)
               != InputEvent.BUTTON1_MASK) || e.isAltDown()) {
             System.gc();
             Runtime currR = Runtime.getRuntime();
             long freeM = currR.freeMemory();
             long totalM = currR.totalMemory();
             long maxM = currR.maxMemory();
             m_logPanel.
             logMessage("[KnowledgeFlow] Memory (free/total/max.) in bytes: " 
                 + String.format("%,d", freeM) + " / " 
                 + String.format("%,d", totalM) + " / " 
                 + String.format("%,d", maxM));
             m_logPanel.statusMessage("[KnowledgeFlow]|Memory (free/total/max.) in bytes: " 
                 + String.format("%,d", freeM) + " / " 
                 + String.format("%,d", totalM) + " / " 
                 + String.format("%,d", maxM)); 
           }
         }
       }
     });
    
     JPanel p1 = new JPanel();
     p1.setLayout(new BorderLayout());
     p1.setBorder(javax.swing.BorderFactory.createCompoundBorder(
			    javax.swing.BorderFactory.
			    createTitledBorder("Knowledge Flow Layout"),
                   javax.swing.BorderFactory.createEmptyBorder(0, 5, 5, 5)
                   ));
     final JScrollPane js = new JScrollPane(m_beanLayout);
     p1.add(js, BorderLayout.CENTER);
     js.getVerticalScrollBar().setUnitIncrement(m_ScrollBarIncrementLayout);
     js.getHorizontalScrollBar().setUnitIncrement(m_ScrollBarIncrementLayout);

     setLayout(new BorderLayout());
     
     add(p1, BorderLayout.CENTER);
     m_beanLayout.setSize(m_FlowWidth, m_FlowHeight);
     Dimension d = m_beanLayout.getPreferredSize();
     m_beanLayout.setMinimumSize(d);
     m_beanLayout.setMaximumSize(d);
     m_beanLayout.setPreferredSize(d);

     Dimension d2 = new Dimension(100, 170);
     m_logPanel.setPreferredSize(d2);
     m_logPanel.setMinimumSize(d2);
     add(m_logPanel, BorderLayout.SOUTH);
     
     setUpToolBars();
     loadUserComponents();
  }
  
  private Image loadImage(String path) {
    Image pic = null;
    // Modified by Zerbetto
    //java.net.URL imageURL = ClassLoader.getSystemResource(path);
    java.net.URL imageURL = this.getClass().getClassLoader().getResource(path);

    // end modifications
    if (imageURL == null) {
      //      System.err.println("Warning: unable to load "+path);
    } else {
      pic = Toolkit.getDefaultToolkit().
	getImage(imageURL);
    }
    return pic;
  }

  /**
   * Describe <code>setUpToolBars</code> method here.
   */
  private void setUpToolBars() {
    JPanel toolBarPanel = new JPanel();
    toolBarPanel.setLayout(new BorderLayout());

    // modifications by Zerbetto
    // first construct the toolbar for saving, loading etc
    if (m_showFileMenu) {
      JToolBar fixedTools = new JToolBar();
      fixedTools.setOrientation(JToolBar.VERTICAL);
      m_saveB = new JButton(new ImageIcon(loadImage(BeanVisual.ICON_PATH +
              "Save24.gif")));
      m_saveB.setToolTipText("Save layout");
      m_loadB = new JButton(new ImageIcon(loadImage(BeanVisual.ICON_PATH +
              "Open24.gif")));
      m_loadB.setToolTipText("Load layout");
      m_newB = new JButton(new ImageIcon(loadImage(BeanVisual.ICON_PATH +
              "New24.gif")));
      m_newB.setToolTipText("Clear the layout");
      fixedTools.add(m_newB);
      fixedTools.add(m_saveB);
      fixedTools.add(m_loadB);

      m_saveB.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent e) {
            saveLayout();
          }
        });

      m_loadB.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent e) {
            m_flowEnvironment = new Environment();
            loadLayout();
          }
        });

      m_newB.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent ae) {
            clearLayout();
          }
        });

      fixedTools.setFloatable(false);
      toolBarPanel.add(fixedTools, BorderLayout.WEST);
    }

    m_stopB = new JButton(new ImageIcon(loadImage(BeanVisual.ICON_PATH +
            "Stop24.gif")));
    m_helpB = new JButton(new ImageIcon(loadImage(BeanVisual.ICON_PATH +
            "Help24.gif")));
    m_stopB.setToolTipText("Stop all execution");
    m_helpB.setToolTipText("Display help");

    Image tempI = loadImage(BeanVisual.ICON_PATH + "Pointer.gif");
    m_pointerB = new JToggleButton(new ImageIcon(tempI));
    m_pointerB.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          m_toolBarBean = null;
          m_mode = NONE;
          setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
        }
      });

    //    Dimension dP = m_saveB.getPreferredSize();
    //    Dimension dM = m_saveB.getMaximumSize();
    //    Dimension dP = m_stopB.getPreferredSize();
    //    Dimension dM = m_stopB.getMaximumSize();
    //    m_pointerB.setPreferredSize(dP);
    //    m_pointerB.setMaximumSize(dM);
    m_toolBarGroup.add(m_pointerB);

    JToolBar fixedTools2 = new JToolBar();
    fixedTools2.setOrientation(JToolBar.VERTICAL);
    fixedTools2.setFloatable(false);
    fixedTools2.add(m_pointerB);
    fixedTools2.add(m_helpB);
    fixedTools2.add(m_stopB);
    //    m_helpB.setPreferredSize(dP);
    //    m_helpB.setMaximumSize(dP);
    m_helpB.setSize(m_pointerB.getSize().width, m_pointerB.getSize().height);
    toolBarPanel.add(fixedTools2, BorderLayout.EAST);
    // end modifications by Zerbetto
    m_stopB.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          m_logPanel.statusMessage("[KnowledgeFlow]|Attempting to stop all components...");
          stopFlow();
          m_logPanel.statusMessage("[KnowledgeFlow]|OK.");
        }
      });

    m_helpB.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent ae) {
          popupHelp();
        }
      });

    final int STANDARD_TOOLBAR = 0;
    final int WEKAWRAPPER_TOOLBAR = 1;

    int toolBarType = STANDARD_TOOLBAR;

    // set up wrapper toolbars
    for (int i = 0; i < TOOLBARS.size(); i++) {
      Vector tempBarSpecs = (Vector) TOOLBARS.elementAt(i);

      // name for the tool bar
      String tempBarName = (String) tempBarSpecs.elementAt(0);

      // Used for weka leaf packages 
      Box singletonHolderPanel = null;

      // name of the bean component to handle this class of weka algorithms
      String tempBeanCompName = (String) tempBarSpecs.elementAt(1);

      // a JPanel holding an instantiated bean + label ready to be added
      // to the current toolbar
      JPanel tempBean;

      // the root package for weka algorithms
      String rootPackage = "";
      weka.gui.HierarchyPropertyParser hpp = null;
      Hashtable hpps = null;

      // Is this a wrapper toolbar?
      if (tempBeanCompName.compareTo("null") != 0) {
        tempBean = null;
        toolBarType = WEKAWRAPPER_TOOLBAR;
        rootPackage = (String) tempBarSpecs.elementAt(2);
        //	hpp = (weka.gui.HierarchyPropertyParser)tempBarSpecs.elementAt(3);
        hpps = (Hashtable) tempBarSpecs.elementAt(3);

        try {
          // modifications by Zerbetto
          // Beans.instantiate(null, tempBeanCompName);
          Beans.instantiate(this.getClass().getClassLoader(), tempBeanCompName);

          // end modifications by Zerbetto
        } catch (Exception ex) {
          // ignore
          System.err.println("[KnowledgeFlow] Failed to instantiate: " + tempBeanCompName);

          break;
        }
      } else {
        toolBarType = STANDARD_TOOLBAR;
      }

      // a toolbar to hold buttons---one for each algorithm
      JToolBar tempToolBar = new JToolBar();

      //      System.err.println(tempToolBar.getLayout());
      //      tempToolBar.setLayout(new FlowLayout());
      int z = 2;

      if (toolBarType == WEKAWRAPPER_TOOLBAR) {
        Enumeration enm = hpps.keys();

        while (enm.hasMoreElements()) {
          String root = (String) enm.nextElement();
          String userPrefix = "";
          hpp = (HierarchyPropertyParser) hpps.get(root);

          if (!hpp.goTo(rootPackage)) {
            System.out.println("[KnowledgeFlow] Processing user package... ");
            //            System.exit(1);
            userPrefix = root + ".";
          }

          String[] primaryPackages = hpp.childrenValues();

          for (int kk = 0; kk < primaryPackages.length; kk++) {
            hpp.goToChild(primaryPackages[kk]);

            // check to see if this is a leaf - if so then there are no
            // sub packages
            if (hpp.isLeafReached()) {
              if (singletonHolderPanel == null) {
                singletonHolderPanel = Box.createHorizontalBox();
                singletonHolderPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(
                    tempBarName));
              }

              String algName = hpp.fullValue();
              tempBean = instantiateToolBarBean(true, tempBeanCompName, algName);

              if (tempBean != null) {
                // tempToolBar.add(tempBean);
                singletonHolderPanel.add(tempBean);
              }

              hpp.goToParent();
            } else {
              // make a titledborder JPanel to hold all the schemes in this
              // package
              //	    JPanel holderPanel = new JPanel();
              Box holderPanel = Box.createHorizontalBox();
              holderPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(userPrefix +
                  primaryPackages[kk]));
              processPackage(holderPanel, tempBeanCompName, hpp);
              tempToolBar.add(holderPanel);
            }
          }

          if (singletonHolderPanel != null) {
            tempToolBar.add(singletonHolderPanel);
            singletonHolderPanel = null;
          }
        }
      } else {
        Box holderPanel = Box.createHorizontalBox();
        holderPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(
            tempBarName));

        for (int j = z; j < tempBarSpecs.size(); j++) {
          tempBean = null;
          tempBeanCompName = (String) tempBarSpecs.elementAt(j);
          tempBean = instantiateToolBarBean((toolBarType == WEKAWRAPPER_TOOLBAR),
              tempBeanCompName, "");

          if (tempBean != null) {
            // set tool tip text (if any)
            // setToolTipText(tempBean)
            holderPanel.add(tempBean);
          }
        }

        tempToolBar.add(holderPanel);
      }

      JScrollPane tempJScrollPane = createScrollPaneForToolBar(tempToolBar);
      // ok, now create tabbed pane to hold this toolbar
      m_toolBars.addTab(tempBarName, null, tempJScrollPane, tempBarName);
    }

    // Any plugin components to process?
    if (BEAN_PLUGINS_PROPERTIES != null && 
        BEAN_PLUGINS_PROPERTIES.size() > 0) {
      for (int i = 0; i < BEAN_PLUGINS_PROPERTIES.size(); i++) {
        Properties tempP = BEAN_PLUGINS_PROPERTIES.get(i);
        JPanel tempBean = null;
        String components = 
        tempP.getProperty("weka.gui.beans.KnowledgeFlow.Plugins");
        StringTokenizer st2 = new StringTokenizer(components, ", ");

        while (st2.hasMoreTokens()) {
          String tempBeanCompName = st2.nextToken().trim();
          tempBean = instantiateToolBarBean(false, tempBeanCompName, "");
          if (m_pluginsToolBar == null) {
            // need to create the plugins tab and toolbar
            setUpPluginsToolBar();
          }
          m_pluginsBoxPanel.add(tempBean);
        }
      }
    }

    toolBarPanel.add(m_toolBars, BorderLayout.CENTER);

    //    add(m_toolBars, BorderLayout.NORTH);
    add(toolBarPanel, BorderLayout.NORTH);
  }
  
  private void stopFlow() {
    Vector components = BeanInstance.getBeanInstances();

    for (int i = 0; i < components.size(); i++) {
      Object temp = ((BeanInstance) components.elementAt(i)).getBean();

      if (temp instanceof BeanCommon) {
        ((BeanCommon) temp).stop();
      }
    }
  }


  private JScrollPane createScrollPaneForToolBar(JToolBar tb) {
    JScrollPane tempJScrollPane = 
      new JScrollPane(tb, 
                      JScrollPane.VERTICAL_SCROLLBAR_NEVER,
                      JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
    
    Dimension d = tb.getPreferredSize();
    tempJScrollPane.setMinimumSize(new Dimension((int)d.getWidth(),
                                                 (int)(d.getHeight()+15)));
    tempJScrollPane.setPreferredSize(new Dimension((int)d.getWidth(),
                                                   (int)(d.getHeight()+15)));
    tempJScrollPane.getHorizontalScrollBar().setUnitIncrement(
        m_ScrollBarIncrementComponents);

    return tempJScrollPane;
  }

  private void processPackage(JComponent holderPanel,
			      String tempBeanCompName,
			      weka.gui.HierarchyPropertyParser hpp) {
    if (hpp.isLeafReached()) {
      // instantiate a bean and add it to the holderPanel
      //      System.err.println("Would add "+hpp.fullValue());
      String algName = hpp.fullValue();
      JPanel tempBean = 
	instantiateToolBarBean(true, tempBeanCompName, algName);
      if (tempBean != null) {
	holderPanel.add(tempBean);
      }
      hpp.goToParent();
      return;
    }
    String [] children = hpp.childrenValues();
    for (int i = 0; i < children.length; i++) {
      hpp.goToChild(children[i]);
      processPackage(holderPanel, tempBeanCompName, hpp);
    }
    hpp.goToParent();
  }

  /**
   * Instantiates a bean for display in the toolbars
   *
   * @param wekawrapper true if the bean to be instantiated is a wekawrapper
   * @param tempBeanCompName the name of the bean to instantiate
   * @param algName holds the name of a weka algorithm to configure the
   * bean with if it is a wekawrapper bean
   * @return a JPanel holding the instantiated (and configured bean)
   */
  private JPanel instantiateToolBarBean(boolean wekawrapper, 
					String tempBeanCompName,
					String algName) {
    Object tempBean;
    if (wekawrapper) {
      try {
        // modifications by Zerbetto
        //tempBean = Beans.instantiate(null, tempBeanCompName);
        tempBean = Beans.instantiate(this.getClass().getClassLoader(),
                                     tempBeanCompName);
        
        // end modifications by Zerbetto
      } catch (Exception ex) {
	System.err.println("[KnowledgeFlow] Failed to instantiate :"+tempBeanCompName
			   +"KnowledgeFlowApp.instantiateToolBarBean()");
	return null;
      }
      if (tempBean instanceof WekaWrapper) {
	//	algName = (String)tempBarSpecs.elementAt(j);
	Class c = null;
	try {
	  c = Class.forName(algName);
	} catch (Exception ex) {
	  System.err.println("[KnowledgeFlow] Can't find class called: "+algName);
	  return null;
	}
	try {
	  Object o = c.newInstance();
	  ((WekaWrapper)tempBean).setWrappedAlgorithm(o);
	} catch (Exception ex) {
	  System.err.println("[KnowledgeFlow] Failed to configure "+tempBeanCompName
			     +" with "+algName);
	  return null;
	}
      }
    } else {
      try {
        // modifications by Zerbetto
        //tempBean = Beans.instantiate(null, tempBeanCompName);
        tempBean = Beans.instantiate(this.getClass().getClassLoader(),
            tempBeanCompName);

        // end modifications
      } catch (Exception ex) {
	ex.printStackTrace();
	System.err.println("[KnowledgeFlow] Failed to instantiate :"+tempBeanCompName
			   +"KnowledgeFlowApp.setUpToolBars()");
	return null;
      }
    }
    
    if (tempBean instanceof BeanContextChild) {
      m_bcSupport.add(tempBean);
    }
    if (tempBean instanceof Visible) {
      ((Visible)tempBean).getVisual().scale(3);
    }

    return makeHolderPanelForToolBarBean(tempBeanCompName, tempBean, 
                                         wekawrapper, algName, false);
  }

  /**
   * Instantiates (by making a serialized copy) the supplied
   * template meta bean for display in the user tool bar
   * 
   * @param bean the prototype MetaBean to display in the toolbar
   */
  private JPanel instantiateToolBarMetaBean(MetaBean bean) {
    // copy the bean via serialization
    ((Visible)bean).getVisual().removePropertyChangeListener(this);
    bean.removePropertyChangeListenersSubFlow(this);
    Object copy = null;
    try {
      SerializedObject so = new SerializedObject(bean);
      copy = (MetaBean)so.getObject();
    } catch (Exception ex) {
      ex.printStackTrace();
      return null;
    }
    ((Visible)bean).getVisual().addPropertyChangeListener(this);
    bean.addPropertyChangeListenersSubFlow(this);

    String displayName ="";
    //
    if (copy instanceof Visible) {
      ((Visible)copy).getVisual().scale(3);
      displayName = ((Visible)copy).getVisual().getText();
    }
    return makeHolderPanelForToolBarBean(displayName,
                                         copy,
                                         false,
                                         null,
                                         true);
  }

  private JPanel makeHolderPanelForToolBarBean(final String tempName,
                                               Object tempBean,
                                               boolean wekawrapper,
                                               String algName,
                                               final boolean metabean) {
    // ---------------------------------------
    JToggleButton tempButton;
    final JPanel tempP = new JPanel();
    JLabel tempL = new JLabel();
    tempL.setFont(new Font(null, Font.PLAIN, 9));

    String labelName = (wekawrapper == true) 
      ? algName 
      : tempName;
    labelName = labelName.substring(labelName.lastIndexOf('.')+1, 
				    labelName.length());
    tempL.setText(" "+labelName+" ");
    tempL.setHorizontalAlignment(JLabel.CENTER);
    tempP.setLayout(new BorderLayout());

    if (tempBean instanceof Visible) {
      BeanVisual bv = ((Visible)tempBean).getVisual();

      tempButton = 
	new JToggleButton(bv.getStaticIcon());
      int width = bv.getStaticIcon().getIconWidth();
      int height = bv.getStaticIcon().getIconHeight();
      
      JPanel labelPanel = 
	multiLineLabelPanel(labelName, width);
      tempP.add(labelPanel, BorderLayout.SOUTH);
    } else {
      tempButton = new JToggleButton();
      tempP.add(tempL, BorderLayout.SOUTH);
    }
    tempP.add(tempButton, BorderLayout.NORTH);
    //    tempP.add(tempL, BorderLayout.SOUTH);
    
    //  holderPanel.add(tempP);
    //    tempToolBar.add(tempP);
    m_toolBarGroup.add(tempButton);
    
    // add an action listener for the button here
    final Object tempBN = tempBean;
    final JToggleButton fButton = tempButton;
    //	  final JToggleButton tempButton2 = tempButton;
    tempButton.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
          boolean changeCursor = true;
	  try {
	    m_toolBarBean = null;
            if (metabean) {
              if ((e.getModifiers() & ActionEvent.SHIFT_MASK) != 0) {
                changeCursor = false;
                m_toolBarGroup.remove(fButton);
                m_userBoxPanel.remove(tempP);
                m_userBoxPanel.revalidate();
                m_userComponents.remove(tempBN);
                if (m_firstUserComponentOpp) {
                  installWindowListenerForSavingUserBeans();
                  m_firstUserComponentOpp = false;
                }
                if (m_userComponents.size() == 0) {
                  m_toolBars.removeTabAt(m_toolBars.getTabCount() - 1);
                  m_userToolBar = null;
                  notifyIsDirty();
                }
              } else {
                SerializedObject so = new SerializedObject(tempBN);
                MetaBean copy = (MetaBean)so.getObject();
                /*((Visible)copy).getVisual().
                  addPropertyChangeListener(KnowledgeFlowApp.this); */
                copy.addPropertyChangeListenersSubFlow(KnowledgeFlowApp.this);
                m_toolBarBean = copy;
              }
            } else {
              // modifications by Zerbetto
              //m_toolBarBean = Beans.instantiate(null, tempName);
              m_toolBarBean = Beans.instantiate(this.getClass().getClassLoader(),
                  tempName);

              // end modifications
            }
            if (m_toolBarBean instanceof WekaWrapper) {
	      Object wrappedAlg = 
		((WekaWrapper)tempBN).getWrappedAlgorithm();
	      
	      ((WekaWrapper)m_toolBarBean).
                setWrappedAlgorithm(wrappedAlg.getClass().newInstance());
	      //		    tempButton2.setSelected(false);
	    }
            if (changeCursor) {
              setCursor(Cursor.
                        getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
              m_mode = ADDING;
            }
	  } catch (Exception ex) {
	    System.err.
	      println("[KnowledgeFlow] Problem adding bean to data flow layout");
            ex.printStackTrace();
	  }
          notifyIsDirty();
	}
      });
    
    if (tempBean instanceof MetaBean) {
      tempButton.setToolTipText("Hold down shift and click to remove");
      m_userComponents.add(tempBean);
    } else {
      // set tool tip text from global info if supplied
      String summary = getGlobalInfo(tempBean);
      if (summary != null) {
        int ci = summary.indexOf('.');
        if (ci != -1) {
          summary = summary.substring(0, ci + 1);
        }
        tempButton.setToolTipText(summary);
      }
    }

    //return tempBean;
    return tempP;
  }

  private JPanel multiLineLabelPanel(String sourceL,
				     int splitWidth) {
    JPanel jp = new JPanel();
    Vector v = new Vector();

    int labelWidth = m_fontM.stringWidth(sourceL);

    if (labelWidth < splitWidth) {
      v.addElement(sourceL);
    } else {
      // find mid point
      int mid = sourceL.length() / 2;
      
      // look for split point closest to the mid
      int closest = sourceL.length();
      int closestI = -1;
      for (int i = 0; i < sourceL.length(); i++) {
	if (sourceL.charAt(i) < 'a') {
	  if (Math.abs(mid - i) < closest) {
	    closest = Math.abs(mid - i);
	    closestI = i;
	  }
	}
      }
      if (closestI != -1) {
	String left = sourceL.substring(0, closestI);
	String right = sourceL.substring(closestI, sourceL.length());
	if (left.length() > 1 && right.length() > 1) {
	  v.addElement(left);
	  v.addElement(right);
	} else {
	  v.addElement(sourceL);
	}
      } else {
	v.addElement(sourceL);
      }
    }

    jp.setLayout(new GridLayout(v.size(), 1));
    for (int i = 0; i < v.size(); i++) {
      JLabel temp = new JLabel();
      temp.setFont(new Font(null, Font.PLAIN, 9));
      temp.setText(" "+((String)v.elementAt(i))+" ");
      temp.setHorizontalAlignment(JLabel.CENTER);
      jp.add(temp);
    }
    return jp;
  }

  private void setUpUserToolBar() {
    m_userBoxPanel = Box.createHorizontalBox();
    m_userBoxPanel.setBorder(javax.swing.BorderFactory.
                             createTitledBorder("User"));
    m_userToolBar = new JToolBar();
    m_userToolBar.add(m_userBoxPanel);
    JScrollPane tempJScrollPane = 
      createScrollPaneForToolBar(m_userToolBar);
    // ok, now create tabbed pane to hold this toolbar
    
    m_toolBars.addTab("User", null, 
                      tempJScrollPane,
                      "User created components");
  }

  private void setUpPluginsToolBar() {
    m_pluginsBoxPanel = Box.createHorizontalBox();
    m_pluginsBoxPanel.setBorder(javax.swing.BorderFactory.
                                createTitledBorder("Plugins"));
    m_pluginsToolBar = new JToolBar();
    m_pluginsToolBar.add(m_pluginsBoxPanel);
    JScrollPane tempJScrollPane = 
      createScrollPaneForToolBar(m_pluginsToolBar);
    // ok, now create tabbed pane to hold this toolbar
    
    m_toolBars.addTab("Plugins", null, 
                      tempJScrollPane,
                      "Plugin components");
  }

  /**
   * Pop up a help window
   */
  private void popupHelp() {
    final JButton tempB = m_helpB;
    try {
      tempB.setEnabled(false);
      // Modified by Zerbetto
      //InputStream inR = 
      //	ClassLoader.
      //        getSystemResourceAsStream("weka/gui/beans/README_KnowledgeFlow");
      InputStream inR = this.getClass().getClassLoader()
                            .getResourceAsStream("weka/gui/beans/README_KnowledgeFlow");

      // end modifications
      StringBuffer helpHolder = new StringBuffer();
      LineNumberReader lnr = new LineNumberReader(new InputStreamReader(inR));
      
      String line;
      
      while ((line = lnr.readLine()) != null) {
	helpHolder.append(line+"\n");
      }
      
      lnr.close();
      final javax.swing.JFrame jf = new javax.swing.JFrame();
      jf.getContentPane().setLayout(new java.awt.BorderLayout());
      final JTextArea ta = new JTextArea(helpHolder.toString());
      ta.setFont(new Font("Monospaced", Font.PLAIN, 12));
      ta.setEditable(false);
      final JScrollPane sp = new JScrollPane(ta);
      jf.getContentPane().add(sp, java.awt.BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
        public void windowClosing(java.awt.event.WindowEvent e) {
	  tempB.setEnabled(true);
          jf.dispose();
        }
      });
      jf.setSize(600,600);
      jf.setVisible(true);
      
    } catch (Exception ex) {
      tempB.setEnabled(true);
    }
  }

  public void clearLayout() {
    stopFlow(); // try and stop any running components
    BeanInstance.reset(m_beanLayout);
    BeanConnection.reset();
    m_beanLayout.revalidate();
    m_beanLayout.repaint();
    m_logPanel.clearStatus();
    m_logPanel.statusMessage("[KnowledgeFlow]|Welcome to the Weka Knowledge Flow");
  }
  
  /**
   * Popup a context sensitive menu for the bean component
   *
   * @param pt holds the panel coordinates for the component
   * @param bi the bean component over which the user right clicked the mouse
   * @param x the x coordinate at which to popup the menu
   * @param y the y coordinate at which to popup the menu
   *
   * Modified by Zerbetto: javax.swing.JPopupMenu transformed into java.awt.PopupMenu
   *
   */
  private void doPopup(Point pt, final BeanInstance bi, int x, int y) {
    final JComponent bc = (JComponent) bi.getBean();
    final int xx = x;
    final int yy = y;
    int menuItemCount = 0;

    // modifications by Zerbetto
    PopupMenu beanContextMenu = new PopupMenu();

    //JPopupMenu beanContextMenu = new JPopupMenu();

    //    beanContextMenu.insert(new JLabel("Edit", 
    //				      SwingConstants.CENTER), 
    //			   menuItemCount);
    MenuItem edit = new MenuItem("Edit:");
    edit.setEnabled(false);
    beanContextMenu.insert(edit, menuItemCount);
    menuItemCount++;

    if (bc instanceof MetaBean) {
      //JMenuItem ungroupItem = new JMenuItem("Ungroup");
      MenuItem ungroupItem = new MenuItem("Ungroup");
      ungroupItem.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent e) {
            // ungroup
            bi.removeBean(m_beanLayout);

            Vector group = ((MetaBean) bc).getBeansInSubFlow();
            Vector associatedConnections = ((MetaBean) bc).getAssociatedConnections();
            ((MetaBean) bc).restoreBeans();

            for (int i = 0; i < group.size(); i++) {
              BeanInstance tbi = (BeanInstance) group.elementAt(i);
              addComponent(tbi, false);
              tbi.addBean(m_beanLayout);
            }

            for (int i = 0; i < associatedConnections.size(); i++) {
              BeanConnection tbc = (BeanConnection) associatedConnections.elementAt(i);
              tbc.setHidden(false);
            }

            m_beanLayout.repaint();
            notifyIsDirty();
          }
        });
      beanContextMenu.add(ungroupItem);
      menuItemCount++;

      // Add to user tab
      //JMenuItem addToUserTabItem = new JMenuItem("Add to user tab");
      MenuItem addToUserTabItem = new MenuItem("Add to user tab");
      addToUserTabItem.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent e) {
            addToUserToolBar((MetaBean) bi.getBean(), true);
            notifyIsDirty();
          }
        });
      beanContextMenu.add(addToUserTabItem);
      menuItemCount++;
    }

    //JMenuItem deleteItem = new JMenuItem("Delete");
    MenuItem deleteItem = new MenuItem("Delete");
    deleteItem.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          BeanConnection.removeConnections(bi);
          bi.removeBean(m_beanLayout);
          if (bc instanceof BeanCommon) {            
            String key = ((BeanCommon)bc).getCustomName()
              + "$" + bc.hashCode();
            m_logPanel.statusMessage(key + "|remove");
          }
          revalidate();
          notifyIsDirty();
        }
      });
    if (bc instanceof BeanCommon) {
      if (((BeanCommon)bc).isBusy()) {
        deleteItem.setEnabled(false);
      }
    }
    beanContextMenu.add(deleteItem);
    menuItemCount++;

    if (bc instanceof BeanCommon) {
      MenuItem nameItem = new MenuItem("Set name");
      nameItem.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent e) {
            String oldName = ((BeanCommon)bc).getCustomName();
            String name = JOptionPane.showInputDialog(KnowledgeFlowApp.this,
                                                      "Enter a name for this component",
                                                      oldName);
            if (name != null) {
              ((BeanCommon)bc).setCustomName(name);
            }
          }
        });
      if (bc instanceof BeanCommon) {
        if (((BeanCommon)bc).isBusy()) {
          nameItem.setEnabled(false);
        }
      }
      beanContextMenu.add(nameItem);
      menuItemCount++;
    }

    try {
      //BeanInfo [] compInfo = null;
      //JComponent [] associatedBeans = null;
      Vector compInfo = new Vector(1);
      Vector associatedBeans = null;
      Vector outputBeans = null;
      Vector compInfoOutputs = null;

      if (bc instanceof MetaBean) {
        compInfo = ((MetaBean) bc).getBeanInfoSubFlow();
        associatedBeans = ((MetaBean) bc).getBeansInSubFlow();

        outputBeans = ((MetaBean) bc).getBeansInOutputs();
        compInfoOutputs = ((MetaBean) bc).getBeanInfoOutputs();
      } else {
        compInfo.add(Introspector.getBeanInfo(bc.getClass()));
        compInfoOutputs = compInfo;
      }

      final Vector tempAssociatedBeans = associatedBeans;

      if (compInfo == null) {
        System.err.println("[KnowledgeFlow] Error in doPopup()");
      } else {
        //	System.err.println("Got bean info");
        for (int zz = 0; zz < compInfo.size(); zz++) {
          final int tt = zz;
          final Class custClass = ((BeanInfo) compInfo.elementAt(zz)).getBeanDescriptor()
            .getCustomizerClass();

          if (custClass != null) {
            //	  System.err.println("Got customizer class");
            //	  popupCustomizer(custClass, bc);
            //JMenuItem custItem = null;
            MenuItem custItem = null;
            boolean customizationEnabled = true;

            if (!(bc instanceof MetaBean)) {
              //custItem = new JMenuItem("Configure...");
              custItem = new MenuItem("Configure...");
              if (bc instanceof BeanCommon) {
                customizationEnabled = 
                  !((BeanCommon)bc).isBusy();
              }
            } else {
              String custName = custClass.getName();
              BeanInstance tbi = (BeanInstance) associatedBeans.elementAt(zz);
              if (tbi.getBean() instanceof BeanCommon) {
                custName = ((BeanCommon)tbi.getBean()).getCustomName();
              } else {
                if (tbi.getBean() instanceof WekaWrapper) {
                  custName = ((WekaWrapper) tbi.getBean()).getWrappedAlgorithm()
                  .getClass().getName();
                } else {
                  custName = custName.substring(0, custName.indexOf("Customizer"));
                }

                custName = custName.substring(custName.lastIndexOf('.') + 1,              
                                            custName.length());
              }
              //custItem = new JMenuItem("Configure: "+ custName);
              custItem = new MenuItem("Configure: " + custName);
              if (tbi.getBean() instanceof BeanCommon) {
                customizationEnabled = 
                  !((BeanCommon)tbi.getBean()).isBusy();
              }
            }

            custItem.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                  if (bc instanceof MetaBean) {
                    popupCustomizer(custClass,
                      (JComponent) ((BeanInstance) tempAssociatedBeans.
                                    elementAt(tt)).getBean());
                  } else {
                    popupCustomizer(custClass, bc);
                  }

                  notifyIsDirty();
                }
              });
            custItem.setEnabled(customizationEnabled);
            beanContextMenu.add(custItem);
            menuItemCount++;
          } else {
            System.err.println("[KnowledgeFlow] No customizer class");
          }
        }

        Vector esdV = new Vector();

        //for (int i = 0; i < compInfoOutputs.size(); i++) {
        for (int i = 0; i < compInfo.size(); i++) {
          EventSetDescriptor[] temp = 
          //  ((BeanInfo) compInfoOutputs.elementAt(i)).getEventSetDescriptors();
          ((BeanInfo) compInfo.elementAt(i)).getEventSetDescriptors();

          if ((temp != null) && (temp.length > 0)) {
            esdV.add(temp);
          }
        }

        //        EventSetDescriptor [] esds = compInfo.getEventSetDescriptors();
        //        if (esds != null && esds.length > 0) {
        if (esdV.size() > 0) {
          //          beanContextMenu.insert(new JLabel("Connections", 
          //                                            SwingConstants.CENTER), 
          //                                 menuItemCount);
          MenuItem connections = new MenuItem("Connections:");
          connections.setEnabled(false);
          beanContextMenu.insert(connections, menuItemCount);
          menuItemCount++;
        }

        //final Vector finalOutputs = outputBeans;
        final Vector finalOutputs = associatedBeans;

        for (int j = 0; j < esdV.size(); j++) {
          final int fj = j;
          String sourceBeanName = "";

          if (bc instanceof MetaBean) {
            //Object sourceBean = ((BeanInstance) outputBeans.elementAt(j)).getBean();
            Object sourceBean = ((BeanInstance) associatedBeans.elementAt(j)).getBean();
            if (sourceBean instanceof BeanCommon) {
              sourceBeanName = ((BeanCommon)sourceBean).getCustomName();
            } else {
              if (sourceBean instanceof WekaWrapper) {
                sourceBeanName = ((WekaWrapper) sourceBean).getWrappedAlgorithm()
                .getClass().getName();
              } else {
                sourceBeanName = sourceBean.getClass().getName();
              }

              sourceBeanName = 
                sourceBeanName.substring(sourceBeanName.lastIndexOf('.') + 1, 
                    sourceBeanName.length());
            }
            sourceBeanName += ": ";
          }

          EventSetDescriptor[] esds = (EventSetDescriptor[]) esdV.elementAt(j);

          for (int i = 0; i < esds.length; i++) {
            //	  System.err.println(esds[i].getName());
            // add each event name to the menu
            //            JMenuItem evntItem = new JMenuItem(sourceBeanName
            //                                               +esds[i].getName());
            MenuItem evntItem = new MenuItem(sourceBeanName +
                                             esds[i].getName());
            final EventSetDescriptor esd = esds[i];

            // Check EventConstraints (if any) here
            boolean ok = true;

            if (bc instanceof EventConstraints) {
              ok = ((EventConstraints) bc).eventGeneratable(esd.getName());
            }

            if (ok) {
              evntItem.addActionListener(new ActionListener() {
                  public void actionPerformed(ActionEvent e) {
                    connectComponents(esd,
                                      (bc instanceof MetaBean)
                                      ? ((BeanInstance) finalOutputs.elementAt(fj)) : bi, xx, yy);
                    notifyIsDirty();
                  }
                });
            } else {
              evntItem.setEnabled(false);
            }

            beanContextMenu.add(evntItem);
            menuItemCount++;
          }
        }
      }
    } catch (IntrospectionException ie) {
      ie.printStackTrace();
    }

    //    System.err.println("Just before look for other options");
    // now look for other options for this bean
    if (bc instanceof UserRequestAcceptor || bc instanceof Startable) {
      Enumeration req = null;
      
      if (bc instanceof UserRequestAcceptor) {
        req = ((UserRequestAcceptor) bc).enumerateRequests();
      }

      if ((bc instanceof Startable) || (req !=null && req.hasMoreElements())) {
        //	beanContextMenu.insert(new JLabel("Actions", 
        //					  SwingConstants.CENTER), 
        //			       menuItemCount);
        MenuItem actions = new MenuItem("Actions:");
        actions.setEnabled(false);
        beanContextMenu.insert(actions, menuItemCount);
        menuItemCount++;
      }

      if (bc instanceof Startable) {
        String tempS = ((Startable)bc).getStartMessage();
        insertUserOrStartableMenuItem(bc, true, tempS, beanContextMenu);
      }
      
      while (req != null && req.hasMoreElements()) {
        String tempS = (String) req.nextElement();
        insertUserOrStartableMenuItem(bc, false, tempS, beanContextMenu);
        menuItemCount++;
      }
    }

    //    System.err.println("Just before showing menu");
    // popup the menu
    if (menuItemCount > 0) {
      //beanContextMenu.show(m_beanLayout, x, y);
      m_beanLayout.add(beanContextMenu);
      beanContextMenu.show(m_beanLayout, x, y);
    }
  }
  
  private void insertUserOrStartableMenuItem(final JComponent bc, 
      final boolean startable, String tempS, PopupMenu beanContextMenu) {

    boolean disabled = false;
    boolean confirmRequest = false;

    // check to see if this item is currently disabled
    if (tempS.charAt(0) == '$') {
      tempS = tempS.substring(1, tempS.length());
      disabled = true;
    }
    
    // check to see if this item requires confirmation
    if (tempS.charAt(0) == '?') {
      tempS = tempS.substring(1, tempS.length());
      confirmRequest = true;
    }

    final String tempS2 = tempS;

    //      JMenuItem custItem = new JMenuItem(tempS2);
    MenuItem custItem = new MenuItem(tempS2);
    if (confirmRequest) {
      custItem.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          // 
          int result = JOptionPane.showConfirmDialog(KnowledgeFlowApp.this,
              tempS2,
              "Confirm action",
              JOptionPane.YES_NO_OPTION);
          if (result == JOptionPane.YES_OPTION) {
            Thread startPointThread = new Thread() {
              public void run() {
                try {
                  if (startable) {
                    ((Startable)bc).start();                    
                  } else if (bc instanceof UserRequestAcceptor) {
                    ((UserRequestAcceptor) bc).performRequest(tempS2);
                  }
                  notifyIsDirty();
                } catch (Exception ex) {
                  ex.printStackTrace();
                }
              }
            };
            startPointThread.setPriority(Thread.MIN_PRIORITY);
            startPointThread.start();
          }
        }
      });
    } else {
      custItem.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          Thread startPointThread = new Thread() {
            public void run() {
              try {
                if (startable) {
                  ((Startable)bc).start();                  
                } else if (bc instanceof UserRequestAcceptor) {
                  ((UserRequestAcceptor) bc).performRequest(tempS2);
                }
                notifyIsDirty();
              } catch (Exception ex) {
                ex.printStackTrace();
              }
            }
          };
          startPointThread.setPriority(Thread.MIN_PRIORITY);
          startPointThread.start();
        }
      });
    }

    if (disabled) {
      custItem.setEnabled(false);
    }

    beanContextMenu.add(custItem); 
  }

  /**
   * Popup the customizer for this bean
   *
   * @param custClass the class of the customizer
   * @param bc the bean to be customized
   */
  private void popupCustomizer(Class custClass, JComponent bc) {
    try {
      // instantiate
      final Object customizer = custClass.newInstance();
      // set environment **before** setting object!!
      if (customizer instanceof EnvironmentHandler) {
        ((EnvironmentHandler)customizer).setEnvironment(m_flowEnvironment);
      }
      ((Customizer)customizer).setObject(bc);
      final javax.swing.JFrame jf = new javax.swing.JFrame();
      jf.getContentPane().setLayout(new BorderLayout());
      jf.getContentPane().add((JComponent)customizer, BorderLayout.CENTER);
      if (customizer instanceof CustomizerCloseRequester) {
	((CustomizerCloseRequester)customizer).setParentFrame(jf);
      }
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	  public void windowClosing(java.awt.event.WindowEvent e) {
	    if (customizer instanceof CustomizerClosingListener) {
	      ((CustomizerClosingListener)customizer).customizerClosing();
	    }
	    jf.dispose();
	  }
	});
      jf.pack();
      jf.setVisible(true);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

                                
  /**
   * Handles adding a custom MetaBean to the user toolbar
   *
   * @param bean the MetaBean
   * @param installListener install a listener for window close
   * events so as to save the user components
   */
  private void addToUserToolBar(MetaBean bean, 
                                boolean installListener) {

    if (m_userToolBar == null) {
      // need to create the user tab and toolbar
      setUpUserToolBar();
    }

    // Disconnect any beans connected to the inputs or outputs
    // of this MetaBean (prevents serialization of the entire
    // KnowledgeFlow!!)
    Vector tempRemovedConnections = new Vector();
    Vector allConnections = BeanConnection.getConnections();
    Vector inputs = bean.getInputs();
    Vector outputs = bean.getOutputs();
    Vector allComps = bean.getSubFlow();
        
    for (int i = 0; i < inputs.size(); i++) {
      BeanInstance temp = (BeanInstance)inputs.elementAt(i);
      // is this input a target for some event?
      for (int j = 0; j < allConnections.size(); j++) {
        BeanConnection tempC = (BeanConnection)allConnections.elementAt(j);
        if (tempC.getTarget() == temp) {
          tempRemovedConnections.add(tempC);
        }
        
        // also check to see if this input is a source for
        // some target that is *not* in the subFlow
        if (tempC.getSource() == temp && !bean.subFlowContains(tempC.getTarget())) {
          tempRemovedConnections.add(tempC);
        }
      }
    }

    for (int i = 0; i < outputs.size(); i++) {
      BeanInstance temp = (BeanInstance)outputs.elementAt(i);
      // is this output a source for some target?
      for (int j = 0; j < allConnections.size(); j++) {
        BeanConnection tempC = (BeanConnection)allConnections.elementAt(j);
        if (tempC.getSource() == temp) {
          tempRemovedConnections.add(tempC);
        }
      }
    }
    
    
    for (int i = 0; i < tempRemovedConnections.size(); i++) {
      BeanConnection temp = 
        (BeanConnection)tempRemovedConnections.elementAt(i);
      temp.remove();
    }
    
    // now add to user tool bar
    JPanel tempUser = instantiateToolBarMetaBean(bean);
    m_userBoxPanel.add(tempUser);
    if (installListener && m_firstUserComponentOpp) {
      try {
        installWindowListenerForSavingUserBeans();
        m_firstUserComponentOpp = false;
      } catch (Exception ex) {
        ex.printStackTrace();
      }
    }

    // Now reinstate any deleted connections to the original MetaBean
    for (int i = 0; i < tempRemovedConnections.size(); i++) {
      BeanConnection temp = 
        (BeanConnection)tempRemovedConnections.elementAt(i);
      BeanConnection newC = 
        new BeanConnection(temp.getSource(), temp.getTarget(),
                           temp.getSourceEventSetDescriptor());
    }    
  }

  /**
   * Popup a menu giving choices for connections to delete (if any)
   *
   * @param closestConnections a vector containing 0 or more BeanConnections
   * @param x the x coordinate at which to popup the menu
   * @param y the y coordinate at which to popup the menu
   *
   * Modified by Zerbetto: javax.swing.JPopupMenu transformed into java.awt.PopupMenu
   */
  private void deleteConnectionPopup(Vector closestConnections, int x, int y) {
    if (closestConnections.size() > 0) {
      int menuItemCount = 0;

      // modifications by Zerbetto
      //JPopupMenu deleteConnectionMenu = new JPopupMenu();
      PopupMenu deleteConnectionMenu = new PopupMenu();

      //      deleteConnectionMenu.insert(new JLabel("Delete Connection", 
      //					     SwingConstants.CENTER), 
      //				  menuItemCount);
      MenuItem deleteConnection = new MenuItem("Delete Connection:");
      deleteConnection.setEnabled(false);
      deleteConnectionMenu.insert(deleteConnection, menuItemCount);
      menuItemCount++;

      for (int i = 0; i < closestConnections.size(); i++) {
        final BeanConnection bc = (BeanConnection) closestConnections.elementAt(i);
        String connName = bc.getSourceEventSetDescriptor().getName();

        //JMenuItem deleteItem = new JMenuItem(connName);
        String targetName = "";
        if (bc.getTarget().getBean() instanceof BeanCommon) {
          targetName = ((BeanCommon)bc.getTarget().getBean()).getCustomName();
        } else {
          targetName = bc.getTarget().getBean().getClass().getName();
          targetName = targetName.substring(targetName.lastIndexOf('.')+1, targetName.length());
        }
        MenuItem deleteItem = new MenuItem(connName + "-->" + targetName);
        deleteItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
              bc.remove();
              m_beanLayout.revalidate();
              m_beanLayout.repaint();
              notifyIsDirty();
            }
          });
        deleteConnectionMenu.add(deleteItem);
        menuItemCount++;
      }

      //deleteConnectionMenu.show(m_beanLayout, x, y);
      m_beanLayout.add(deleteConnectionMenu);
      deleteConnectionMenu.show(m_beanLayout, x, y);
    }
  }

  /**
   * Initiates the connection process for two beans
   *
   * @param esd the EventSetDescriptor for the source bean
   * @param bi the source bean
   * @param x the x coordinate to start connecting from
   * @param y the y coordinate to start connecting from
   */
  private void connectComponents(EventSetDescriptor esd, 
				 BeanInstance bi,
				 int x,
				 int y) {
    // record the event set descriptior for this event
    m_sourceEventSetDescriptor = esd;

    Class listenerClass = esd.getListenerType(); // class of the listener
    JComponent source = (JComponent)bi.getBean();
    // now determine which (if any) of the other beans implement this
    // listener
    int targetCount = 0;
    Vector beanInstances = BeanInstance.getBeanInstances();
    for (int i = 0; i < beanInstances.size(); i++) {
      JComponent bean = 
	(JComponent)((BeanInstance)beanInstances.elementAt(i)).getBean();
      boolean connectable = false;
      boolean canContinue = false;
      if (bean != source) {
        if (bean instanceof MetaBean) {
          if (((MetaBean)bean).canAcceptConnection(listenerClass)) {
            canContinue = true;
          }
        } else if (listenerClass.isInstance(bean) && bean != source) {
          canContinue = true;
        }
      }
      if (canContinue) {
	if (!(bean instanceof BeanCommon)) {
	  connectable = true; // assume this bean is happy to receive a connection
	} else {
	  // give this bean a chance to veto any proposed connection via
	  // the listener interface
	  if (((BeanCommon)bean).
	      //connectionAllowed(esd.getName())) {
	      connectionAllowed(esd)) {
	    connectable = true;
	  }
	}
	if (connectable) {
	  if (bean instanceof Visible) {
	    targetCount++;
	    ((Visible)bean).getVisual().setDisplayConnectors(true);
	  }
	}
      }
    }
    
    // have some possible beans to connect to?
    if (targetCount > 0) {
      //      System.err.println("target count "+targetCount);
      if (source instanceof Visible) {
	((Visible)source).getVisual().setDisplayConnectors(true);
      }

      m_editElement = bi;
      Point closest = ((Visible)source).getVisual().
	getClosestConnectorPoint(new Point(x, y));

      m_startX = (int)closest.getX();
      m_startY = (int)closest.getY();
      m_oldX = m_startX;
      m_oldY = m_startY;

      Graphics2D gx = (Graphics2D)m_beanLayout.getGraphics();
      gx.setXORMode(java.awt.Color.white);
      gx.drawLine(m_startX, m_startY, m_startX, m_startY);
      gx.dispose();
      m_mode = CONNECTING;
    }
  }

  private void addComponent(BeanInstance comp, boolean repaint) {
    if (comp.getBean() instanceof Visible) {
      ((Visible)comp.getBean()).getVisual().addPropertyChangeListener(this);
    }
    if (comp.getBean() instanceof BeanCommon) {
      ((BeanCommon)comp.getBean()).setLog(m_logPanel);
    }
    if (comp.getBean() instanceof MetaBean) {
      // re-align sub-beans
      Vector list;
      
      list = ((MetaBean) comp.getBean()).getInputs();
      for (int i = 0; i < list.size(); i++) {
        ((BeanInstance) list.get(i)).setX(comp.getX());
        ((BeanInstance) list.get(i)).setY(comp.getY());
      }

      list = ((MetaBean) comp.getBean()).getOutputs();
      for (int i = 0; i < list.size(); i++) {
        ((BeanInstance) list.get(i)).setX(comp.getX());
        ((BeanInstance) list.get(i)).setY(comp.getY());
      }
    }
    setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
    if (repaint) {
      m_beanLayout.repaint();
    }
    m_pointerB.setSelected(true);
    m_mode = NONE;
  }

  private void addComponent(int x, int y) {
    if (m_toolBarBean instanceof MetaBean) {
      // need to add the MetaBean's internal connections
      // to BeanConnection's vector
      Vector associatedConnections = 
        ((MetaBean)m_toolBarBean).getAssociatedConnections();
      BeanConnection.getConnections().addAll(associatedConnections);
    }

    if (m_toolBarBean instanceof BeanContextChild) {
      m_bcSupport.add(m_toolBarBean);
    }
    BeanInstance bi = new BeanInstance(m_beanLayout, m_toolBarBean, x, y);
    //    addBean((JComponent)bi.getBean());
    m_toolBarBean = null;
    addComponent(bi, true);
  }

  /**
   * Handles the checking of a selected set of components
   * for suitability for grouping. If suitable the user
   * is prompted for a name and then a MetaBean is used
   * group the components.
   */
  private void checkSubFlow(int startX, int startY,
                            int endX, int endY) {

    java.awt.Rectangle r = 
      new java.awt.Rectangle((startX < endX) ? startX : endX,
                             (startY < endY) ? startY: endY,
                             Math.abs(startX - endX),
                             Math.abs(startY - endY));
    //    System.err.println(r);
    Vector selected = BeanInstance.findInstances(r);
    //    System.err.println(r);
    // check if sub flow is valid
    Vector inputs = BeanConnection.inputs(selected);
    Vector outputs = BeanConnection.outputs(selected);
    
    // screen the inputs and outputs
    if (inputs.size() == 0 || outputs.size() == 0) {
      return;
    }

    // dissallow MetaBeans in the selected set (for the
    // time being).
    for (int i = 0; i < selected.size(); i++) {
      BeanInstance temp = (BeanInstance)selected.elementAt(i);
      if (temp.getBean() instanceof MetaBean) {
        return;
      }
    }

    // show connector dots for selected beans
    for (int i = 0; i < selected.size(); i++) {
      BeanInstance temp = (BeanInstance)selected.elementAt(i);
      if (temp.getBean() instanceof Visible) {
        ((Visible)temp.getBean()).getVisual().setDisplayConnectors(true);
      }
    }

    // show connector dots for input beans
    for (int i = 0; i < inputs.size(); i++) {
      BeanInstance temp = (BeanInstance)inputs.elementAt(i);
      if (temp.getBean() instanceof Visible) {
        ((Visible)temp.getBean()).getVisual().
          setDisplayConnectors(true, java.awt.Color.red);
      }
    }

    // show connector dots for output beans
    for (int i = 0; i < outputs.size(); i++) {
      BeanInstance temp = (BeanInstance)outputs.elementAt(i);
      if (temp.getBean() instanceof Visible) {
        ((Visible)temp.getBean()).getVisual().
          setDisplayConnectors(true, java.awt.Color.green);
      }
    }
    
    BufferedImage subFlowPreview = null; 
    try {
      	subFlowPreview = createImage(m_beanLayout, r);              
    } catch (IOException ex) {
      ex.printStackTrace();
      // drop through quietly
    }

    // Confirmation pop-up
    int result = JOptionPane.showConfirmDialog(KnowledgeFlowApp.this,
                                               "Group this sub-flow?",
                                               "Group Components",
                                               JOptionPane.YES_NO_OPTION);
    if (result == JOptionPane.YES_OPTION) {
      Vector associatedConnections = 
        BeanConnection.associatedConnections(selected);

      String name = JOptionPane.showInputDialog(KnowledgeFlowApp.this,
                                                "Enter a name for this group",
                                                "MyGroup");
      if (name != null) {       
        MetaBean group = new MetaBean();
        group.setSubFlow(selected);
        group.setAssociatedConnections(associatedConnections);
        group.setInputs(inputs);
        group.setOutputs(outputs);
        group.setSubFlowPreview(new ImageIcon(subFlowPreview));
        if (name.length() > 0) {
          //          group.getVisual().setText(name);
          group.setCustomName(name);
        }
        
        if (group instanceof BeanContextChild) {
          m_bcSupport.add(group);
        }
        BeanInstance bi = new BeanInstance(m_beanLayout, group, 
                                           (int)r.getX()+(int)(r.getWidth()/2),
                                           (int)r.getY()+(int)(r.getHeight()/2));
        for (int i = 0; i < selected.size(); i++) {
          BeanInstance temp = (BeanInstance)selected.elementAt(i);
          temp.removeBean(m_beanLayout);
          if (temp.getBean() instanceof Visible) {
            ((Visible)temp.getBean()).getVisual().removePropertyChangeListener(this);
          }
        }
        for (int i = 0; i < associatedConnections.size(); i++) {
          BeanConnection temp = (BeanConnection)associatedConnections.elementAt(i);
          temp.setHidden(true);
        }
        group.shiftBeans(bi, true);
        
        addComponent(bi, true);
      }
    }

    // hide connector dots
    for (int i = 0; i < selected.size(); i++) {
      BeanInstance temp = (BeanInstance)selected.elementAt(i);
      if (temp.getBean() instanceof Visible) {
        ((Visible)temp.getBean()).getVisual().setDisplayConnectors(false);
      }
    }    
  }

  /**
   * Accept property change events
   *
   * @param e a <code>PropertyChangeEvent</code> value
   */
  public void propertyChange(PropertyChangeEvent e) {
    revalidate();
    m_beanLayout.repaint();
  }
  
  /**
   * Load a pre-saved layout
   */
  private void loadLayout() {
    m_loadB.setEnabled(false);
    m_saveB.setEnabled(false);
    int returnVal = m_FileChooser.showOpenDialog(this);
    if (returnVal == JFileChooser.APPROVE_OPTION) {
      stopFlow();

      // determine filename
      File oFile = m_FileChooser.getSelectedFile();
      // set internal flow directory environment variable
      m_flowEnvironment.addVariable("Internal.knowledgeflow.directory", oFile.getParent());

      // add extension if necessary
      if (m_FileChooser.getFileFilter() == m_KfFilter) {
        if (!oFile.getName().toLowerCase().endsWith(FILE_EXTENSION)) {
          oFile = new File(oFile.getParent(), 
                           oFile.getName() + FILE_EXTENSION);
        }
      } else if (m_FileChooser.getFileFilter() == m_KOMLFilter) {
        if (!oFile.getName().toLowerCase().endsWith(KOML.FILE_EXTENSION + "kf")) {
          oFile = new File(oFile.getParent(), 
                           oFile.getName() + KOML.FILE_EXTENSION + "kf");
        }
      } else if (m_FileChooser.getFileFilter() == m_XMLFilter) {
        if (!oFile.getName().toLowerCase().endsWith(FILE_EXTENSION_XML)) {
          oFile = new File(oFile.getParent(), 
                           oFile.getName() + FILE_EXTENSION_XML);
        }
      } else if (m_FileChooser.getFileFilter() == m_XStreamFilter) {
        if (!oFile.getName().toLowerCase().endsWith(XStream.FILE_EXTENSION +"kf")) {
          oFile = new File(oFile.getParent(), 
                           oFile.getName() + XStream.FILE_EXTENSION + "kf");
        }
      }
    
      try {
        Vector beans       = new Vector();
        Vector connections = new Vector();

        // KOML?
        if ( (KOML.isPresent()) && 
             (oFile.getAbsolutePath().toLowerCase().
              endsWith(KOML.FILE_EXTENSION + "kf")) ) {
          Vector v     = (Vector) KOML.read(oFile.getAbsolutePath());
          beans        = (Vector) v.get(XMLBeans.INDEX_BEANINSTANCES);
          connections  = (Vector) v.get(XMLBeans.INDEX_BEANCONNECTIONS);
        } /* XStream */ else if ( (XStream.isPresent()) && 
             (oFile.getAbsolutePath().toLowerCase().
              endsWith(XStream.FILE_EXTENSION + "kf")) ) {
          Vector v     = (Vector) XStream.read(oFile.getAbsolutePath());
          beans        = (Vector) v.get(XMLBeans.INDEX_BEANINSTANCES);
          connections  = (Vector) v.get(XMLBeans.INDEX_BEANCONNECTIONS);
        } /* XML? */ else if (oFile.getAbsolutePath().toLowerCase().
                              endsWith(FILE_EXTENSION_XML)) {
          XMLBeans xml = new XMLBeans(m_beanLayout, m_bcSupport); 
          Vector v     = (Vector) xml.read(oFile);
          beans        = (Vector) v.get(XMLBeans.INDEX_BEANINSTANCES);
          connections  = (Vector) v.get(XMLBeans.INDEX_BEANCONNECTIONS);
          //connections  = new Vector();
        } /* binary */ else {
          InputStream is = new FileInputStream(oFile);
          ObjectInputStream ois = new ObjectInputStream(is);
          beans = (Vector) ois.readObject();
          connections = (Vector) ois.readObject();
          ois.close();
        }

        integrateFlow(beans, connections);
        setEnvironment();
        m_logPanel.clearStatus();
        m_logPanel.statusMessage("[KnowledgeFlow]|Flow loaded.");
      } catch (Exception ex) {
        m_logPanel.statusMessage("[KnowledgeFlow]|Unable to load flow (see log).");
        m_logPanel.logMessage("[KnowledgeFlow] Unable to load flow ("
            + ex.getMessage() + ").");
	ex.printStackTrace();
      }
    }
    m_loadB.setEnabled(true);
    m_saveB.setEnabled(true);
  }

  // Link the supplied beans into the KnowledgeFlow gui
  private void integrateFlow(Vector beans, Vector connections) {
    java.awt.Color bckC = getBackground();
    m_bcSupport = new BeanContextSupport();
    m_bcSupport.setDesignTime(true);

    // register this panel as a property change listener with each
    // bean
    for (int i = 0; i < beans.size(); i++) {
      BeanInstance tempB = (BeanInstance)beans.elementAt(i);
      if (tempB.getBean() instanceof Visible) {
        ((Visible)(tempB.getBean())).getVisual().
          addPropertyChangeListener(this);

        // A workaround to account for JPanel's with their default
        // background colour not being serializable in Apple's JRE
        ((Visible)(tempB.getBean())).getVisual().
          setBackground(bckC);
        ((JComponent)(tempB.getBean())).setBackground(bckC);
      }
      if (tempB.getBean() instanceof BeanCommon) {
        ((BeanCommon)(tempB.getBean())).setLog(m_logPanel);
      }
      if (tempB.getBean() instanceof BeanContextChild) {
        m_bcSupport.add(tempB.getBean());
      }
    }
    BeanInstance.setBeanInstances(beans, m_beanLayout);
    BeanConnection.setConnections(connections);
    m_beanLayout.revalidate();
    m_beanLayout.repaint();
  }

  /**
   * Set the flow for the KnowledgeFlow to edit. Assumes that client
   * has loaded a Vector of beans and a Vector of connections. the supplied
   * beans and connections are deep-copied via serialization before being
   * set in the layout.
   *
   * @param v a Vector containing a Vector of beans and a Vector of connections
   * @exception Exception if something goes wrong
   */
  public void setFlow(Vector v) throws Exception {
    //    Vector beansCopy = null, connectionsCopy = null;
    clearLayout();
    SerializedObject so = new SerializedObject(v);
    Vector copy = (Vector)so.getObject();
    
    Vector beans = (Vector)copy.elementAt(0);
    Vector connections = (Vector)copy.elementAt(1);
    
    // reset environment variables
    m_flowEnvironment = new Environment();
    integrateFlow(beans, connections);
  }

  /**
   * Gets the current flow being edited. The flow is returned as a single
   * Vector containing two other Vectors: the beans and the connections.
   * These two vectors are deep-copied via serialization before being
   * returned.
   *
   * @return the current flow being edited
   */
  public Vector getFlow() throws Exception {
    Vector v = new Vector();
    Vector beans = BeanInstance.getBeanInstances();
    Vector connections = BeanConnection.getConnections();
    detachFromLayout(beans);
    v.add(beans);
    v.add(connections);

    SerializedObject so = new SerializedObject(v);
    Vector copy = (Vector)so.getObject();

    //    tempWrite(beans, connections);
    
    integrateFlow(beans, connections);
    return copy;
  }
  
  /**
   * Utility method to create an image of a region of the given component
   * @param component the component to create an image of
   * @param region the region of the component to put into the image
   * @return the image
   * @throws IOException
   */
  protected static BufferedImage createImage(JComponent component, Rectangle region)
  throws IOException {
    boolean opaqueValue = component.isOpaque();
    component.setOpaque( true );
    BufferedImage image = new BufferedImage(region.width, 
	region.height, BufferedImage.TYPE_INT_RGB);
    Graphics2D g2d = image.createGraphics();
    g2d.translate(-region.getX(), -region.getY());
    //g2d.setClip( region );
    component.paint( g2d );
    g2d.dispose();
    component.setOpaque( opaqueValue );
    
    return image;
  }

  // Remove this panel as a property changle listener from
  // each bean
  private void detachFromLayout(Vector beans) {
    for (int i = 0; i < beans.size(); i++) {
      BeanInstance tempB = (BeanInstance)beans.elementAt(i);
      if (tempB.getBean() instanceof Visible) {
        ((Visible)(tempB.getBean())).getVisual().
          removePropertyChangeListener(this);
          
        if (tempB.getBean() instanceof MetaBean) {
          ((MetaBean)tempB.getBean()).
            removePropertyChangeListenersSubFlow(this);
        }

        // A workaround to account for JPanel's with their default
        // background colour not being serializable in Apple's JRE.
        // JComponents are rendered with a funky stripy background
        // under OS X using java.awt.TexturePaint - unfortunately
        // TexturePaint doesn't implement Serializable.
        ((Visible)(tempB.getBean())).getVisual().
          setBackground(java.awt.Color.white);
        ((JComponent)(tempB.getBean())).setBackground(java.awt.Color.white);
      }
    }
  }

  /**
   * Serialize the layout to a file
   */
  private void saveLayout() {
    //    m_loadB.setEnabled(false);
    //    m_saveB.setEnabled(false);
    int returnVal = m_FileChooser.showSaveDialog(this);
    java.awt.Color bckC = getBackground();
    if (returnVal == JFileChooser.APPROVE_OPTION) {
      // temporarily remove this panel as a property changle listener from
      // each bean

      Vector beans = BeanInstance.getBeanInstances();
      detachFromLayout(beans);

      // determine filename
      File sFile = m_FileChooser.getSelectedFile();

      // add extension if necessary
      if (m_FileChooser.getFileFilter() == m_KfFilter) {
        if (!sFile.getName().toLowerCase().endsWith(FILE_EXTENSION)) {
          sFile = new File(sFile.getParent(), 
                           sFile.getName() + FILE_EXTENSION);
        }
      } else if (m_FileChooser.getFileFilter() == m_KOMLFilter) {
        if (!sFile.getName().toLowerCase().endsWith(KOML.FILE_EXTENSION + "kf")) {
          sFile = new File(sFile.getParent(), 
                           sFile.getName() + KOML.FILE_EXTENSION + "kf");
        }
      } else if (m_FileChooser.getFileFilter() == m_XStreamFilter) {
        if (!sFile.getName().toLowerCase().endsWith(XStream.FILE_EXTENSION + "kf")) {
          sFile = new File(sFile.getParent(), 
                           sFile.getName() + XStream.FILE_EXTENSION + "kf");
        }
      } else if (m_FileChooser.getFileFilter() == m_XMLFilter) {
        if (!sFile.getName().toLowerCase().endsWith(FILE_EXTENSION_XML)) {
          sFile = new File(sFile.getParent(), 
                           sFile.getName() + FILE_EXTENSION_XML);
        }
      }
    
      // now serialize components vector and connections vector
      try {
        // KOML?
        if ((KOML.isPresent()) && 
            (sFile.getAbsolutePath().toLowerCase().
             endsWith(KOML.FILE_EXTENSION + "kf")) ) {
          Vector v = new Vector();
          v.setSize(2);
          v.set(XMLBeans.INDEX_BEANINSTANCES, beans);
          v.set(XMLBeans.INDEX_BEANCONNECTIONS, BeanConnection.getConnections());
          KOML.write(sFile.getAbsolutePath(), v);
        } /* XStream */ else if ((XStream.isPresent()) && 
            (sFile.getAbsolutePath().toLowerCase().
             endsWith(XStream.FILE_EXTENSION + "kf")) ) {
          Vector v = new Vector();
          v.setSize(2);
          v.set(XMLBeans.INDEX_BEANINSTANCES, beans);
          v.set(XMLBeans.INDEX_BEANCONNECTIONS, BeanConnection.getConnections());
          XStream.write(sFile.getAbsolutePath(), v);
        } /* XML? */ else if (sFile.getAbsolutePath().
                              toLowerCase().endsWith(FILE_EXTENSION_XML)) {
          Vector v = new Vector();
          v.setSize(2);
          v.set(XMLBeans.INDEX_BEANINSTANCES, beans);
          v.set(XMLBeans.INDEX_BEANCONNECTIONS, BeanConnection.getConnections());
          XMLBeans xml = new XMLBeans(m_beanLayout, m_bcSupport); 
          xml.write(sFile, v);
        } /* binary */ else {
          OutputStream os = new FileOutputStream(sFile);
          ObjectOutputStream oos = new ObjectOutputStream(os);
          oos.writeObject(beans);
          oos.writeObject(BeanConnection.getConnections());
          oos.flush();
          oos.close();
        }
        m_logPanel.statusMessage("[KnowledgeFlow]|Flow saved.");
        
        // set the internal knowledgeflow directory environment var for this flow
        m_flowEnvironment.addVariable("Internal.knowledgeflow.directory", sFile.getParent());
        setEnvironment();
      } catch (Exception ex) {
        m_logPanel.statusMessage("[KnowledgeFlow]|Unable to save flow (see log).");
        m_logPanel.logMessage("[KnowledgeFlow] Unable to save flow ("
            + ex.getMessage() + ").");
	ex.printStackTrace();
      } finally {
	// restore this panel as a property change listener in the beans
	for (int i = 0; i < beans.size(); i++) {
	  BeanInstance tempB = (BeanInstance)beans.elementAt(i);
	  if (tempB.getBean() instanceof Visible) {
	    ((Visible)(tempB.getBean())).getVisual().
	      addPropertyChangeListener(this);

            if (tempB.getBean() instanceof MetaBean) {
              ((MetaBean)tempB.getBean()).
                addPropertyChangeListenersSubFlow(this);
            }
	    // Restore the default background colour
	    ((Visible)(tempB.getBean())).getVisual().
	      setBackground(bckC);
	    ((JComponent)(tempB.getBean())).setBackground(bckC);
	  }
	}
      }
    }
    //    m_saveB.setEnabled(true);
    //    m_loadB.setEnabled(true);
  }

  /**
   * Save the knowledge flow into the OutputStream passed at input. Only
   * supports saving the layout data (no trained models) to XML.
   *
   * @param out		the output stream to save the layout in
   */
  public void saveLayout(OutputStream out) {
    // temporarily remove this panel as a property changle listener from
    // each bean
    Vector beans = BeanInstance.getBeanInstances();

    for (int i = 0; i < beans.size(); i++) {
      BeanInstance tempB = (BeanInstance) beans.elementAt(i);

      if (tempB.getBean() instanceof Visible) {
        ((Visible) (tempB.getBean())).getVisual()
         .removePropertyChangeListener(this);

        if (tempB.getBean() instanceof MetaBean) {
          ((MetaBean) tempB.getBean()).removePropertyChangeListenersSubFlow(this);
        }
      }
    }

    // now serialize components vector and connections vector
    try {
      Vector v = new Vector();
      v.setSize(2);
      v.set(XMLBeans.INDEX_BEANINSTANCES, beans);
      v.set(XMLBeans.INDEX_BEANCONNECTIONS, BeanConnection.getConnections());

      XMLBeans xml = new XMLBeans(m_beanLayout, m_bcSupport);
      xml.write(out, v);
    } catch (Exception ex) {
      ex.printStackTrace();
    } finally {
      // restore this panel as a property change listener in the beans
      for (int i = 0; i < beans.size(); i++) {
        BeanInstance tempB = (BeanInstance) beans.elementAt(i);

        if (tempB.getBean() instanceof Visible) {
          ((Visible) (tempB.getBean())).getVisual()
           .addPropertyChangeListener(this);

          if (tempB.getBean() instanceof MetaBean) {
            ((MetaBean) tempB.getBean()).addPropertyChangeListenersSubFlow(this);
          }
        }
      }
    }
  }

  private void loadUserComponents() {
    Vector tempV = null;
    String ext = "";
    if (m_UserComponentsInXML)
      ext = USERCOMPONENTS_XML_EXTENSION;
    File sFile = 
      new File(System.getProperty("user.home")
               +File.separator + ".knowledgeFlow"
               +File.separator + "userComponents"
               +ext);
    if (sFile.exists()) {
      try {
        if (m_UserComponentsInXML) {
          XMLBeans xml = new XMLBeans(m_beanLayout, m_bcSupport, XMLBeans.DATATYPE_USERCOMPONENTS);
          tempV = (Vector) xml.read(sFile);
        }
        else {
          InputStream is = new FileInputStream(sFile);
          ObjectInputStream ois = new ObjectInputStream(is);
          tempV = (Vector)ois.readObject();
          ois.close();
        }
      } catch (Exception ex) {
        System.err.println("[KnowledgeFlow] Problem reading user components.");
        ex.printStackTrace();
        return;
      }
      if (tempV.size() > 0) {
        // create the user tab and add the components
        for (int i = 0; i < tempV.size(); i++) {
          MetaBean tempB = (MetaBean)tempV.elementAt(i);
          addToUserToolBar(tempB, false);
        }
      }
    }
  }

  private void installWindowListenerForSavingUserBeans() {
    ((java.awt.Window)getTopLevelAncestor()).
      addWindowListener(new java.awt.event.WindowAdapter() {
          public void windowClosing(java.awt.event.WindowEvent e) {
            System.out.println("[KnowledgeFlow] Saving user components....");
            File sFile = 
              new File(System.getProperty("user.home")
                       +File.separator+".knowledgeFlow");
            if (!sFile.exists()) {
              if (!sFile.mkdir()) {
                System.err.println("[KnowledgeFlow] Unable to create .knowledgeFlow "
                                   +"directory in your HOME.");
              } else {
                // make the plugins subdirectory for the user
                sFile = new File(sFile.toString() + File.separator 
                    + "plugins");
                sFile.mkdir();
              }
            }
            try {
              String ext = "";
              if (m_UserComponentsInXML)
                ext = USERCOMPONENTS_XML_EXTENSION;
              File sFile2 = new File(sFile.getAbsolutePath()
                                     +File.separator
                                     +"userComponents"
                                     +ext);
                
              if (m_UserComponentsInXML) {
                XMLBeans xml = new XMLBeans(m_beanLayout, m_bcSupport, XMLBeans.DATATYPE_USERCOMPONENTS);
                xml.write(sFile2, m_userComponents);
              }
              else {
                OutputStream os = new FileOutputStream(sFile2);
                ObjectOutputStream oos = new ObjectOutputStream(os);
                oos.writeObject(m_userComponents);
                oos.flush();
                oos.close();
              }
            } catch (Exception ex) {
              System.err.println("[KnowledgeFlow] Unable to save user components");
              ex.printStackTrace();
            } 

          }
        });
  }
  
  /**
   * Utility method for grabbing the global info help (if it exists) from
   * an arbitrary object
   *
   * @param tempBean the object to grab global info from 
   * @return the global help info or null if global info does not exist
   */
  public static String getGlobalInfo(Object tempBean) {
    // set tool tip text from global info if supplied
    String gi = null;
    try {
      BeanInfo bi = Introspector.getBeanInfo(tempBean.getClass());
      MethodDescriptor [] methods = bi.getMethodDescriptors();
      for (int i = 0; i < methods.length; i++) {
	String name = methods[i].getDisplayName();
	Method meth = methods[i].getMethod();
	if (name.equals("globalInfo")) {
	  if (meth.getReturnType().equals(String.class)) {
	    Object args[] = { };
	    String globalInfo = (String)(meth.invoke(tempBean, args));
	    gi = globalInfo;
	    break;
	  }
	}
      }
    } catch (Exception ex) {
      
    }
    return gi;
  }

  /** variable for the KnowLedgeFlow class which would be set to null by the 
      memory monitoring thread to free up some memory if we running out of 
      memory.
   */
  private static KnowledgeFlowApp m_knowledgeFlow;

  /** for monitoring the Memory consumption */
  private static Memory m_Memory = new Memory(true);

  // list of things to be notified when the startup process of
  // the KnowledgeFlow is complete
  public static Vector s_startupListeners = new Vector();

  // modifications by Zerbetto
  // If showFileMenu is true, the file menu (open file, new file, save file buttons) is showed
  private boolean m_showFileMenu = true;
  
  /**
   * Create the singleton instance of the KnowledgeFlow
   * @param args can contain a file argument for loading a flow layout 
   * (format: "file=[path to layout file]")
   * Modified by Zerbetto: you can specify the path of a knowledge flow layout file at input
   */
  public static void createSingleton(String[] args) {
    //modifications by Zerbetto 05-12-2007
    String fileName = null;
    boolean showFileMenu = true;

    if ((args != null) && (args.length > 0)) {
      for (int i = 0; i < args.length; i++) {
        String arg = args[i];

        if (arg.startsWith("file=")) {
          fileName = arg.substring("file=".length());
        } else if (arg.startsWith("showFileMenu=")) {
          showFileMenu = Boolean.parseBoolean(arg.substring(
                "showFileMenu=".length()));
        }
      }
    }

    if (m_knowledgeFlow == null) {
      m_knowledgeFlow = new KnowledgeFlowApp(showFileMenu);
    }

    // end modifications by Zerbetto

    // notify listeners (if any)
    for (int i = 0; i < s_startupListeners.size(); i++) {
      ((StartUpListener) s_startupListeners.elementAt(i)).startUpComplete();
    }

    //modifications by Zerbetto 05-12-2007
    if (fileName != null) {
      m_knowledgeFlow.loadInitialLayout(fileName);
    }

    // end modifications 
  }
  
  public static void disposeSingleton() {
    m_knowledgeFlow = null;
  }

  /**
   * Return the singleton instance of the KnowledgeFlow
   *
   * @return the singleton instance
   */
  public static KnowledgeFlowApp getSingleton() {
    return m_knowledgeFlow;
  }

  /**
   * Add a listener to be notified when startup is complete
   * 
   * @param s a listener to add
   */
  public static void addStartupListener(StartUpListener s) {
    s_startupListeners.add(s);
  }

  /**
   * Loads the specified file at input
   *
   * Added by Zerbetto
   */
  //modifications by Zerbetto 05-12-2007
  private void loadInitialLayout(String fileName) {
    File oFile = new File(fileName);

    if (oFile.exists() && oFile.isFile()) {
      m_FileChooser.setSelectedFile(oFile);

      int index = fileName.lastIndexOf('.');

      if (index != -1) {
        String extension = fileName.substring(index);

        if (FILE_EXTENSION_XML.equalsIgnoreCase(extension)) {
          m_FileChooser.setFileFilter(m_knowledgeFlow.m_XMLFilter);
        } else if (FILE_EXTENSION.equalsIgnoreCase(extension)) {
          m_FileChooser.setFileFilter(m_knowledgeFlow.m_KfFilter);
        }
      }
    } else {
      System.err.println("[KnowledgeFlow] File '" + fileName + "' does not exists.");
    }

    try {
      Vector beans = new Vector();
      Vector connections = new Vector();

      // KOML?
      if ((KOML.isPresent()) &&
            (oFile.getAbsolutePath().toLowerCase().endsWith(KOML.FILE_EXTENSION))) {
        Vector v = (Vector) KOML.read(oFile.getAbsolutePath());
        beans = (Vector) v.get(XMLBeans.INDEX_BEANINSTANCES);
        connections = (Vector) v.get(XMLBeans.INDEX_BEANCONNECTIONS);
      } /* XML? */ else if (oFile.getAbsolutePath().toLowerCase()
                                     .endsWith(FILE_EXTENSION_XML)) {
        XMLBeans xml = new XMLBeans(m_beanLayout, m_bcSupport);
        Vector v = (Vector) xml.read(oFile);
        beans = (Vector) v.get(XMLBeans.INDEX_BEANINSTANCES);
        connections = (Vector) v.get(XMLBeans.INDEX_BEANCONNECTIONS);

        //connections  = new Vector();
      } /* binary */ else {
        InputStream is = new FileInputStream(oFile);
        ObjectInputStream ois = new ObjectInputStream(is);
        beans = (Vector) ois.readObject();
        connections = (Vector) ois.readObject();
        ois.close();
      }

      java.awt.Color bckC = getBackground();
      m_bcSupport = new BeanContextSupport();
      m_bcSupport.setDesignTime(true);

      // register this panel as a property change listener with each
      // bean
      for (int i = 0; i < beans.size(); i++) {
        BeanInstance tempB = (BeanInstance) beans.elementAt(i);

        if (tempB.getBean() instanceof Visible) {
          ((Visible) (tempB.getBean())).getVisual()
           .addPropertyChangeListener(this);

          // A workaround to account for JPanel's with their default
          // background colour not being serializable in Apple's JRE
          ((Visible) (tempB.getBean())).getVisual().setBackground(bckC);
          ((JComponent) (tempB.getBean())).setBackground(bckC);
        }

        if (tempB.getBean() instanceof BeanCommon) {
          ((BeanCommon) (tempB.getBean())).setLog(m_logPanel);
        }

        if (tempB.getBean() instanceof BeanContextChild) {
          m_bcSupport.add(tempB.getBean());
        }
      }

      BeanInstance.setBeanInstances(beans, m_beanLayout);
      BeanConnection.setConnections(connections);
      m_beanLayout.revalidate();
      m_beanLayout.repaint();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  //end modifications

  /**
   * Notifies to the parent swt that the layout is dirty
   *
   * Added by Zerbetto
   */
  private void notifyIsDirty() {
    //this.firePropertyChange(new Integer(IEditorPart.PROP_DIRTY).toString(), null, null);
    this.firePropertyChange("PROP_DIRTY", null, null);
  }

  /**
   * Main method.
   *
   * @param args a <code>String[]</code> value
   */
  public static void main(String [] args) {

    LookAndFeel.setLookAndFeel();
    
    try {
      // uncomment to disable the memory management:
      //m_Memory.setEnabled(false);

      final javax.swing.JFrame jf = new javax.swing.JFrame();
      jf.getContentPane().setLayout(new java.awt.BorderLayout());
      //final KnowledgeFlowApp tm = new KnowledgeFlowApp();
      m_knowledgeFlow = new KnowledgeFlowApp(true);

      jf.getContentPane().add(m_knowledgeFlow, java.awt.BorderLayout.CENTER);
      jf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

      jf.setSize(1000,750);
      jf.setVisible(true);     

      
      Thread memMonitor = new Thread() {
        public void run() {
          while(true) {
            try {
              //System.out.println("Before sleeping");
              this.sleep(4000);
              
              System.gc();

              if (m_Memory.isOutOfMemory()) {
                // clean up
                jf.dispose();
                m_knowledgeFlow = null;
                System.gc();

                // stop threads
                m_Memory.stopThreads();
               
                // display error
                System.err.println("\n[KnowledgeFlow] displayed message:");
                m_Memory.showOutOfMemory();
                System.err.println("\nexiting");
                System.exit(-1);
              }

            } catch(InterruptedException ex) { ex.printStackTrace(); }
          }
        }
      };

      memMonitor.setPriority(Thread.NORM_PRIORITY);
      memMonitor.start();
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
