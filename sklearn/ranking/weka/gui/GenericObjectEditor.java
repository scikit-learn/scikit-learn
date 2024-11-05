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
 *    GenericObjectEditor.java
 *    Copyright (C) 2002-2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import weka.core.Capabilities;
import weka.core.CapabilitiesHandler;
import weka.core.ClassDiscovery;
import weka.core.OptionHandler;
import weka.core.WekaPackageManager;
import weka.core.SerializedObject;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.gui.CheckBoxList.CheckBoxListModel;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.beans.PropertyEditor;
import java.beans.PropertyEditorManager;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.lang.reflect.Array;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Properties;
import java.util.StringTokenizer;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import weka.core.Capabilities;
import weka.core.CapabilitiesHandler;
import weka.core.ClassDiscovery;
import weka.core.CustomDisplayStringProvider;
import weka.core.OptionHandler;
import weka.core.SerializedObject;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.gui.CheckBoxList.CheckBoxListModel;

/**
 * A PropertyEditor for objects. It can be used either in a static or a dynamic
 * way. <br>
 * <br>
 * In the <b>static</b> way (<code>USE_DYNAMIC</code> is <code>false</code>) the
 * objects have been defined as editable in the GenericObjectEditor
 * configuration file, which lists possible values that can be selected from,
 * and themselves configured. The configuration file is called
 * "GenericObjectEditor.props" and may live in either the location given by
 * "user.home" or the current directory (this last will take precedence), and a
 * default properties file is read from the Weka distribution. For speed, the
 * properties file is read only once when the class is first loaded -- this may
 * need to be changed if we ever end up running in a Java OS ;-). <br>
 * <br>
 * If it is used in a <b>dynamic</b> way (the <code>UseDynamic</code> property 
 * of the GenericPropertiesCreator props file is set to <code>true</code>) 
 * then the classes to list are discovered by the 
 * <code>GenericPropertiesCreator</code> class (it checks the complete classpath). 
 * 
 * @see GenericPropertiesCreator
 * @see GenericPropertiesCreator#useDynamic()
 * @see GenericPropertiesCreator#CREATOR_FILE
 * @see weka.core.ClassDiscovery
 * 
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Xin Xu (xx5@cs.waikato.ac.nz)
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6624 $
 */
public class GenericObjectEditor implements PropertyEditor, CustomPanelSupplier {
  
  /** The object being configured. */
  protected Object m_Object;
  
  /** Holds a copy of the current object that can be reverted to
      if the user decides to cancel. */
  protected Object m_Backup;
    
  /** Handles property change notification. */
  protected PropertyChangeSupport m_Support = new PropertyChangeSupport(this);
    
  /** The Class of objects being edited. */
  protected Class m_ClassType;
    
  /** The model containing the list of names to select from. */
  protected Hashtable m_ObjectNames;

  /** The GUI component for editing values, created when needed. */
  protected GOEPanel m_EditorComponent;
    
  /** True if the GUI component is needed. */
  protected boolean m_Enabled = true;
    
  /** The name of the properties file. */
  protected static String PROPERTY_FILE = "weka/gui/GenericObjectEditor.props";
    
  /** Contains the editor properties. */
  protected static Properties EDITOR_PROPERTIES;

  /** the properties files containing the class/editor mappings. */
  public static final String GUIEDITORS_PROPERTY_FILE = "weka/gui/GUIEditors.props";

  /** The tree node of the current object so we can re-select it for the user. */
  protected GOETreeNode m_treeNodeOfCurrentObject;

  /** The property panel created for the objects. */
  protected PropertyPanel m_ObjectPropertyPanel;
    
  /** whether the class can be changed. */
  protected boolean m_canChangeClassInDialog;
  
  /** whether the Weka Editors were already registered. */
  protected static boolean m_EditorsRegistered;

  /** for filtering the tree based on the Capabilities of the leaves. */
  protected Capabilities m_CapabilitiesFilter = null;
  
  public static void determineClasses() {
    try {
      // make sure we load all packages first!!!
      WekaPackageManager.loadPackages(false);

      GenericPropertiesCreator creator = new GenericPropertiesCreator();

      // dynamic approach?
      if (creator.useDynamic()) {
        try {
          creator.execute(false);
          EDITOR_PROPERTIES = creator.getOutputProperties();
        }
        catch (Exception e) {
          JOptionPane.showMessageDialog(
              null,
              "Could not determine the properties for the generic object\n"
              + "editor. This exception was produced:\n"
              + e.toString(),
              "GenericObjectEditor",
              JOptionPane.ERROR_MESSAGE);
        }
      }
      else {
        // Allow a properties file in the current directory to override
        try {
          EDITOR_PROPERTIES = Utils.readProperties(PROPERTY_FILE);
          java.util.Enumeration keys = 
            (java.util.Enumeration)EDITOR_PROPERTIES.propertyNames();
          if (!keys.hasMoreElements()) {
            throw new Exception("Failed to read a property file for the "
                +"generic object editor");
          }
        }
        catch (Exception ex) {
          JOptionPane.showMessageDialog(
              null,
              "Could not read a configuration file for the generic object\n"
              +"editor. An example file is included with the Weka distribution.\n"
              +"This file should be named \"" + PROPERTY_FILE + "\" and\n"
              +"should be placed either in your user home (which is set\n"
              + "to \"" + System.getProperties().getProperty("user.home") + "\")\n"
              + "or the directory that java was started from\n",
              "GenericObjectEditor",
              JOptionPane.ERROR_MESSAGE);
        }
      }
    }
    catch (Exception e) {
      JOptionPane.showMessageDialog(
          null,
          "Could not initialize the GenericPropertiesCreator. "
          + "This exception was produced:\n"
          + e.toString(),
          "GenericObjectEditor",
          JOptionPane.ERROR_MESSAGE);
    }
  }
  
  /** 
   * Loads the configuration property file (USE_DYNAMIC is FALSE) or determines
   * the classes dynamically (USE_DYNAMIC is TRUE)
   * @see #USE_DYNAMIC
   * @see GenericPropertiesCreator
   */
  static {
    determineClasses();
  }

  /**
   * A specialized TreeNode for supporting filtering via Capabilities.
   */
  public class GOETreeNode
    extends DefaultMutableTreeNode {
    
    /** for serialization. */
    static final long serialVersionUID = -1707872446682150133L;
    
    /** color for "no support". */
    public final static String NO_SUPPORT = "silver";
    
    /** color for "maybe support". */
    public final static String MAYBE_SUPPORT = "blue";
    
    /** the Capabilities object to use for filtering. */
    protected Capabilities m_Capabilities = null;
    
    /**
     * Creates a tree node that has no parent and no children, but which 
     * allows children.
     */
    public GOETreeNode() {
      super();
    }
    
    /**
     * Creates a tree node with no parent, no children, but which allows 
     * children, and initializes it with the specified user object.
     * 
     * @param userObject	an Object provided by the user that constitutes 
     * 				the node's data
     */
    public GOETreeNode(Object userObject) {
      super(userObject);
    }
    
    /**
     * Creates a tree node with no parent, no children, initialized with the 
     * specified user object, and that allows children only if specified.

     * @param userObject	an Object provided by the user that constitutes 
     * 				the node's data
     * @param allowsChildren	if true, the node is allowed to have child nodes 
     * 				-- otherwise, it is always a leaf node
     */
    public GOETreeNode(Object userObject, boolean allowsChildren) {
      super(userObject, allowsChildren);
    }
    
    /**
     * generates if necessary a Capabilities object for the given leaf.
     */
    protected void initCapabilities() {
      String 	classname;
      Class	cls;
      Object	obj;
      
      if (m_Capabilities != null)
	return;
      if (!isLeaf())
	return;
      
      classname = getClassnameFromPath(new TreePath(getPath()));
      try {
	cls = Class.forName(classname);
	if (!ClassDiscovery.hasInterface(CapabilitiesHandler.class, cls))
	  return;
	
	obj = cls.newInstance();
	m_Capabilities = ((CapabilitiesHandler) obj).getCapabilities();
      }
      catch (Exception e) {
	// ignore it
      }
    }
    
    /**
     * returns a string representation of this treenode.
     * 
     * @return 		the text to display 
     */
    public String toString() {
      String	result;
      
      result = super.toString();
      
      if (m_CapabilitiesFilter != null) {
	initCapabilities();
	if (m_Capabilities != null) {
	  if (m_Capabilities.supportsMaybe(m_CapabilitiesFilter) && !m_Capabilities.supports(m_CapabilitiesFilter))
	    result = "<html><font color=\"" + MAYBE_SUPPORT + "\">" + result + "</font></i><html>";
	  else if (!m_Capabilities.supports(m_CapabilitiesFilter))
	    result = "<html><font color=\"" + NO_SUPPORT + "\">" + result + "</font></i><html>";
	}
      }
      
      return result;
    }
  }
  
  /**
   * A dialog for selecting Capabilities to look for in the GOE tree.
   */
  public class CapabilitiesFilterDialog 
    extends JDialog {
    
    /** for serialization. */
    static final long serialVersionUID = -7845503345689646266L;
    
    /** the dialog itself. */
    protected JDialog m_Self;
    
    /** the popup to display again. */
    protected JPopupMenu m_Popup = null;
    
    /** the capabilities used for initializing the dialog. */
    protected Capabilities m_Capabilities = new Capabilities(null);

    /** the label, listing the name of the superclass. */
    protected JLabel m_InfoLabel = new JLabel();
    
    /** the list with all the capabilities. */
    protected CheckBoxList m_List = new CheckBoxList();
    
    /** the OK button. */
    protected JButton m_OkButton = new JButton("OK");
    
    /** the Cancel button. */
    protected JButton m_CancelButton = new JButton("Cancel");
    
    /**
     * creates a dialog to choose Capabilities from.
     */
    public CapabilitiesFilterDialog() {
      super();

      m_Self = this;
      
      initGUI();
    }
    
    /**
     * sets up the GUI.
     */
    protected void initGUI() {
      JPanel			panel;
      CheckBoxListModel		model;

      setTitle("Filtering Capabilities...");
      setLayout(new BorderLayout());
      
      panel = new JPanel(new BorderLayout());
      panel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
      getContentPane().add(panel, BorderLayout.NORTH);
      m_InfoLabel.setText(
	    "<html>"
	  + m_ClassType.getName().replaceAll(".*\\.", "") + "s"
	  + " have to support <i>at least</i> the following capabilities <br>"
	  + "(the ones highlighted <font color=\"" + GOETreeNode.NO_SUPPORT + "\">" + GOETreeNode.NO_SUPPORT + "</font> don't meet these requirements <br>"
	  + "the ones highlighted  <font color=\"" + GOETreeNode.MAYBE_SUPPORT + "\">" + GOETreeNode.MAYBE_SUPPORT + "</font> possibly meet them):"
	  + "</html>");
      panel.add(m_InfoLabel, BorderLayout.CENTER);
      
      // list
      getContentPane().add(new JScrollPane(m_List), BorderLayout.CENTER);
      model = (CheckBoxListModel) m_List.getModel();
      for (Capability cap: Capability.values())
	model.addElement(cap);
      
      // buttons
      panel = new JPanel(new FlowLayout(FlowLayout.CENTER));
      getContentPane().add(panel, BorderLayout.SOUTH);
      
      m_OkButton.setMnemonic('O');
      m_OkButton.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          updateCapabilities();
          if (m_CapabilitiesFilter == null)
            m_CapabilitiesFilter = new Capabilities(null);
          m_CapabilitiesFilter.assign(m_Capabilities);
          m_Self.setVisible(false);
          showPopup();
        }
      });
      panel.add(m_OkButton);
      
      m_CancelButton.setMnemonic('C');
      m_CancelButton.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          m_Self.setVisible(false);
          showPopup();
        }
      });
      panel.add(m_CancelButton);
      pack();
    }

    /**
     * transfers the Capabilities object to the JList.
     * 
     * @see #m_Capabilities
     * @see #m_List
     */
    protected void updateList() {
      CheckBoxListModel		model;
      
      model = (CheckBoxListModel) m_List.getModel();

      for (Capability cap: Capability.values())
	model.setChecked(model.indexOf(cap), m_Capabilities.handles(cap));
    }
    
    /**
     * transfers the selected Capabilities from the JList to the 
     * Capabilities object.
     * 
     * @see #m_Capabilities
     * @see #m_List
     */
    protected void updateCapabilities() {
      CheckBoxListModel		model;
      
      model = (CheckBoxListModel) m_List.getModel();

      for (Capability cap: Capability.values()) {
	if (model.getChecked(model.indexOf(cap)))
          m_Capabilities.enable(cap);
	else
	  m_Capabilities.disable(cap);
      }
    }
    
    /**
     * sets the initial capabilities.
     * 
     * @param value the capabilities to use
     */
    public void setCapabilities(Capabilities value) {
      if (value != null)
	m_Capabilities.assign(value);
      else
	m_Capabilities = new Capabilities(null);
      
      updateList();
    }
    
    /**
     * returns the currently selected capabilities.
     * 
     * @return the currently selected capabilities
     */
    public Capabilities getCapabilities() {
      return m_Capabilities;
    }
    
    /**
     * sets the JPopupMenu to display again after closing the dialog.
     * 
     * @param value the JPopupMenu to display again
     */
    public void setPopup(JPopupMenu value) {
      m_Popup = value;
    }
    
    /**
     * returns the currently set JPopupMenu.
     * 
     * @return the current JPopupMenu
     */
    public JPopupMenu getPopup() {
      return m_Popup;
    }
    
    /**
     * if a JPopupMenu is set, it is displayed again. Displaying this dialog
     * closes any JPopupMenu automatically.
     */
    public void showPopup() {
      if (getPopup() != null)
	getPopup().setVisible(true);
    }
  }
  
  /**
   * Creates a popup menu containing a tree that is aware
   * of the screen dimensions.
   */
  public class JTreePopupMenu 
    extends JPopupMenu {
    
    /** for serialization. */
    static final long serialVersionUID = -3404546329655057387L;

    /** the popup itself. */
    private JPopupMenu m_Self;
    
    /** The tree. */
    private JTree m_tree;

    /** The scroller. */
    private JScrollPane m_scroller;

    /** The filter button in case of CapabilitiesHandlers. */
    private JButton m_FilterButton = new JButton("Filter...");

    /** The remove filter button in case of CapabilitiesHandlers. */
    private JButton m_RemoveFilterButton = new JButton("Remove filter");
    
    /** The button for closing the popup again. */
    private JButton m_CloseButton = new JButton("Close");
    
    /**
     * Constructs a new popup menu.
     *
     * @param tree the tree to put in the menu
     */
    public JTreePopupMenu(JTree tree) {

      m_Self = this;
      
      setLayout(new BorderLayout());
      JPanel panel = new JPanel(new FlowLayout(FlowLayout.RIGHT));
      add(panel, BorderLayout.SOUTH);

      if (ClassDiscovery.hasInterface(CapabilitiesHandler.class, m_ClassType)) {
	// filter
	m_FilterButton.setMnemonic('F');
	m_FilterButton.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    if (e.getSource() == m_FilterButton) {
	      CapabilitiesFilterDialog dialog = new CapabilitiesFilterDialog();
	      dialog.setCapabilities(m_CapabilitiesFilter);
	      dialog.setPopup(m_Self);
	      dialog.setVisible(true);
	      m_Support.firePropertyChange("", null, null);
	      repaint();
	    }
	  }
	});
	panel.add(m_FilterButton);
	
	// remove
	m_RemoveFilterButton.setMnemonic('R');
	m_RemoveFilterButton.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    if (e.getSource() == m_RemoveFilterButton) {
	      m_CapabilitiesFilter = null;
	      m_Support.firePropertyChange("", null, null);
	      repaint();
	    }
	  }
	});
	panel.add(m_RemoveFilterButton);
      }

      // close
      m_CloseButton.setMnemonic('C');
      m_CloseButton.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          if (e.getSource() == m_CloseButton) {
            m_Self.setVisible(false);
          }
        }
      });
      panel.add(m_CloseButton);
      
      m_tree = tree;
      
      JPanel treeView = new JPanel();
      treeView.setLayout(new BorderLayout());
      treeView.add(m_tree, BorderLayout.NORTH);
      
      // make backgrounds look the same
      treeView.setBackground(m_tree.getBackground());

      m_scroller = new JScrollPane(treeView);
      
      m_scroller.setPreferredSize(new Dimension(300, 400));
      m_scroller.getVerticalScrollBar().setUnitIncrement(20);

      add(m_scroller);
    }

    /**
     * Displays the menu, making sure it will fit on the screen.
     *
     * @param invoker the component thast invoked the menu
     * @param x the x location of the popup
     * @param y the y location of the popup
     */
    public void show(Component invoker, int x, int y) {

      super.show(invoker, x, y);

      // calculate available screen area for popup
      java.awt.Point location = getLocationOnScreen();
      java.awt.Dimension screenSize = getToolkit().getScreenSize();
      int maxWidth = (int) (screenSize.getWidth() - location.getX());
      int maxHeight = (int) (screenSize.getHeight() - location.getY());

      // if the part of the popup goes off the screen then resize it
      Dimension scrollerSize = m_scroller.getPreferredSize();
      int height = (int) scrollerSize.getHeight();
      int width = (int) scrollerSize.getWidth();
      if (width > maxWidth) width = maxWidth;
      if (height > maxHeight) height = maxHeight;
      
      // commit any size changes
      m_scroller.setPreferredSize(new Dimension(width, height));
      revalidate();
      pack();
    }
  }

  /**
   * Handles the GUI side of editing values.
   */
  public class GOEPanel 
    extends JPanel {
    
    /** for serialization. */
    static final long serialVersionUID = 3656028520876011335L;
    
    /** The component that performs classifier customization. */
    protected PropertySheetPanel m_ChildPropertySheet;
    
    /** The name of the current class. */
    protected JLabel m_ClassNameLabel;

    /** Open object from disk. */
    protected JButton m_OpenBut;
    
    /** Save object to disk. */
    protected JButton m_SaveBut;
    
    /** ok button. */
    protected JButton m_okBut;
    
    /** cancel button. */
    protected JButton m_cancelBut;
    
    /** The filechooser for opening and saving object files. */
    protected JFileChooser m_FileChooser;
    
    /** Creates the GUI editor component. */
    public GOEPanel() {
    
      m_Backup = copyObject(m_Object);
      
      m_ClassNameLabel = new JLabel("None");
      m_ClassNameLabel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));

      m_ChildPropertySheet = new PropertySheetPanel();
      m_ChildPropertySheet.addPropertyChangeListener
	(new PropertyChangeListener() {
	    public void propertyChange(PropertyChangeEvent evt) {
	      m_Support.firePropertyChange("", null, null);
	    }
	  });
      
      m_OpenBut = new JButton("Open...");
      m_OpenBut.setToolTipText("Load a configured object");
      m_OpenBut.setEnabled(true);
      m_OpenBut.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    Object object = openObject();
	    if (object != null) {
	      // setValue takes care of: Making sure obj is of right type,
	      // and firing property change.
	      setValue(object);
	      // Need a second setValue to get property values filled in OK.
	      // Not sure why.
	      setValue(object);
	    }
	  }
	});
      
      m_SaveBut = new JButton("Save...");
      m_SaveBut.setToolTipText("Save the current configured object");
      m_SaveBut.setEnabled(true);
      m_SaveBut.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    saveObject(m_Object);
	  }
	});
      
      m_okBut = new JButton("OK");
      m_okBut.setEnabled(true);
      m_okBut.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {

	    m_Backup = copyObject(m_Object);
	    if ((getTopLevelAncestor() != null)
		&& (getTopLevelAncestor() instanceof Window)) {
	      Window w = (Window) getTopLevelAncestor();
	      w.dispose();
	    }
	  }
	});
      
      m_cancelBut = new JButton("Cancel");
      m_cancelBut.setEnabled(true);
      m_cancelBut.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {		 
	    if (m_Backup != null) {
	
	      m_Object = copyObject(m_Backup);
	      
	      // To fire property change
	      m_Support.firePropertyChange("", null, null);
	      m_ObjectNames = getClassesFromProperties();
	      updateObjectNames();
	      updateChildPropertySheet();
	    }
	    if ((getTopLevelAncestor() != null)
		&& (getTopLevelAncestor() instanceof Window)) {
	      Window w = (Window) getTopLevelAncestor();
	      w.dispose();
	    }
	  }
	});
      
      setLayout(new BorderLayout());

      if (m_canChangeClassInDialog) {
	JButton chooseButton = createChooseClassButton();
	JPanel top = new JPanel();
	top.setLayout(new BorderLayout());
	top.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
	top.add(chooseButton, BorderLayout.WEST);
	top.add(m_ClassNameLabel, BorderLayout.CENTER);
	add(top, BorderLayout.NORTH);
      } else {
	add(m_ClassNameLabel, BorderLayout.NORTH);
      }

      add(m_ChildPropertySheet, BorderLayout.CENTER);
      // Since we resize to the size of the property sheet, a scrollpane isn't
      // typically needed
      // add(new JScrollPane(m_ChildPropertySheet), BorderLayout.CENTER);
      
      JPanel okcButs = new JPanel();
      okcButs.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
      okcButs.setLayout(new GridLayout(1, 4, 5, 5));
      okcButs.add(m_OpenBut);
      okcButs.add(m_SaveBut);
      okcButs.add(m_okBut);
      okcButs.add(m_cancelBut);
      add(okcButs, BorderLayout.SOUTH);

      if (m_ClassType != null) {
	m_ObjectNames = getClassesFromProperties();
	if (m_Object != null) {
	  updateObjectNames();
	  updateChildPropertySheet();
	}
      }
    }
    
    /**
     * Enables/disables the cancel button.
     *
     * @param flag true to enable cancel button, false
     * to disable it
     */
    protected void setCancelButton(boolean flag) {

      if(m_cancelBut != null)
	m_cancelBut.setEnabled(flag);
    }
    
    /**
     * Opens an object from a file selected by the user.
     * 
     * @return the loaded object, or null if the operation was cancelled
     */
    protected Object openObject() {
      
      if (m_FileChooser == null) {
	createFileChooser();
      }
      int returnVal = m_FileChooser.showOpenDialog(this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
	File selected = m_FileChooser.getSelectedFile();
	try {
	  ObjectInputStream oi = new ObjectInputStream(new BufferedInputStream(new FileInputStream(selected)));
	  Object obj = oi.readObject();
	  oi.close();
	  if (!m_ClassType.isAssignableFrom(obj.getClass())) {
	    throw new Exception("Object not of type: " + m_ClassType.getName());
	  }
	  return obj;
	} catch (Exception ex) {
	  JOptionPane.showMessageDialog(this,
					"Couldn't read object: "
					+ selected.getName() 
					+ "\n" + ex.getMessage(),
					"Open object file",
					JOptionPane.ERROR_MESSAGE);
	}
      }
      return null;
    }
    
    /**
     * Saves an object to a file selected by the user.
     * 
     * @param object the object to save
     */
    protected void saveObject(Object object) {
      
      if (m_FileChooser == null) {
	createFileChooser();
      }
      int returnVal = m_FileChooser.showSaveDialog(this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
	File sFile = m_FileChooser.getSelectedFile();
	try {
	  ObjectOutputStream oo = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(sFile)));
	  oo.writeObject(object);
	  oo.close();
	} catch (Exception ex) {
	  JOptionPane.showMessageDialog(this,
					"Couldn't write to file: "
					+ sFile.getName() 
					+ "\n" + ex.getMessage(),
					"Save object",
					JOptionPane.ERROR_MESSAGE);
	}
      }
    }

    /**
     * Creates the file chooser the user will use to save/load files with.
     */
    protected void createFileChooser() {
      
      m_FileChooser = new JFileChooser(new File(System.getProperty("user.dir")));
      m_FileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
    }
    
    /**
     * Makes a copy of an object using serialization.
     * 
     * @param source the object to copy
     * @return a copy of the source object
     */
    protected Object copyObject(Object source) {

      Object result = null;
      try {
        result = GenericObjectEditor.this.makeCopy(source);
	setCancelButton(true);
	
      } catch (Exception ex) {
	setCancelButton(false);
	System.err.println("GenericObjectEditor: Problem making backup object");
	System.err.println(ex);
      }
      return result;
    }
    
    /**
     * Allows customization of the action label on the dialog.
     * 
     * @param newLabel the new string for the ok button
     */
    public void setOkButtonText(String newLabel) {

      m_okBut.setText(newLabel);
    }

    /** 
     * This is used to hook an action listener to the ok button.
     * 
     * @param a The action listener.
     */
    public void addOkListener(ActionListener a) {

      m_okBut.addActionListener(a);
    }
    
    /**
     * This is used to hook an action listener to the cancel button.
     * 
     * @param a The action listener.
     */
    public void addCancelListener(ActionListener a) {

      m_cancelBut.addActionListener(a);
    }
	
    /**
     * This is used to remove an action listener from the ok button.
     * 
     * @param a The action listener
     */
    public void removeOkListener(ActionListener a) {

      m_okBut.removeActionListener(a);
    }
    
    /**
     * This is used to remove an action listener from the cancel button.
     * 
     * @param a The action listener
     */
    public void removeCancelListener(ActionListener a) {

      m_cancelBut.removeActionListener(a);
    }
    
    /**
     * Updates the child property sheet, and creates if needed.
     */
    public void updateChildPropertySheet() {
      
      // Update the object name displayed
      String className = "None";
      if (m_Object != null) {
	className = m_Object.getClass().getName();
      }
      m_ClassNameLabel.setText(className);

      // Set the object as the target of the propertysheet
      m_ChildPropertySheet.setTarget(m_Object);

      // Adjust size of containing window if possible
      if ((getTopLevelAncestor() != null)
	  && (getTopLevelAncestor() instanceof Window)) {
	((Window) getTopLevelAncestor()).pack();
      }
    }	
  }

  /**
   * Default constructor.
   */
  public GenericObjectEditor() {

    this(false);
  }

  /**
   * Constructor that allows specifying whether it is possible
   * to change the class within the editor dialog.
   *
   * @param canChangeClassInDialog whether the user can change the class
   */
  public GenericObjectEditor(boolean canChangeClassInDialog) {

    m_canChangeClassInDialog = canChangeClassInDialog;
  }
  
  /**
   * registers all the editors in Weka.
   */
  public static void registerEditors() {
    Properties 		props;
    Enumeration 	enm;
    String 		name;
    String 		value;
    Class 		baseCls;
    Class		cls;

    if (m_EditorsRegistered)
      return;
    
    System.err.println("---Registering Weka Editors---");
    m_EditorsRegistered = true;

    // load properties
    try {
      props = Utils.readProperties(GUIEDITORS_PROPERTY_FILE);
    }
    catch (Exception e) {
      props = new Properties();
      e.printStackTrace();
    }
    
    enm = props.propertyNames();
    while (enm.hasMoreElements()) {
      name  = enm.nextElement().toString();
      value = props.getProperty(name, "");
      
      registerEditor(name, value);
    }
  }
  
  public static void registerEditor(String name, String value) {
    Class              baseCls;
    Class               cls;
    
    try {
      // array class?
      if (name.endsWith("[]")) {
        baseCls = Class.forName(name.substring(0, name.indexOf("[]")));
        cls = Array.newInstance(baseCls, 1).getClass();
      }
      else {
        cls = Class.forName(name);
      }
      // register
      PropertyEditorManager.registerEditor(cls, Class.forName(value));
    }
    catch (Exception e) {
      System.err.println("Problem registering " + name + "/" + value + ": " + e);
    }
  }

  /**
   * Sets whether the user can change the class in the dialog.
   * 
   * @param value	if true then the user can change the class
   */
  public void setCanChangeClassInDialog(boolean value) {
    m_canChangeClassInDialog = value;
  }
  
  /**
   * Returns whether the user can change the class in the dialog.
   * 
   * @return		true if the user can change the class
   */
  public boolean getCanChangeClassInDialog() {
    return m_canChangeClassInDialog;
  }
  
  /**
   * Returns the backup object (may be null if there is no
   * backup.
   *
   * @return the backup object
   */
  public Object getBackup() {
    return m_Backup;
  }
  
  /**
   * returns the name of the root element of the given class name, 
   * <code>null</code> if it doesn't contain the separator.
   * 
   * @param clsname the full classname
   * @param separator the separator 
   * @return string the root element
   */
  protected static String getRootFromClass(String clsname, String separator) {
    if (clsname.indexOf(separator) > -1)
      return clsname.substring(0, clsname.indexOf(separator));
    else
      return null;
  }
  
  /**
   * parses the given string of classes separated by ", " and returns the
   * a hashtable with as many entries as there are different root elements in 
   * the class names (the key is the root element). E.g. if there's only 
   * "weka." as the prefix for all classes the a hashtable of size 1 is returned. 
   * if NULL is the input, then NULL is also returned.
   * 
   * @param classes the classnames to work on
   * @return for each distinct root element in the classnames, one entry in
   * the hashtable (with the root element as key)
   */
  public static Hashtable sortClassesByRoot(String classes) {
    Hashtable                 roots;
    Hashtable                 result;
    Enumeration               enm;
    int                       i;
    StringTokenizer           tok;
    String                    clsname;
    Vector                    list;
    HierarchyPropertyParser   hpp;
    String                    separator;
    String                    root;
    String                    tmpStr;
    
    if (classes == null)
      return null;
    
    roots     = new Hashtable();
    hpp       = new HierarchyPropertyParser();
    separator = hpp.getSeperator();
    
    // go over all classnames and store them in the hashtable, with the
    // root element as the key
    tok   = new StringTokenizer(classes, ", ");
    while (tok.hasMoreElements()) {
      clsname = tok.nextToken();
      root    = getRootFromClass(clsname, separator);
      if (root == null)
        continue;
      
      // already stored?
      if (!roots.containsKey(root)) {
        list = new Vector();
        roots.put(root, list);
      }
      else {
        list = (Vector) roots.get(root);
      }
      
      list.add(clsname);
    }
    
    // build result
    result = new Hashtable();
    enm    = roots.keys();
    while (enm.hasMoreElements()) {
      root = (String) enm.nextElement();
      list = (Vector) roots.get(root);
      tmpStr = "";
      for (i = 0; i < list.size(); i++) {
        if (i > 0)
          tmpStr += ",";
        tmpStr += (String) list.get(i);
      }
      result.put(root, tmpStr);
    }
      
    return result;
  }

  /** 
   * Called when the class of object being edited changes. 
   * 
   * @return the hashtable containing the HierarchyPropertyParsers for the root
   *         elements
   */
  protected Hashtable getClassesFromProperties() {	    

    Hashtable hpps = new Hashtable();
    String className = m_ClassType.getName();
    Hashtable typeOptions = sortClassesByRoot(EDITOR_PROPERTIES.getProperty(className));
    if (typeOptions == null) {
      /*
      System.err.println("Warning: No configuration property found in\n"
			 + PROPERTY_FILE + "\n"
			 + "for " + className);
      */
    } else {		    
      try {
        Enumeration enm = typeOptions.keys();
        while (enm.hasMoreElements()) {
          String root = (String) enm.nextElement();
          String typeOption = (String) typeOptions.get(root);
          HierarchyPropertyParser hpp = new HierarchyPropertyParser();
          hpp.build(typeOption, ", ");
	  hpps.put(root, hpp);
        }
      } catch (Exception ex) {
	System.err.println("Invalid property: " + typeOptions);
      }	    
    }
    return hpps;
  }
  
  /**
   * Updates the list of selectable object names, adding any new names to the list.
   */
  protected void updateObjectNames() {
    
    if (m_ObjectNames == null) {
      m_ObjectNames = getClassesFromProperties();
    }
    
    if (m_Object != null) {
      String className = m_Object.getClass().getName();
      String root = getRootFromClass(className, new HierarchyPropertyParser().getSeperator());
      HierarchyPropertyParser hpp = (HierarchyPropertyParser) m_ObjectNames.get(root);
      if (hpp != null) {
        if(!hpp.contains(className)){
          hpp.add(className);
        }
      }
    }
  }
  
  /**
   * Sets whether the editor is "enabled", meaning that the current
   * values will be painted.
   *
   * @param newVal a value of type 'boolean'
   */
  public void setEnabled(boolean newVal) {
    
    if (newVal != m_Enabled) {
      m_Enabled = newVal;
    }
  }
  
  /**
   * Sets the class of values that can be edited.
   *
   * @param type a value of type 'Class'
   */
  public void setClassType(Class type) {
    
    m_ClassType = type;
    m_ObjectNames = getClassesFromProperties();
  }
  
  /**
   * Sets the current object to be the default, taken as the first item in
   * the chooser.
   */
  public void setDefaultValue() {
    
    if (m_ClassType == null) {
      System.err.println("No ClassType set up for GenericObjectEditor!!");
      return;
    }	
    
    Hashtable hpps = getClassesFromProperties();
    HierarchyPropertyParser hpp = null;
    Enumeration enm = hpps.elements();
    
    try{
      while (enm.hasMoreElements()) {
        hpp = (HierarchyPropertyParser) enm.nextElement(); 
        if(hpp.depth() > 0) {		
          hpp.goToRoot();
          while(!hpp.isLeafReached())
            hpp.goToChild(0);
          
          String defaultValue = hpp.fullValue();
          setValue(Class.forName(defaultValue).newInstance());
        }
      }
    }catch(Exception ex){
      System.err.println("Problem loading the first class: "+
			 hpp.fullValue());
      ex.printStackTrace();
    }
  }
  
  /**
   * Sets the current Object. If the Object is in the
   * Object chooser, this becomes the selected item (and added
   * to the chooser if necessary).
   *
   * @param o an object that must be a Object.
   */
  public void setValue(Object o) {
    
    if (m_ClassType == null) {
      System.err.println("No ClassType set up for GenericObjectEditor!!");
      return;
    }
    if (!m_ClassType.isAssignableFrom(o.getClass())) {
      System.err.println("setValue object not of correct type!");
      return;
    }
    
    setObject(o);

    if (m_EditorComponent != null) m_EditorComponent.repaint();

    updateObjectNames();
  }
  
  /**
   * Sets the current Object.
   *
   * @param c a value of type 'Object'
   */
  protected void setObject(Object c) {
    
    // This should really call equals() for comparison.
    boolean trueChange ;
    if (getValue() != null) {
      trueChange = (!c.equals(getValue()));
    }
    else
      trueChange = true;
    
    m_Backup = m_Object;
    
    m_Object = c;
    
    if (m_EditorComponent != null) {
      m_EditorComponent.updateChildPropertySheet();
    }
    if (trueChange) {
      m_Support.firePropertyChange("", null, null);
    }
  }
  
  /**
   * Gets the current Object.
   *
   * @return the current Object
   */
  public Object getValue() {
    
    Object result = null;
    try {
      result = makeCopy(m_Object);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    return result;
  }
  
  /**
   * Supposedly returns an initialization string to create a Object
   * identical to the current one, including it's state, but this doesn't
   * appear possible given that the initialization string isn't supposed to
   * contain multiple statements.
   *
   * @return the java source code initialisation string
   */
  public String getJavaInitializationString() {

    return "new " + m_Object.getClass().getName() + "()";
  }

  /**
   * Returns true to indicate that we can paint a representation of the
   * Object.
   *
   * @return true
   */
  public boolean isPaintable() {

    return true;
  }

  /**
   * Paints a representation of the current Object.
   *
   * @param gfx the graphics context to use
   * @param box the area we are allowed to paint into
   */
  public void paintValue(java.awt.Graphics gfx, java.awt.Rectangle box) {

    if (m_Enabled) {
      String rep;
      if (m_Object != null) {
	if (m_Object instanceof CustomDisplayStringProvider) {
	  rep = ((CustomDisplayStringProvider) m_Object).toDisplay();
	}
	else {
	  rep = m_Object.getClass().getName();
	  int dotPos = rep.lastIndexOf('.');
	  if (dotPos != -1) {
	    rep = rep.substring(dotPos + 1);
	  }
	}
      } else {
	rep = "None";
      }
      java.awt.Font originalFont = gfx.getFont();
      gfx.setFont(originalFont.deriveFont(java.awt.Font.BOLD));

      FontMetrics fm = gfx.getFontMetrics();
      int vpad = (box.height - fm.getHeight());
      gfx.drawString(rep, 2, fm.getAscent() + vpad);
      int repwidth = fm.stringWidth(rep);

      gfx.setFont(originalFont);
      if ((m_Object instanceof OptionHandler) && !(m_Object instanceof CustomDisplayStringProvider)) {
	gfx.drawString(" " + Utils.joinOptions(((OptionHandler)m_Object).getOptions()),
					       repwidth + 2, fm.getAscent() + vpad);
      }
    }
  }

  /**
   * Returns null as we don't support getting/setting values as text.
   *
   * @return null
   */
  public String getAsText() {

    return null;
  }

  /**
   * Returns null as we don't support getting/setting values as text. 
   *
   * @param text the text value
   * @throws IllegalArgumentException as we don't support
   * getting/setting values as text.
   */
  public void setAsText(String text) {

    throw new IllegalArgumentException(text);
  }

  /**
   * Returns null as we don't support getting values as tags.
   *
   * @return null
   */
  public String[] getTags() {

    return null;
  }

  /**
   * Returns true because we do support a custom editor.
   *
   * @return true
   */
  public boolean supportsCustomEditor() {

    return true;
  }
  
  /**
   * Returns the array editing component.
   *
   * @return a value of type 'java.awt.Component'
   */
  public java.awt.Component getCustomEditor() {

    if (m_EditorComponent == null) {
      m_EditorComponent = new GOEPanel();
    }
    return m_EditorComponent;
  }

  /**
   * Adds a PropertyChangeListener who will be notified of value changes.
   *
   * @param l a value of type 'PropertyChangeListener'
   */
  public void addPropertyChangeListener(PropertyChangeListener l) {

    m_Support.addPropertyChangeListener(l);
  }

  /**
   * Removes a PropertyChangeListener.
   *
   * @param l a value of type 'PropertyChangeListener'
   */
  public void removePropertyChangeListener(PropertyChangeListener l) {

    m_Support.removePropertyChangeListener(l);
  }

  /**
   * Gets the custom panel used for editing the object.
   *
   * @return the panel
   */
  public JPanel getCustomPanel() {
    final JButton chooseButton = createChooseClassButton();    
    m_ObjectPropertyPanel = new PropertyPanel(this, true);
    
    JPanel customPanel = new JPanel() {
      public void setEnabled(boolean enabled) {
        super.setEnabled(enabled);
        chooseButton.setEnabled(enabled);
      }
    };
    customPanel.setLayout(new BorderLayout());
    customPanel.add(chooseButton, BorderLayout.WEST);
    customPanel.add(m_ObjectPropertyPanel, BorderLayout.CENTER);
    return customPanel;
  }

  /**
   * Creates a button that when clicked will enable the user to change
   * the class of the object being edited.
   * 
   * @return the choose button
   */
  protected JButton createChooseClassButton() {

    JButton setButton = new JButton("Choose");

    // anonymous action listener shows a JTree popup and allows the user
    // to choose the class they want
    setButton.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {

	  JPopupMenu popup = getChooseClassPopupMenu();

	  // show the popup where the source component is
	  if (e.getSource() instanceof Component) {
	    Component comp = (Component) e.getSource();
	    popup.show(comp, comp.getX(), comp.getY());
	    popup.pack();
	    popup.repaint();
	  }
	}
      });

    return setButton;
  }

  /**
   * creates a classname from the given path.
   * 
   * @param path	the path to generate the classname from
   * @return		the generated classname
   */
  protected String getClassnameFromPath(TreePath path) {
    StringBuffer classname = new StringBuffer();
    
    // recreate class name from path
    int start = 0;
    if (m_ObjectNames.size() > 1)
      start = 1;

    for (int i = start; i < path.getPathCount(); i++) {
      if (i>start) classname.append(".");
      classname.append(
	  (String) ((GOETreeNode) path.getPathComponent(i)).getUserObject());
    }
    
    return classname.toString();
  }
  
  /**
   * Returns a popup menu that allows the user to change
   * the class of object.
   *
   * @return a JPopupMenu that when shown will let the user choose the class
   */
  public JPopupMenu getChooseClassPopupMenu() {

    updateObjectNames();

    // create the tree, and find the path to the current class
    m_treeNodeOfCurrentObject = null;
    final JTree tree = createTree(m_ObjectNames);
    if (m_treeNodeOfCurrentObject != null) {
      tree.setSelectionPath(new TreePath(m_treeNodeOfCurrentObject.getPath()));
    }
    tree.getSelectionModel().setSelectionMode
      (TreeSelectionModel.SINGLE_TREE_SELECTION);

    // create the popup
    final JPopupMenu popup = new JTreePopupMenu(tree);

    // respond when the user chooses a class
    tree.addTreeSelectionListener(new TreeSelectionListener() {
	public void valueChanged(TreeSelectionEvent e) {
	  GOETreeNode node = (GOETreeNode) tree.getLastSelectedPathComponent();
	  
	  if (node == null) 
	    return;
	  
	  if (node.isLeaf()) {
	    classSelected(getClassnameFromPath(tree.getSelectionPath()));
	    popup.setVisible(false);
	  }
	}
      });
    
    return popup;
  }

  /**
   * Creates a JTree from an object heirarchy.
   *
   * @param hpps the hierarchy of objects to mirror in the tree
   * @return a JTree representation of the hierarchy
   */
  protected JTree createTree(Hashtable hpps) {
    GOETreeNode             superRoot;
    Enumeration             enm;
    HierarchyPropertyParser hpp;
    
    if (hpps.size() > 1)
      superRoot = new GOETreeNode("root");
    else
      superRoot = null;

    enm = hpps.elements();
    while (enm.hasMoreElements()) {
      hpp = (HierarchyPropertyParser) enm.nextElement();
      hpp.goToRoot();
      GOETreeNode root = new GOETreeNode(hpp.getValue());
      addChildrenToTree(root, hpp);
      
      if (superRoot == null)
        superRoot = root;
      else
        superRoot.add(root);
    }
    
    JTree tree = new JTree(superRoot);
    
    return tree;
  }

  /**
   * Recursively builds a JTree from an object heirarchy.
   * Also updates m_treeNodeOfCurrentObject if the current object
   * is discovered during creation.
   *
   * @param tree the root of the tree to add children to
   * @param hpp the hierarchy of objects to mirror in the tree
   */
  protected void addChildrenToTree(GOETreeNode tree,
				   HierarchyPropertyParser hpp) {

    try {
      for (int i=0; i<hpp.numChildren(); i++) {
	hpp.goToChild(i);
	GOETreeNode child = new GOETreeNode(hpp.getValue());
	if ((m_Object != null) &&
	    m_Object.getClass().getName().equals(hpp.fullValue())) {
	  m_treeNodeOfCurrentObject = child;
	}
	tree.add(child);
	addChildrenToTree(child, hpp);
	hpp.goToParent();
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Called when the user selects an class type to change to.
   *
   * @param className the name of the class that was selected
   */
  protected void classSelected(String className) {

    try {		    
      if ((m_Object != null) && m_Object.getClass().getName().equals(className)) {
	return;
      }
      
      setValue(Class.forName(className).newInstance());
      //m_ObjectPropertyPanel.showPropertyDialog();
      if (m_EditorComponent != null) {
	m_EditorComponent.updateChildPropertySheet();
      }
    } catch (Exception ex) {
      JOptionPane.showMessageDialog(null,
				    "Could not create an example of\n"
				    + className + "\n"
				    + "from the current classpath",
				    "Class load failed",
				    JOptionPane.ERROR_MESSAGE);
      ex.printStackTrace();
      try {
	if(m_Backup != null)
	  setValue(m_Backup);
	else
	  setDefaultValue();			
      } catch(Exception e) {
	System.err.println(ex.getMessage());
	ex.printStackTrace();
      }
    }
  }
  
  /**
   * Sets the capabilities to use for filtering.
   * 
   * @param value	the object to get the filter capabilities from
   */
  public void setCapabilitiesFilter(Capabilities value) {
    m_CapabilitiesFilter = new Capabilities(null);
    m_CapabilitiesFilter.assign(value);
  }
  
  /**
   * Returns the current Capabilities filter, can be null.
   * 
   * @return		the current Capabiliities used for filtering
   */
  public Capabilities getCapabilitiesFilter() {
    return m_CapabilitiesFilter;
  }
  
  /**
   * Removes the current Capabilities filter.
   */
  public void removeCapabilitiesFilter() {
    m_CapabilitiesFilter = null;
  }

  /**
   * Makes a copy of an object using serialization.
   * 
   * @param source the object to copy
   * @return a copy of the source object
   * @exception Exception if the copy fails
   */
  public static Object makeCopy(Object source) throws Exception {
    SerializedObject so = new SerializedObject(source);
    Object result = so.getObject();
    return result;
  }
  
  /**
   * Returns the available classnames for a certain property in the 
   * props file.
   * 
   * @param property	the property to get the classnames for
   * @return		the classnames
   */
  public static Vector<String> getClassnames(String property) {
    Vector<String>	result;
    String		value;
    String[]		items;
    int			i;
    
    result = new Vector<String>();
    
    value = EDITOR_PROPERTIES.getProperty(property, "").replaceAll(" ", "").trim();
    if (value.length() > 0) {
      items = value.split(",");
      for (i = 0; i < items.length; i++)
	result.add(items[i]);
    }
    
    return result;
  }

  /**
   * Tests out the Object editor from the command line.
   *
   * @param args may contain the class name of a Object to edit
   */
  public static void main(String [] args) {

    try {
      GenericObjectEditor.registerEditors();
      GenericObjectEditor ce = new GenericObjectEditor(true);
      ce.setClassType(weka.classifiers.Classifier.class);
      Object initial = new weka.classifiers.rules.ZeroR();
      if (args.length > 0){
	ce.setClassType(Class.forName(args[0]));
	if(args.length > 1){
	  initial = (Object)Class.forName(args[1]).newInstance();
	  ce.setValue(initial);
	}
	else
	  ce.setDefaultValue();
      }
      else	  
	ce.setValue(initial);
      
      PropertyDialog pd = new PropertyDialog((Frame) null, ce, 100, 100);
      pd.addWindowListener(new WindowAdapter() {
	  public void windowClosing(WindowEvent e) {
	    PropertyEditor pe = ((PropertyDialog)e.getSource()).getEditor();
	    Object c = (Object)pe.getValue();
	    String options = "";
	    if (c instanceof OptionHandler) {
	      options = Utils.joinOptions(((OptionHandler)c).getOptions());
	    }
	    System.out.println(c.getClass().getName() + " " + options);
	    System.exit(0);
	  }
	});
      pd.setVisible(true);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
