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
 *    PropertySelectorDialog.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import weka.experiment.PropertyNode;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.BeanInfo;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.beans.PropertyEditor;
import java.beans.PropertyEditorManager;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

/** 
 * Allows the user to select any (supported) property of an object, including
 * properties that any of it's property values may have.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.9 $
 */
public class PropertySelectorDialog
  extends JDialog {

  /** for serialization */
  private static final long serialVersionUID = -3155058124137930518L;
  
  /** Click to choose the currently selected property */
  protected JButton m_SelectBut = new JButton("Select");

  /** Click to cancel the property selection */
  protected JButton m_CancelBut = new JButton("Cancel");

  /** The root of the property tree */
  protected DefaultMutableTreeNode m_Root;

  /** The object at the root of the tree */
  protected Object m_RootObject;

  /** Whether the selection was made or cancelled */
  protected int m_Result;

  /** Stores the path to the selected property */
  protected Object [] m_ResultPath;

  /** The component displaying the property tree */
  protected JTree m_Tree;

  /** Signifies an OK property selection */
  public static final int APPROVE_OPTION = 0;

  /** Signifies a cancelled property selection */
  public static final int CANCEL_OPTION = 1;
  
  /**
   * Create the property selection dialog.
   *
   * @param parentFrame the parent frame of the dialog
   * @param rootObject the object containing properties to select from
   */
  public PropertySelectorDialog(Frame parentFrame, Object rootObject) {
    
    super(parentFrame, "Select a property", true);
    m_CancelBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	m_Result = CANCEL_OPTION;
	setVisible(false);
      }
    });
    m_SelectBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	// value = path from root to selected;
	TreePath tPath = m_Tree.getSelectionPath();
	if (tPath == null) {
	  m_Result = CANCEL_OPTION;
	} else {
	  m_ResultPath = tPath.getPath();
	  if ((m_ResultPath == null) || (m_ResultPath.length < 2)) {
	    m_Result = CANCEL_OPTION;
	  } else {
	    m_Result = APPROVE_OPTION;
	  }
	} 
	setVisible(false);
      }
    });
    m_RootObject = rootObject;
    m_Root = new DefaultMutableTreeNode(
	     new PropertyNode(m_RootObject));
    createNodes(m_Root);
    
    Container c = getContentPane();
    c.setLayout(new BorderLayout());
    //    setBorder(BorderFactory.createTitledBorder("Select a property"));
    Box b1 = new Box(BoxLayout.X_AXIS);
    b1.add(m_SelectBut);
    b1.add(Box.createHorizontalStrut(10));
    b1.add(m_CancelBut);
    c.add(b1, BorderLayout.SOUTH);
    m_Tree = new JTree(m_Root);
    m_Tree.getSelectionModel()
      .setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
    c.add(new JScrollPane(m_Tree), BorderLayout.CENTER);
    pack();
  }

  /**
   * Pops up the modal dialog and waits for cancel or a selection.
   *
   * @return either APPROVE_OPTION, or CANCEL_OPTION
   */
  public int showDialog() {

    m_Result = CANCEL_OPTION;
    setVisible(true);
    return m_Result;
  }

  /**
   * Gets the path of property nodes to the selected property.
   *
   * @return an array of PropertyNodes
   */
  public PropertyNode [] getPath() {

    PropertyNode [] result = new PropertyNode [m_ResultPath.length - 1];
    for (int i = 0; i < result.length; i++) {
      result[i] = (PropertyNode) ((DefaultMutableTreeNode) m_ResultPath[i + 1])
	.getUserObject();
    }
    return result;
  }

  /**
   * Creates the property tree below the current node.
   *
   * @param localNode a value of type 'DefaultMutableTreeNode'
   */
  protected void createNodes(DefaultMutableTreeNode localNode) {

    PropertyNode pNode = (PropertyNode)localNode.getUserObject();
    Object localObject = pNode.value;
    // Find all the properties of the object in the root node
    PropertyDescriptor localProperties[];
    try {
      BeanInfo bi = Introspector.getBeanInfo(localObject.getClass());
      localProperties = bi.getPropertyDescriptors();
    } catch (IntrospectionException ex) {
      System.err.println("PropertySelectorDialog: Couldn't introspect");
      return;
    }

    // Put their values into child nodes.
    for (int i = 0; i < localProperties.length; i++) {
      // Don't display hidden or expert properties.
      if (localProperties[i].isHidden() || localProperties[i].isExpert()) {
	continue;
      }
      String name = localProperties[i].getDisplayName();
      Class type = localProperties[i].getPropertyType();
      Method getter = localProperties[i].getReadMethod();
      Method setter = localProperties[i].getWriteMethod();
      Object value = null;
      // Only display read/write properties.
      if (getter == null || setter == null) {
	continue;
      }
      try {
	Object args[] = { };
	value = getter.invoke(localObject, args);
	PropertyEditor editor = null;
	Class pec = localProperties[i].getPropertyEditorClass();
	if (pec != null) {
	  try {
	    editor = (PropertyEditor)pec.newInstance();
	  } catch (Exception ex) {
	  }
	}
	if (editor == null) {
	  editor = PropertyEditorManager.findEditor(type);
	}
	if ((editor == null) || (value == null)) {
	  continue;
	}
      } catch (InvocationTargetException ex) {
	System.err.println("Skipping property " + name
			   + " ; exception on target: "
			   + ex.getTargetException());
	ex.getTargetException().printStackTrace();
	continue;
      } catch (Exception ex) {
	System.err.println("Skipping property " + name
			   + " ; exception: " + ex);
	ex.printStackTrace();
	continue;
      }
      // Make a child node
      DefaultMutableTreeNode child = new DefaultMutableTreeNode(
				     new PropertyNode(value,
						      localProperties[i],
						      localObject.getClass()));
      localNode.add(child);
      createNodes(child);
    }
  }

  
  /**
   * Tests out the property selector from the command line.
   *
   * @param args ignored
   */
  public static void main(String [] args) {

    try {
      GenericObjectEditor.registerEditors();

      Object rp
	= new weka.experiment.AveragingResultProducer();
      final PropertySelectorDialog jd = new PropertySelectorDialog(null, rp);
      int result = jd.showDialog();
      if (result == PropertySelectorDialog.APPROVE_OPTION) {
	System.err.println("Property Selected");
	PropertyNode [] path = jd.getPath();
	for (int i = 0; i < path.length; i++) {
	  PropertyNode pn = path[i];
	  System.err.println("" + (i + 1) + "  " + pn.toString()
			     + " " + pn.value.toString());
	}
      } else {
	System.err.println("Cancelled");
      }
      System.exit(0);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
