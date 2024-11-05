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
 *    PropertyDialog.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dialog;
import java.awt.Frame;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyEditor;

import javax.swing.JDialog;
import javax.swing.JInternalFrame;

/** 
 * Support for PropertyEditors with custom editors: puts the editor into
 * a separate frame.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 6596 $
 */
public class PropertyDialog
  extends JDialog {

  /** for serialization. */
  private static final long serialVersionUID = -2314850859392433539L;

  /** The property editor. */
  private PropertyEditor m_Editor;

  /** The custom editor component. */
  private Component m_EditorComponent;
  
  /**
   * Creates the editor frame - only kept for backward-compatibility.
   *
   * @param pe 		the PropertyEditor
   * @param x 		initial x coord for the frame
   * @param y 		initial y coord for the frame
   * @deprecated 	instead of this constructor, one should use the constructors
   * 			with an explicit owner (either derived from 
   * 			<code>java.awt.Dialog</code> or from 
   * 			<code>java.awt.Frame</code>) or, if none available,
   * 			using <code>(Frame) null</code> as owner.
   */
  public PropertyDialog(PropertyEditor pe, int x, int y) {
    this((Frame) null, pe, x, y);
    setVisible(true);
  }
  
  /**
   * Creates the (screen-centered) editor dialog. The dialog is automatically
   * modal in case the owner is non-null.
   *
   * @param owner	the dialog that opens this dialog
   * @param pe 		the PropertyEditor
   */
  public PropertyDialog(Dialog owner, PropertyEditor pe) {
    this(owner, pe, -1, -1);
  }
  
  /**
   * Creates the editor dialog at the given position. The dialog is automatically
   * modal in case the owner is non-null.
   *
   * @param owner	the dialog that opens this dialog
   * @param pe 		the PropertyEditor
   * @param x 		initial x coord for the dialog
   * @param y 		initial y coord for the dialog
   */
  public PropertyDialog(Dialog owner, PropertyEditor pe, int x, int y) {
    super(owner, pe.getClass().getName(), true);
    initialize(pe, x, y);
  }
  
  /**
   * Creates the (screen-centered) editor dialog. The dialog is automatically
   * modal in case the owner is non-null.
   *
   * @param owner	the frame that opens this dialog
   * @param pe 		the PropertyEditor
   */
  public PropertyDialog(Frame owner, PropertyEditor pe) {
    this(owner, pe, -1, -1);
  }
  
  /**
   * Creates the editor dialog at the given position. The dialog is automatically
   * modal in case the owner is non-null.
   *
   * @param owner	the frame that opens this dialog
   * @param pe 		the PropertyEditor
   * @param x 		initial x coord for the dialog
   * @param y 		initial y coord for the dialog
   */
  public PropertyDialog(Frame owner, PropertyEditor pe, int x, int y) {
    super(owner, pe.getClass().getName(), true);
    
    initialize(pe, x, y);
  }
  
  /**
   * Initializes the dialog.
   *
   * @param pe 		the PropertyEditor
   * @param x 		initial x coord for the dialog
   * @param y 		initial y coord for the dialog
   */
  protected void initialize(PropertyEditor pe, int x, int y) {
    addWindowListener(new WindowAdapter() {
      public void windowClosing(WindowEvent e) {
	e.getWindow().dispose();
      }
    });
    getContentPane().setLayout(new BorderLayout());

    m_Editor = pe;
    m_EditorComponent = pe.getCustomEditor();
    getContentPane().add(m_EditorComponent, BorderLayout.CENTER);

    pack();
    
    int screenWidth = getGraphicsConfiguration().getBounds().width;
    int screenHeight = getGraphicsConfiguration().getBounds().height;

    // adjust height to a maximum of 95% of screen height
    if (getHeight() > (double) screenHeight * 0.95)
      setSize(getWidth(), (int) ((double) screenHeight * 0.95));
    
    if ((x == -1) && (y == -1)) {
      setLocationRelativeTo(null);
    }
    else {
      // adjust position if necessary
      if (x + getWidth() > screenWidth)
	x = screenWidth - getWidth();
      if (y + getHeight() > screenHeight)
	y = screenHeight - getHeight();
      setLocation(x, y);
    }
  }

  /**
   * Gets the current property editor.
   *
   * @return a value of type 'PropertyEditor'
   */
  public PropertyEditor getEditor() {
    return m_Editor;
  }

  /**
   * Tries to determine the frame this panel is part of.
   * 
   * @param c		the container to start with
   * @return		the parent frame if one exists or null if not
   */
  public static Frame getParentFrame(Container c) {
    Frame	result;
    Container	parent;
    
    result = null;
    
    parent = c;
    while (parent != null) {
      if (parent instanceof Frame) {
	result = (Frame) parent;
	break;
      }
      else {
	parent = parent.getParent();
      }
    }
    
    return result;
  }

  /**
   * Tries to determine the internal frame this panel is part of.
   * 
   * @param c		the container to start with
   * @return		the parent internal frame if one exists or null if not
   */
  public static JInternalFrame getParentInternalFrame(Container c) {
    JInternalFrame	result;
    Container		parent;
    
    result = null;
    
    parent = c;
    while (parent != null) {
      if (parent instanceof JInternalFrame) {
	result = (JInternalFrame) parent;
	break;
      }
      else {
	parent = parent.getParent();
      }
    }
    
    return result;
  }

  /**
   * Tries to determine the dialog this panel is part of.
   * 
   * @param c		the container to start with
   * @return		the parent dialog if one exists or null if not
   */
  public static Dialog getParentDialog(Container c) {
    Dialog	result;
    Container	parent;
    
    result = null;
    
    parent = c;
    while (parent != null) {
      if (parent instanceof Dialog) {
	result = (Dialog) parent;
	break;
      }
      else {
	parent = parent.getParent();
      }
    }
    
    return result;
  }
}

