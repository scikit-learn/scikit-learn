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
 *    PropertyText.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.beans.PropertyEditor;

import javax.swing.JTextField;

/** 
 * Support for a PropertyEditor that uses text.
 * Isn't going to work well if the property gets changed
 * somewhere other than this field simultaneously
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
class PropertyText
  extends JTextField {

  /** for serialization */
  private static final long serialVersionUID = -3915342928825822730L;

  /** The property editor */
  private PropertyEditor m_Editor;

  /**
   * Sets up the editing component with the supplied editor.
   *
   * @param pe the PropertyEditor
   */
  PropertyText(PropertyEditor pe) {
 
    //super(pe.getAsText());
    super((pe.getAsText().equals("null"))?"":pe.getAsText());
    m_Editor = pe;
    
    /*    m_Editor.addPropertyChangeListener(new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent evt) {
	updateUs();
      }
      }); */
    addKeyListener(new KeyAdapter() {
      public void keyReleased(KeyEvent e) {
	//	if (e.getKeyCode() == KeyEvent.VK_ENTER) {
	updateEditor();
	//	}
      }
    });
    addFocusListener(new FocusAdapter() {
      public void focusLost(FocusEvent e) {
	updateEditor();
      }
    });
  }

  /**
   * Attempts to update the textfield value from the editor.
   */
  protected void updateUs() {
    try {
      setText(m_Editor.getAsText());
    } catch (IllegalArgumentException ex) {
      // Quietly ignore.
    }
  }

  /**
   * Attempts to update the editor value from the textfield.
   */
  protected void updateEditor() {
    try {
      m_Editor.setAsText(getText());
    } catch (IllegalArgumentException ex) {
      // Quietly ignore.
    }
  }
}
