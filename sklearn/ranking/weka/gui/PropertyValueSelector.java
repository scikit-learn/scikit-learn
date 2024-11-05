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
 *    PropertyValueSelector.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.beans.PropertyEditor;

import javax.swing.ComboBoxModel;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;

/** 
 * Support for any PropertyEditor that uses tags.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.7 $
 */
class PropertyValueSelector
  extends JComboBox {

  /** for serialization */
  private static final long serialVersionUID = 128041237745933212L;

  /** The property editor */
  PropertyEditor m_Editor;
  
  /**
   * Sets up the editing component with the supplied editor.
   *
   * @param pe the PropertyEditor
   */
  public PropertyValueSelector(PropertyEditor pe) {
      
    m_Editor = pe;
    Object value = m_Editor.getAsText();
    String tags[] = m_Editor.getTags();
    ComboBoxModel model = new DefaultComboBoxModel(tags) {
      private static final long serialVersionUID = 7942587653040180213L;
      
      public Object getSelectedItem() {
	return m_Editor.getAsText();
      }
      public void setSelectedItem(Object o) {
	m_Editor.setAsText((String)o);
      }
    };
    setModel(model);
    setSelectedItem(value);
  }
}


