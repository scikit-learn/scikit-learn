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
 * InfoPanelCellRenderer.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.sql;

import java.awt.Component;

import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.ListCellRenderer;

/**
 * A specialized renderer that takes care of JLabels in a JList.
 *
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.2 $
 */

public class InfoPanelCellRenderer 
  extends JLabel 
  implements ListCellRenderer {

  /** for serialization */
  private static final long serialVersionUID = -533380118807178531L;
  
  /**
   * the constructor
   */
  public InfoPanelCellRenderer() {
    super();
    setOpaque(true);
  }
  
  /**
   * Return a component that has been configured to display the specified value.
   * @param list The JList we're painting.
   * @param value The value returned by list.getModel().getElementAt(index).
   * @param index The cells index.
   * @param isSelected True if the specified cell was selected.
   * @param cellHasFocus True if the specified cell has the focus.
   */
  public Component getListCellRendererComponent(
      JList list, Object value,
      int index, boolean isSelected, boolean cellHasFocus) {

    if (value instanceof JLabel) {
      setIcon(((JLabel) value).getIcon());
      setText(((JLabel) value).getText());
    }
    else {
      setIcon(null);
      setText(value.toString());
    }

    if (isSelected) {
      setBackground(list.getSelectionBackground());
      setForeground(list.getSelectionForeground());
    }
    else {
      setBackground(list.getBackground());
      setForeground(list.getForeground());
    }
    setEnabled(list.isEnabled());
    setFont(list.getFont());

    return this;
  }
}

