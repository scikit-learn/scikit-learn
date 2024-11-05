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
 *    PropertyPanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui;

import weka.core.OptionHandler;
import weka.core.Utils;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyEditor;

import javax.swing.BorderFactory;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;

/** 
 * Support for drawing a property value in a component.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 5841 $
 */
public class PropertyPanel 
  extends JPanel {
  
  /** for serialization */
  static final long serialVersionUID = 5370025273466728904L;

  /** The property editor */
  private PropertyEditor m_Editor;

  /** The currently displayed property dialog, if any */
  private PropertyDialog m_PD;

  /** Whether the editor has provided its own panel */
  private boolean m_HasCustomPanel = false;
  
  /** The custom panel (if any) */
  private JPanel m_CustomPanel;
  
  /**
   * Create the panel with the supplied property editor.
   *
   * @param pe the PropertyEditor
   */
  public PropertyPanel(PropertyEditor pe) {

    this(pe, false);
  }

  /**
   * Create the panel with the supplied property editor,
   * optionally ignoring any custom panel the editor can provide.
   *
   * @param pe the PropertyEditor
   * @param ignoreCustomPanel whether to make use of any available custom panel
   */
  public PropertyPanel(PropertyEditor pe, boolean ignoreCustomPanel) {

    m_Editor = pe;
    
    if (!ignoreCustomPanel && m_Editor instanceof CustomPanelSupplier) {
      setLayout(new BorderLayout());
      m_CustomPanel = ((CustomPanelSupplier)m_Editor).getCustomPanel();
      add(m_CustomPanel, BorderLayout.CENTER);
      m_HasCustomPanel = true;
    } else {
      createDefaultPanel();
    }
  }

  /**
   * Creates the default style of panel for editors that do not
   * supply their own.
   */
  protected void createDefaultPanel() {

    setBorder(BorderFactory.createEtchedBorder());
    setToolTipText("Left-click to edit properties for this object, right-click/Alt+Shift+left-click for menu");
    setOpaque(true);
    final Component comp = this;
    addMouseListener(new MouseAdapter() {
      public void mouseClicked(MouseEvent evt) {
        if (evt.getClickCount() == 1) {
          if (    (evt.getButton() == MouseEvent.BUTTON1) && !evt.isAltDown() && !evt.isShiftDown() ) {
            showPropertyDialog();
          }
          else if (    (evt.getButton() == MouseEvent.BUTTON3) 
              	    || ((evt.getButton() == MouseEvent.BUTTON1) && evt.isAltDown() && evt.isShiftDown()) ) {
            JPopupMenu menu = new JPopupMenu();
            JMenuItem item;

            if (m_Editor.getValue() != null) {
              item = new JMenuItem("Show properties...");
              item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                  showPropertyDialog();
                }
              });
              menu.add(item);

              item = new JMenuItem("Copy configuration to clipboard");
              item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                  String str = m_Editor.getValue().getClass().getName();
                  if (m_Editor.getValue() instanceof OptionHandler)
                    str += " " + Utils.joinOptions(((OptionHandler) m_Editor.getValue()).getOptions());
                  StringSelection selection = new StringSelection(str.trim());
                  Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                  clipboard.setContents(selection, selection);
                }
              });
              menu.add(item);
            }
            
            item = new JMenuItem("Enter configuration...");
            item.addActionListener(new ActionListener() {
              public void actionPerformed(ActionEvent e) {
        	String str = JOptionPane.showInputDialog(
        	                 comp, 
        	                 "Configuration (<classname> [<options>])");
        	if (str != null) {
        	  try {
        	    String[] options = Utils.splitOptions(str);
        	    String classname = options[0];
        	    options[0] = "";
        	    m_Editor.setValue(
        		Utils.forName(
        		    Object.class, classname, options));
        	  }
        	  catch (Exception ex) {
        	    ex.printStackTrace();
        	    JOptionPane.showMessageDialog(
        		comp, 
        		"Error parsing commandline:\n" + ex, 
        		"Error...",
        		JOptionPane.ERROR_MESSAGE);
        	  }
        	}
              }
            });
            menu.add(item);
            
            menu.show(comp, evt.getX(), evt.getY());
          }
        }
      }
    });
    Dimension newPref = getPreferredSize();
    newPref.height = getFontMetrics(getFont()).getHeight() * 5 / 4;
    newPref.width = newPref.height * 5;
    setPreferredSize(newPref);

    m_Editor.addPropertyChangeListener(new PropertyChangeListener () {
	public void propertyChange(PropertyChangeEvent evt) {
	  repaint();
	}
      });
  }

  /**
   * Displays the property edit dialog for the panel.
   */
  public void showPropertyDialog() {

    if (m_Editor.getValue() != null) {
      if (m_PD == null) {
	int x = getLocationOnScreen().x;
	int y = getLocationOnScreen().y;
	if (PropertyDialog.getParentDialog(this) != null)
	  m_PD = new PropertyDialog(PropertyDialog.getParentDialog(this), m_Editor, x, y);
	else
	  m_PD = new PropertyDialog(PropertyDialog.getParentFrame(this), m_Editor, x, y);
	m_PD.setVisible(true);
      } else {
	m_PD.setVisible(true);
      }
      // make sure that m_Backup is correctly initialized!
      m_Editor.setValue(m_Editor.getValue());
    }
  }

  /**
   * Cleans up when the panel is destroyed.
   */
  public void removeNotify() {

    super.removeNotify();
    if (m_PD != null) {
      m_PD.dispose();
      m_PD = null;
    }
  }
  
  /**
   * Passes on enabled/disabled status to the custom
   * panel (if one is set).
   * 
   * @param enabled true if this panel (and the custom panel is enabled)
   */
  public void setEnabled(boolean enabled) {
    super.setEnabled(enabled);
    if (m_HasCustomPanel) {
      m_CustomPanel.setEnabled(enabled);
    }
    
  }

  /**
   * Paints the component, using the property editor's paint method.
   *
   * @param g the current graphics context
   */
  public void paintComponent(Graphics g) {

    if (!m_HasCustomPanel) {
      Insets i = getInsets();
      Rectangle box = new Rectangle(i.left, i.top,
				    getSize().width - i.left - i.right - 1,
				    getSize().height - i.top - i.bottom - 1);
      
      g.clearRect(i.left, i.top,
		  getSize().width - i.right - i.left,
		  getSize().height - i.bottom - i.top);
      m_Editor.paintValue(g, box);
    }
  }

}
