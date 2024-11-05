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
 * InfoPanel.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.sql;

import weka.gui.ComponentHelper;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

/**
 * A simple panel for displaying information, e.g. progress information etc.
 *
 * @author      FracPete (fracpete at waikato dot ac dot nz)
 * @version     $Revision: 1.3 $
 */

public class InfoPanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -7701133696481997973L;
  
  /** the parent of this panel */
  protected JFrame m_Parent;
  
  /** the list the contains the messages */
  protected JList m_Info;

  /** the model for the list */
  protected DefaultListModel m_Model;

  /** the button to clear the area */
  protected JButton m_ButtonClear;

  /** the button to copy the selected message */
  protected JButton m_ButtonCopy;
  
  /**
   * creates the panel
   * @param parent      the parent of this panel
   */
  public InfoPanel(JFrame parent) {
    super();
    m_Parent = parent;
    createPanel();
  }

  /**
   * inserts the components into the panel
   */
  protected void createPanel() {
    JPanel          panel;
    JPanel          panel2;
    
    setLayout(new BorderLayout());
    setPreferredSize(new Dimension(0, 80));

    // text
    m_Model = new DefaultListModel();
    m_Info  = new JList(m_Model);
    m_Info.setCellRenderer(new InfoPanelCellRenderer());
    m_Info.addListSelectionListener(new ListSelectionListener() {
      public void valueChanged(ListSelectionEvent e) {
        setButtons(e);
      }
    });
    add(new JScrollPane(m_Info), BorderLayout.CENTER);

    // clear button
    panel = new JPanel(new BorderLayout());
    add(panel, BorderLayout.EAST);
    m_ButtonClear = new JButton("Clear");
    m_ButtonClear.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  clear();
	}
      });
    panel.add(m_ButtonClear, BorderLayout.NORTH);

    // clear button
    panel2 = new JPanel(new BorderLayout());
    panel.add(panel2, BorderLayout.CENTER);
    m_ButtonCopy = new JButton("Copy");
    m_ButtonCopy.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  copyToClipboard();
	}
      });
    panel2.add(m_ButtonCopy, BorderLayout.NORTH);
  }
  
  /**
   * sets the state of the buttons according to the selection state of the
   * JList
   */
  protected void setButtons(ListSelectionEvent e) {
    if ( (e == null) || (e.getSource() == m_Info) ) {
      m_ButtonClear.setEnabled(m_Model.getSize() > 0);
      m_ButtonCopy.setEnabled(m_Info.getSelectedIndices().length == 1);
    }
  }

  /**
   * sets the focus in a designated control
   */
  public void setFocus() {
    m_Info.requestFocus();
  }

  /**
   * clears the content of the panel
   */
  public void clear() {
    m_Model.clear();
    setButtons(null);
  }

  /**
   * copies the currently selected error message to the clipboard
   * 
   * @return		true if the content was copied
   */
  public boolean copyToClipboard() {
    StringSelection      selection;
    Clipboard            clipboard;
    
    if (m_Info.getSelectedIndices().length != 1)
      return false;
    
    selection = new StringSelection(((JLabel) m_Info.getSelectedValue()).getText());
    clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
    clipboard.setContents(selection, selection);
    return true;
  }
  
  /**
   * adds the given message to the end of the list (with the associated icon
   * at the beginning)
   * @param msg       the message to append to the list
   * @param icon      the filename of the icon
   */
  public void append(String msg, String icon) {
    append(new JLabel(msg, ComponentHelper.getImageIcon(icon), JLabel.LEFT));
  }

  /**
   * adds the given message to the end of the list
   * @param msg       the message to append to the list
   */
  public void append(Object msg) {
    if (msg instanceof String) {
      append(msg.toString(), "empty_small.gif");
      return;
    }

    m_Model.addElement(msg);
    m_Info.setSelectedIndex(m_Model.getSize() - 1);
    m_Info.ensureIndexIsVisible(m_Info.getSelectedIndex());
    
    setButtons(null);
  }
}

