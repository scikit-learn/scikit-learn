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
 *    FilterCustomizer.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.filters.Filter;
import weka.gui.GenericObjectEditor;
import weka.gui.PropertySheetPanel;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.Customizer;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * GUI customizer for the filter bean
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 5340 $
 */
public class FilterCustomizer
  extends JPanel
  implements Customizer, CustomizerCloseRequester {

  /** for serialization */
  private static final long serialVersionUID = 2049895469240109738L;
  
  static {
     GenericObjectEditor.registerEditors();
  }

  private PropertyChangeSupport m_pcSupport = 
    new PropertyChangeSupport(this);

  private weka.gui.beans.Filter m_filter;
/*  private GenericObjectEditor m_filterEditor = 
    new GenericObjectEditor(true); */
  
  /** Backup if user presses cancel */
  private weka.filters.Filter m_backup;
  
  private PropertySheetPanel m_filterEditor = 
    new PropertySheetPanel();
  
  private JFrame m_parentFrame;
 
  public FilterCustomizer() {
    m_filterEditor.
    setBorder(BorderFactory.createTitledBorder("Filter options"));



    setLayout(new BorderLayout());
    add(m_filterEditor, BorderLayout.CENTER);

    JPanel butHolder = new JPanel();
    butHolder.setLayout(new GridLayout(1,2));
    JButton OKBut = new JButton("OK");
    OKBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        m_parentFrame.dispose();
      }
    });

    JButton CancelBut = new JButton("Cancel");
    CancelBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        // cancel requested, so revert to backup and then
        // close the dialog
        if (m_backup != null) {
          m_filter.setFilter(m_backup);
        }
        m_parentFrame.dispose();
      }
    });
    
    butHolder.add(OKBut);
    butHolder.add(CancelBut);
    add(butHolder, BorderLayout.SOUTH);
  }
  
  /**
   * Set the filter bean to be edited
   *
   * @param object a Filter bean
   */
  public void setObject(Object object) {
    m_filter = (weka.gui.beans.Filter)object;
    try {
      m_backup = 
        (weka.filters.Filter)GenericObjectEditor.makeCopy(m_filter.getFilter());
    } catch (Exception ex) {
      // ignore
    }
    m_filterEditor.setTarget(m_filter.getFilter());
  }

  /**
   * Add a property change listener
   *
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void addPropertyChangeListener(PropertyChangeListener pcl) {
    m_pcSupport.addPropertyChangeListener(pcl);
  }

  /**
   * Remove a property change listener
   *
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void removePropertyChangeListener(PropertyChangeListener pcl) {
    m_pcSupport.removePropertyChangeListener(pcl);
  }

  public void setParentFrame(JFrame parent) {
    m_parentFrame = parent;
  }
}

