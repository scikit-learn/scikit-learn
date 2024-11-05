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
 *    ClassValuePickerCustomizer.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instances;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.Customizer;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;

import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * @author Mark Hall
 * @version $Revision: 1.6 $
 */
public class ClassValuePickerCustomizer
  extends JPanel
  implements Customizer, CustomizerClosingListener, DataFormatListener {

  /** for serialization */
  private static final long serialVersionUID = 8213423053861600469L;

  private boolean m_displayValNames = false;

  private ClassValuePicker m_classValuePicker;
  
  private PropertyChangeSupport m_pcSupport = 
    new PropertyChangeSupport(this);

  private JComboBox m_ClassValueCombo = new JComboBox();
  private JPanel m_holderP = new JPanel();

  private JLabel m_messageLabel = new JLabel("No customization possible at present.");

  public ClassValuePickerCustomizer() {
    setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 5, 5, 5));
    
    setLayout(new BorderLayout());
    add(new javax.swing.JLabel("ClassValuePickerCustomizer"), 
	BorderLayout.NORTH);
    m_holderP.setLayout(new BorderLayout());
    m_holderP.setBorder(BorderFactory.createTitledBorder("Choose class value"));
    m_holderP.add(m_ClassValueCombo, BorderLayout.CENTER);
    m_ClassValueCombo.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  if (m_classValuePicker != null) {
	    m_classValuePicker.setClassValueIndex(m_ClassValueCombo.getSelectedIndex());
	  }
	}
      });

    add(m_messageLabel, BorderLayout.CENTER);
  }

  private void setUpNoCustPossible() {
    if (m_displayValNames == true) {
      remove(m_holderP);
      add(m_messageLabel, BorderLayout.CENTER);
      m_displayValNames = false;
    }
    validate(); repaint();
  }

  private void setUpValueSelection(Instances format) {
    if (format.classIndex() < 0 || format.classAttribute().isNumeric()) {
      // cant do anything in this case
      return;
    }
    if (m_displayValNames == false) {
      remove(m_messageLabel);
    }
    
    int existingClassVal = m_classValuePicker.getClassValueIndex();
    String [] attribValNames = new String [format.classAttribute().numValues()];
    for (int i = 0; i < attribValNames.length; i++) {
      attribValNames[i] = format.classAttribute().value(i);
    }
    m_ClassValueCombo.setModel(new DefaultComboBoxModel(attribValNames));
    if (attribValNames.length > 0) {
      if (existingClassVal < attribValNames.length) {
	m_ClassValueCombo.setSelectedIndex(existingClassVal);
      }
    }
    if (m_displayValNames == false) {
      add(m_holderP, BorderLayout.CENTER);
      m_displayValNames = true;
    }
    validate(); repaint();
  }

  /**
   * Set the bean to be edited
   *
   * @param object an <code>Object</code> value
   */
  public void setObject(Object object) {
    if (m_classValuePicker != (ClassValuePicker)object) {
      // remove ourselves as a listener from the old ClassvaluePicker (if necessary)
      if (m_classValuePicker != null) {
	m_classValuePicker.removeDataFormatListener(this);
      }
      m_classValuePicker = (ClassValuePicker)object;
      // add ourselves as a data format listener
      m_classValuePicker.addDataFormatListener(this);
      if (m_classValuePicker.getConnectedFormat() != null) {
	setUpValueSelection(m_classValuePicker.getConnectedFormat());	
      }
    }
  }
  
  public void customizerClosing() {
    // remove ourselves as a listener from the ClassValuePicker (if necessary)
    if (m_classValuePicker != null) {
      System.err.println("Customizer deregistering with class value picker");
      m_classValuePicker.removeDataFormatListener(this);
    }
  }

  public void newDataFormat(DataSetEvent dse) {
    if (dse.getDataSet() != null) {
      setUpValueSelection(m_classValuePicker.getConnectedFormat());
    } else {
      setUpNoCustPossible();
    }
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
}
