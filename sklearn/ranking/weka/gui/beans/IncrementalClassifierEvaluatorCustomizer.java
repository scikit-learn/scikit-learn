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
 *    IncrementalClassifierEvaluatorCustomizer.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.gui.PropertySheetPanel;

import java.awt.BorderLayout;
import java.beans.Customizer;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;

import javax.swing.JPanel;

/**
 * GUI Customizer for the incremental classifier evaluator bean
 *
 * @author Mark Hall (mhall{[at}]pentaho{[dot]}org
 * @version $Revision: 5563 $
 */
public class IncrementalClassifierEvaluatorCustomizer
  extends JPanel
  implements Customizer {

  /** for serialization */
  //  private static final long serialVersionUID = 1229878140258668581L;

  private PropertyChangeSupport m_pcSupport = 
    new PropertyChangeSupport(this);

  private PropertySheetPanel m_ieEditor = 
    new PropertySheetPanel();

  public IncrementalClassifierEvaluatorCustomizer() {
    setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 5, 5, 5));

    setLayout(new BorderLayout());
    add(m_ieEditor, BorderLayout.CENTER);
    add(new javax.swing.JLabel("IncrementalClassifierEvaluatorCustomizer"), 
	BorderLayout.NORTH);
  }
  
  /**
   * Set the object to be edited
   *
   * @param object a IncrementalClassifierEvaluator object
   */
  public void setObject(Object object) {
    m_ieEditor.setTarget((IncrementalClassifierEvaluator)object);
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
