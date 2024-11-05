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
 *    DataSetEvent.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instances;

import java.util.EventObject;

/**
 * Event encapsulating a data set
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.4 $
 * @see EventObject
 */
public class DataSetEvent
  extends EventObject {

  /** for serialization */
  private static final long serialVersionUID = -5111218447577318057L;

  private Instances m_dataSet;
  private boolean m_structureOnly;

  public DataSetEvent(Object source, Instances dataSet) {
    super(source);
    m_dataSet = dataSet;
    if (m_dataSet != null && m_dataSet.numInstances() == 0) {
      m_structureOnly = true;
    }
  }
  
  /**
   * Return the instances of the data set
   *
   * @return an <code>Instances</code> value
   */
  public Instances getDataSet() {
    return m_dataSet;
  }

  /**
   * Returns true if the encapsulated instances
   * contain just header information
   *
   * @return true if only header information is
   * available in this DataSetEvent
   */
  public boolean isStructureOnly() {
    return m_structureOnly;
  }
}
