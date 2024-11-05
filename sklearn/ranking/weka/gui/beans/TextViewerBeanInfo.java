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
 *    TextViewerBeanInfo.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.beans.*;

/**
 * Bean info class for the text viewer
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.3 $
 */
public class TextViewerBeanInfo extends SimpleBeanInfo {
  
  /**
   * Get the event set descriptors for this bean
   *
   * @return an <code>EventSetDescriptor[]</code> value
   */
  public EventSetDescriptor [] getEventSetDescriptors() {
    try {
      EventSetDescriptor [] esds = { 
        new EventSetDescriptor(TextViewer.class,
                               "text",
                               TextListener.class,
                               "acceptText")
      };      
      return esds;
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    return null;
  }
}
