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
 *    ClustererBeanInfo.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.beans.*;

/**
 * BeanInfo class for the Clusterer wrapper bean
 *
 * @author <a href="mailto:mutter@cs.waikato.ac.nz">Stefan Mutter</a>
 * @version $Revision: 5247 $
 */
public class ClustererBeanInfo extends SimpleBeanInfo {
 
  public EventSetDescriptor [] getEventSetDescriptors() {
    try {
      EventSetDescriptor [] esds = { 
	new EventSetDescriptor(Clusterer.class, 
			       "batchClusterer", 
			       BatchClustererListener.class, 
			       "acceptClusterer"),
        new EventSetDescriptor(Clusterer.class,
			       "graph",
			       GraphListener.class,
			       "acceptGraph"),                       
	new EventSetDescriptor(Clusterer.class,
			       "text",
			       TextListener.class,
			       "acceptText"),
	new EventSetDescriptor(Clusterer.class,
	                       "configuration",
	                       ConfigurationListener.class,
	                       "acceptConfiguration")
	};
      return esds;
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    return null;
  }

  /**
   * Get the bean descriptor for this bean
   *
   * @return a <code>BeanDescriptor</code> value
   */
  public BeanDescriptor getBeanDescriptor() {
    return new BeanDescriptor(weka.gui.beans.Clusterer.class, 
			      ClustererCustomizer.class);
  }
}

