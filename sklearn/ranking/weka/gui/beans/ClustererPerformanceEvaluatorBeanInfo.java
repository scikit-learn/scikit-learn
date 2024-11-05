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
 *    ClustererPerformanceEvaluatorBeanInfo.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.beans.*;

/**
 * Bean info class for the clusterer performance evaluator
 *
 * @author <a href="mailto:mutter@cs.waikato.ac.nz">Stefan Mutter</a>
 * @version $Revision: 1.2 $
 */
public class ClustererPerformanceEvaluatorBeanInfo extends SimpleBeanInfo {
  
  public EventSetDescriptor [] getEventSetDescriptors() {
    try {
      EventSetDescriptor [] esds = { 
	new EventSetDescriptor(ClustererPerformanceEvaluator.class, 
			       "text", 
			       TextListener.class, 
			       "acceptText")};
	
      return esds;
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    return null;
  }
}
