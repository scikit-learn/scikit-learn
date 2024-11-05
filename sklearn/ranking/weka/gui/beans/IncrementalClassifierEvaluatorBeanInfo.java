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
 *    IncrementalClassifierEvaluatorBeanInfo.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.beans.*;

/**
 * Bean info class for the incremental classifier evaluator bean
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.4 $
 */
public class IncrementalClassifierEvaluatorBeanInfo extends SimpleBeanInfo {
  
  /**
   * Return the property descriptors for this bean
   *
   * @return a <code>PropertyDescriptor[]</code> value
   */
  public PropertyDescriptor[] getPropertyDescriptors() {
    try {
      PropertyDescriptor p1;
      PropertyDescriptor p2;
      p1 = new PropertyDescriptor("statusFrequency", IncrementalClassifierEvaluator.class);
      p2 = new PropertyDescriptor("outputPerClassInfoRetrievalStats", 
                                  IncrementalClassifierEvaluator.class);
      PropertyDescriptor [] pds = { p1, p2 };
      return pds;
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    return null;
  }

  /**
   * Get the event set descriptors for this bean
   *
   * @return an <code>EventSetDescriptor[]</code> value
   */
  public EventSetDescriptor [] getEventSetDescriptors() {
    try {
      EventSetDescriptor [] esds = { 
	new EventSetDescriptor(IncrementalClassifierEvaluator.class, 
			       "chart", 
			       ChartListener.class, 
			       "acceptDataPoint"),
	new EventSetDescriptor(IncrementalClassifierEvaluator.class,
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

  /**
   * Return the bean descriptor for this bean
   *
   * @return a <code>BeanDescriptor</code> value
   */
  public BeanDescriptor getBeanDescriptor() {
    return new BeanDescriptor(weka.gui.beans.IncrementalClassifierEvaluator.class,
			      IncrementalClassifierEvaluatorCustomizer.class);
  }
}
