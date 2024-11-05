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
 *    ClassValuePickerBeanInfo.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.beans.*;

/**
 * BeanInfo class for the class value picker bean
 *
 * @author Mark Hall
 * @version $Revision: 1.3 $
 */
public class ClassValuePickerBeanInfo 
  extends SimpleBeanInfo {

  /**
   * Returns the event set descriptors
   *
   * @return an <code>EventSetDescriptor[]</code> value
   */
  public EventSetDescriptor [] getEventSetDescriptors() {
    try {
      EventSetDescriptor [] esds = 
      { new EventSetDescriptor(ClassValuePicker.class, 
			       "dataSet", 
			       DataSourceListener.class, 
			       "acceptDataSet") };
      return esds;
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    return null;
  }
  
  /**
   * Returns the property descriptors
   *
   * @return a <code>PropertyDescriptor[]</code> value
   */
  public PropertyDescriptor[] getPropertyDescriptors() {
    try {
      PropertyDescriptor p1;
      p1 = new PropertyDescriptor("classValueIndex", ClassValuePicker.class);
      PropertyDescriptor [] pds = { p1 };
      return pds;
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    return null;
  }

  public BeanDescriptor getBeanDescriptor() {
    return new BeanDescriptor(weka.gui.beans.ClassValuePicker.class,
			      ClassValuePickerCustomizer.class);
  }
}
