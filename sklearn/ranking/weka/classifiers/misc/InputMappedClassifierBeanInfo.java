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
 * Copyright 2010 University of Waikato
 */

package weka.classifiers.misc;

import java.beans.PropertyDescriptor;
import java.beans.SimpleBeanInfo;
import java.util.ArrayList;

/**
 * Bean info class for the InputMappedClassifier. Ensures that the 
 * FileEnvironmentField class is used as the custom property editor
 * in the GOE for the modelPath property.
 * 
 * @author Mar Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6803 $
 *
 */
public class InputMappedClassifierBeanInfo extends SimpleBeanInfo {
  
  /**
   * Get an array of PropertyDescriptors for the InputMappedClassifier's
   * public properties.
   * 
   * @return an array of PropertyDescriptors
   */
  public PropertyDescriptor[] getPropertyDescriptors() {
    try {
      PropertyDescriptor p1;
      ArrayList<PropertyDescriptor> pds = new ArrayList<PropertyDescriptor>();
      
      p1 = new PropertyDescriptor("modelPath", InputMappedClassifier.class);
      p1.setPropertyEditorClass(weka.gui.beans.FileEnvironmentField.class);
      pds.add(p1);
      
      pds.add(new PropertyDescriptor("ignoreCaseForNames", InputMappedClassifier.class));
      pds.add(new PropertyDescriptor("suppressMappingReport", InputMappedClassifier.class));
      pds.add(new PropertyDescriptor("trim", InputMappedClassifier.class));
      pds.add(new PropertyDescriptor("classifier", InputMappedClassifier.class));
      
      // this one is only really needed for XMLSerialization
      pds.add(new PropertyDescriptor("options", InputMappedClassifier.class));
      
      return pds.toArray(new PropertyDescriptor[1]);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    return null;
  }
}
