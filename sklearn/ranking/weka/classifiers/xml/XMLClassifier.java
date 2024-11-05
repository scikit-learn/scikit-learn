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
 * XMLClassifier.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.xml;

import weka.core.RevisionUtils;
import weka.core.xml.XMLBasicSerialization;

/**
 * This class serializes and deserializes a Classifier instance to and
 * fro XML.<br>
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.6 $ 
 */
public class XMLClassifier
   extends XMLBasicSerialization {

   /**
    * initializes the serialization
    * 
    * @throws Exception if initialization fails
    */
   public XMLClassifier() throws Exception {
      super();
   }

   /**
    * generates internally a new XML document and clears also the IgnoreList and
    * the mappings for the Read/Write-Methods
    */
   public void clear() throws Exception {
      super.clear();

      // allow
      m_Properties.addAllowed(weka.classifiers.Classifier.class, "debug");
      m_Properties.addAllowed(weka.classifiers.Classifier.class, "options");
   }
   
   /**
    * Returns the revision string.
    * 
    * @return		the revision
    */
   public String getRevision() {
     return RevisionUtils.extract("$Revision: 1.6 $");
   }
}
