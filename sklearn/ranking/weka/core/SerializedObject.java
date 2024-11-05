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
 *    SerializedObject.java
 *    Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import weka.core.scripting.Jython;
import weka.core.scripting.JythonSerializableObject;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Class for storing an object in serialized form in memory. It can be used 
 * to make deep copies of objects, and also allows compression to conserve
 * memory. <p>
 *
 * @author Richard Kirkby (rbk1@cs.waikato.ac.nz)
 * @version $Revision: 5018 $ 
 */
public class SerializedObject
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 6635502953928860434L;

  /** The array storing the object. */
  private byte[] m_storedObjectArray;

  /** Whether or not the object is compressed. */
  private boolean m_isCompressed;

  /** Whether it is a Jython object or not */
  private boolean m_isJython;
  
  /**
   * Creates a new serialized object (without compression).
   *
   * @param toStore the object to store
   * @exception Exception if the object couldn't be serialized
   */ 
  public SerializedObject(Object toStore) throws Exception {

    this(toStore, false);
  }

  /**
   * Creates a new serialized object.
   *
   * @param toStore the object to store
   * @param compress whether or not to use compression
   * @exception Exception if the object couldn't be serialized
   */ 
  public SerializedObject(Object toStore, boolean compress) throws Exception {

    ByteArrayOutputStream ostream = new ByteArrayOutputStream();
    OutputStream os = ostream;
    ObjectOutputStream p;
    if (!compress)
      p = new ObjectOutputStream(new BufferedOutputStream(os));
    else
      p = new ObjectOutputStream(new BufferedOutputStream(new GZIPOutputStream(os)));
    p.writeObject(toStore);
    p.flush();
    p.close(); // used to be ostream.close() !
    m_storedObjectArray = ostream.toByteArray();

    m_isCompressed = compress;
    m_isJython     = (toStore instanceof JythonSerializableObject);
  }

  /*
   * Checks to see whether this object is equal to another.
   *
   * @param compareTo the object to compare to
   * @return whether or not the objects are equal
   */
  public final boolean equals(Object compareTo) {

    if (compareTo == null) return false;
    if (!compareTo.getClass().equals(this.getClass())) return false;
    byte[] compareArray = ((SerializedObject)compareTo).m_storedObjectArray;
    if (compareArray.length != m_storedObjectArray.length) return false;
    for (int i=0; i<compareArray.length; i++) {
      if (compareArray[i] != m_storedObjectArray[i]) return false;
    }
    return true;
  }

  /**
   * Returns a hashcode for this object.
   *
   * @return the hashcode
   */
  public int hashCode() {

    return m_storedObjectArray.length;
  }

  /**
   * Returns a serialized object. Uses org.python.util.PythonObjectInputStream 
   * for Jython objects (read 
   * <a href="http://aspn.activestate.com/ASPN/Mail/Message/Jython-users/1001401">here</a>
   * for more details).
   *
   * @return the restored object
   * @exception Exception if the object couldn't be restored
   */ 
  public Object getObject() {

    try {
      ByteArrayInputStream istream = new ByteArrayInputStream(m_storedObjectArray);
      ObjectInputStream p;
      Object toReturn = null;
      if (m_isJython) {
	if (!m_isCompressed)
	  toReturn = Jython.deserialize(new BufferedInputStream(istream));
	else 
	  toReturn = Jython.deserialize(new BufferedInputStream(new GZIPInputStream(istream)));
      }
      else {
	if (!m_isCompressed)
	  p = new ObjectInputStream(new BufferedInputStream(istream));
	else 
	  p = new ObjectInputStream(new BufferedInputStream(new GZIPInputStream(istream)));
	toReturn = p.readObject();
      }
      istream.close();
      return toReturn;
    } catch (Exception e) {
      e.printStackTrace();
      return null;
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5018 $");
  }
}
