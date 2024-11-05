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
 * KOML.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 */

package weka.core.xml;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * This class is a helper class for XML serialization using <a href="http://koala.ilog.fr/XML/serialization/" target="_blank">KOML</a> .
 * KOML does not need to be present, since the class-calls are done generically via Reflection.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision 1.0$
 */
public class KOML
   implements RevisionHandler {
  
   /**
    * indicates whether <a href="http://koala.ilog.fr/XML/serialization/" target="_blank">KOML</a> 
    * (Koala Object Markup Language) is present
    */
   protected static boolean m_Present = false;

   /** the extension for KOML files (including '.') */
   public final static String FILE_EXTENSION = ".koml";
   
   /** check for KOML statically (needs only to be done once) */
   static {
      checkForKOML();
   }

   /**
    * checks whether the KOML is present in the class path
    */
   private static void checkForKOML() {
      try {
         Class.forName("fr.dyade.koala.xml.koml.KOMLSerializer");
         m_Present = true;
      }
      catch (Exception e) {
         m_Present = false;
      }
   }
  
   /**
    * returns whether KOML is present or not, i.e. whether the classes are in the
    * classpath or not
    *
    * @return whether KOML is available
    */
   public static boolean isPresent() {
      return m_Present;
   }
   
   /**
    * reads the XML-serialized object from the given file
    * @param filename the file to deserialize the object from
    * @return the deserialized object
    * @throws Exception if something goes wrong while reading from the file
    */
   public static Object read(String filename) throws Exception {
      return read(new FileInputStream(filename));
   }
   
   /**
    * reads the XML-serialized object from the given file
    * @param file the file to deserialize the object from
    * @return the deserialized object
    * @throws Exception if something goes wrong while reading from the file
    */
   public static Object read(File file) throws Exception {
      return read(new FileInputStream(file));
   }
   
   /**
    * reads the XML-serialized object from a stream
    * @param stream the stream to deserialize the object from
    * @return the deserialized object
    * @throws Exception if something goes wrong while reading from the stream
    */
   public static Object read(InputStream stream) throws Exception {
      Class<?>                            komlClass;
      Class[]                          komlClassArgs;
      Object[]                         komlArgs;
      java.lang.reflect.Constructor    constructor;
      Object                           koml;
      java.lang.reflect.Method         methodRead;
      java.lang.reflect.Method         methodClose;
      Class[]                          readArgsClasses;
      Class[]                          closeArgsClasses;
      Object[]                         readArgs;
      Object[]                         closeArgs;
      Object                           result;

      result = null;
      
      // get Deserializer
      komlClass        = Class.forName("fr.dyade.koala.xml.koml.KOMLDeserializer");
      komlClassArgs    = new Class[2];
      komlClassArgs[0] = java.io.InputStream.class;
      komlClassArgs[1] = Boolean.TYPE;
      komlArgs         = new Object[2];
      komlArgs[0]      = stream;
      komlArgs[1]      = new Boolean(false);
      constructor      = komlClass.getConstructor(komlClassArgs);
      koml             = constructor.newInstance(komlArgs);
      readArgsClasses  = new Class[0];
      methodRead       = komlClass.getMethod("readObject", readArgsClasses);
      readArgs         = new Object[0];
      closeArgsClasses = new Class[0];
      methodClose      = komlClass.getMethod("close", closeArgsClasses);
      closeArgs        = new Object[0];

      // execute it
      try {
         result = methodRead.invoke(koml, readArgs);
      }
      catch (Exception e) {
         result = null;
      } 
      finally {
         methodClose.invoke(koml, closeArgs);
      }
      
      return result;
   }
   
   /**
    * writes the XML-serialized object to the given file
    * @param filename the file to serialize the object to
    * @param o the object to write to the file
    * @return whether writing was successful or not
    * @throws Exception if something goes wrong while writing to the file
    */
   public static boolean write(String filename, Object o) throws Exception {
      return write(new FileOutputStream(filename), o);
   }
   
   /**
    * write the XML-serialized object to the given file
    * @param file the file to serialize the object to
    * @param o the object to write to the file
    * @return whether writing was successful or not
    * @throws Exception if something goes wrong while writing to the file
    */
   public static boolean write(File file, Object o) throws Exception {
      return write(new FileOutputStream(file), o);
   }
   
   /**
    * writes the XML-serialized object to a stream
    * @param stream the stream to serialize the object to
    * @param o the object to write to the stream
    * @return whether writing was successful or not
    * @throws Exception if something goes wrong while writing to the stream
    */
   public static boolean write(OutputStream stream, Object o) throws Exception {
      Class<?>                            komlClass;
      Class[]                          komlClassArgs;
      Object[]                         komlArgs;
      java.lang.reflect.Constructor    constructor;
      Object                           koml;
      java.lang.reflect.Method         methodAdd;
      java.lang.reflect.Method         methodClose;
      Class[]                          addArgsClasses;
      Class[]                          closeArgsClasses;
      Object[]                         addArgs;
      Object[]                         closeArgs;
      boolean                          result;
      
      result = false;

      // get Deserializer
      komlClass        = Class.forName("fr.dyade.koala.xml.koml.KOMLSerializer");
      komlClassArgs    = new Class[2];
      komlClassArgs[0] = java.io.OutputStream.class;
      komlClassArgs[1] = Boolean.TYPE;
      komlArgs         = new Object[2];
      komlArgs[0]      = stream;
      komlArgs[1]      = new Boolean(false);
      constructor      = komlClass.getConstructor(komlClassArgs);
      koml             = constructor.newInstance(komlArgs);
      addArgsClasses   = new Class[1];
      addArgsClasses[0] = Object.class;
      methodAdd        = komlClass.getMethod("addObject", addArgsClasses);
      addArgs          = new Object[1];
      addArgs[0]       = o;
      closeArgsClasses = new Class[0];
      methodClose      = komlClass.getMethod("close", closeArgsClasses);
      closeArgs        = new Object[0];

      // execute it
      try {
         methodAdd.invoke(koml, addArgs);
         result = true;
      }
      catch (Exception e) {
         result = false;
      } 
      finally {
         methodClose.invoke(koml, closeArgs);
      }
      
      return result;
   }
   
   /**
    * Returns the revision string.
    * 
    * @return		the revision
    */
   public String getRevision() {
     return RevisionUtils.extract("$Revision: 5953 $");
   }
}
