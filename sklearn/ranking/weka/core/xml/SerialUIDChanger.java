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
 * SerialUIDChanger.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 */

package weka.core.xml;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

/**
 * This class enables one to change the UID of a serialized object and therefore
 * not losing the data stored in the binary format.
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $ 
 */
public class SerialUIDChanger
   implements RevisionHandler {
  
   /**
    * checks whether KOML is present
    * 
    * @return returns <code>true</code> if KOML is present
    * @throws Exception if KOML is not present 
    */
   protected static boolean checkKOML() throws Exception {
      if (!KOML.isPresent())
         throw new Exception("KOML is not present!");
      else 
         return true; 
   }
   
   /**
    * checks whether the given filename ends with ".koml"
    * 
    * @param filename the filename to check
    * @return whether it is a KOML file or not
    * @see KOML#FILE_EXTENSION
    */
   public static boolean isKOML(String filename) {
      return filename.toLowerCase().endsWith(KOML.FILE_EXTENSION);
   }
   
   /**
    * loads a serialized object and returns it
    * 
    * @param binary the filename that points to the file containing the
    *        serialized object
    * @return the object from the file
    * @throws Exception if reading fails
    */
   protected static Object readBinary(String binary) throws Exception {
      FileInputStream         fi;
      ObjectInputStream       oi;
      Object                  o;
      
      fi = new FileInputStream(binary);
      oi = new ObjectInputStream(new BufferedInputStream(fi));
      o  = oi.readObject();
      oi.close();
      
      return o;
   }
   
   /**
    * serializes the given object into the given file
    * 
    * @param binary the file to store the object in
    * @param o the object to serialize
    * @throws Exception if saving fails 
    */
   protected static void writeBinary(String binary, Object o) throws Exception {
      FileOutputStream        fo;
      ObjectOutputStream      oo;

      fo = new FileOutputStream(binary);
      oo = new ObjectOutputStream(new BufferedOutputStream(fo));
      oo.writeObject(o);
      oo.close();
   }
   
   /**
    * converts a binary file into a KOML XML file
    * 
    * @param binary the binary file to convert
    * @param koml where to store the XML output
    * @throws Exception if conversion fails
    */
   public static void binaryToKOML(String binary, String koml) throws Exception {
      Object            o;
      
      // can we use KOML?
      checkKOML();

      // read binary
      o = readBinary(binary);
      if (o == null)
         throw new Exception("Failed to deserialize object from binary file '" + binary + "'!");
      
      // save as KOML
      KOML.write(koml, o);
   }
   
   /**
    * converts a KOML file into a binary one
    * 
    * @param koml the filename with the XML data
    * @param binary the name of the 
    */
   public static void komlToBinary(String koml, String binary) throws Exception {
      Object         o;
      
      // can we use KOML? 
      checkKOML();

      // read KOML
      o = KOML.read(koml);
      if (o == null)
         throw new Exception("Failed to deserialize object from XML file '" + koml + "'!");
      
      // write binary
      writeBinary(binary, o);
   }
   
   /**
    * changes the oldUID into newUID from the given file (binary/KOML) into the
    * other one (binary/KOML). it basically does a replace in the XML, i.e. it
    * looks for " uid='oldUID'" and replaces it with " uid='newUID'".
    * 
    * @param oldUID the old UID to change
    * @param newUID the new UID to use
    * @param fromFile the original file with the old UID
    * @param toFile the new file where to store the modified UID
    * @throws Exception if conversion fails
    */
   public static void changeUID(long oldUID, long newUID, String fromFile, String toFile) throws Exception {
      String            inputFile;
      String            tempFile;
      File              file;
      String            content;
      String            line;
      BufferedReader    reader;
      BufferedWriter    writer;
      
      // input
      if (!isKOML(fromFile)) {
         inputFile = fromFile + ".koml";
         binaryToKOML(fromFile, inputFile);
      }
      else {
         inputFile = fromFile;
      }
      
      // load KOML
      reader = new BufferedReader(new FileReader(inputFile));
      content = "";
      while ((line = reader.readLine()) != null) {
         if (!content.equals(""))
            content += "\n";
         content += line;
      }
      reader.close();
      
      // transform UID
      content = content.replaceAll(" uid='" + Long.toString(oldUID) + "'", " uid='" + Long.toString(newUID) + "'");
      
      // save to tempFile
      tempFile = inputFile + ".temp";
      writer = new BufferedWriter(new FileWriter(tempFile));
      writer.write(content);
      writer.flush();
      writer.close();
      
      // output
      if (!isKOML(toFile)) {
         komlToBinary(tempFile, toFile);
      }
      else {
         writer = new BufferedWriter(new FileWriter(toFile));
         writer.write(content);
         writer.flush();
         writer.close();
      }
      
      // remove tempFile
      file = new File(tempFile);
      file.delete();
   }
   
   /**
    * Returns the revision string.
    * 
    * @return		the revision
    */
   public String getRevision() {
     return RevisionUtils.extract("$Revision: 5953 $");
   }
   
   /**
    * exchanges an old UID for a new one. a file that doesn't end with ".koml"
    * is considered being binary.
    * takes four arguments: oldUID newUID oldFilename newFilename
    * 
    * @param args the command line parameters
    * @see KOML#FILE_EXTENSION
    */
   public static void main(String[] args) throws Exception {
      if (args.length != 4) {
         System.out.println();
         System.out.println("Usage: " + SerialUIDChanger.class.getName() + " <oldUID> <newUID> <oldFilename> <newFilename>");
         System.out.println("       <oldFilename> and <newFilename> have to be different");
         System.out.println();
      }
      else {
         if (args[2].equals(args[3]))
            throw new Exception("Filenames have to be different!");
         
         changeUID( Long.parseLong(args[0]),
                    Long.parseLong(args[1]),
                    args[2],
                    args[3] );
      }
   }
}
