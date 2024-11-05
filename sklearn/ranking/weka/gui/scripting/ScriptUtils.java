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
 * ScriptUtils.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.scripting;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * A helper class for Script related stuff.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5928 $
 */
public class ScriptUtils {
  
  /**
   * Copies or moves files and directories (recursively).
   * If targetLocation does not exist, it will be created.
   * <p/>
   * Original code from <a href="http://www.java-tips.org/java-se-tips/java.io/how-to-copy-a-directory-from-one-location-to-another-loc.html" target="_blank">Java-Tips.org</a>.
   * 
   * @param sourceLocation	the source file/dir
   * @param targetLocation	the target file/dir
   * @param move		if true then the source files/dirs get deleted
   * 				as soon as copying finished
   * @throws IOException	if copying/moving fails
   */
  protected static void copyOrMove(File sourceLocation, File targetLocation, boolean move) throws IOException {
    String[] 		children;
    int 		i;
    InputStream 	in;
    OutputStream 	out;
    byte[] 		buf;
    int 		len;
    
    if (sourceLocation.isDirectory()) {
      if (!targetLocation.exists())
	targetLocation.mkdir();

      children = sourceLocation.list();
      for (i = 0; i < children.length; i++) {
	copyOrMove(
	    new File(sourceLocation, children[i]),
	    new File(targetLocation, children[i]),
	    move);
      }
      
      if (move)
	sourceLocation.delete();
    }
    else {
      in = new FileInputStream(sourceLocation);
      // do we need to append the filename?
      if (targetLocation.isDirectory())
	out = new FileOutputStream(targetLocation.getAbsolutePath() + File.separator + sourceLocation.getName());
      else
	out = new FileOutputStream(targetLocation);

      // Copy the content from instream to outstream
      buf = new byte[1024];
      while ((len = in.read(buf)) > 0)
	out.write(buf, 0, len);
      
      in.close();
      out.close();
      
      if (move)
	sourceLocation.delete();
    }
  }
  
  /**
   * Copies the file/directory (recursively).
   * 
   * @param sourceLocation	the source file/dir
   * @param targetLocation	the target file/dir
   * @throws IOException	if copying fails
   */
  public static void copy(File sourceLocation, File targetLocation) throws IOException {
    copyOrMove(sourceLocation, targetLocation, false);
  }
  
  /**
   * Moves the file/directory (recursively).
   * 
   * @param sourceLocation	the source file/dir
   * @param targetLocation	the target file/dir
   * @throws IOException	if moving fails
   */
  public static void move(File sourceLocation, File targetLocation) throws IOException {
    copyOrMove(sourceLocation, targetLocation, true);
  }
  
  /**
   * Saves the content to a file.
   * 
   * @param file		the file to save to
   * @param content		the content to save
   * @return			true if successfully saved
   */
  public static boolean save(File file, String content) {
    boolean		result;
    BufferedWriter	writer;
    
    writer = null;
    try {
      writer = new BufferedWriter(new FileWriter(file));
      writer.write(content);
      writer.flush();
      result = true;
    }
    catch (Exception e) {
      e.printStackTrace();
      result = false;
    }
    finally {
      if (writer != null) {
	try {
	  writer.close();
	}
	catch (Exception e) {
	  // ignored
	}
      }
    }
    
    return result;
  }
  
  /**
   * Tries to load the file and return its content.
   * 
   * @param file	the file to open
   * @return		the content, otherwise null
   */
  public static String load(File file) {
    StringBuffer	result;
    BufferedReader	reader;
    String		line;
    String		newLine;
    
    result  = new StringBuffer();
    newLine = System.getProperty("line.separator");
    reader  = null;
    try {
      // add new content
      reader = new BufferedReader(new FileReader(file));
      while ((line = reader.readLine()) != null) {
	result.append(line);
	result.append(newLine);
      }
    }
    catch (Exception e) {
      e.printStackTrace();
      result = null;
    }
    finally {
      if (reader != null) {
	try {
	  reader.close();
	}
	catch (Exception e) {
	  // ignored
	}
      }
    }
    
    return ((result != null) ? result.toString() : null);
  }
}
