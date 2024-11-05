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
 * Loader.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.io.Reader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;

/**
 * This class is for loading resources from a JAR archive.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.2 $ 
 */

public class Loader {
  /** the dir to use as prefix if filenames are w/o it, must have a slash
   * at the end (the path separator is a slash!) */
  private String          dir;
  
  /**
   * initializes the object
   */
  public Loader(String dir) {
    this.dir = dir;
  }
  
  /**
   * returns the dir prefix
   */
  public String getDir() {
    return dir;
  }
  
  /**
   * returns the processed filename, i.e. with the dir-prefix if it's
   * missing
   */
  public String processFilename(String filename) {
    if (!filename.startsWith(getDir()))
      filename = getDir() + filename;
    
    return filename;
  }
  
  /**
   * returns a URL for the given filename, can be NULL if it fails
   */
  public static URL getURL(String dir, String filename) {
    Loader         loader;
    
    loader = new Loader(dir);
    return loader.getURL(filename);
  }
  
  /**
   * returns a URL for the given filename, can be NULL if it fails
   */
  public URL getURL(String filename) {
    filename = processFilename(filename);
    return Loader.class.getClassLoader().getResource(filename);
  }
  
  /**
   * returns an InputStream for the given dir and filename, can be NULL if it 
   * fails
   */
  public static InputStream getInputStream(String dir, String filename) {
    Loader         loader;
    
    loader = new Loader(dir);
    return loader.getInputStream(filename);
  }
  
  /**
   * returns an InputStream for the given filename, can be NULL if it fails
   */
  public InputStream getInputStream(String filename) {
    filename = processFilename(filename);
    return Loader.class.getResourceAsStream(filename);
  }
  
  /**
   * returns a Reader for the given filename and dir, can be NULL if it fails
   */
  public static Reader getReader(String dir, String filename) {
    Loader            loader;
    
    loader = new Loader(dir);
    return loader.getReader(filename);
  }
  
  /**
   * returns a Reader for the given filename, can be NULL if it fails
   */
  public Reader getReader(String filename) {
    InputStream          in;
    
    in = getInputStream(filename);
    
    if (in == null)
      return null;
    else
      return new InputStreamReader(in);
  }
}
