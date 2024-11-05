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
 *    ExtensionFileFilter.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui;

import java.io.File;
import java.io.Serializable;
import java.io.FilenameFilter;
import javax.swing.filechooser.FileFilter;

/**
 * Provides a file filter for FileChoosers that accepts or rejects files
 * based on their extension. Compatible with both java.io.FilenameFilter and
 * javax.swing.filechooser.FileFilter (why there are two I have no idea).
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public class ExtensionFileFilter
  extends FileFilter
  implements FilenameFilter, Serializable {

  /** The text description of the types of files accepted */
  protected String m_Description;

  /** The filename extensions of accepted files */
  protected String[] m_Extension;

  /**
   * Creates the ExtensionFileFilter
   *
   * @param extension the extension of accepted files.
   * @param description a text description of accepted files.
   */
  public ExtensionFileFilter(String extension, String description) {
    m_Extension = new String [1];
    m_Extension[0] = extension;
    m_Description = description;
  }

  /**
   * Creates an ExtensionFileFilter that accepts files that have any of
   * the extensions contained in the supplied array.
   *
   * @param extensions an array of acceptable file extensions (as Strings).
   * @param description a text description of accepted files.
   */
  public ExtensionFileFilter(String [] extensions, String description) {
    m_Extension = extensions;
    m_Description = description;
  }
  
  /**
   * Gets the description of accepted files.
   *
   * @return the description.
   */
  public String getDescription() {
    
    return m_Description;
  }
  
  /**
   * Returns a copy of the acceptable extensions.
   * 
   * @return the accepted extensions
   */
  public String[] getExtensions() {
    return (String[]) m_Extension.clone();
  }
  
  /**
   * Returns true if the supplied file should be accepted (i.e.: if it
   * has the required extension or is a directory).
   *
   * @param file the file of interest.
   * @return true if the file is accepted by the filter.
   */
  public boolean accept(File file) {
    
    String name = file.getName().toLowerCase();
    if (file.isDirectory()) {
      return true;
    }
    for (int i = 0; i < m_Extension.length; i++) {
      if (name.endsWith(m_Extension[i])) {
	return true;
      }
    }
    return false;
  }
  
  /**
   * Returns true if the file in the given directory with the given name
   * should be accepted.
   *
   * @param dir the directory where the file resides.
   * @param name the name of the file.
   * @return true if the file is accepted.
   */
  public boolean accept(File dir, String name) {
    return accept(new File(dir, name));
  }
}
