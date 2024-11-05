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
 *    FileSourcedConverter.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.core.converters;

import java.io.File;
import java.io.IOException;

/**
 * Interface to a loader/saver that loads/saves from a file source.
 *
 * @author Mark Hall
 * @version $Revision: 5953 $
 */
public interface FileSourcedConverter {

  /**
   * Get the file extension used for this type of file
   *
   * @return the file extension
   */
  public String getFileExtension();

  /**
   * Gets all the file extensions used for this type of file
   *
   * @return the file extensions
   */
  public String[] getFileExtensions();

  /**
   * Get a one line description of the type of file
   *
   * @return a description of the file type
   */
  public String getFileDescription();

  /**
   * Set the file to load from/ to save in
   *
   * @param file the file to load from
   * @exception IOException if an error occurs
   */
  public void setFile(File file) throws IOException;

  /**
   * Return the current source file/ destination file
   *
   * @return a <code>File</code> value
   */
  public File retrieveFile();

  /**
   * Set whether to use relative rather than absolute paths
   *
   * @param rp true if relative paths are to be used
   */
  public void setUseRelativePath(boolean rp);

  /**
   * Gets whether relative paths are to be used
   *
   * @return true if relative paths are to be used
   */
  public boolean getUseRelativePath();

}
