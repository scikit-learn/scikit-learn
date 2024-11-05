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
 *    OutputZipper.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 * OutputZipper writes output to either gzipped files or to a
 * multi entry zip file. If the destination file is a directory
 * each output string will be written to an individually named
 * gzip file. If the destination file is a file, then each
 * output string is appended as a named entry to the zip file until
 * finished() is called to close the file.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public class OutputZipper
  implements RevisionHandler {
  
  File m_destination;
  DataOutputStream m_zipOut = null;
  ZipOutputStream m_zs = null;

  /**
   * Constructor.
   * 
   * @param destination a destination file or directory
   * @throws Exception if something goes wrong.
   */
  public OutputZipper(File destination) throws Exception { 

    m_destination = destination;

    // if a directory is specified then use gzip format, otherwise
    // use zip
    if (!m_destination.isDirectory()) {
      m_zs = new ZipOutputStream(new FileOutputStream(m_destination));
      m_zipOut = new DataOutputStream(m_zs);
    }    
  }

  /**
   * Saves a string to either an individual gzipped file or as
   * an entry in a zip file.
   * 
   * @param outString the output string to save
   * @param name the name of the file/entry to save it to
   * @throws Exception if something goes wrong
   */
  public void zipit(String outString, String name) throws Exception {
    File saveFile;
    ZipEntry ze;
    
    if (m_zipOut == null) {
      saveFile = new File(m_destination, name+".gz");
      DataOutputStream dout = 
	new DataOutputStream(new GZIPOutputStream(
			     new FileOutputStream(saveFile)));
      
      dout.writeBytes(outString);
      dout.close();
    } else {
      ze = new ZipEntry(name);
      m_zs.putNextEntry(ze);
      m_zipOut.writeBytes(outString);
      m_zs.closeEntry();
    }
  }

  /**
   * Closes the zip file.
   * 
   * @throws Exception if something goes wrong
   */
  public void finished() throws Exception {
    if (m_zipOut != null) {
      m_zipOut.close();
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.8 $");
  }

  /**
   * Main method for testing this class
   */
  public static void main(String [] args) {
    
    try {
      File testF = new File(new File(System.getProperty("user.dir")), 
			    "testOut.zip");
      OutputZipper oz = new OutputZipper(testF);
      
      /*      OutputZipper oz = new OutputZipper(
	      new File(System.getProperty("user.dir"))); */
      oz.zipit("Here is some test text to be zipped","testzip");
      oz.zipit("Here is a second entry to be zipped","testzip2");
      oz.finished();
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
