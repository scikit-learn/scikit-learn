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
  *    PostscriptWriter.java
  *    Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
  *
  */

package weka.gui.visualize;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;

import javax.swing.JComponent;

/** 
 * This class takes any Component and outputs it to a Postscript file.<p>
 * <b>Note:</b><br>
 * This writer does not work with Components that rely on clipping, like e.g.
 * scroll lists. Here the complete list is printed, instead of only in the
 * borders of the scroll list (may overlap other components!). This is due to
 * the way, clipping is handled in Postscript. There was no easy way around 
 * this issue. :-(
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.3 $
 * @see PostscriptGraphics
 */
public class PostscriptWriter
  extends JComponentWriter {
  
  /**
   * initializes the object 
   */
  public PostscriptWriter() {
    super(null);
  }
  
  /**
   * initializes the object with the given Component
   * 
   * @param c         the component to print in the output format
   */
  public PostscriptWriter(JComponent c) {
    super(c);
  }
  
  /**
   * initializes the object with the given Component and filename
   * 
   * @param c         the component to print in the output format
   * @param f         the file to store the output in
   */
  public PostscriptWriter(JComponent c, File f) {
    super(c, f);
  }
  
  /**
   * returns the name of the writer, to display in the FileChooser.
   * must be overridden in the derived class.
   */
  public String getDescription() {
    return "Postscript-File";
  }
  
  /**
   * returns the extension (incl. ".") of the output format, to use in the
   * FileChooser. 
   * must be overridden in the derived class.
   */
  public String getExtension() {
    return ".eps";
  }
  
  /**
   * generates the actual output
   * 
   * @throws Exception	if something goes wrong
   */
  public void generateOutput() throws Exception {
    BufferedOutputStream      ostrm;
    PostscriptGraphics        psg;

    ostrm = null;
    
    try { 
      ostrm = new BufferedOutputStream(new FileOutputStream(getFile()));
      psg = new PostscriptGraphics(getComponent().getHeight(), getComponent().getWidth(), ostrm);
      psg.setFont(getComponent().getFont());
      psg.scale(getXScale(), getYScale());
      getComponent().printAll(psg);
      psg.finished();
    } 
    catch (Exception e) {
      System.err.println(e); 
    } 
    finally { 
      if (ostrm != null) {
        try {
          ostrm.close();
        } catch (Exception e) {
          // Nothing to really do for error on close
        }
      }
    }
  }
  
  /**
   * for testing only
   */
  public static void main(String[] args) throws Exception {
    System.out.println("building TreeVisualizer...");
    weka.gui.treevisualizer.TreeBuild builder = new weka.gui.treevisualizer.TreeBuild();
    weka.gui.treevisualizer.NodePlace arrange = new weka.gui.treevisualizer.PlaceNode2();
    weka.gui.treevisualizer.Node top = builder.create(new java.io.StringReader("digraph atree { top [label=\"the top\"] a [label=\"the first node\"] b [label=\"the second nodes\"] c [label=\"comes off of first\"] top->a top->b b->c }"));
    weka.gui.treevisualizer.TreeVisualizer tv = new weka.gui.treevisualizer.TreeVisualizer(null, top, arrange);
    tv.setSize(800 ,600);
    
    String filename = System.getProperty("java.io.tmpdir") + "test.eps";
    System.out.println("outputting to '" + filename + "'...");
    toOutput(new PostscriptWriter(), tv, new File(filename));

    System.out.println("done!");
  }
}
