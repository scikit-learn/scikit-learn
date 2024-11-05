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
 * BMPWriter.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.visualize;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.imageio.ImageIO;
import javax.swing.JComponent;

/**
 * This class takes any JComponent and outputs it to a BMP-file.
 * Scaling is by default disabled, since we always take a screenshot.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5920 $
 */
public class BMPWriter
  extends JComponentWriter {

  /** the background color. */
  protected Color m_Background;
  
  /**
   * initializes the object.
   */
  public BMPWriter() {
    super();
  }

  /**
   * initializes the object with the given Component.
   * 
   * @param c		the component to print in the output format
   */
  public BMPWriter(JComponent c) {
    super(c);
  }

  /**
   * initializes the object with the given Component and filename.
   * 
   * @param c		the component to print in the output format
   * @param f		the file to store the output in
   */
  public BMPWriter(JComponent c, File f) {
    super(c, f);
  }
  
  /**
   * further initialization.
   */
  public void initialize() {
    super.initialize();
    
    setScalingEnabled(false);
  }

  /**
   * returns the name of the writer, to display in the FileChooser.
   * must be overridden in the derived class.
   * 
   * @return 		the name of the writer
   */
  public String getDescription() {
    return "BMP-Image";
  }
  
  /**
   * returns the extension (incl. ".") of the output format, to use in the
   * FileChooser. 
   * 
   * @return 		the file extension
   */
  public String getExtension() {
    return ".bmp";
  }
  
  /**
   * returns the current background color.
   * 
   * @return		the current background color
   */
  public Color getBackground() {
    return m_Background;
  }
  
  /**
   * sets the background color to use in creating the BMP.
   * 
   * @param c 		the color to use for background
   */
  public void setBackground(Color c) {
    m_Background = c;
  }
  
  /**
   * generates the actual output.
   * 
   * @throws Exception	if something goes wrong
   */
  public void generateOutput() throws Exception {
    BufferedImage	bi;
    Graphics		g;

    bi = new BufferedImage(getComponent().getWidth(), getComponent().getHeight(), BufferedImage.TYPE_INT_RGB);
    g  = bi.getGraphics();
    g.setPaintMode();
    g.setColor(getBackground());
    if (g instanceof Graphics2D)
      ((Graphics2D) g).scale(getXScale(), getYScale());
    g.fillRect(0, 0, getComponent().getWidth(), getComponent().getHeight());
    getComponent().printAll(g);
    ImageIO.write(bi, "bmp", getFile());
  }
  
  /**
   * for testing only.
   * 
   * @param args 	the commandline arguments
   * @throws Exception 	if something goes wrong
   */
  public static void main(String[] args) throws Exception {
    System.out.println("building TreeVisualizer...");
    weka.gui.treevisualizer.TreeBuild builder = new weka.gui.treevisualizer.TreeBuild();
    weka.gui.treevisualizer.NodePlace arrange = new weka.gui.treevisualizer.PlaceNode2();
    weka.gui.treevisualizer.Node top = builder.create(new java.io.StringReader("digraph atree { top [label=\"the top\"] a [label=\"the first node\"] b [label=\"the second nodes\"] c [label=\"comes off of first\"] top->a top->b b->c }"));
    weka.gui.treevisualizer.TreeVisualizer tv = new weka.gui.treevisualizer.TreeVisualizer(null, top, arrange);
    tv.setSize(800 ,600);
    
    String filename = System.getProperty("java.io.tmpdir") + File.separator + "test.bmp";
    System.out.println("outputting to '" + filename + "'...");
    toOutput(new BMPWriter(), tv, new File(filename));

    System.out.println("done!");
  }
}
