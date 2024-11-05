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
 * MemoryUsagePanel.java
 * Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 */

package weka.gui;

import weka.core.Memory;
import weka.core.Utils;
import weka.gui.visualize.VisualizeUtils;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Properties;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

/**
 * A panel for displaying the memory usage.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.1 $
 */
public class MemoryUsagePanel
  extends JPanel {

  /** for serialization. */
  private static final long serialVersionUID = -4812319791687471721L;

  /**
   * Specialized thread for monitoring the memory usage.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 1.1 $
   */
  protected class MemoryMonitor
    extends Thread {
    
    /** the refresh interval in msecs. */
    protected int m_Interval;
    
    /** whether the thread is still running. */
    protected boolean m_Monitoring;
    
    /**
     * default constructor.
     */
    public MemoryMonitor() {
      super();

      setInterval(1000);  // TODO: via props file
    }
    
    /**
     * Returns the refresh interval in msecs.
     * 
     * @return		returns the refresh interval
     */
    public int getInterval() {
      return m_Interval;
    }
    
    /**
     * Sets the refresh interval in msecs.
     * 
     * @param value	the refresh interval
     */
    public void setInterval(int value) {
      m_Interval = value;
    }
    
    /**
     * Returns whether the thread is still running.
     * 
     * @return		true if the thread is still running
     */
    public boolean isMonitoring() {
      return m_Monitoring;
    }
    
    /**
     * stops the monitoring thread.
     */
    public void stopMonitoring() {
      m_Monitoring = false;
    }

    /**
     * The run method.
     */
    public void run() {
      m_Monitoring = true;
      
      while (m_Monitoring) {
	try {
	  Thread.sleep(m_Interval);
	  
	  // update GUI
	  if (m_Monitoring) {
	    Runnable doUpdate = new Runnable() {
	      public void run() {
		update();
	      }
	    };
	    SwingUtilities.invokeLater(doUpdate);
	  }
	} 
	catch(InterruptedException ex) { 
	  ex.printStackTrace();
	}
      }
    }
    
    /**
     * Updates the GUI.
     */
    protected void update() {
      double		perc;
      Dimension		size;

      // current usage
      perc = (double) m_Memory.getCurrent() / (double) m_Memory.getMax();
      perc = Math.round(perc * 1000) / 10;
      
      // tool tip
      setToolTipText("" + perc + "% used");
      
      // update history
      m_History.insertElementAt(perc, 0);
      size = getSize();
      while (m_History.size() > size.getWidth())
	m_History.remove(m_History.size() - 1);
      
      // display history
      repaint();
    }
  }
  
  /** The name of the properties file. */
  protected static String PROPERTY_FILE = "weka/gui/MemoryUsage.props";
    
  /** Contains the properties. */
  protected static Properties PROPERTIES;
  
  /** the memory usage over time. */
  protected Vector<Double> m_History;

  /** for monitoring the memory usage. */
  protected Memory m_Memory;

  /** the thread for monitoring the memory usage. */
  protected MemoryMonitor m_Monitor;
  
  /** the button for running the garbage collector. */
  protected JButton m_ButtonGC;

  /** the threshold percentages to change color. */
  protected Vector<Double> m_Percentages;
  
  /** the corresponding colors for the thresholds. */
  protected Hashtable<Double,Color> m_Colors;
  
  /** the default color. */
  protected Color m_DefaultColor;
  
  /** the background color. */
  protected Color m_BackgroundColor;

  /** the position for the dialog. */
  protected Point m_FrameLocation;
  
  /** 
   * Loads the configuration property file (USE_DYNAMIC is FALSE) or determines
   * the classes dynamically (USE_DYNAMIC is TRUE)
   * @see #USE_DYNAMIC
   * @see GenericPropertiesCreator
   */
  static {
    // Allow a properties file in the current directory to override
    try {
      PROPERTIES = Utils.readProperties(PROPERTY_FILE);
      Enumeration keys = PROPERTIES.propertyNames();
      if (!keys.hasMoreElements())
	throw new Exception("Failed to read a property file for the "
	    +"memory usage panel");
    }
    catch (Exception ex) {
      JOptionPane.showMessageDialog(
	  null,
	  "Could not read a configuration file for the memory usage\n"
	  +"panel. An example file is included with the Weka distribution.\n"
	  +"This file should be named \"" + PROPERTY_FILE + "\" and\n"
	  +"should be placed either in your user home (which is set\n"
	  + "to \"" + System.getProperties().getProperty("user.home") + "\")\n"
	  + "or the directory that java was started from\n",
	  "MemoryUsagePanel",
	  JOptionPane.ERROR_MESSAGE);
    }
  }
  
  /**
   * default constructor.
   */
  public MemoryUsagePanel() {
    super();

    // initializes members
    m_Memory      = new Memory();
    m_History     = new Vector<Double>();
    m_Percentages = new Vector<Double>();
    m_Colors      = new Hashtable<Double,Color>();

    // colors and percentages
    m_BackgroundColor = parseColor("BackgroundColor", Color.WHITE);
    m_DefaultColor    = parseColor("DefaultColor", Color.GREEN);
    String[] percs    = PROPERTIES.getProperty("Percentages", "70,80,90").split(",");
    for (int i = 0; i < percs.length; i++) {
      // do we have a color associated with percentage?
      if (PROPERTIES.getProperty(percs[i]) != null) {
	double perc;
	Color color;
	
	// try parsing the number
	try {
	  perc = Double.parseDouble(percs[i]);
	}
	catch (Exception e) {
	  System.err.println(
	      "MemoryUsagePanel: cannot parse percentage '" 
	      + percs[i] + "' - ignored!");
	  continue;
	}

	// try parsing the color
	color = parseColor(percs[i], null);
	if (color == null)
	  continue;
	
	// store color and percentage
	m_Percentages.add(perc);
	m_Colors.put(perc, color);
      }
      else {
	System.err.println(
	    "MemoryUsagePanel: cannot find color for percentage '" 
	    + percs[i] + "' - ignored!");
      }
    }
    Collections.sort(m_Percentages);
    
    // layout
    setLayout(new BorderLayout());

    JPanel panel = new JPanel(new BorderLayout());
    add(panel, BorderLayout.EAST);
    
    m_ButtonGC = new JButton("GC");
    m_ButtonGC.setToolTipText("Runs the garbage collector.");
    m_ButtonGC.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent evt) {
	System.gc();
      }
    });
    panel.add(m_ButtonGC, BorderLayout.NORTH);

    // dimensions
    int height;
    int width;
    try {
      height = Integer.parseInt(PROPERTIES.getProperty("Height", "" + (int) m_ButtonGC.getPreferredSize().getHeight()));
      width  = Integer.parseInt(PROPERTIES.getProperty("Width", "400"));
    }
    catch (Exception e) {
      System.err.println("MemoryUsagePanel: Problem parsing the dimensions - " + e);
      height = (int) m_ButtonGC.getPreferredSize().getHeight();
      width = 400;
    }
    setPreferredSize(new Dimension(width, height));

    // position
    int top;
    int left;
    try {
      top  = Integer.parseInt(PROPERTIES.getProperty("Top", "0"));
      left = Integer.parseInt(PROPERTIES.getProperty("Left", "0"));
    }
    catch (Exception e) {
      System.err.println("MemoryUsagePanel: Problem parsing the position - " + e);
      top  = 0;
      left = 0;
    }
    m_FrameLocation = new Point(left, top);
    
    // monitoring thread
    int interval;
    try {
      interval = Integer.parseInt(PROPERTIES.getProperty("Interval", "1000"));
    }
    catch (Exception e) {
      System.err.println("MemoryUsagePanel: Problem parsing the refresh interval - " + e);
      interval = 1000;
    }
    m_Monitor = new MemoryMonitor();
    m_Monitor.setInterval(interval);
    m_Monitor.setPriority(Thread.MAX_PRIORITY);
    m_Monitor.start();
  }

  /**
   * parses the color and returns the corresponding Color object.
   * 
   * @param prop	the color property to read and parse
   * @param defValue	the default color
   * @return		the parsed color or the default color of the 
   */
  protected Color parseColor(String prop, Color defValue) {
    Color	result;
    Color	color;
    String 	colorStr;
    
    result = defValue;
    
    try {
      colorStr = PROPERTIES.getProperty(prop);
      color    = VisualizeUtils.processColour(colorStr, result);
      if (color == null)
	throw new Exception(colorStr);
      result = color;
    }
    catch (Exception e) {
      System.err.println(
	  "MemoryUsagePanel: cannot parse color '" 
	  + e.getMessage() + "' - ignored!");
    }
    
    return result;
  }
  
  /**
   * Returns whether the thread is still running.
   * 
   * @return		true if the thread is still running
   */
  public boolean isMonitoring() {
    return m_Monitor.isMonitoring();
  }
  
  /**
   * stops the monitoring thread.
   */
  public void stopMonitoring() {
    m_Monitor.stopMonitoring();
  }
  
  /**
   * Returns the default position for the dialog.
   * 
   * @return		the default position
   */
  public Point getFrameLocation() {
    return m_FrameLocation;
  }
  
  /**
   * draws the background image.
   * 
   * @param g		the graphics context
   */
  public void paintComponent(Graphics g) {
    int		i;
    int		n;
    int		len;
    double	scale;
    double	perc;
    Color	color;
    
    super.paintComponent(g);
    
    g.setColor(m_BackgroundColor);
    g.fillRect(0, 0, getWidth(), getHeight());
    scale = (double) getHeight() / 100.0;
    for (i = 0; i < m_History.size(); i++) {
      perc = m_History.get(i);
      
      // determine color
      color = m_DefaultColor;
      for (n = m_Percentages.size() - 1; n >= 0; n--) {
	if (perc >= m_Percentages.get(n)) {
	  color = m_Colors.get(m_Percentages.get(n));
	  break;
	}
      }
      
      // paint line
      g.setColor(color);
      len = (int) Math.round(perc * scale);
      g.drawLine(i, getHeight() - 1, i, getHeight() - len);
    }
  }
}
