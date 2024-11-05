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
 *    Plot2D.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.visualize;

import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.Random;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * This class plots datasets in two dimensions. It can also plot
 * classifier errors and clusterer predictions.
 * 
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5987 $
 */
public class Plot2D
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -1673162410856660442L;  

  /* constants for shape types */
  public static final int MAX_SHAPES = 5;
  public static final int ERROR_SHAPE = 1000;
  public static final int MISSING_SHAPE = 2000;
  public static final int CONST_AUTOMATIC_SHAPE = -1;
  public static final int X_SHAPE = 0;
  public static final int PLUS_SHAPE = 1;
  public static final int DIAMOND_SHAPE = 2;
  public static final int TRIANGLEUP_SHAPE = 3;
  public static final int TRIANGLEDOWN_SHAPE = 4;
  public static final int DEFAULT_SHAPE_SIZE = 2;

  /** Default colour for the axis */
  protected Color m_axisColour = Color.green;

  /** Default colour for the plot background */
  protected Color m_backgroundColour = Color.black;

  /** The plots to display */
  protected FastVector m_plots = new FastVector();

  /** The master plot */
  protected PlotData2D m_masterPlot = null;

  /** The name of the master plot */
  protected String m_masterName = "master plot";

  /** The instances to be plotted */
  protected Instances m_plotInstances=null;

  /** An optional "compainion" of the panel. If specified, this
      class will get to do its thing with our graphics context
      before we do any drawing. Eg. the visualize panel may need
      to draw polygons etc. before we draw plot axis and data points */
  protected Plot2DCompanion m_plotCompanion=null;

  /** the class for displaying instance info. */
  protected Class m_InstanceInfoFrameClass = null;
  
  /** For popping up text info on data points */
  protected JFrame m_InstanceInfo = null;

  /** The list of the colors used */
  protected FastVector m_colorList;

  /** default colours for colouring discrete class */
  protected Color [] m_DefaultColors = {Color.blue,
					Color.red,
					Color.green,
					Color.cyan,
					Color.pink,
					new Color(255, 0, 255),
					Color.orange,
					new Color(255, 0, 0),
					new Color(0, 255, 0),
					Color.white};

  /** Indexes of the attributes to go on the x and y axis and the attribute
      to use for colouring and the current shape for drawing */
  protected int m_xIndex=0;
  protected int m_yIndex=0;
  protected int m_cIndex=0;
  protected int m_sIndex=0;

  /** Holds the min and max values of the x, y and colouring attributes 
   over all plots */
  protected double m_maxX;
  protected double m_minX;
  protected double m_maxY;
  protected double m_minY;
  protected double m_maxC;
  protected double m_minC;
    
  /** Axis padding */
  protected final int m_axisPad = 5;

  /** Tick size */
  protected final int m_tickSize = 5;

  /**the offsets of the axes once label metrics are calculated */
  protected int m_XaxisStart=0;
  protected int m_YaxisStart=0;
  protected int m_XaxisEnd=0;
  protected int m_YaxisEnd=0;

  /** if the user resizes the window, or the attributes selected for
      the attributes change, then the lookup table for points needs
      to be recalculated */
  protected boolean m_plotResize = true;
  
  /** if the user changes attribute assigned to an axis */
  protected boolean m_axisChanged = false;

  /** An array used to show if a point is hidden or not.
   * This is used for speeding up the drawing of the plot panel
   * although I am not sure how much performance this grants over
   * not having it.
   */
  protected int[][] m_drawnPoints;

  /** Font for labels */
  protected Font m_labelFont;
  protected FontMetrics m_labelMetrics=null; 

  /** the level of jitter */
  protected int m_JitterVal=0;

  /** random values for perterbing the data points */
  protected Random m_JRand = new Random(0);

  /** lookup table for plotted points */
  protected double [][] m_pointLookup=null;

  /** Constructor */
  public Plot2D() {
    super();
    setProperties();
    this.setBackground(m_backgroundColour);

    m_drawnPoints = new int[this.getWidth()][this.getHeight()];

    /** Set up some default colours */
    m_colorList = new FastVector(10);
    for (int noa = m_colorList.size(); noa < 10; noa++) {
      Color pc = m_DefaultColors[noa % 10];
      int ija =  noa / 10;
      ija *= 2; 
      for (int j=0;j<ija;j++) {
	pc = pc.darker();
      }
      
      m_colorList.addElement(pc);
    }
  }

  /**
   * Set the properties for Plot2D
   */
  private void setProperties() {
    if (VisualizeUtils.VISUALIZE_PROPERTIES != null) {
      String thisClass = this.getClass().getName();
      String axisKey = thisClass+".axisColour";
      String backgroundKey = thisClass+".backgroundColour";

      String axisColour = VisualizeUtils.VISUALIZE_PROPERTIES.
	getProperty(axisKey);
      if (axisColour == null) {
	/*
	System.err.println("Warning: no configuration property found in "
			   +VisualizeUtils.PROPERTY_FILE
			   +" for "+axisKey);*/
      } else {
	//System.err.println("Setting axis colour to: "+axisColour);
	m_axisColour = VisualizeUtils.processColour(axisColour, m_axisColour);
      }

      String backgroundColour = 
	VisualizeUtils.VISUALIZE_PROPERTIES.getProperty(backgroundKey);
      if (backgroundColour == null) {
	/*
	System.err.println("Warning: no configuration property found in "
			   +VisualizeUtils.PROPERTY_FILE
			   +" for "+backgroundKey);*/
      } else {
	//System.err.println("Setting background colour to: "+backgroundColour);
	m_backgroundColour = VisualizeUtils.processColour(backgroundColour, 
							  m_backgroundColour);
      }
      
      try {
	m_InstanceInfoFrameClass = Class.forName(VisualizeUtils.VISUALIZE_PROPERTIES.getProperty(thisClass + ".instanceInfoFrame", "weka.gui.visualize.InstanceInfoFrame"));
      }
      catch (Exception e) {
	e.printStackTrace();
	m_InstanceInfoFrameClass = InstanceInfoFrame.class;
      }
    }
  }

  /** 
   * This will check the values of the screen points passed and make sure 
   * that they land on the screen
   * @param x1 The x coord.
   * @param y1 The y coord.
   */
  private boolean checkPoints(double x1, double y1) {
    if (x1 < 0 || x1 > this.getSize().width || y1 < 0 
	|| y1 > this.getSize().height) {
      return false;
    }
    return true;
  }

  /**
   * Set a companion class. This is a class that might want
   * to render something on the plot before we do our thing. Eg,
   * Malcolm's shape drawing stuff needs to happen before we plot
   * axis and points
   * @param p a companion class
   */
  public void setPlotCompanion(Plot2DCompanion p) {
    m_plotCompanion = p;
  }

  /**
   * Set level of jitter and repaint the plot using the new jitter value
   * @param j the level of jitter
   */
  public void setJitter(int j) {
    if (m_plotInstances.numAttributes() > 0 
	&& m_plotInstances.numInstances() > 0) {
      if (j >= 0) {
	m_JitterVal = j;
	m_JRand = new Random(m_JitterVal);
	//      if (m_pointLookup != null) {
	m_drawnPoints = new int[m_XaxisEnd - m_XaxisStart + 1]
	  [m_YaxisEnd - m_YaxisStart + 1];
	updatePturb();
	//      }
	this.repaint();
      }
    }
  }

  /**
   * Set a list of colours to use when colouring points according
   * to class values or cluster numbers
   * @param cols the list of colours to use
   */
  public void setColours (FastVector cols) {
    m_colorList = cols;
  }

  /**
   * Set the index of the attribute to go on the x axis
   * @param x the index of the attribute to use on the x axis
   */
  public void setXindex(int x) {
    m_xIndex = x;
    for (int i=0;i<m_plots.size();i++) {
      ((PlotData2D)m_plots.elementAt(i)).setXindex(m_xIndex);
    }
    determineBounds();
    if (m_JitterVal != 0) {
      updatePturb();
    }
    m_axisChanged = true;
    this.repaint();
  }
    
  /**
   * Set the index of the attribute to go on the y axis
   * @param y the index of the attribute to use on the y axis
   */
  public void setYindex(int y) {
    m_yIndex = y;
    for (int i=0;i<m_plots.size();i++) {
      ((PlotData2D)m_plots.elementAt(i)).setYindex(m_yIndex);
    }
    determineBounds();
    if (m_JitterVal != 0) {
      updatePturb();
    }
    m_axisChanged = true;
    this.repaint();
  }

  /**
   * Set the index of the attribute to use for colouring
   * @param c the index of the attribute to use for colouring
   */
  public void setCindex(int c) {
    m_cIndex = c;
    for (int i=0;i<m_plots.size();i++) {
      ((PlotData2D)m_plots.elementAt(i)).setCindex(m_cIndex);
    }
    determineBounds();
    m_axisChanged = true;
    this.repaint();
  }

  /**
   * Return the list of plots
   * @return the list of plots
   */
  public FastVector getPlots() {
    return m_plots;
  }

  /**
   * Get the master plot
   * @return the master plot
   */
  public PlotData2D getMasterPlot() {
    return m_masterPlot;
  }

  /** 
   * Return the current max value of the attribute plotted on the x axis
   * @return the max x value
   */
  public double getMaxX() {
    return m_maxX;
  }

  /** 
   * Return the current max value of the attribute plotted on the y axis
   * @return the max y value
   */
  public double getMaxY() {
    return m_maxY;
  }

  /** 
   * Return the current min value of the attribute plotted on the x axis
   * @return the min x value
   */
  public double getMinX() {
    return m_minX;
  }
  
  /** 
   * Return the current min value of the attribute plotted on the y axis
   * @return the min y value
   */
  public double getMinY() {
    return m_minY;
  }

  /** 
   * Return the current max value of the colouring attribute
   * @return the max colour value
   */
  public double getMaxC() {
    return m_maxC;
  }
  
  /** 
   * Return the current min value of the colouring attribute
   * @return the min colour value
   */
  public double getMinC() {
    return m_minC;
  }
    
  /**
   * Sets the master plot from a set of instances
   * @param inst the instances
   * @exception Exception if instances could not be set
   */
  public void setInstances(Instances inst) throws Exception {
    //System.err.println("Setting Instances");
    PlotData2D tempPlot = new PlotData2D(inst);
    tempPlot.setPlotName("master plot");
    setMasterPlot(tempPlot);
  }

  /**
   * Set the master plot.
   * @param master the plot to make the master plot
   * @exception Exception if the plot could not be set.
   */
  public void setMasterPlot(PlotData2D master) throws Exception {
    if (master.m_plotInstances == null) {
      throw new Exception("No instances in plot data!");
    }
    removeAllPlots();
    m_masterPlot = master;
    m_plots.addElement(m_masterPlot);
    m_plotInstances = m_masterPlot.m_plotInstances;
    
    m_xIndex=0;
    m_yIndex=0;
    m_cIndex=0;
    
    determineBounds();
  }

  /**
   * Clears all plots
   */
  public void removeAllPlots() {
    m_masterPlot = null;
    m_plotInstances = null;
    m_plots = new FastVector();
    m_xIndex = 0; m_yIndex = 0; m_cIndex = 0;
  }

  /**
   * Add a plot to the list of plots to display
   * @param newPlot the new plot to add
   * @exception Exception if the plot could not be added
   */
  public void addPlot(PlotData2D newPlot) throws Exception {
    if (newPlot.m_plotInstances == null) {
      throw new Exception("No instances in plot data!");
    }

    if (m_masterPlot != null) {
      if (m_masterPlot.m_plotInstances.
	  equalHeaders(newPlot.m_plotInstances) == false) {
	throw new Exception("Plot2D :Plot data's instances are incompatable "
			    +" with master plot");
      }
    } else {
      m_masterPlot = newPlot;
      m_plotInstances = m_masterPlot.m_plotInstances;
    }
    m_plots.addElement(newPlot);
    setXindex(m_xIndex);
    setYindex(m_yIndex);
    setCindex(m_cIndex);
  }

  /**
   * Set up fonts and font metrics
   * @param gx the graphics context
   */
  private void setFonts(Graphics gx) {
    if (m_labelMetrics == null) {
      m_labelFont = new Font("Monospaced", Font.PLAIN, 12);
      m_labelMetrics = gx.getFontMetrics(m_labelFont);
    }
    gx.setFont(m_labelFont);
  }

  /**
   * Pops up a window displaying attribute information on any instances
   * at a point+-plotting_point_size (in panel coordinates)
   *
   * @param x the x value of the clicked point
   * @param y the y value of the clicked point
   * @param newFrame true if instance info is to be displayed in a
   * new frame.
   */
  public void searchPoints(int x, int y, final boolean newFrame) {
    if (m_masterPlot.m_plotInstances != null) {
      int longest=0;
      for (int j=0;j<m_masterPlot.m_plotInstances.numAttributes();j++) {
	if (m_masterPlot.m_plotInstances.attribute(j).name().length() > 
	    longest) {
	  longest = m_masterPlot.m_plotInstances.attribute(j).name().length();
	}
      }

      StringBuffer insts = new StringBuffer(); 
      Vector<Instances> data = new Vector<Instances>();
      for (int jj=0;jj<m_plots.size();jj++) {
	PlotData2D temp_plot = (PlotData2D)(m_plots.elementAt(jj));
	data.add(new Instances(temp_plot.m_plotInstances, 0));
	
	for (int i=0;i<temp_plot.m_plotInstances.numInstances();i++) {
	  if (temp_plot.m_pointLookup[i][0] != Double.NEGATIVE_INFINITY) {
	    double px = temp_plot.m_pointLookup[i][0] + 
	      temp_plot.m_pointLookup[i][2];
	    double py = temp_plot.m_pointLookup[i][1] + 
	      temp_plot.m_pointLookup[i][3];
	    //	    double size = temp_plot.m_pointLookup[i][2];
	    double size = temp_plot.m_shapeSize[i];
	    if ((x >= px-size) && (x <= px+size) &&
		(y >= py-size) && (y <= py+size)) {
	      {
		data.get(jj).add((Instance)temp_plot.m_plotInstances.instance(i).copy());
		insts.append("\nPlot : "+temp_plot.m_plotName
			     + "\nInstance: " + (i + 1 ) + "\n");
		for (int j=0;j<temp_plot.m_plotInstances.numAttributes();j++) {
		  for (int k = 0;k < 
			 (longest-temp_plot.m_plotInstances.
			  attribute(j).name().length()); k++) {
		    insts.append(" ");
		  }
		  insts.append(temp_plot.m_plotInstances.attribute(j).name());  
		  insts.append(" : ");
		  
		  if (temp_plot.m_plotInstances.instance(i).isMissing(j)) {
		    insts.append("Missing");
		  } else if (temp_plot.m_plotInstances.attribute(j).
			     isNominal() || temp_plot.m_plotInstances.attribute(j).
			     isRanking()) {
		    insts.append(temp_plot.m_plotInstances.
				 attribute(j).
				 value((int)temp_plot.m_plotInstances.
				       instance(i).value(j)));
		  } else {
		    insts.append(temp_plot.m_plotInstances.
				 instance(i).value(j));
		  }
		  insts.append("\n");
		}
	      }
	    }
	  }
	}
      }
      // remove datasets that contain no instances
      int i = 0;
      while (data.size() > i) {
	if (data.get(i).numInstances() == 0)
	  data.remove(i);
	else
	  i++;
      }

      if (insts.length() > 0) {
	// Pop up a new frame
	if (newFrame || m_InstanceInfo == null) {
	  try {
	    final JFrame jf = (JFrame) m_InstanceInfoFrameClass.newInstance();
	    ((InstanceInfo) jf).setInfoText(insts.toString());
	    ((InstanceInfo) jf).setInfoData(data);
	    final JFrame testf = m_InstanceInfo;
	    jf.addWindowListener(new WindowAdapter() {
	      public void windowClosing(WindowEvent e) {
		if (!newFrame || testf == null) {
		  m_InstanceInfo = null;
		}
		jf.dispose();
	      }
	    });
	    jf.setVisible(true);
	    if (m_InstanceInfo == null)
	      m_InstanceInfo = jf;
	  }
	  catch (Exception e) {
	    e.printStackTrace();
	  }
	}
	else {
	  // Overwrite info in existing frame	  
	  ((InstanceInfo) m_InstanceInfo).setInfoText(insts.toString());
	  ((InstanceInfo) m_InstanceInfo).setInfoData(data);
	}
      }
    }
  }
  

  /**
   * Determine the min and max values for axis and colouring attributes
   */
  public void determineBounds() {
    double value,min,max;
    
    // find maximums minimums over all plots
    m_minX = ((PlotData2D)m_plots.elementAt(0)).m_minX;
    m_maxX = ((PlotData2D)m_plots.elementAt(0)).m_maxX;
    m_minY = ((PlotData2D)m_plots.elementAt(0)).m_minY;
    m_maxY = ((PlotData2D)m_plots.elementAt(0)).m_maxY;
    m_minC = ((PlotData2D)m_plots.elementAt(0)).m_minC;
    m_maxC = ((PlotData2D)m_plots.elementAt(0)).m_maxC;
    for (int i=1;i<m_plots.size();i++) {
      value = ((PlotData2D)m_plots.elementAt(i)).m_minX;
      if (value < m_minX) {
	m_minX = value;
      }
      value = ((PlotData2D)m_plots.elementAt(i)).m_maxX;
      if (value > m_maxX) {
	m_maxX = value;
      }
      value = ((PlotData2D)m_plots.elementAt(i)).m_minY;
      if (value < m_minY) {
	m_minY= value;
      }
      value = ((PlotData2D)m_plots.elementAt(i)).m_maxY;
      if (value > m_maxY) {
	m_maxY = value;
      }
      value = ((PlotData2D)m_plots.elementAt(i)).m_minC;
      if (value < m_minC) {
	m_minC = value;
      }
      value = ((PlotData2D)m_plots.elementAt(i)).m_maxC;
      if (value > m_maxC) {
	m_maxC = value;
      }
    }

    fillLookup();
    this.repaint();
  }

    
  //to convert screen coords to attrib values
  // note that I use a double to avoid accuracy 
  //headaches with ints
  /**
   * convert a Panel x coordinate to a raw x value.
   * @param scx The Panel x coordinate
   * @return A raw x value.
   */
  public double convertToAttribX(double scx) {
    double temp = m_XaxisEnd - m_XaxisStart;
    double temp2 = ((scx - m_XaxisStart) * (m_maxX - m_minX)) / temp;
      
    temp2 = temp2 + m_minX;
      
    return temp2;
  }
    
  /**
   * convert a Panel y coordinate to a raw y value.
   * @param scy The Panel y coordinate
   * @return A raw y value.
   */
  public double convertToAttribY(double scy) {
    double temp = m_YaxisEnd - m_YaxisStart;
    double temp2 = ((scy - m_YaxisEnd) * (m_maxY - m_minY)) / temp;
      
    temp2 = -(temp2 - m_minY);
      
    return temp2;
  }
  //////
    
  /**
   * returns a value by which an x value can be peturbed. Makes sure
   * that the x value+pturb stays within the plot bounds
   * @param xvalP the x coordinate to be peturbed
   * @param xj a random number to use in calculating a peturb value
   * @return a peturb value
   */
  int pturbX(double xvalP, double xj) {
    int xpturb = 0;
    if (m_JitterVal > 0) {
      xpturb = (int)((double)m_JitterVal * (xj / 2.0));
      if (((xvalP + xpturb) < m_XaxisStart) || 
	  ((xvalP + xpturb) > m_XaxisEnd)) {
	xpturb *= -1;
      }
    }
    return xpturb;
  }

  /**
   * Convert an raw x value to Panel x coordinate.
   * @param xval the raw x value
   * @return an x value for plotting in the panel.
   */
  public double convertToPanelX(double xval) {
    double temp = (xval - m_minX)/(m_maxX - m_minX);
    double temp2 = temp * (m_XaxisEnd - m_XaxisStart);
      
    temp2 = temp2 + m_XaxisStart;
	
    return temp2;
  }

  /**
   * returns a value by which a y value can be peturbed. Makes sure
   * that the y value+pturb stays within the plot bounds
   * @param yvalP the y coordinate to be peturbed
   * @param yj a random number to use in calculating a peturb value
   * @return a peturb value
   */
  int pturbY(double yvalP, double yj) {
    int ypturb = 0;
    if (m_JitterVal > 0) {
      ypturb = (int)((double)m_JitterVal * (yj / 2.0));
      if (((yvalP + ypturb) < m_YaxisStart) || 
	  ((yvalP + ypturb) > m_YaxisEnd)) {
	ypturb *= -1;
      }
    }
    return ypturb;
  }

  /**
   * Convert an raw y value to Panel y coordinate.
   * @param yval the raw y value
   * @return an y value for plotting in the panel.
   */
  public double convertToPanelY(double yval) {
    double temp = (yval - m_minY)/(m_maxY - m_minY);
    double temp2 = temp * (m_YaxisEnd - m_YaxisStart);
      
    temp2 = m_YaxisEnd - temp2;

    return temp2;
  }

  /**
   * Draws an X.
   * @param gx the graphics context
   * @param x the x coord
   * @param y the y coord
   * @param size the size of the shape
   */
  private static void drawX(Graphics gx, double x, double y, int size) {
     gx.drawLine((int)(x-size),(int)(y-size),
		  (int)(x+size),(int)(y+size));
     
      gx.drawLine((int)(x+size),(int)(y-size),
		  (int)(x-size),(int)(y+size));     
  }

  /**
   * Draws a plus.
   * @param gx the graphics context
   * @param x the x coord
   * @param y the y coord
   * @param size the size of the shape
   */
  private static void drawPlus(Graphics gx, double x, double y, int size) {
     gx.drawLine((int)(x-size),(int)(y),
		  (int)(x+size),(int)(y));
     
      gx.drawLine((int)(x),(int)(y-size),
		  (int)(x),(int)(y+size));     
  }

  /**
   * Draws a diamond.
   * @param gx the graphics context
   * @param x the x coord
   * @param y the y coord
   * @param size the size of the shape
   */
  private static void drawDiamond(Graphics gx, double x, double y, int size) {
    gx.drawLine((int)(x-size),(int)(y),
		(int)(x),(int)(y-size));
    
    gx.drawLine((int)(x),(int)(y-size),
		  (int)(x+size),(int)(y));

    gx.drawLine((int)(x+size),(int)(y),
		  (int)(x),(int)(y+size));

     gx.drawLine((int)(x),(int)(y+size),
		  (int)(x-size),(int)(y));
  }

  /**
   * Draws an triangle (point at top).
   * @param gx the graphics context
   * @param x the x coord
   * @param y the y coord
   * @param size the size of the shape
   */
  private static void drawTriangleUp(Graphics gx, double x, 
				     double y, int size) {
    gx.drawLine((int)(x),(int)(y-size),
		(int)(x-size),(int)(y+size));

    gx.drawLine((int)(x-size),(int)(y+size),
		(int)(x+size),(int)(y+size));

    gx.drawLine((int)(x+size),(int)(y+size),
		(int)(x),(int)(y-size));

  }

  /**
   * Draws an triangle (point at bottom).
   * @param gx the graphics context
   * @param x the x coord
   * @param y the y coord
   * @param size the size of the shape
   */
  private static void drawTriangleDown(Graphics gx, double x, 
				       double y, int size) {
    gx.drawLine((int)(x),(int)(y+size),
		(int)(x-size),(int)(y-size));

    gx.drawLine((int)(x-size),(int)(y-size),
		(int)(x+size),(int)(y-size));

    gx.drawLine((int)(x+size),(int)(y-size),
		(int)(x),(int)(y+size));

  }

  /**
   * Draws a data point at a given set of panel coordinates at a given
   * size and connects a line to the previous point.
   * @param x the x coord
   * @param y the y coord
   * @param xprev the x coord of the previous point
   * @param yprev the y coord of the previous point
   * @param size the size of the point
   * @param shape the shape of the data point (square is reserved for nominal
   * error data points). Shapes: 0=x, 1=plus, 2=diamond, 3=triangle(up),
   * 4 = triangle (down).
   * @param gx the graphics context
   */
  protected static void drawDataPoint(double x, 
			       double y,
			       double xprev,
			       double yprev,
			       int size,
			       int shape,
			       Graphics gx) {

    drawDataPoint(x,y,size,shape,gx);

    // connect a line to the previous point
    gx.drawLine((int)x, (int)y, (int)xprev, (int)yprev);
  }

  /**
   * Draws a data point at a given set of panel coordinates at a given
   * size.
   * @param x the x coord
   * @param y the y coord
   * @param size the size of the point
   * @param shape the shape of the data point (square is reserved for nominal
   * error data points). Shapes: 0=x, 1=plus, 2=diamond, 3=triangle(up),
   * 4 = triangle (down).
   * @param gx the graphics context
   */
  protected static void drawDataPoint(double x, 
				      double y,
				      int size,
				      int shape,
				      Graphics gx) {

    Font lf = new Font("Monospaced", Font.PLAIN, 12);
    FontMetrics fm = gx.getFontMetrics(lf);

    if (size == 0) {
      size = 1;
    }

    if (shape != ERROR_SHAPE && shape != MISSING_SHAPE) {
      shape = shape % 5;
    }

    switch (shape) {
    case X_SHAPE:
      drawX(gx, x, y, size);
      break;      
    case PLUS_SHAPE:
      drawPlus(gx, x, y, size);
      break;
    case DIAMOND_SHAPE:
      drawDiamond(gx, x, y, size);
      break;
    case TRIANGLEUP_SHAPE:
      drawTriangleUp(gx, x, y, size);
      break;
    case TRIANGLEDOWN_SHAPE:
      drawTriangleDown(gx, x, y, size);
      break;
    case ERROR_SHAPE: // draws the nominal error shape 
      gx.drawRect((int)(x-size),(int)(y-size),(size*2),(size*2));
      break;
    case MISSING_SHAPE:
      int hf = fm.getAscent();
      int width = fm.stringWidth("M");
      gx.drawString("M",(int)(x-(width / 2)), (int)(y+(hf / 2)));
      break;
    }
  }

  /**
   * Updates the perturbed values for the plots when the jitter value is
   * changed
   */
  private void updatePturb() {
    double xj=0;
    double yj=0;
    for (int j=0;j<m_plots.size();j++) {
      PlotData2D temp_plot = (PlotData2D)(m_plots.elementAt(j));
      for (int i=0;i<temp_plot.m_plotInstances.numInstances();i++) {
	if (temp_plot.m_plotInstances.instance(i).isMissing(m_xIndex) ||
	    temp_plot.m_plotInstances.instance(i).isMissing(m_yIndex)) {
	} else {
	  if (m_JitterVal > 0) {
	    xj = m_JRand.nextGaussian();
	    yj = m_JRand.nextGaussian();
	  }
	  temp_plot.m_pointLookup[i][2] = 
	    pturbX(temp_plot.m_pointLookup[i][0],xj);
	  temp_plot.m_pointLookup[i][3] = 
	    pturbY(temp_plot.m_pointLookup[i][1],yj);
	}
      }
    }
  }

  /**
   * Fills the lookup caches for the plots. Also calculates errors for
   * numeric predictions (if any) in plots
   */
  private void fillLookup() {

    for (int j=0;j<m_plots.size();j++) {
      PlotData2D temp_plot = (PlotData2D)(m_plots.elementAt(j));

      if (temp_plot.m_plotInstances.numInstances() > 0 &&
	  temp_plot.m_plotInstances.numAttributes() > 0) {
	for (int i=0;i<temp_plot.m_plotInstances.numInstances();i++) {
	  if (temp_plot.m_plotInstances.instance(i).isMissing(m_xIndex) ||
	      temp_plot.m_plotInstances.instance(i).isMissing(m_yIndex)) {
	    temp_plot.m_pointLookup[i][0] = Double.NEGATIVE_INFINITY;
	    temp_plot.m_pointLookup[i][1] = Double.NEGATIVE_INFINITY;
	  } else {
	    double x = convertToPanelX(temp_plot.m_plotInstances.
				       instance(i).value(m_xIndex));
	    double y = convertToPanelY(temp_plot.m_plotInstances.
				       instance(i).value(m_yIndex));
	    temp_plot.m_pointLookup[i][0] = x;
	    temp_plot.m_pointLookup[i][1] = y;
	  }
	}
      }
    }
  }
    
  /**
   * Draws the data points and predictions (if provided).
   * @param gx the graphics context
   */
  private void paintData(Graphics gx) {

    for (int j=0;j<m_plots.size();j++) {
      PlotData2D temp_plot = (PlotData2D)(m_plots.elementAt(j));

      for (int i=0;i<temp_plot.m_plotInstances.numInstances();i++) {
	if (temp_plot.m_plotInstances.instance(i).isMissing(m_xIndex) ||
	    temp_plot.m_plotInstances.instance(i).isMissing(m_yIndex)) {
	} else {
	  double x = (temp_plot.m_pointLookup[i][0] + 
		      temp_plot.m_pointLookup[i][2]);
	  double y = (temp_plot.m_pointLookup[i][1] + 
		      temp_plot.m_pointLookup[i][3]);

	  double prevx = 0;
	  double prevy = 0;
	  if (i > 0) {
	    prevx = (temp_plot.m_pointLookup[i - 1][0] + 
			temp_plot.m_pointLookup[i - 1][2]);
	    prevy = (temp_plot.m_pointLookup[i - 1][1] + 
			temp_plot.m_pointLookup[i - 1][3]);
	  }

	  int x_range = (int)x - m_XaxisStart;
	  int y_range = (int)y - m_YaxisStart;

	  if (x_range >= 0 && y_range >= 0) {
	    if (m_drawnPoints[x_range][y_range] == i 
		|| m_drawnPoints[x_range][y_range] == 0
		|| temp_plot.m_shapeSize[i] == temp_plot.m_alwaysDisplayPointsOfThisSize 
		|| temp_plot.m_displayAllPoints == true) {
	      m_drawnPoints[x_range][y_range] = i;
	      if (temp_plot.m_plotInstances.attribute(m_cIndex).isNominal() || temp_plot.m_plotInstances.attribute(m_cIndex).isRanking()) {
		if (temp_plot.m_plotInstances.attribute(m_cIndex).numValues() >
		    m_colorList.size() && 
		    !temp_plot.m_useCustomColour) {
		  extendColourMap(temp_plot.m_plotInstances.
				  attribute(m_cIndex).numValues());
		}

		Color ci;
		if (temp_plot.m_plotInstances.instance(i).
		    isMissing(m_cIndex)) {
		  ci = Color.gray;
		} else {
		  int ind = (int)temp_plot.m_plotInstances.instance(i).
		    value(m_cIndex);
		  ci = (Color)m_colorList.elementAt(ind);
		}

		if (!temp_plot.m_useCustomColour) {
		  gx.setColor(ci);	    
		} else {
		  gx.setColor(temp_plot.m_customColour);
		}

		if (temp_plot.m_plotInstances.instance(i).
		    isMissing(m_cIndex)) {
		  if (temp_plot.m_connectPoints[i] == true) {
		    drawDataPoint(x,y,prevx,prevy,temp_plot.m_shapeSize[i],
				  MISSING_SHAPE,gx);
		  } else {
		    drawDataPoint(x,y,temp_plot.m_shapeSize[i],
				  MISSING_SHAPE,gx);
		  }
		} else {
		  if (temp_plot.m_shapeType[i] == CONST_AUTOMATIC_SHAPE) {
		    if (temp_plot.m_connectPoints[i] == true) {
		      drawDataPoint(x,y,prevx,prevy,
				    temp_plot.m_shapeSize[i],j,gx);
		    } else {
		      drawDataPoint(x,y,temp_plot.m_shapeSize[i],j,gx);
		    }
		  } else {
		    if (temp_plot.m_connectPoints[i] == true) {
		       drawDataPoint(x,y,prevx,prevy,temp_plot.m_shapeSize[i],
				     temp_plot.m_shapeType[i],gx);
		    } else {
		      drawDataPoint(x,y,temp_plot.m_shapeSize[i],
				    temp_plot.m_shapeType[i],gx);
		    }
		  }
		}
	      } else {
		double r;
		Color ci = null;
		if (!temp_plot.m_plotInstances.instance(i).
		    isMissing(m_cIndex)) {
		  r = (temp_plot.m_plotInstances.instance(i).
		       value(m_cIndex) - m_minC) / (m_maxC - m_minC);
		  r = (r * 240) + 15;
		  ci = new Color((int)r,150,(int)(255-r));
		} else {
		  ci = Color.gray;
		}
		if (!temp_plot.m_useCustomColour) {
		  gx.setColor(ci);
		} else {
		  gx.setColor(temp_plot.m_customColour);
		}
		if (temp_plot.m_plotInstances.instance(i).
		    isMissing(m_cIndex)) {
		  if (temp_plot.m_connectPoints[i] == true) {
		    drawDataPoint(x,y,prevx,prevy,temp_plot.m_shapeSize[i],
				  MISSING_SHAPE,gx);
		  } else {
		    drawDataPoint(x,y,temp_plot.m_shapeSize[i],
				  MISSING_SHAPE,gx);
		  }
		} else {
		  if (temp_plot.m_shapeType[i] == CONST_AUTOMATIC_SHAPE) {
		    if (temp_plot.m_connectPoints[i] == true) {
		      drawDataPoint(x,y,prevx,prevy,
				    temp_plot.m_shapeSize[i],j,gx);
		    } else {
		      drawDataPoint(x,y,temp_plot.m_shapeSize[i],j,gx);
		    }
		  } else {
		    if (temp_plot.m_connectPoints[i] == true) {
		      drawDataPoint(x,y,prevx,prevy,temp_plot.m_shapeSize[i],
				    temp_plot.m_shapeType[i],gx);
		    } else {
		      drawDataPoint(x,y,temp_plot.m_shapeSize[i],
				    temp_plot.m_shapeType[i],gx);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  /*
  public void determineAxisPositions(Graphics gx) {
    setFonts(gx);
    int mxs = m_XaxisStart;
    int mxe = m_XaxisEnd;
    int mys = m_YaxisStart;
    int mye = m_YaxisEnd;
    m_axisChanged = false;

    int h = this.getHeight();
    int w = this.getWidth();
    int hf = m_labelMetrics.getAscent();
    int mswx=0;
    int mswy=0;

    //      determineBounds();
    int fieldWidthX = (int)((Math.log(m_maxX)/Math.log(10)))+1;
    int precisionX = 1;
    if ((Math.abs(m_maxX-m_minX) < 1) && ((m_maxY-m_minX) != 0)) {
      precisionX = (int)Math.abs(((Math.log(Math.abs(m_maxX-m_minX)) / 
				   Math.log(10))))+1;
    }
    String maxStringX = Utils.doubleToString(m_maxX,
					     fieldWidthX+1+precisionX
					     ,precisionX);
    mswx = m_labelMetrics.stringWidth(maxStringX);
    int fieldWidthY = (int)((Math.log(m_maxY)/Math.log(10)))+1;
    int precisionY = 1;
    if (Math.abs((m_maxY-m_minY)) < 1 && ((m_maxY-m_minY) != 0)) {
      precisionY = (int)Math.abs(((Math.log(Math.abs(m_maxY-m_minY)) / 
				   Math.log(10))))+1;
    }
    String maxStringY = Utils.doubleToString(m_maxY,
					     fieldWidthY+1+precisionY
					     ,precisionY);
    String minStringY = Utils.doubleToString(m_minY,
					     fieldWidthY+1+precisionY
					     ,precisionY);

    if (m_plotInstances.attribute(m_yIndex).isNumeric()) {
      mswy = (m_labelMetrics.stringWidth(maxStringY) > 
	      m_labelMetrics.stringWidth(minStringY))
	? m_labelMetrics.stringWidth(maxStringY)
	: m_labelMetrics.stringWidth(minStringY);
    } else {
      mswy = m_labelMetrics.stringWidth("MM");
    }

    m_YaxisStart = m_axisPad;
    m_XaxisStart = 0+m_axisPad+m_tickSize+mswy;

    m_XaxisEnd = w-m_axisPad-(mswx/2);
      
    m_YaxisEnd = h-m_axisPad-(2 * hf)-m_tickSize;
    } */

  /**
   * Draws the axis and a spectrum if the colouring attribute is numeric
   * @param gx the graphics context
   */
  private void paintAxis(Graphics gx) {
    setFonts(gx);
    int mxs = m_XaxisStart;
    int mxe = m_XaxisEnd;
    int mys = m_YaxisStart;
    int mye = m_YaxisEnd;
    m_plotResize = false;

    int h = this.getHeight();
    int w = this.getWidth();
    int hf = m_labelMetrics.getAscent();
    int mswx=0;
    int mswy=0;

    //      determineBounds();
    int precisionXmax = 1;
    int precisionXmin = 1;
    int precisionXmid = 1;
    /*if ((Math.abs(m_maxX-m_minX) < 1) && ((m_maxY-m_minX) != 0)) {
      precisionX = (int)Math.abs(((Math.log(Math.abs(m_maxX-m_minX)) / 
				   Math.log(10))))+1;
				   } */
    int whole = (int)Math.abs(m_maxX);
    double decimal = Math.abs(m_maxX) - whole;
    int nondecimal;
    nondecimal = (whole > 0) 
      ? (int)(Math.log(whole) / Math.log(10))
      : 1;
    
    precisionXmax = (decimal > 0) 
      ? (int)Math.abs(((Math.log(Math.abs(m_maxX)) / 
				      Math.log(10))))+2
      : 1;
    if (precisionXmax > VisualizeUtils.MAX_PRECISION) {
      precisionXmax = 1;
    }

    String maxStringX = Utils.doubleToString(m_maxX,
					     nondecimal+1+precisionXmax
					     ,precisionXmax);

    whole = (int)Math.abs(m_minX);
    decimal = Math.abs(m_minX) - whole;
    nondecimal = (whole > 0) 
      ? (int)(Math.log(whole) / Math.log(10))
      : 1;
    precisionXmin = (decimal > 0) 
      ? (int)Math.abs(((Math.log(Math.abs(m_minX)) / 
				      Math.log(10))))+2
      : 1;
    if (precisionXmin > VisualizeUtils.MAX_PRECISION) {
      precisionXmin = 1;
    }
   
    String minStringX = Utils.doubleToString(m_minX,
					     nondecimal+1+precisionXmin,
					     precisionXmin);

    mswx = m_labelMetrics.stringWidth(maxStringX);

    int precisionYmax = 1;
    int precisionYmin = 1;
    int precisionYmid = 1;
    whole = (int)Math.abs(m_maxY);
    decimal = Math.abs(m_maxY) - whole;
    nondecimal = (whole > 0) 
      ? (int)(Math.log(whole) / Math.log(10))
      : 1;
    precisionYmax = (decimal > 0) 
      ? (int)Math.abs(((Math.log(Math.abs(m_maxY)) / 
				      Math.log(10))))+2
      : 1;
    if (precisionYmax > VisualizeUtils.MAX_PRECISION) {
      precisionYmax = 1;
    }
    
    String maxStringY = Utils.doubleToString(m_maxY,
					     nondecimal+1+precisionYmax
					     ,precisionYmax);


    whole = (int)Math.abs(m_minY);
    decimal = Math.abs(m_minY) - whole;
    nondecimal = (whole > 0) 
      ? (int)(Math.log(whole) / Math.log(10))
      : 1;
    precisionYmin = (decimal > 0) 
      ? (int)Math.abs(((Math.log(Math.abs(m_minY)) / 
				      Math.log(10))))+2
      : 1;
    if (precisionYmin > VisualizeUtils.MAX_PRECISION) {
      precisionYmin = 1;
    }
   
    String minStringY = Utils.doubleToString(m_minY,
					     nondecimal+1+precisionYmin
					     ,precisionYmin);

    if (m_plotInstances.attribute(m_yIndex).isNumeric()) {
      mswy = (m_labelMetrics.stringWidth(maxStringY) > 
	      m_labelMetrics.stringWidth(minStringY))
	? m_labelMetrics.stringWidth(maxStringY)
	: m_labelMetrics.stringWidth(minStringY);
      mswy += m_labelMetrics.stringWidth("M");
    } else {
      mswy = m_labelMetrics.stringWidth("MM");
    }

    m_YaxisStart = m_axisPad;
    m_XaxisStart = 0+m_axisPad+m_tickSize+mswy;

    m_XaxisEnd = w-m_axisPad-(mswx/2);
      
    m_YaxisEnd = h-m_axisPad-(2 * hf)-m_tickSize;

    // draw axis
    gx.setColor(m_axisColour);
    if (m_plotInstances.attribute(m_xIndex).isNumeric()) {
      if (w > (2 * mswx)) {
	
	gx.drawString(maxStringX, 
		      m_XaxisEnd-(mswx/2),
		      m_YaxisEnd+hf+m_tickSize);
	
	mswx = m_labelMetrics.stringWidth(minStringX);
	gx.drawString(minStringX,
		      (m_XaxisStart-(mswx/2)),
		      m_YaxisEnd+hf+m_tickSize);

	// draw the middle value
	if (w > (3 * mswx) && 
	    (m_plotInstances.attribute(m_xIndex).isNumeric())) {
	  double mid = m_minX+((m_maxX-m_minX)/2.0);
	   whole = (int)Math.abs(mid);
	   decimal = Math.abs(mid) - whole;
	   nondecimal = (whole > 0) 
	     ? (int)(Math.log(whole) / Math.log(10))
	     : 1;
	   precisionXmid = (decimal > 0) 
	     ? (int)Math.abs(((Math.log(Math.abs(mid)) / 
			       Math.log(10))))+2
	     : 1;
	   if (precisionXmid > VisualizeUtils.MAX_PRECISION) {
	     precisionXmid = 1;
	   }
	  
	  String maxString = Utils.doubleToString(mid,
						  nondecimal+1+precisionXmid,
						  precisionXmid);
	  int sw = m_labelMetrics.stringWidth(maxString);
	  double mx = m_XaxisStart+((double)(m_XaxisEnd-m_XaxisStart)/2.0);
	  gx.drawString(maxString,
			(int)(mx-(((double)sw)/2.0)),
			m_YaxisEnd+hf+m_tickSize);
	  gx.drawLine((int)mx,m_YaxisEnd,(int)mx,m_YaxisEnd+m_tickSize);
	}
      }
    } else {
      int numValues = m_plotInstances.attribute(m_xIndex).numValues();
      int div = (numValues % 2 > 0) ? (numValues/2)+1 : (numValues/2);
      int maxXStringWidth = (m_XaxisEnd - m_XaxisStart) / numValues;

      for (int i=0;i<numValues;i++) {
	String val = m_plotInstances.attribute(m_xIndex).value(i);
	int sw = m_labelMetrics.stringWidth(val);
	int rm;
	// truncate string if necessary
	if (sw > maxXStringWidth) {
	  int incr = (sw / val.length());
	  rm = (sw - maxXStringWidth) / incr;
	  if (rm == 0) {
	    rm = 1;
	  }
	  val = val.substring(0,val.length()-rm);
	  sw = m_labelMetrics.stringWidth(val);
	}
	if (i == 0) {
	  gx.drawString(val,(int)convertToPanelX(i),
			m_YaxisEnd+hf+m_tickSize);
	} else if (i == numValues -1) {
	  if ((i % 2) == 0) {
	    gx.drawString(val,
			  m_XaxisEnd-sw,
			  m_YaxisEnd+hf+m_tickSize);
	  } else {
	    gx.drawString(val,
			  m_XaxisEnd-sw,
			  m_YaxisEnd+(2*hf)+m_tickSize);
	  }
	} else {
	  if ((i % 2) == 0) {
	    gx.drawString(val,
			  (int)convertToPanelX(i)-(sw/2),
			  m_YaxisEnd+hf+m_tickSize);
	  } else {
	    gx.drawString(val,
			  (int)convertToPanelX(i)-(sw/2),
			  m_YaxisEnd+(2*hf)+m_tickSize);
	  }
	}
	gx.drawLine((int)convertToPanelX(i),
		    m_YaxisEnd,
		    (int)convertToPanelX(i),
		    m_YaxisEnd+m_tickSize);
      }
	
    }

    // draw the y axis
    if (m_plotInstances.attribute(m_yIndex).isNumeric()) {
      if (h > (2 * hf)) {
	gx.drawString(maxStringY, 
		      m_XaxisStart-mswy-m_tickSize,
		      m_YaxisStart+(hf));

	gx.drawString(minStringY,
		      (m_XaxisStart-mswy-m_tickSize),
		      m_YaxisEnd);

	// draw the middle value
	if (w > (3 * hf) && 
	    (m_plotInstances.attribute(m_yIndex).isNumeric())) {
	  double mid = m_minY+((m_maxY-m_minY)/2.0);
	  whole = (int)Math.abs(mid);
	  decimal = Math.abs(mid) - whole;
	  nondecimal = (whole > 0) 
	    ? (int)(Math.log(whole) / Math.log(10))
	    : 1;
	  precisionYmid = (decimal > 0) 
	    ? (int)Math.abs(((Math.log(Math.abs(mid)) / 
			      Math.log(10))))+2
	    : 1;
	  if (precisionYmid > VisualizeUtils.MAX_PRECISION) {
	    precisionYmid = 1;
	  }
	 
	  String maxString = Utils.doubleToString(mid,
						  nondecimal+1+precisionYmid,
						  precisionYmid);
	  int sw = m_labelMetrics.stringWidth(maxString);
	  double mx = m_YaxisStart+((double)(m_YaxisEnd-m_YaxisStart)/2.0);
	  gx.drawString(maxString,
			m_XaxisStart-sw-m_tickSize-1,
			(int)(mx+(((double)hf)/2.0)));
	  gx.drawLine(m_XaxisStart-m_tickSize,(int)mx,m_XaxisStart,(int)mx);
	}
      }
    } else {
      int numValues = m_plotInstances.attribute(m_yIndex).numValues();
      int div = ((numValues % 2) == 0) ? (numValues/2) : (numValues/2+1);
      int maxYStringHeight = (m_YaxisEnd - m_XaxisStart) / div;
      int sw = m_labelMetrics.stringWidth("M");
      for (int i=0;i<numValues;i++) {
	// can we at least print 2 characters
	if (maxYStringHeight >= (2*hf)) {
	  String val = m_plotInstances.attribute(m_yIndex).value(i);
	  int numPrint = ((maxYStringHeight/hf) > val.length()) ?
	    val.length() :
	    (maxYStringHeight/hf);
	    
	  for (int j=0;j<numPrint;j++) {
	    String ll = val.substring(j,j+1);
	    if (val.charAt(j) == '_' || val.charAt(j) == '-') {
	      ll = "|";
	    }
	    if (i == 0) {
	      gx.drawString(ll,m_XaxisStart-sw-m_tickSize-1,
			    (int)convertToPanelY(i)
			    -((numPrint-1)*hf)
			    +(j*hf)
			    +(hf/2));
	    } else if (i == (numValues-1)) {
	      if ((i % 2) == 0) {
		gx.drawString(ll,m_XaxisStart-sw-m_tickSize-1,
			      (int)convertToPanelY(i)
			      +(j*hf)
			      +(hf/2));
	      } else {
		gx.drawString(ll,m_XaxisStart-(2*sw)-m_tickSize-1,
			      (int)convertToPanelY(i)
			      +(j*hf)
			      +(hf/2));
	      }
	    } else {
	      if ((i % 2) == 0) {
		gx.drawString(ll,m_XaxisStart-sw-m_tickSize-1,
			      (int)convertToPanelY(i)
			      -(((numPrint-1)*hf)/2)
			      +(j*hf)
			      +(hf/2));
	      } else {
		gx.drawString(ll,m_XaxisStart-(2*sw)-m_tickSize-1,
			      (int)convertToPanelY(i)
			      -(((numPrint-1)*hf)/2)
			      +(j*hf)
			      +(hf/2));
	      }
	    }
	  }
	}
	gx.drawLine(m_XaxisStart-m_tickSize,
		    (int)convertToPanelY(i),
		    m_XaxisStart,
		    (int)convertToPanelY(i));
      }
    }

    gx.drawLine(m_XaxisStart,
		m_YaxisStart,
		m_XaxisStart,
		m_YaxisEnd);
    gx.drawLine(m_XaxisStart,
		m_YaxisEnd,
		m_XaxisEnd,
		m_YaxisEnd);

    if (m_XaxisStart != mxs || m_XaxisEnd != mxe ||
	m_YaxisStart != mys || m_YaxisEnd != mye) {
      m_plotResize = true;
    }
  }

  /**
   * Add more colours to the colour map
   */
  private void extendColourMap(int highest) {
    //System.err.println("Extending colour map");
    for (int i = m_colorList.size(); i < highest; i++) {
      Color pc = m_DefaultColors[i % 10];
      int ija =  i / 10;
      ija *= 2; 
      for (int j=0;j<ija;j++) {
	pc = pc.brighter();
      }
      
      m_colorList.addElement(pc);
    }
  }

  /**
   * Renders this component
   * @param gx the graphics context
   */
  public void paintComponent(Graphics gx) {
    super.paintComponent(gx);
    if (m_plotInstances != null 
	&& m_plotInstances.numInstances() > 0
	&& m_plotInstances.numAttributes() > 0) {
      if (m_plotCompanion != null) {
	m_plotCompanion.prePlot(gx);
      }

      m_JRand = new Random(m_JitterVal);
      paintAxis(gx);
      if (m_axisChanged || m_plotResize) {
	int x_range = m_XaxisEnd - m_XaxisStart;
	int y_range = m_YaxisEnd - m_YaxisStart;
	if (x_range < 10) {
	  x_range = 10;
	}
	if (y_range < 10) {
	  y_range = 10;
	}

	m_drawnPoints = new int[x_range + 1][y_range + 1];
	fillLookup();
	m_plotResize = false;
	m_axisChanged = false;
      }
      paintData(gx);
    }
  }
  
  protected static Color checkAgainstBackground(Color c, Color background) {
    if (background == null) {
      return c;
    }
    
    if (c.equals(background)) {
      int red = c.getRed();
      int blue = c.getBlue();
      int green = c.getGreen();
      red += (red < 128) ? (255 - red) / 2 : -(red / 2);
      blue += (blue < 128) ? (blue - red) / 2 : -(blue / 2);
      green += (green< 128) ? (255 - green) / 2 : -(green / 2);
      c = new Color(red, green, blue);
    }
    return c;
  }

  /**
   * Main method for testing this class
   * @param args arguments
   */
  public static void main(String [] args) {
    try {
      if (args.length < 1) {
	System.err.println("Usage : weka.gui.visualize.Plot2D "
			   +"<dataset> [<dataset> <dataset>...]");
	System.exit(1);
      }

      final javax.swing.JFrame jf = 
	new javax.swing.JFrame("Weka Explorer: Visualize");
      jf.setSize(500,400);
      jf.getContentPane().setLayout(new BorderLayout());
      final Plot2D p2 = new Plot2D();
      jf.getContentPane().add(p2, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	  public void windowClosing(java.awt.event.WindowEvent e) {
	    jf.dispose();
	    System.exit(0);
	  }
	});
      
      p2.addMouseListener(new MouseAdapter() {
	  public void mouseClicked(MouseEvent e) {
	    if ((e.getModifiers() & InputEvent.BUTTON1_MASK) ==
		InputEvent.BUTTON1_MASK) {
	      p2.searchPoints(e.getX(), e.getY(), false);
	    } else {
	      p2.searchPoints(e.getX(), e.getY(), true);
	    }
	  }
	});

      jf.setVisible(true);
      if (args.length >= 1) {
	for (int j = 0; j < args.length; j++) {
	  System.err.println("Loading instances from " + args[j]);
	  java.io.Reader r = new java.io.BufferedReader(
			     new java.io.FileReader(args[j]));
	  Instances i = new Instances(r);
	  i.setClassIndex(i.numAttributes()-1);
	  PlotData2D pd1 = new PlotData2D(i);

	  if (j == 0) {
	    pd1.setPlotName("Master plot");
	    p2.setMasterPlot(pd1);
	    p2.setXindex(2);
	    p2.setYindex(3);
	    p2.setCindex(i.classIndex());
	  } else {
	    pd1.setPlotName("Plot "+(j+1));
	    pd1.m_useCustomColour = true;
	    pd1.m_customColour = (j % 2 == 0) ? Color.red : Color.blue; 
	    p2.addPlot(pd1);
	  }
	}
      }
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
