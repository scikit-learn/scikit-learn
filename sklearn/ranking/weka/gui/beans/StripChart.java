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
 *    StripChart.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instance;
import weka.core.Instances;
import weka.gui.visualize.PrintableComponent;
import weka.gui.visualize.VisualizeUtils;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Image;
import java.beans.EventSetDescriptor;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.Enumeration;
import java.util.LinkedList;
import java.util.Random;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

/**
 * Bean that can display a horizontally scrolling strip chart. Can
 * display multiple plots simultaneously
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 4806 $
 */
public class StripChart
  extends JPanel
  implements ChartListener, InstanceListener, Visible,
	     BeanCommon, UserRequestAcceptor {

  /** for serialization */
  private static final long serialVersionUID = 1483649041577695019L;

  /** default colours for colouring lines */
  protected Color [] m_colorList = {Color.green,
				    Color.red,
				    Color.blue,
				    Color.cyan,
				    Color.pink,
				    new Color(255, 0, 255),
				    Color.orange,
				    new Color(255, 0, 0),
				    new Color(0, 255, 0),
				    Color.white};

  /** the background color. */
  protected Color m_BackgroundColor;
  
  /** the color of the legend panel's border. */
  protected Color m_LegendPanelBorderColor;
  
  /**
   * Class providing a panel for the plot.
   */
  private class StripPlotter
    extends JPanel {

    /** for serialization. */
    private static final long serialVersionUID = -7056271598761675879L;

    public void paintComponent(Graphics g) {
      super.paintComponent(g);
      if (m_osi != null) {
	g.drawImage(m_osi,0,0,this);
      }
    }
  }

  private transient JFrame m_outputFrame = null;
  private transient StripPlotter m_plotPanel = null;

  /**
   * The off screen image for rendering to.
   */
  private transient Image m_osi = null;

  /**
   * Width and height of the off screen image.
   */
  private int m_iheight;
  private int m_iwidth;

  /**
   * Max value for the y axis.
   */
  private double m_max = 1;

  /**
   * Min value for the y axis.
   */
  private double m_min = 0;

  /**
   * Scale update requested.
   */
  private boolean m_yScaleUpdate = false;
  private double m_oldMax;
  private double m_oldMin;

  private Font m_labelFont = new Font("Monospaced", Font.PLAIN, 10);
  private FontMetrics m_labelMetrics;

  //  private int m_plotCount = 0;

  private Vector m_legendText = new Vector();

  /**
   * Class providing a panel for displaying the y axis.
   */
  private class ScalePanel 
    extends JPanel {

    /** for serialization. */
    private static final long serialVersionUID = 6416998474984829434L;

    public void paintComponent(Graphics gx) {
      super.paintComponent(gx);
      if (m_labelMetrics == null) {
	m_labelMetrics = gx.getFontMetrics(m_labelFont);
      }
      gx.setFont(m_labelFont);
      int hf = m_labelMetrics.getAscent();
      String temp = ""+m_max;
      gx.setColor(m_colorList[m_colorList.length-1]);
      gx.drawString(temp, 1, hf-2);
      temp = ""+(m_min + ((m_max - m_min)/2.0));
      gx.drawString(temp, 1, (this.getHeight() / 2)+(hf / 2));
      temp = ""+m_min;
      gx.drawString(temp, 1, this.getHeight()-1);
    }
  };
  
  /** the scale. */
  private ScalePanel m_scalePanel = new ScalePanel();

  /**
   * Class providing a panel for the legend.
   */
  private class LegendPanel 
    extends JPanel {

    /** for serialization. */
    private static final long serialVersionUID = 7713986576833797583L;

    public void paintComponent(Graphics gx) {
      super.paintComponent(gx);

      if (m_labelMetrics == null) {
	m_labelMetrics = gx.getFontMetrics(m_labelFont);
      }
      int hf = m_labelMetrics.getAscent();
      int x = 10; int y = hf+15;
      gx.setFont(m_labelFont);
      for (int i = 0; i < m_legendText.size(); i++) {
	String temp = (String)m_legendText.elementAt(i);
	gx.setColor(m_colorList[(i % m_colorList.length)]);
	gx.drawString(temp,x,y);
	y+=hf;
      }
      StripChart.this.revalidate();
    }
  };
  
  /** the legend. */
  private LegendPanel m_legendPanel = new LegendPanel();

  /**
   * Holds the data to be plotted. Entries in the list are arrays of
   * y points.
   */
  private LinkedList m_dataList = new LinkedList();
  //  private double [] m_dataPoint = new double[1];
  private double [] m_previousY = new double[1];

  private transient Thread m_updateHandler;

  protected BeanVisual m_visual =
    new BeanVisual("StripChart",
		   BeanVisual.ICON_PATH+"StripChart.gif",
		   BeanVisual.ICON_PATH+"StripChart_animated.gif");

  private Object m_listenee = null;
  private transient weka.gui.Logger m_log = null;

  /**
   * Print x axis labels every m_xValFreq points
   */
  private int m_xValFreq = 500;
  private int m_xCount = 0;

  /**
   * Shift the plot by this many pixels every time a point is plotted
   */
  private int m_refreshWidth = 1;

  /**
   * Plot every m_refreshFrequency'th point
   */
  private int m_refreshFrequency = 5;

  /** the class responsible for printing */
  protected PrintableComponent m_Printer = null;

  public StripChart() {

    //    m_plotPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));

    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);

    // start a thread to handle plot updates
    initPlot();
  }

  /**
   * Set a custom (descriptive) name for this bean
   *
   * @param name the name to use
   */
  public void setCustomName(String name) {
    m_visual.setText(name);
  }

  /**
   * Get the custom (descriptive) name for this bean (if one has been set)
   *
   * @return the custom name (or the default name)
   */
  public String getCustomName() {
    return m_visual.getText();
  }

  /**
   * Global info for this bean
   *
   * @return a <code>String</code> value
   */
  public String globalInfo() {
    return "Visualize incremental classifier performance as a scrolling plot.";
  }

  /**
   * GUI Tip text
   *
   * @return a <code>String</code> value
   */
  public String xLabelFreqTipText() {
    return "Show x axis labels this often";
  }

  /**
   * Set the frequency for printing x label values
   *
   * @param freq an <code>int</code> value
   */
  public void setXLabelFreq(int freq) {
    m_xValFreq = freq;
    if (getGraphics() != null)
      setRefreshWidth();
  }

  /**
   * Get the frequency by which x axis values are printed
   *
   * @return an <code>int</code> value
   */
  public int getXLabelFreq() {
    return m_xValFreq;
  }

  /**
   * GUI Tip text
   *
   * @return a <code>String</code> value
   */
  public String refreshFreqTipText() {
    return "Plot every x'th data point";
  }

  /**
   * Set how often (in x axis points) to refresh the display
   *
   * @param freq an <code>int</code> value
   */
  public void setRefreshFreq(int freq) {
    m_refreshFrequency = freq;
    if (getGraphics() != null)
      setRefreshWidth();
  }

  /**
   * Get the refresh frequency
   *
   * @return an <code>int</code> value
   */
  public int getRefreshFreq() {
    return m_refreshFrequency;
  }

  private void setRefreshWidth() {
    m_refreshWidth = 1;
    if (m_labelMetrics == null) {
      getGraphics().setFont(m_labelFont);
      m_labelMetrics = getGraphics().getFontMetrics(m_labelFont);
    }

    int refWidth = m_labelMetrics.stringWidth("99000");
    // compute how often x label will be rendered
    int z = (getXLabelFreq() / getRefreshFreq());
    if (z < 1) {
      z = 1;
    }

    if (z * m_refreshWidth < refWidth+5) {
      m_refreshWidth *= (((refWidth+5) / z) + 1) ;
    }
  }

  /**
   * Provide some necessary initialization after object has
   * been deserialized.
   *
   * @param ois an <code>ObjectInputStream</code> value
   * @exception IOException if an error occurs
   * @exception ClassNotFoundException if an error occurs
   */
  private void readObject(ObjectInputStream ois)
    throws IOException, ClassNotFoundException {
    try {
      ois.defaultReadObject();
      initPlot();
      //      startHandler();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  /**
   * Loads properties from properties file.
   * 
   * @see	KnowledgeFlowApp#BEAN_PROPERTIES
   */
  private void setProperties() {
    String 	key;
    String 	color;

    // background color
    key   = this.getClass().getName() + ".backgroundColour";
    color = KnowledgeFlowApp.BEAN_PROPERTIES.getProperty(key);
    m_BackgroundColor = Color.BLACK;
    if (color != null)
      m_BackgroundColor = VisualizeUtils.processColour(color, m_BackgroundColor);

    // legend color (border)
    key   = m_legendPanel.getClass().getName() + ".borderColour";
    color = KnowledgeFlowApp.BEAN_PROPERTIES.getProperty(key);
    m_LegendPanelBorderColor = Color.BLUE;
    if (color != null)
      m_LegendPanelBorderColor = VisualizeUtils.processColour(color, m_LegendPanelBorderColor);
  }

  private void initPlot() {
    setProperties();
    m_plotPanel = new StripPlotter();
    m_plotPanel.setBackground(m_BackgroundColor);
    m_scalePanel.setBackground(m_BackgroundColor);
    m_legendPanel.setBackground(m_BackgroundColor);
    m_xCount = 0;
  }

  private void startHandler() {
    if (m_updateHandler == null) {
      m_updateHandler = new Thread() {
	  private double [] dataPoint;
	  public void run() {
	    while (true) {
	      if (m_outputFrame != null) {
		synchronized(m_dataList) {
		  while(m_dataList.isEmpty()) {
		  //		  while (m_dataList.empty()) {
		    try {
		      m_dataList.wait();
		    } catch (InterruptedException ex) {
		      return;
		    }
		  }
		  dataPoint = (double [])m_dataList.remove(0);
		  //dataPoint  = (double [])m_dataList.pop();
		}
		if (m_outputFrame != null) {
		  StripChart.this.updateChart(dataPoint);
		}
	      }
	    }
	  }
	};
      //      m_updateHandler.setPriority(Thread.MIN_PRIORITY);
      m_updateHandler.start();
    }
  }

  /**
   * Popup the chart panel
   */
  public void showChart() {
    if (m_outputFrame == null) {
      m_outputFrame = new JFrame("Strip Chart");
      m_outputFrame.getContentPane().setLayout(new BorderLayout());
      JPanel panel = new JPanel(new BorderLayout());
      new PrintableComponent(panel);
      m_outputFrame.getContentPane().add(panel, BorderLayout.CENTER);
      panel.add(m_legendPanel, BorderLayout.WEST);
      panel.add(m_plotPanel, BorderLayout.CENTER);
      panel.add(m_scalePanel, BorderLayout.EAST);
      m_legendPanel.setMinimumSize(new Dimension(100,getHeight()));
      m_legendPanel.setPreferredSize(new Dimension(100,getHeight()));
      m_scalePanel.setMinimumSize(new Dimension(30, getHeight()));
      m_scalePanel.setPreferredSize(new Dimension(30, getHeight()));
      Font lf = new Font("Monospaced", Font.PLAIN, 12);
      m_legendPanel.setBorder(BorderFactory.
			      createTitledBorder(BorderFactory.
				    createEtchedBorder(Color.gray,
						       Color.darkGray),
				    "Legend" ,
				    TitledBorder.CENTER,
				    TitledBorder.DEFAULT_POSITION, lf,
				    m_LegendPanelBorderColor));
      m_outputFrame.addWindowListener(new java.awt.event.WindowAdapter() {
	  public void windowClosing(java.awt.event.WindowEvent e) {
	    if (m_updateHandler != null) {
	      System.err.println("Interrupting");
	      m_updateHandler.interrupt();
	      m_updateHandler = null;
	    }
	    synchronized (m_dataList) {
	      m_dataList = new LinkedList();
	    }
	    m_outputFrame.dispose();
	    m_outputFrame = null;
	  }
	});
      m_outputFrame.pack();
      m_outputFrame.setSize(600,150);
      m_outputFrame.setResizable(false);
      m_outputFrame.setVisible(true);
      int iwidth = m_plotPanel.getWidth();
      int iheight = m_plotPanel.getHeight();
      m_osi = m_plotPanel.createImage(iwidth, iheight);
      Graphics m = m_osi.getGraphics();
      m.setColor(m_BackgroundColor);
      m.fillRect(0,0,iwidth,iheight);
      m_previousY[0] = -1;
      setRefreshWidth();
      if (m_updateHandler == null) {
	System.err.println("Starting handler");
	startHandler();
      }
    } else {
      m_outputFrame.toFront();
    }
  }

  private int convertToPanelY(double yval) {
    int height = m_plotPanel.getHeight();
    double temp = (yval - m_min) / (m_max - m_min);
    temp = temp * height;
    temp = height - temp;
    return (int)temp;
  }

  /**
   * Update the plot
   *
   * @param dataPoint contains y values to plot
   */
  protected void updateChart(double [] dataPoint) {
    if (m_previousY[0] == -1) {
      int iw = m_plotPanel.getWidth();
      int ih = m_plotPanel.getHeight();
      m_osi = m_plotPanel.createImage(iw, ih);
      Graphics m = m_osi.getGraphics();
      m.setColor(m_BackgroundColor);
      m.fillRect(0,0,iw,ih);
      m_previousY[0] = convertToPanelY(0);
      m_iheight = ih; m_iwidth = iw;
    }

    if (dataPoint.length-1 != m_previousY.length) {
      m_previousY = new double [dataPoint.length-1];
      //      m_plotCount = 0;
      for (int i = 0; i < dataPoint.length-1; i++) {
	m_previousY[i] = convertToPanelY(0);
      }
    }

    Graphics osg = m_osi.getGraphics();
    Graphics g = m_plotPanel.getGraphics();

    osg.copyArea(m_refreshWidth,0,m_iwidth-m_refreshWidth,
		 m_iheight,-m_refreshWidth,0);
    osg.setColor(m_BackgroundColor);
    osg.fillRect(m_iwidth-m_refreshWidth,0, m_iwidth, m_iheight);

    // paint the old scale onto the plot if a scale update has occured
    if (m_yScaleUpdate) {
      String maxVal = numToString(m_oldMax);
      String minVal = numToString(m_oldMin);
      String midVal = numToString((m_oldMax - m_oldMin) / 2.0);
      if (m_labelMetrics == null) {
	m_labelMetrics = g.getFontMetrics(m_labelFont);
      }
      osg.setFont(m_labelFont);
      int wmx = m_labelMetrics.stringWidth(maxVal);
      int wmn = m_labelMetrics.stringWidth(minVal);
      int wmd = m_labelMetrics.stringWidth(midVal);

      int hf = m_labelMetrics.getAscent();
      osg.setColor(m_colorList[m_colorList.length-1]);
      osg.drawString(maxVal, m_iwidth-wmx, hf-2);
      osg.drawString(midVal, m_iwidth-wmd, (m_iheight / 2)+(hf / 2));
      osg.drawString(minVal, m_iwidth-wmn, m_iheight-1);
      m_yScaleUpdate = false;
    }

    double pos;
    for (int i = 0; i < dataPoint.length-1; i++) {
      osg.setColor(m_colorList[(i % m_colorList.length)]);
      pos = convertToPanelY(dataPoint[i]);
      osg.drawLine(m_iwidth-m_refreshWidth, (int)m_previousY[i],
		   m_iwidth-1, (int)pos);
      m_previousY[i] = pos;
      if (dataPoint[dataPoint.length-1] % m_xValFreq == 0) {
	// draw the actual y value onto the plot for this curve
	String val = numToString(dataPoint[i]);
	if (m_labelMetrics == null) {
	  m_labelMetrics = g.getFontMetrics(m_labelFont);
	}
	int hf = m_labelMetrics.getAscent();
	if (pos - hf < 0) {
	  pos += hf;
	}
	int w = m_labelMetrics.stringWidth(val);
	osg.setFont(m_labelFont);
	osg.drawString(val, m_iwidth-w, (int)pos);
      }
    }

    // last element in the data point array contains the data point number
    if (dataPoint[dataPoint.length-1] % m_xValFreq == 0) {

      String xVal = ""+(int)dataPoint[dataPoint.length-1];
      osg.setColor(m_colorList[m_colorList.length-1]);
      int w = m_labelMetrics.stringWidth(xVal);
      osg.setFont(m_labelFont);
      osg.drawString(xVal, m_iwidth-w, m_iheight - 1);
    }
    g.drawImage(m_osi,0,0,m_plotPanel);
    //    System.err.println("Finished");
    //    m_plotCount++;
  }

  private static String numToString(double num) {
    int precision = 1;
    int whole = (int)Math.abs(num);
    double decimal = Math.abs(num) - whole;
    int nondecimal;
    nondecimal = (whole > 0)
      ? (int)(Math.log(whole) / Math.log(10))
      : 1;

    precision = (decimal > 0)
      ? (int)Math.abs(((Math.log(Math.abs(num)) /
				      Math.log(10))))+2
      : 1;
    if (precision > 5) {
      precision = 1;
    }

    String numString = weka.core.Utils.doubleToString(num,
						      nondecimal+1+precision
						      ,precision);

    return numString;
  }

  ChartEvent m_ce = new ChartEvent(this);
  double [] m_dataPoint = null;
  public void acceptInstance(InstanceEvent e) {
    if (e.getStatus() == InstanceEvent.FORMAT_AVAILABLE) {
      Instances structure = e.getStructure();
      m_legendText = new Vector();
      m_max = 1.0;
      m_min = 0;
      int i = 0;
      for (i =0; i < structure.numAttributes(); i++) {
	if (i > 10) {
	  i--;
	  break;
	}
	m_legendText.addElement(structure.attribute(i).name());
	m_legendPanel.repaint();
	m_scalePanel.repaint();
      }
      m_dataPoint = new double[i];
      m_xCount = 0;
      return;
    }

    // process data point
    Instance inst = e.getInstance();
    for (int i = 0; i < m_dataPoint.length; i++) {
      if (!inst.isMissing(i)) {
	m_dataPoint[i] = inst.value(i);
      }
    }
    acceptDataPoint(m_dataPoint);
    m_xCount++;
  }

  /**
   * Accept a data point (encapsulated in a chart event) to plot
   *
   * @param e a <code>ChartEvent</code> value
   */
  public void acceptDataPoint(ChartEvent e) {
    if (e.getReset()) {
      m_xCount = 0;
      m_max = 1;
      m_min = 0;
    }
    if (m_outputFrame != null) {
      boolean refresh = false;
      if (e.getLegendText() != null & e.getLegendText() != m_legendText) {
	m_legendText = e.getLegendText();
	refresh = true;
      }

      if (e.getMin() != m_min || e.getMax() != m_max) {
	m_oldMax = m_max; m_oldMin = m_min;
	m_max = e.getMax();
	m_min = e.getMin();
	refresh = true;
	m_yScaleUpdate = true;
      }

      if (refresh) {
	m_legendPanel.repaint();
	m_scalePanel.repaint();
      }

      acceptDataPoint(e.getDataPoint());
    }
    m_xCount++;
  }

  /**
   * Accept a data point to plot
   *
   * @param dataPoint a <code>double[]</code> value
   */
  public void acceptDataPoint(double [] dataPoint) {

    if (m_outputFrame != null && (m_xCount % m_refreshFrequency == 0 )) {
      double [] dp = new double[dataPoint.length+1];
      dp[dp.length-1] = m_xCount;
      System.arraycopy(dataPoint, 0, dp, 0, dataPoint.length);
      // check for out of scale values
      for (int i = 0; i < dataPoint.length; i++) {
	if (dataPoint[i] < m_min) {
	  m_oldMin = m_min; m_min = dataPoint[i];
	  m_yScaleUpdate = true;
	}

	if (dataPoint[i] > m_max) {
	  m_oldMax = m_max; m_max = dataPoint[i];
	  m_yScaleUpdate = true;
	}
      }
      if (m_yScaleUpdate) {
	m_scalePanel.repaint();
	m_yScaleUpdate = false;
      }
      synchronized(m_dataList) {
	m_dataList.add(m_dataList.size(), dp);
	//	m_dataList.push(dp);
	m_dataList.notifyAll();
	/*	if (m_dataList.size() != 0) {
	  System.err.println("***** "+m_dataList.size());
	  } */

	//      System.err.println(m_xCount);
      }
    }
  }

  /**
   * Set the visual appearance of this bean
   *
   * @param newVisual a <code>BeanVisual</code> value
   */
  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  /**
   * Get the visual appearance of this bean
   */
  public BeanVisual getVisual() {
    return m_visual;
  }

  /**
   * Use the default visual appearance for this bean
   */
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"StripChart.gif",
		       BeanVisual.ICON_PATH+"StripChart_animated.gif");
  }

  /**
   * Stop any processing that the bean might be doing.
   */
  public void stop() {
    // tell the listenee (upstream bean) to stop
    if (m_listenee instanceof BeanCommon) {
      ((BeanCommon)m_listenee).stop();
    }
  }
  
  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return (m_updateHandler != null);
  }

  /**
   * Set a logger
   *
   * @param logger a <code>weka.gui.Logger</code> value
   */
  public void setLog(weka.gui.Logger logger) {
    m_log = logger;
  }

  /**
   * Returns true if, at this time,
   * the object will accept a connection via the named event
   *
   * @param eventName the name of the event
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {
    if (m_listenee == null) {
      return true;
    }
    return false;
  }

  /**
   * Returns true if, at this time,
   * the object will accept a connection according to the supplied
   * EventSetDescriptor
   *
   * @param esd the EventSetDescriptor
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(EventSetDescriptor esd) {
    return connectionAllowed(esd.getName());
  }

  /**
   * Notify this object that it has been registered as a listener with
   * a source for recieving events described by the named event
   * This object is responsible for recording this fact.
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public void connectionNotification(String eventName, Object source) {
    if (connectionAllowed(eventName)) {
      m_listenee = source;
    }
  }

  /**
   * Notify this object that it has been deregistered as a listener with
   * a source for named event. This object is responsible
   * for recording this fact.
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public void disconnectionNotification(String eventName, Object source) {
    m_listenee = null;
  }

  /**
   * Describe <code>enumerateRequests</code> method here.
   *
   * @return an <code>Enumeration</code> value
   */
  public Enumeration enumerateRequests() {
    Vector newVector = new Vector(0);
    newVector.addElement("Show chart");
    return newVector.elements();
  }

  /**
   * Describe <code>performRequest</code> method here.
   *
   * @param request a <code>String</code> value
   * @exception IllegalArgumentException if an error occurs
   */
  public void performRequest(String request) {
    if (request.compareTo("Show chart") == 0) {
      showChart();
    } else {
      throw new
	IllegalArgumentException(request
				 + " not supported (StripChart)");
    }
  }

  /**
   * Tests out the StripChart from the command line
   *
   * @param args ignored
   */
  public static void main(String [] args) {

    try {
      final javax.swing.JFrame jf =
	new javax.swing.JFrame("Weka Knowledge Flow : StipChart");
      jf.getContentPane().setLayout(new BorderLayout());
      final StripChart jd = new StripChart();
      jf.getContentPane().add(jd, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	public void windowClosing(java.awt.event.WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
      jd.showChart();
      Random r = new Random(1);
      for (int i = 0; i < 1020; i++) {
	double [] pos = new double[1];
	pos[0] = r.nextDouble();
	jd.acceptDataPoint(pos);
      }
      System.err.println("Done sending data");
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
