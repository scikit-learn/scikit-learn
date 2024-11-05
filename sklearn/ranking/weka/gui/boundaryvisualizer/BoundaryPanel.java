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
 *   BoundaryPanel.java
 *   Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.boundaryvisualizer;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Utils;
import weka.gui.visualize.JPEGWriter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.Iterator;
import java.util.Locale;
import java.util.Random;
import java.util.Vector;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.plugins.jpeg.JPEGImageWriteParam;
import javax.imageio.stream.ImageOutputStream;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.ToolTipManager;

/**
 * BoundaryPanel. A class to handle the plotting operations
 * associated with generating a 2D picture of a classifier's decision
 * boundaries.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 5987 $
 * @since 1.0
 * @see JPanel
 */
public class BoundaryPanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -8499445518744770458L;
  
  /** default colours for classes */
  public static final Color [] DEFAULT_COLORS = {
    Color.red,
    Color.green,
    Color.blue,
    new Color(0, 255, 255), // cyan
    new Color(255, 0, 255), // pink
    new Color(255, 255, 0), // yellow
    new Color(255, 255, 255), //white
    new Color(0, 0, 0)};
    
  /** The distance we can click away from a point in the GUI and still remove it. */
  public static final double REMOVE_POINT_RADIUS = 7.0;

  protected FastVector m_Colors = new FastVector();

  /** training data */
  protected Instances m_trainingData;

  /** distribution classifier to use */
  protected Classifier m_classifier;

  /** data generator to use */
  protected DataGenerator m_dataGenerator;

  /** index of the class attribute */
  private int m_classIndex = -1;

  // attributes for visualizing on
  protected int m_xAttribute;
  protected int m_yAttribute;

  // min, max and ranges of these attributes
  protected double m_minX;
  protected double m_minY;
  protected double m_maxX;
  protected double m_maxY;
  private double m_rangeX;
  private double m_rangeY;

  // pixel width and height in terms of attribute values
  protected double m_pixHeight;
  protected double m_pixWidth;

  /** used for offscreen drawing */
  protected Image m_osi = null;

  // width and height of the display area
  protected int m_panelWidth;
  protected int m_panelHeight;

  // number of samples to take from each region in the fixed dimensions
  protected int m_numOfSamplesPerRegion = 2;

  // number of samples per kernel = base ^ (# non-fixed dimensions)
  protected int m_numOfSamplesPerGenerator;
  protected double m_samplesBase = 2.0;

  /** listeners to be notified when plot is complete */
  private Vector m_listeners = new Vector();

  /**
   * small inner class for rendering the bitmap on to
   */
  private class PlotPanel
    extends JPanel {

    /** for serialization */
    private static final long serialVersionUID = 743629498352235060L;
    
    public PlotPanel() {
      this.setToolTipText("");
    }
    
    public void paintComponent(Graphics g) {
      super.paintComponent(g);
      if (m_osi != null) {
	g.drawImage(m_osi,0,0,this);
      }
    }
    
    public String getToolTipText(MouseEvent event) {
      if (m_probabilityCache == null) {
	return null;
      }
      
      if (m_probabilityCache[event.getY()][event.getX()] == null) {
	return null;
      }
      
      String pVec = "(X: "
	+Utils.doubleToString(convertFromPanelX((double)event.getX()), 2)
	+" Y: "
	+Utils.doubleToString(convertFromPanelY((double)event.getY()), 2)+") ";
      // construct a string holding the probability vector
      for (int i = 0; i < m_trainingData.classAttribute().numValues(); i++) {
	pVec += 
	  Utils.
	  doubleToString(m_probabilityCache[event.getY()][event.getX()][i],
			 3)+" ";
      }
      return pVec;
    }
  }

  /** the actual plotting area */
  private PlotPanel m_plotPanel = new PlotPanel();

  /** thread for running the plotting operation in */
  private Thread m_plotThread = null;

  /** Stop the plotting thread */
  protected boolean m_stopPlotting = false;

  /** Stop any replotting threads */
  protected boolean m_stopReplotting = false;

  // Used by replotting threads to pause and resume the main plot thread
  private Double m_dummy = new Double(1.0);
  private boolean m_pausePlotting = false;
  /** what size of tile is currently being plotted */
  private int m_size = 1;
  /** is the main plot thread performing the initial coarse tiling */
  private boolean m_initialTiling;

  /** A random number generator  */
  private Random m_random = null;

  /** cache of probabilities for fast replotting */
  protected double [][][] m_probabilityCache;

  /** plot the training data */
  protected boolean m_plotTrainingData = true;

  /**
   * Creates a new <code>BoundaryPanel</code> instance.
   *
   * @param panelWidth the width in pixels of the panel
   * @param panelHeight the height in pixels of the panel
   */
  public BoundaryPanel(int panelWidth, int panelHeight) {
    ToolTipManager.sharedInstance().setDismissDelay(Integer.MAX_VALUE);
    m_panelWidth = panelWidth;
    m_panelHeight = panelHeight;
    setLayout(new BorderLayout());
    m_plotPanel.setMinimumSize(new Dimension(m_panelWidth, m_panelHeight));
    m_plotPanel.setPreferredSize(new Dimension(m_panelWidth, m_panelHeight));
    m_plotPanel.setMaximumSize(new Dimension(m_panelWidth, m_panelHeight));
    add(m_plotPanel, BorderLayout.CENTER);
    setPreferredSize(m_plotPanel.getPreferredSize());
    setMaximumSize(m_plotPanel.getMaximumSize());
    setMinimumSize(m_plotPanel.getMinimumSize());

    m_random = new Random(1);
    for (int i = 0; i < DEFAULT_COLORS.length; i++) {
      m_Colors.addElement(new Color(DEFAULT_COLORS[i].getRed(),
				    DEFAULT_COLORS[i].getGreen(),
				    DEFAULT_COLORS[i].getBlue()));
    }
    m_probabilityCache = new double[m_panelHeight][m_panelWidth][];
    
  }

  /**
   * Set the number of points to uniformly sample from a region (fixed
   * dimensions).
   *
   * @param num an <code>int</code> value
   */
  public void setNumSamplesPerRegion(int num) {
    m_numOfSamplesPerRegion = num;
  }

  /**
   * Get the number of points to sample from a region (fixed dimensions).
   *
   * @return an <code>int</code> value
   */
  public int getNumSamplesPerRegion() {
    return m_numOfSamplesPerRegion;
  }

  /**
   * Set the base for computing the number of samples to obtain from each
   * generator. number of samples = base ^ (# non fixed dimensions)
   *
   * @param ksb a <code>double</code> value
   */
  public void setGeneratorSamplesBase(double ksb) {
    m_samplesBase = ksb;
  }

  /**
   * Get the base used for computing the number of samples to obtain from
   * each generator
   *
   * @return a <code>double</code> value
   */
  public double getGeneratorSamplesBase() {
    return m_samplesBase;
  }

  /**
   * Set up the off screen bitmap for rendering to
   */
  protected void initialize() {
    int iwidth = m_plotPanel.getWidth();
    int iheight = m_plotPanel.getHeight();
    //    System.err.println(iwidth+" "+iheight);
    m_osi = m_plotPanel.createImage(iwidth, iheight);
    Graphics m = m_osi.getGraphics();
    m.fillRect(0,0,iwidth,iheight);
  }

  /**
   * Stop the plotting thread
   */
  public void stopPlotting() {
    m_stopPlotting = true;
    try {
    	m_plotThread.join(100);
    } catch (Exception e){};
  }

  /** Set up the bounds of our graphic based by finding the smallest reasonable
      area in the instance space to surround our data points.
  */
  public void computeMinMaxAtts() {
    m_minX = Double.MAX_VALUE;
    m_minY = Double.MAX_VALUE;
    m_maxX = Double.MIN_VALUE;
    m_maxY = Double.MIN_VALUE;
    
    boolean allPointsLessThanOne = true;
    
    if (m_trainingData.numInstances() == 0) {
      m_minX = m_minY = 0.0;
      m_maxX = m_maxY = 1.0;
    }
    else
    {
	for (int i = 0; i < m_trainingData.numInstances(); i++) {
		Instance inst = m_trainingData.instance(i);
		double x = inst.value(m_xAttribute);
		double y = inst.value(m_yAttribute);
		if (!Utils.isMissingValue(x) && !Utils.isMissingValue(y)) {
			if (x < m_minX) {
			m_minX = x;
			}
			if (x > m_maxX) {
			m_maxX = x;
			}
		
			if (y < m_minY) {
			m_minY = y;
			}
			if (y > m_maxY) {
			m_maxY = y;
			}
			if (x > 1.0 || y > 1.0)
				allPointsLessThanOne = false;
		}
	}
    }
    
    if (m_minX == m_maxX)
    	m_minX = 0;
    if (m_minY == m_maxY)
    	m_minY = 0;
    if (m_minX == Double.MAX_VALUE)
    	m_minX = 0;
    if (m_minY == Double.MAX_VALUE)
    	m_minY = 0;
    if (m_maxX == Double.MIN_VALUE)
    	m_maxX = 1;
    if (m_maxY == Double.MIN_VALUE)
    	m_maxY = 1;
    if (allPointsLessThanOne) {
    	m_minX = m_minY = 0.0;
    	m_maxX = m_maxY = 1.0;
    }
    
    
    
    m_rangeX = (m_maxX - m_minX);
    m_rangeY = (m_maxY - m_minY);
    
    m_pixWidth = m_rangeX / (double)m_panelWidth;
    m_pixHeight = m_rangeY / (double) m_panelHeight;
  }

  /**
   * Return a random x attribute value contained within
   * the pix'th horizontal pixel
   *
   * @param pix the horizontal pixel number
   * @return a value in attribute space
   */
  private double getRandomX(int pix) {

    double minPix =  m_minX + (pix * m_pixWidth);

    return minPix + m_random.nextDouble() * m_pixWidth;
  }

  /**
   * Return a random y attribute value contained within
   * the pix'th vertical pixel
   *
   * @param pix the vertical pixel number
   * @return a value in attribute space
   */
  private double getRandomY(int pix) {
    
    double minPix = m_minY + (pix * m_pixHeight);
    
    return minPix +  m_random.nextDouble() * m_pixHeight;
  }
  
  /**
   * Start the plotting thread
   *
   * @exception Exception if an error occurs
   */
  public void start() throws Exception {
    m_numOfSamplesPerGenerator = 
      (int)Math.pow(m_samplesBase, m_trainingData.numAttributes()-3);

    m_stopReplotting = true;
    if (m_trainingData == null) {
      throw new Exception("No training data set (BoundaryPanel)");
    }
    if (m_classifier == null) {
      throw new Exception("No classifier set (BoundaryPanel)");
    }
    if (m_dataGenerator == null) {
      throw new Exception("No data generator set (BoundaryPanel)");
    }
    if (m_trainingData.attribute(m_xAttribute).isNominal() || m_trainingData.attribute(m_xAttribute).isRanking() ||
	m_trainingData.attribute(m_yAttribute).isNominal() || m_trainingData.attribute(m_yAttribute).isRanking()) {
      throw new Exception("Visualization dimensions must be numeric "
			  +"(BoundaryPanel)");
    }
    
    computeMinMaxAtts();
    
    startPlotThread();
    /*if (m_plotThread == null) {
      m_plotThread = new PlotThread();
      m_plotThread.setPriority(Thread.MIN_PRIORITY);
      m_plotThread.start();
    }*/
  }
  
  // Thread for main plotting operation
  protected class PlotThread extends Thread {
    double [] m_weightingAttsValues;
    boolean [] m_attsToWeightOn;
    double [] m_vals;
    double [] m_dist;
    Instance m_predInst;
    public void run() {

      m_stopPlotting = false;
      try {
        initialize();
        repaint();
        
        // train the classifier
        m_probabilityCache = new double[m_panelHeight][m_panelWidth][];
        m_classifier.buildClassifier(m_trainingData);
        
        // build DataGenerator
        m_attsToWeightOn = new boolean[m_trainingData.numAttributes()];
        m_attsToWeightOn[m_xAttribute] = true;
        m_attsToWeightOn[m_yAttribute] = true;
	      
        m_dataGenerator.setWeightingDimensions(m_attsToWeightOn);
	      
        m_dataGenerator.buildGenerator(m_trainingData);

        // generate samples
        m_weightingAttsValues = new double [m_attsToWeightOn.length];
        m_vals = new double[m_trainingData.numAttributes()];
        m_predInst = new DenseInstance(1.0, m_vals);
        m_predInst.setDataset(m_trainingData);

	
        m_size = 1 << 4;  // Current sample region size
	
	m_initialTiling = true;
        // Display the initial coarse image tiling.
      abortInitial:
        for (int i = 0; i <= m_panelHeight; i += m_size) {   
          for (int j = 0; j <= m_panelWidth; j += m_size) {   
	    if (m_stopPlotting) {
	      break abortInitial;
	    }
	    if (m_pausePlotting) {
	      synchronized (m_dummy) {
		try {
		  m_dummy.wait();
		} catch (InterruptedException ex) {
		  m_pausePlotting = false;
		}
	      }
	    }
            plotPoint(j, i, m_size, m_size, 
		      calculateRegionProbs(j, i), (j == 0));
          }
        }
	if (!m_stopPlotting) {
	  m_initialTiling = false;
	}
        
        // Sampling and gridding loop
        int size2 = m_size / 2;
        abortPlot: 
        while (m_size > 1) { // Subdivide down to the pixel level
          for (int i = 0; i <= m_panelHeight; i += m_size) {
            for (int j = 0; j <= m_panelWidth; j += m_size) {
              if (m_stopPlotting) {
                break abortPlot;
              }
	      if (m_pausePlotting) {
		synchronized (m_dummy) {
		  try {
		    m_dummy.wait();
		  } catch (InterruptedException ex) {
		    m_pausePlotting = false;
		  }
		}
	      }
              boolean update = (j == 0 && i % 2 == 0);
              // Draw the three new subpixel regions
              plotPoint(j, i + size2, size2, size2, 
			calculateRegionProbs(j, i + size2), update);
              plotPoint(j + size2, i + size2, size2, size2, 
			calculateRegionProbs(j + size2, i + size2), update);
              plotPoint(j + size2, i, size2, size2, 
			calculateRegionProbs(j + size2, i), update);
            }
          }
          // The new region edge length is half the old edge length
          m_size = size2;
          size2 = size2 / 2;
        }
	update();
	

	/*
        // Old method without sampling.
        abortPlot: 
        for (int i = 0; i < m_panelHeight; i++) {
          for (int j = 0; j < m_panelWidth; j++) {
            if (m_stopPlotting) {
              break abortPlot;
            }
            plotPoint(j, i, calculateRegionProbs(j, i), (j == 0));
          }
        }
        */


        if (m_plotTrainingData) {
          plotTrainingData();
        }
	      
      } catch (Exception ex) {
        ex.printStackTrace();
	JOptionPane.showMessageDialog(null,"Error while plotting: \"" + ex.getMessage() + "\"");
      } finally {
        m_plotThread = null;
        // notify any listeners that we are finished
        Vector l;
        ActionEvent e = new ActionEvent(this, 0, "");
        synchronized(this) {
          l = (Vector)m_listeners.clone();
        }
        for (int i = 0; i < l.size(); i++) {
          ActionListener al = (ActionListener)l.elementAt(i);
          al.actionPerformed(e);
        }
      }
    }
    
    private double [] calculateRegionProbs(int j, int i) throws Exception {
      double [] sumOfProbsForRegion = 
	new double [m_trainingData.classAttribute().numValues()];

      for (int u = 0; u < m_numOfSamplesPerRegion; u++) {
      
        double [] sumOfProbsForLocation = 
	  new double [m_trainingData.classAttribute().numValues()];
      
        m_weightingAttsValues[m_xAttribute] = getRandomX(j);
        m_weightingAttsValues[m_yAttribute] = getRandomY(m_panelHeight-i-1);
      
        m_dataGenerator.setWeightingValues(m_weightingAttsValues);
      
        double [] weights = m_dataGenerator.getWeights();
        double sumOfWeights = Utils.sum(weights);
        int [] indices = Utils.sort(weights);
      
        // Prune 1% of weight mass
        int [] newIndices = new int[indices.length];
        double sumSoFar = 0; 
        double criticalMass = 0.99 * sumOfWeights;
        int index = weights.length - 1; int counter = 0;
        for (int z = weights.length - 1; z >= 0; z--) {
          newIndices[index--] = indices[z];
          sumSoFar += weights[indices[z]];
          counter++;
          if (sumSoFar > criticalMass) {
            break;
          }
        }
        indices = new int[counter];
        System.arraycopy(newIndices, index + 1, indices, 0, counter);
      
        for (int z = 0; z < m_numOfSamplesPerGenerator; z++) {
        
          m_dataGenerator.setWeightingValues(m_weightingAttsValues);
          double [][] values = m_dataGenerator.generateInstances(indices);
        
          for (int q = 0; q < values.length; q++) {
            if (values[q] != null) {
              System.arraycopy(values[q], 0, m_vals, 0, m_vals.length);
              m_vals[m_xAttribute] = m_weightingAttsValues[m_xAttribute];
              m_vals[m_yAttribute] = m_weightingAttsValues[m_yAttribute];
            
              // classify the instance
              m_dist = m_classifier.distributionForInstance(m_predInst);
              for (int k = 0; k < sumOfProbsForLocation.length; k++) {
                sumOfProbsForLocation[k] += (m_dist[k] * weights[q]); 
              }
            }
          }
        }
      
        for (int k = 0; k < sumOfProbsForRegion.length; k++) {
          sumOfProbsForRegion[k] += (sumOfProbsForLocation[k] * 
				     sumOfWeights); 
        }
      }
    
      // average
      Utils.normalize(sumOfProbsForRegion);

      // cache
      if ((i < m_panelHeight) && (j < m_panelWidth)) {
        m_probabilityCache[i][j] = new double[sumOfProbsForRegion.length];
        System.arraycopy(sumOfProbsForRegion, 0, m_probabilityCache[i][j], 
			 0, sumOfProbsForRegion.length);
      }
		
      return sumOfProbsForRegion;
    }
  }

  /** Render the training points on-screen.
  */
  public void plotTrainingData() {
    
    Graphics2D osg = (Graphics2D)m_osi.getGraphics();
    Graphics g = m_plotPanel.getGraphics();
    osg.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                         RenderingHints.VALUE_ANTIALIAS_ON);
    double xval = 0; double yval = 0;
    
    for (int i = 0; i < m_trainingData.numInstances(); i++) {
      if (!m_trainingData.instance(i).isMissing(m_xAttribute) &&
          !m_trainingData.instance(i).isMissing(m_yAttribute)) {
	  
	if (m_trainingData.instance(i).isMissing(m_classIndex)) //jimmy.
		continue; //don't plot if class is missing. TODO could we plot it differently instead?
	
        xval = m_trainingData.instance(i).value(m_xAttribute);
        yval = m_trainingData.instance(i).value(m_yAttribute);
       
        int panelX = convertToPanelX(xval);
        int panelY = convertToPanelY(yval);
        Color ColorToPlotWith = 
          ((Color)m_Colors.elementAt((int)m_trainingData.instance(i).
                                     value(m_classIndex) % m_Colors.size()));
	
        if (ColorToPlotWith.equals(Color.white)) {
          osg.setColor(Color.black);
        } else {
          osg.setColor(Color.white);
        }
        osg.fillOval(panelX-3, panelY-3, 7, 7);
        osg.setColor(ColorToPlotWith);
        osg.fillOval(panelX-2, panelY-2, 5, 5);
      }
    }
    g.drawImage(m_osi,0,0,m_plotPanel);
  }
  
  /** Convert an X coordinate from the instance space to the panel space.
  */
  private int convertToPanelX(double xval) {
    double temp = (xval - m_minX) / m_rangeX;
    temp = temp * (double) m_panelWidth;

    return (int)temp;
  }

  /** Convert a Y coordinate from the instance space to the panel space.
  */
  private int convertToPanelY(double yval) {
    double temp = (yval - m_minY) / m_rangeY;
    temp = temp * (double) m_panelHeight;
    temp = m_panelHeight - temp;
    
    return (int)temp;
  }
  
  /** Convert an X coordinate from the panel space to the instance space.
  */
  private double convertFromPanelX(double pX) {
    pX /= (double) m_panelWidth;
    pX *= m_rangeX;
    return pX + m_minX;
  }

  /** Convert a Y coordinate from the panel space to the instance space.
  */
  private double convertFromPanelY(double pY) {
    pY  = m_panelHeight - pY;
    pY /= (double) m_panelHeight;
    pY *= m_rangeY;
    
    return pY + m_minY;
  }


  /** Plot a point in our visualization on-screen.
  */
  protected  void plotPoint(int x, int y, double [] probs, boolean update) {
    plotPoint(x, y, 1, 1, probs, update);
  }
  
  /** Plot a point in our visualization on-screen.
  */
  private void plotPoint(int x, int y, int width, int height, 
			 double [] probs, boolean update) {

    // draw a progress line
    Graphics osg = m_osi.getGraphics();
    if (update) {
      osg.setXORMode(Color.white);
      osg.drawLine(0, y, m_panelWidth-1, y);
      update();
      osg.drawLine(0, y, m_panelWidth-1, y);
    }

    // plot the point
    osg.setPaintMode();
    float [] colVal = new float[3];
    
    float [] tempCols = new float[3];
    for (int k = 0; k < probs.length; k++) {
      Color curr = (Color)m_Colors.elementAt(k % m_Colors.size());

      curr.getRGBColorComponents(tempCols);
      for (int z = 0 ; z < 3; z++) {
        colVal[z] += probs[k] * tempCols[z];
      }
    }

    for (int z = 0; z < 3; z++) {
      if (colVal[z] < 0) {
	colVal[z] = 0;
      } else if (colVal[z] > 1) {
	colVal[z] = 1;
      }
    }
    
    osg.setColor(new Color(colVal[0], 
                           colVal[1], 
                           colVal[2]));
    osg.fillRect(x, y, width, height);
  }
  
  /** Update the rendered image.
  */
  private void update() {
    Graphics g = m_plotPanel.getGraphics();
    g.drawImage(m_osi, 0, 0, m_plotPanel);
  }

  /**
   * Set the training data to use
   *
   * @param trainingData the training data
   * @exception Exception if an error occurs
   */
  public void setTrainingData(Instances trainingData) throws Exception {

    m_trainingData = trainingData;
    if (m_trainingData.classIndex() < 0) {
      throw new Exception("No class attribute set (BoundaryPanel)");
    }
    m_classIndex = m_trainingData.classIndex();
  }
  
  /** Adds a training instance to the visualization dataset.
  */
  public void addTrainingInstance(Instance instance) {
  	
  	if (m_trainingData == null) {
		//TODO
		System.err.println("Trying to add to a null training set (BoundaryPanel)");
	}
  	
  	m_trainingData.add(instance);
  }

  /**
   * Register a listener to be notified when plotting completes
   *
   * @param newListener the listener to add
   */
  public void addActionListener(ActionListener newListener) {
    m_listeners.add(newListener);
  }
  
  /**
   * Remove a listener
   *
   * @param removeListener the listener to remove
   */
  public void removeActionListener(ActionListener removeListener) {
    m_listeners.removeElement(removeListener);
  }
  
  /**
   * Set the classifier to use.
   *
   * @param classifier the classifier to use
   */
  public void setClassifier(Classifier classifier) {
    m_classifier = classifier;
  }
  
  /**
   * Set the data generator to use for generating new instances
   *
   * @param dataGenerator the data generator to use
   */
  public void setDataGenerator(DataGenerator dataGenerator) {
    m_dataGenerator = dataGenerator;
  }
  
  /**
   * Set the x attribute index
   *
   * @param xatt index of the attribute to use on the x axis
   * @exception Exception if an error occurs
   */
  public void setXAttribute(int xatt) throws Exception {
    if (m_trainingData == null) {
      throw new Exception("No training data set (BoundaryPanel)");
    }
    if (xatt < 0 || 
        xatt > m_trainingData.numAttributes()) {
      throw new Exception("X attribute out of range (BoundaryPanel)");
    }
    if (m_trainingData.attribute(xatt).isNominal() || m_trainingData.attribute(xatt).isRanking()) {
      throw new Exception("Visualization dimensions must be numeric "
                          +"(BoundaryPanel)");
    }
    /*if (m_trainingData.numDistinctValues(xatt) < 2) {
      throw new Exception("Too few distinct values for X attribute "
                          +"(BoundaryPanel)");
    }*/ //removed by jimmy. TESTING!
    m_xAttribute = xatt;
  }

  /**
   * Set the y attribute index
   *
   * @param yatt index of the attribute to use on the y axis
   * @exception Exception if an error occurs
   */
  public void setYAttribute(int yatt) throws Exception {
    if (m_trainingData == null) {
      throw new Exception("No training data set (BoundaryPanel)");
    }
    if (yatt < 0 || 
        yatt > m_trainingData.numAttributes()) {
      throw new Exception("X attribute out of range (BoundaryPanel)");
    }
    if (m_trainingData.attribute(yatt).isNominal() || m_trainingData.attribute(yatt).isRanking()) {
      throw new Exception("Visualization dimensions must be numeric "
                          +"(BoundaryPanel)");
    }
    /*if (m_trainingData.numDistinctValues(yatt) < 2) {
      throw new Exception("Too few distinct values for Y attribute "
                          +"(BoundaryPanel)");
    }*/ //removed by jimmy. TESTING!
    m_yAttribute = yatt;
  }
  
  /**
   * Set a vector of Color objects for the classes
   *
   * @param colors a <code>FastVector</code> value
   */
  public void setColors(FastVector colors) {
    synchronized (m_Colors) {
      m_Colors = colors;
    }
    //replot(); //commented by jimmy
    update(); //added by jimmy
  }

  /**
   * Set whether to superimpose the training data
   * plot
   *
   * @param pg a <code>boolean</code> value
   */
  public void setPlotTrainingData(boolean pg) {
    m_plotTrainingData = pg;
  }

  /**
   * Returns true if training data is to be superimposed
   *
   * @return a <code>boolean</code> value
   */
  public boolean getPlotTrainingData() {
    return m_plotTrainingData;
  }
  
  /**
   * Get the current vector of Color objects used for the classes
   *
   * @return a <code>FastVector</code> value
   */
  public FastVector getColors() {
    return m_Colors;
  }
  
  /**
   * Quickly replot the display using cached probability estimates
   */
  public void replot() {
    if (m_probabilityCache[0][0] == null) {
      return;
    }
    m_stopReplotting = true;
    m_pausePlotting = true;
    // wait 300 ms to give any other replot threads a chance to halt
    try {
      Thread.sleep(300);
    } catch (Exception ex) {}

    final Thread replotThread = new Thread() {
        public void run() {
          m_stopReplotting = false;
	  int size2 = m_size / 2;
          finishedReplot: for (int i = 0; i < m_panelHeight; i += m_size) {
            for (int j = 0; j < m_panelWidth; j += m_size) {
              if (m_probabilityCache[i][j] == null || m_stopReplotting) {
                break finishedReplot;
              }

	      boolean update = (j == 0 && i % 2 == 0);
	      if (i < m_panelHeight && j < m_panelWidth) {
		// Draw the three new subpixel regions or single course tiling
		if (m_initialTiling || m_size == 1) {
		  if (m_probabilityCache[i][j] == null) {
		    break finishedReplot;
		  }
		  plotPoint(j, i, m_size, m_size, 
			    m_probabilityCache[i][j], update);
		} else {
		  if (m_probabilityCache[i+size2][j] == null) {
		    break finishedReplot;
		  }
		  plotPoint(j, i + size2, size2, size2, 
			    m_probabilityCache[i + size2][j], update);
		  if (m_probabilityCache[i+size2][j+size2] == null) {
		    break finishedReplot;
		  }
		  plotPoint(j + size2, i + size2, size2, size2, 
			    m_probabilityCache[i + size2][j + size2], update);
		  if (m_probabilityCache[i][j+size2] == null) {
		    break finishedReplot;
		  }
		  plotPoint(j + size2, i, size2, size2, 
			    m_probabilityCache[i + size2][j], update);
		}
	      }
	    }
          }
	  update();
          if (m_plotTrainingData) {
            plotTrainingData();
          }
	  m_pausePlotting = false;
	  if (!m_stopPlotting) {
	    synchronized (m_dummy) {
	      m_dummy.notifyAll();
	    }
	  }
        }
      };
    
    replotThread.start();      
  }

  protected void saveImage(String fileName) {
    BufferedImage	bi;
    Graphics2D 		gr2;
    ImageWriter 	writer;
    Iterator 		iter;
    ImageOutputStream 	ios;
    ImageWriteParam 	param;

    try {
      // render image
      bi  = new BufferedImage(m_panelWidth, m_panelHeight, BufferedImage.TYPE_INT_RGB);
      gr2 = bi.createGraphics();
      gr2.drawImage(m_osi, 0, 0, m_panelWidth, m_panelHeight, null);

      // get jpeg writer
      writer = null;
      iter   = ImageIO.getImageWritersByFormatName("jpg");
      if (iter.hasNext())
	writer = (ImageWriter) iter.next();
      else
	throw new Exception("No JPEG writer available!");

      // prepare output file
      ios = ImageIO.createImageOutputStream(new File(fileName));
      writer.setOutput(ios);

      // set the quality
      param = new JPEGImageWriteParam(Locale.getDefault());
      param.setCompressionMode(ImageWriteParam.MODE_EXPLICIT) ;
      param.setCompressionQuality(1.0f);

      // write the image
      writer.write(null, new IIOImage(bi, null, null), param);

      // cleanup
      ios.flush();
      writer.dispose();
      ios.close();    
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  /** Adds a training instance to our dataset, based on the coordinates of the mouse on the panel.
      This method sets the x and y attributes and the class (as defined by classAttIndex), and sets
      all other values as Missing.
   *  @param mouseX the x coordinate of the mouse, in pixels.
   *  @param mouseY the y coordinate of the mouse, in pixels.
   *  @param classAttIndex the index of the attribute that is currently selected as the class attribute.
   *  @param classValue the value to set the class to in our new point.
   */
  public void addTrainingInstanceFromMouseLocation(int mouseX, int mouseY, int classAttIndex, double classValue) {
  	//convert to coordinates in the training instance space.
	double x = convertFromPanelX(mouseX);
	double y = convertFromPanelY(mouseY);
	
	//build the training instance
	Instance newInstance = new DenseInstance(m_trainingData.numAttributes());
	for (int i = 0; i < newInstance.numAttributes(); i++) {
		if (i == classAttIndex) {
			newInstance.setValue(i,classValue);
		}
		else if (i == m_xAttribute)
			newInstance.setValue(i,x);
		else if (i == m_yAttribute)
			newInstance.setValue(i,y);
		else newInstance.setMissing(i);
	}
	
	//add it to our data set.
	addTrainingInstance(newInstance);
  }
  
  /** Deletes all training instances from our dataset.
  */
  public void removeAllInstances() {
  	if (m_trainingData != null)
	{
  		m_trainingData.delete();
		try { initialize();} catch (Exception e) {};
	}
	
  }
  
  /** Removes a single training instance from our dataset, if there is one that is close enough
      to the specified mouse location.
  */
  public void removeTrainingInstanceFromMouseLocation(int mouseX, int mouseY) {
  	
	//convert to coordinates in the training instance space.
	double x = convertFromPanelX(mouseX);
	double y = convertFromPanelY(mouseY);
	
	int bestIndex = -1;
	double bestDistanceBetween = Integer.MAX_VALUE;
	
	//find the closest point.
	for (int i = 0; i < m_trainingData.numInstances(); i++) {
		Instance current = m_trainingData.instance(i);
		double distanceBetween = (current.value(m_xAttribute) - x) * (current.value(m_xAttribute) - x) + (current.value(m_yAttribute) - y) * (current.value(m_yAttribute) - y); // won't bother to sqrt, just used square values.
		
		if (distanceBetween < bestDistanceBetween)
		{
			bestIndex = i;
			bestDistanceBetween = distanceBetween;
		}
	}
	if (bestIndex == -1)
		return;
	Instance best = m_trainingData.instance(bestIndex);
	double panelDistance = (convertToPanelX(best.value(m_xAttribute)) - mouseX) * (convertToPanelX(best.value(m_xAttribute)) - mouseX)
		+ (convertToPanelY(best.value(m_yAttribute)) - mouseY) * (convertToPanelY(best.value(m_yAttribute)) - mouseY);
	if (panelDistance < REMOVE_POINT_RADIUS * REMOVE_POINT_RADIUS) {//the best point is close enough. (using squared distances)
		m_trainingData.delete(bestIndex);
	}
  }
  
  /** Starts the plotting thread.  Will also create it if necessary.
  */
  public void startPlotThread() {
  	if (m_plotThread == null) { //jimmy
      		m_plotThread = new PlotThread();
      		m_plotThread.setPriority(Thread.MIN_PRIORITY);
      		m_plotThread.start();
    	}
  }
  
  /** Adds a mouse listener.
  */
  public void addMouseListener(MouseListener l) {
  	m_plotPanel.addMouseListener(l);
  }
  
  /** Gets the minimum x-coordinate bound, in training-instance units (not mouse coordinates).
  */
  public double getMinXBound() {
  	return m_minX;
  }
  
  /** Gets the minimum y-coordinate bound, in training-instance units (not mouse coordinates).
  */
  public double getMinYBound() {
  	return m_minY;
  }
  
  /** Gets the maximum x-coordinate bound, in training-instance units (not mouse coordinates).
  */
  public double getMaxXBound() {
  	return m_maxX;
  }
  
  /** Gets the maximum x-coordinate bound, in training-instance units (not mouse coordinates).
  */
  public double getMaxYBound() {
  	return m_maxY;
  }

  /**
   * Main method for testing this class
   *
   * @param args a <code>String[]</code> value
   */
  public static void main (String [] args) {
    try {
      if (args.length < 8) {
	System.err.println("Usage : BoundaryPanel <dataset> "
			   +"<class col> <xAtt> <yAtt> "
			   +"<base> <# loc/pixel> <kernel bandwidth> "
			   +"<display width> "
			   +"<display height> <classifier "
			   +"[classifier options]>");
	System.exit(1);
      }
      final javax.swing.JFrame jf = 
	new javax.swing.JFrame("Weka classification boundary visualizer");
      jf.getContentPane().setLayout(new BorderLayout());

      System.err.println("Loading instances from : "+args[0]);
      java.io.Reader r = new java.io.BufferedReader(
			 new java.io.FileReader(args[0]));
      final Instances i = new Instances(r);
      i.setClassIndex(Integer.parseInt(args[1]));

      //      bv.setClassifier(new Logistic());
      final int xatt = Integer.parseInt(args[2]);
      final int yatt = Integer.parseInt(args[3]);
      int base = Integer.parseInt(args[4]);
      int loc = Integer.parseInt(args[5]);

      int bandWidth = Integer.parseInt(args[6]);
      int panelWidth = Integer.parseInt(args[7]);
      int panelHeight = Integer.parseInt(args[8]);

      final String classifierName = args[9];
      final BoundaryPanel bv = new BoundaryPanel(panelWidth,panelHeight);
      bv.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    String classifierNameNew = 
	      classifierName.substring(classifierName.lastIndexOf('.')+1, 
				       classifierName.length());
	    bv.saveImage(classifierNameNew+"_"+i.relationName()
			 +"_X"+xatt+"_Y"+yatt+".jpg");
	  }
	});

      jf.getContentPane().add(bv, BorderLayout.CENTER);
      jf.setSize(bv.getMinimumSize());
      //      jf.setSize(200,200);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	  public void windowClosing(java.awt.event.WindowEvent e) {
	    jf.dispose();
	    System.exit(0);
	  }
	});

      jf.pack();
      jf.setVisible(true);
      //      bv.initialize();
      bv.repaint();
      

      String [] argsR = null;
      if (args.length > 10) {
	argsR = new String [args.length-10];
	for (int j = 10; j < args.length; j++) {
	  argsR[j-10] = args[j];
	}
      }
      Classifier c = AbstractClassifier.forName(args[9], argsR);
      KDDataGenerator dataGen = new KDDataGenerator();
      dataGen.setKernelBandwidth(bandWidth);
      bv.setDataGenerator(dataGen);
      bv.setNumSamplesPerRegion(loc);
      bv.setGeneratorSamplesBase(base);
      bv.setClassifier(c);
      bv.setTrainingData(i);
      bv.setXAttribute(xatt);
      bv.setYAttribute(yatt);

      try {
	// try and load a color map if one exists
	FileInputStream fis = new FileInputStream("colors.ser");
	ObjectInputStream ois = new ObjectInputStream(fis);
	FastVector colors = (FastVector)ois.readObject();
	bv.setColors(colors);	
      } catch (Exception ex) {
	System.err.println("No color map file");
      }
      bv.start();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}

