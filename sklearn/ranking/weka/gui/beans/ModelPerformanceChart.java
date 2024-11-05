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
 *    ModelPerformanceChart.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instances;
import weka.gui.visualize.PlotData2D;
import weka.gui.visualize.VisualizePanel;

import java.awt.BorderLayout;
import java.beans.PropertyChangeListener;
import java.beans.VetoableChangeListener;
import java.beans.beancontext.BeanContext;
import java.beans.beancontext.BeanContextChild;
import java.beans.beancontext.BeanContextChildSupport;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * Bean that can be used for displaying threshold curves (e.g. ROC
 * curves) and scheme error plots
 *
 * @author Mark Hall
 * @version $Revision: 4762 $
 */
public class ModelPerformanceChart
  extends JPanel
  implements ThresholdDataListener, VisualizableErrorListener, 
             Visible, UserRequestAcceptor,
	     Serializable, BeanContextChild {

  /** for serialization */
  private static final long serialVersionUID = -4602034200071195924L;
  
  protected BeanVisual m_visual;

  protected transient PlotData2D m_masterPlot;
  
  protected transient JFrame m_popupFrame;

  protected boolean m_framePoppedUp = false;

  /**
   * True if this bean's appearance is the design mode appearance
   */
  protected boolean m_design;

  /**
   * BeanContex that this bean might be contained within
   */
  protected transient BeanContext m_beanContext = null;

  private transient VisualizePanel m_visPanel;
  
  /**
   * BeanContextChild support
   */
  protected BeanContextChildSupport m_bcSupport = 
    new BeanContextChildSupport(this);

  public ModelPerformanceChart() {
    java.awt.GraphicsEnvironment ge = 
      java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
    if (!ge.isHeadless()) {
      appearanceFinal();
    }
  }

    /**
   * Global info for this bean
   *
   * @return a <code>String</code> value
   */
  public String globalInfo() {
    return "Visualize performance charts (such as ROC).";
  }

  protected void appearanceDesign() {
    removeAll();
    m_visual = new BeanVisual("ModelPerformanceChart", 
			      BeanVisual.ICON_PATH+"ModelPerformanceChart.gif",
			      BeanVisual.ICON_PATH
			      +"ModelPerformanceChart_animated.gif");
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
  }

  protected void appearanceFinal() {
    removeAll();
    setLayout(new BorderLayout());
    setUpFinal();
  }

  protected void setUpFinal() {
    if (m_visPanel == null) {
      m_visPanel = new VisualizePanel();
    }
    add(m_visPanel, BorderLayout.CENTER);
  }

  /**
   * Display a threshold curve.
   *
   * @param e a ThresholdDataEvent
   */
  public synchronized void acceptDataSet(ThresholdDataEvent e) {
    if (m_visPanel == null) {
      m_visPanel = new VisualizePanel();
    }
    if (m_masterPlot == null) {
      m_masterPlot = e.getDataSet();
    }
    try {
    // check for compatable data sets
      if (!m_masterPlot.getPlotInstances().relationName().
	  equals(e.getDataSet().getPlotInstances().relationName())) {
	
	// if not equal then remove all plots and set as new master plot
	m_masterPlot = e.getDataSet();
	m_visPanel.setMasterPlot(m_masterPlot);
	m_visPanel.validate(); m_visPanel.repaint();
      } else {
	// add as new plot
	m_visPanel.addPlot(e.getDataSet());
	m_visPanel.validate(); m_visPanel.repaint();
      }
      m_visPanel.setXIndex(4); m_visPanel.setYIndex(5);
    } catch (Exception ex) {
      System.err.println("Problem setting up visualization (ModelPerformanceChart)");
      ex.printStackTrace();
    }
  }

  /**
   * Display a scheme error plot.
   *
   * @param e a VisualizableErrorEvent
   */
  public synchronized void acceptDataSet(VisualizableErrorEvent e) {
    if (m_visPanel == null) {
      m_visPanel = new VisualizePanel();
    }
    if (m_masterPlot == null) {
      m_masterPlot = e.getDataSet();
    }
    try {
      m_visPanel.setMasterPlot(m_masterPlot);
    } catch (Exception ex) {
      System.err.println("Problem setting up visualization (ModelPerformanceChart)");
      ex.printStackTrace();
    }
    m_visPanel.validate();
    m_visPanel.repaint();
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
   * Return the visual appearance of this bean
   */
  public BeanVisual getVisual() {
    return m_visual;
  }

  /**
   * Use the default appearance for this bean
   */
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"DefaultDataVisualizer.gif",
		       BeanVisual.ICON_PATH+"DefaultDataVisualizer_animated.gif");
  }

  /**
   * Describe <code>enumerateRequests</code> method here.
   *
   * @return an <code>Enumeration</code> value
   */
  public Enumeration enumerateRequests() {
    Vector newVector = new Vector(0);
    if (m_masterPlot != null) {
      newVector.addElement("Show chart");
      newVector.addElement("?Clear all plots");
    }
    return newVector.elements();
  }

  /**
   * Add a property change listener to this bean
   *
   * @param name the name of the property of interest
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void addPropertyChangeListener(String name,
					PropertyChangeListener pcl) {
    m_bcSupport.addPropertyChangeListener(name, pcl);
  }

  /**
   * Remove a property change listener from this bean
   *
   * @param name the name of the property of interest
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void removePropertyChangeListener(String name,
					   PropertyChangeListener pcl) {
    m_bcSupport.removePropertyChangeListener(name, pcl);
  }

  /**
   * Add a vetoable change listener to this bean
   *
   * @param name the name of the property of interest
   * @param vcl a <code>VetoableChangeListener</code> value
   */
  public void addVetoableChangeListener(String name,
				       VetoableChangeListener vcl) {
    m_bcSupport.addVetoableChangeListener(name, vcl);
  }
  
  /**
   * Remove a vetoable change listener from this bean
   *
   * @param name the name of the property of interest
   * @param vcl a <code>VetoableChangeListener</code> value
   */
  public void removeVetoableChangeListener(String name,
					   VetoableChangeListener vcl) {
    m_bcSupport.removeVetoableChangeListener(name, vcl);
  }

  /**
   * Set a bean context for this bean
   *
   * @param bc a <code>BeanContext</code> value
   */
  public void setBeanContext(BeanContext bc) {
    m_beanContext = bc;
    m_design = m_beanContext.isDesignTime();
    if (m_design) {
      appearanceDesign();
    } else {
      java.awt.GraphicsEnvironment ge = 
        java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment(); 
      if (!ge.isHeadless()) {
        appearanceFinal();
      }
    }
  }

  /**
   * Return the bean context (if any) that this bean is embedded in
   *
   * @return a <code>BeanContext</code> value
   */
  public BeanContext getBeanContext() {
    return m_beanContext;
  }

  /**
   * Describe <code>performRequest</code> method here.
   *
   * @param request a <code>String</code> value
   * @exception IllegalArgumentException if an error occurs
   */
  public void performRequest(String request) {
    if (request.compareTo("Show chart") == 0) {
      try {
	// popup visualize panel
	if (!m_framePoppedUp) {
	  m_framePoppedUp = true;

	  final javax.swing.JFrame jf = 
	    new javax.swing.JFrame("Model Performance Chart");
	  jf.setSize(800,600);
	  jf.getContentPane().setLayout(new BorderLayout());
	  jf.getContentPane().add(m_visPanel, BorderLayout.CENTER);
	  jf.addWindowListener(new java.awt.event.WindowAdapter() {
	      public void windowClosing(java.awt.event.WindowEvent e) {
		jf.dispose();
		m_framePoppedUp = false;
	      }
	    });
	  jf.setVisible(true);
	  m_popupFrame = jf;
	} else {
	  m_popupFrame.toFront();
	}
      } catch (Exception ex) {
	ex.printStackTrace();
	m_framePoppedUp = false;
      }
    } else if (request.equals("Clear all plots")) {
        m_visPanel.removeAllPlots();
        m_visPanel.validate(); m_visPanel.repaint();
        m_visPanel = null;
        m_masterPlot = null;
    } else {
      throw new IllegalArgumentException(request
					 + " not supported (Model Performance Chart)");
    }
  }
  
  public static void main(String [] args) {
    try {
      if (args.length != 1) {
	System.err.println("Usage: ModelPerformanceChart <dataset>");
	System.exit(1);
      }
      java.io.Reader r = new java.io.BufferedReader(
			 new java.io.FileReader(args[0]));
      Instances inst = new Instances(r);
      final javax.swing.JFrame jf = new javax.swing.JFrame();
      jf.getContentPane().setLayout(new java.awt.BorderLayout());
      final ModelPerformanceChart as = new ModelPerformanceChart();
      PlotData2D pd = new PlotData2D(inst);
      pd.setPlotName(inst.relationName());
      ThresholdDataEvent roc = new ThresholdDataEvent(as, pd);
      as.acceptDataSet(roc);      

      jf.getContentPane().add(as, java.awt.BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
        public void windowClosing(java.awt.event.WindowEvent e) {
          jf.dispose();
          System.exit(0);
        }
      });
      jf.setSize(800,600);
      jf.setVisible(true);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
