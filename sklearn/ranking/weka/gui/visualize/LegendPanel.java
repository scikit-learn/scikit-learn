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
 *    LegendPanel.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.visualize;

import weka.core.FastVector;
import weka.core.Instances;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JColorChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

/**
 * This panel displays legends for a list of plots. If a given plot
 * has a custom colour defined then this panel allows the colour to
 * be changed.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 4752 $
 */
public class LegendPanel
  extends JScrollPane {

  /** for serialization */
  private static final long serialVersionUID = -1262384440543001505L;

  /** the list of plot elements */
  protected FastVector m_plots;

  /** the panel that contains the legend entries */
  protected JPanel m_span=null;

  /** a list of components that need to be repainted when a colour is
      changed */
  protected FastVector m_Repainters = new FastVector();

  /**
   * Inner class for handling legend entries
   */
  protected class LegendEntry
    extends JPanel {

    /** for serialization */
    private static final long serialVersionUID = 3879990289042935670L;

    /** the data for this legend entry */
    private PlotData2D m_plotData=null;

    /** the index (in the list of plots) of the data for this legend---
	used to draw the correct shape for this data */
    private int m_dataIndex;

    /** the text part of this legend */
    private JLabel m_legendText;

    /** displays the point shape associated with this legend entry */
    private JPanel m_pointShape;

    public LegendEntry(PlotData2D data, int dataIndex) {
      javax.swing.ToolTipManager.sharedInstance().setDismissDelay(5000);
      m_plotData = data;
      m_dataIndex = dataIndex;
      //      this.setBackground(Color.black);
      /*      this.setPreferredSize(new Dimension(0, 20));
	      this.setMinimumSize(new Dimension(0, 20)); */

      if (m_plotData.m_useCustomColour) {
	this.addMouseListener(new MouseAdapter() {
	    public void mouseClicked(MouseEvent e) {
	      
	      if ((e.getModifiers() & e.BUTTON1_MASK) == e.BUTTON1_MASK) {
		Color tmp = JColorChooser.showDialog
		  (LegendPanel.this, "Select new Color", 
		   m_plotData.m_customColour);
		
		if (tmp != null) {
		  m_plotData.m_customColour = tmp;
		  m_legendText.setForeground(tmp);

		  if (m_Repainters.size() > 0) {
		    for (int i=0;i<m_Repainters.size();i++) {
		      ((Component)(m_Repainters.elementAt(i))).repaint();
		    }
		  }
		  LegendPanel.this.repaint();
		}
	      }
	    }
	  });
      }

      m_legendText = new JLabel(m_plotData.m_plotName);
      m_legendText.setToolTipText(m_plotData.getPlotNameHTML());
      if (m_plotData.m_useCustomColour) {
	m_legendText.setForeground(m_plotData.m_customColour);
      }
      this.setLayout(new BorderLayout());
      this.add(m_legendText, BorderLayout.CENTER);
      /*      GridBagLayout gb = new GridBagLayout();
      GridBagConstraints constraints = new GridBagConstraints();
      constraints.fill = GridBagConstraints.HORIZONTAL;
      constraints.gridx=0;constraints.gridy=0;constraints.weightx=5; */
      m_pointShape = new JPanel() {
	private static final long serialVersionUID = -7048435221580488238L;
	
	public void paintComponent(Graphics gx) {
	  super.paintComponent(gx);
	  if (!m_plotData.m_useCustomColour) {
	    gx.setColor(Color.black);
	  } else {
	    gx.setColor(m_plotData.m_customColour);
	  }
	  Plot2D.drawDataPoint(10,10,3,m_dataIndex,gx);
	}
      };
      //      m_pointShape.setBackground(Color.black);
      m_pointShape.setPreferredSize(new Dimension(20, 20));
      m_pointShape.setMinimumSize(new Dimension(20, 20));
      this.add(m_pointShape, BorderLayout.WEST);
    }
  }

  /**
   * Constructor
   */
  public LegendPanel() {
    this.setBackground(Color.blue);
    setVerticalScrollBarPolicy(VERTICAL_SCROLLBAR_ALWAYS);
  }

  /**
   * Set the list of plots to generate legend entries for
   * @param pl a list of plots
   */
  public void setPlotList(FastVector pl) {
    m_plots = pl;
    updateLegends();
  }

  /**
   * Adds a component that will need to be repainted if the user
   * changes the colour of a label.
   * @param c the component to be repainted
   */
  public void addRepaintNotify(Component c) {
    m_Repainters.addElement(c);
  }

  /**
   * Redraw the panel with the legend entries
   */
  private void updateLegends() {
    if (m_span == null) {
      m_span = new JPanel();
    }
      
    JPanel padder = new JPanel();
    JPanel padd2 = new JPanel();

    m_span.setPreferredSize(new Dimension(m_span.getPreferredSize().width, 
					  (m_plots.size() + 1) * 20));
    m_span.setMaximumSize(new Dimension(m_span.getPreferredSize().width, 
					  (m_plots.size() + 1) * 20));

    LegendEntry tmp;

    GridBagLayout gb = new GridBagLayout();
    GridBagLayout gb2 = new GridBagLayout();
    GridBagConstraints constraints = new GridBagConstraints();
      
    m_span.removeAll();

    padder.setLayout(gb);
    m_span.setLayout(gb2);
    constraints.anchor = GridBagConstraints.CENTER;
    constraints.gridx=0;constraints.gridy=0;constraints.weightx=5;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridwidth=1;constraints.gridheight=1;
    constraints.insets = new Insets(0, 0, 0, 0);
    padder.add(m_span, constraints);

    constraints.gridx=0;constraints.gridy=1;constraints.weightx=5;
    constraints.fill = GridBagConstraints.BOTH;
    constraints.gridwidth=1;constraints.gridheight=1;constraints.weighty=5;
    constraints.insets = new Insets(0, 0, 0, 0);
    padder.add(padd2, constraints);

    constraints.weighty=0;
    setViewportView(padder);

    constraints.anchor = GridBagConstraints.CENTER;
    constraints.gridx=0;constraints.gridy=0;constraints.weightx=5;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridwidth=1;constraints.gridheight=1;constraints.weighty=5;
    constraints.insets = new Insets(2,4,2,4);
    //int numLines = ((PlotData2D)m_plots.elementAt(0)).getPlotName().split("<br>").length;
     for (int i=0;i<m_plots.size();i++) {
       tmp = new LegendEntry((PlotData2D)m_plots.elementAt(i),i);
       constraints.gridy = i;
/*       constraints.gridheight = 1;
       if (numLines > 0) {
         constraints.gridheight = (numLines + 2);
       } */
       m_span.add(tmp, constraints);
     }
  }

  /**
   * Main method for testing this class
   * @param args a list of arff files
   */
  public static void main(String [] args) {
    try {
      if (args.length < 1) {
	System.err.println("Usage : weka.gui.visualize.LegendPanel "
			   +"<dataset> [dataset2], [dataset3],...");
	System.exit(1);
      }

      final javax.swing.JFrame jf = 
	new javax.swing.JFrame("Weka Explorer: Legend");
      jf.setSize(100,100);
      jf.getContentPane().setLayout(new BorderLayout());
      final LegendPanel p2 = new LegendPanel();
      jf.getContentPane().add(p2, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	public void windowClosing(java.awt.event.WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });

      FastVector plotList = new FastVector();
      for (int j=0;j<args.length;j++) {
	System.err.println("Loading instances from " + args[j]);
	java.io.Reader r = new java.io.BufferedReader(
			   new java.io.FileReader(args[j]));
	Instances i = new Instances(r);
	PlotData2D tmp = new PlotData2D(i);
	if (j != 1) {
	  tmp.m_useCustomColour = true;
	  tmp.m_customColour = Color.red;
	}
	tmp.setPlotName(i.relationName());
	plotList.addElement(tmp);
      }
      
      p2.setPlotList(plotList);
      jf.setVisible(true);
    } catch (Exception ex) {
      System.err.println(ex.getMessage());
      ex.printStackTrace();
    }
  }
}
