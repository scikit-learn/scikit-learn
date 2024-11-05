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
 *    CostBenefitAnalysis.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.beans.EventSetDescriptor;
import java.beans.PropertyVetoException;
import java.beans.VetoableChangeListener;
import java.beans.beancontext.BeanContext;
import java.beans.beancontext.BeanContextChild;
import java.beans.beancontext.BeanContextChildSupport;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import weka.classifiers.evaluation.ThresholdCurve;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Utils;
import weka.gui.Logger;
import weka.gui.visualize.VisualizePanel;
import weka.gui.visualize.Plot2D;
import weka.gui.visualize.PlotData2D;


/**
 * Bean that aids in analyzing cost/benefit tradeoffs.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6137 $
 */
public class CostBenefitAnalysis extends JPanel 
  implements BeanCommon, ThresholdDataListener, Visible, UserRequestAcceptor,
  Serializable, BeanContextChild {
  
  /** For serialization */
  private static final long serialVersionUID = 8647471654613320469L;

  protected BeanVisual m_visual;
  
  protected transient JFrame m_popupFrame;

  protected boolean m_framePoppedUp = false;
  
  private transient AnalysisPanel m_analysisPanel;
  
  /**
   * True if this bean's appearance is the design mode appearance
   */
  protected boolean m_design;

  /**
   * BeanContex that this bean might be contained within
   */
  protected transient BeanContext m_beanContext = null;
  
  /**
   * BeanContextChild support
   */
  protected BeanContextChildSupport m_bcSupport = 
    new BeanContextChildSupport(this);
  
  /**
   * The object sending us data (we allow only one connection at any one time)
   */
  protected Object m_listenee;
  
  /**
   * Inner class for displaying the plots and all control widgets.
   * 
   * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
   */
  protected static class AnalysisPanel extends JPanel {
    
    /** For serialization */
    private static final long serialVersionUID = 5364871945448769003L;

    /** Displays the performance graphs(s) */
    protected VisualizePanel m_performancePanel = new VisualizePanel();
    
    /** Displays the cost/benefit (profit/loss) graph */
    protected VisualizePanel m_costBenefitPanel = new VisualizePanel();
    
    /** 
     * The class attribute from the data that was used to generate
     * the threshold curve
     */
    protected Attribute m_classAttribute;
    
    /** Data for the threshold curve */
    protected PlotData2D m_masterPlot;
    
    /** Data for the cost/benefit curve */
    protected PlotData2D m_costBenefit;
    
    /** The size of the points being plotted */
    protected int[] m_shapeSizes;
    
    /** The index of the previous plotted point that was highlighted */
    protected int m_previousShapeIndex = -1;
        
    /** The slider for adjusting the threshold */
    protected JSlider m_thresholdSlider = new JSlider(0,100,0);
    
    protected JRadioButton m_percPop = new JRadioButton("% of Population");
    protected JRadioButton m_percOfTarget = new JRadioButton("% of Target (recall)");
    protected JRadioButton m_threshold = new JRadioButton("Score Threshold");
    
    protected JLabel m_percPopLab = new JLabel();
    protected JLabel m_percOfTargetLab = new JLabel();
    protected JLabel m_thresholdLab = new JLabel();
    
    // Confusion matrix stuff
    protected JLabel m_conf_predictedA = new JLabel("Predicted (a)", SwingConstants.RIGHT);
    protected JLabel m_conf_predictedB = new JLabel("Predicted (b)", SwingConstants.RIGHT);
    protected JLabel m_conf_actualA = new JLabel(" Actual (a):");
    protected JLabel m_conf_actualB = new JLabel(" Actual (b):");
    protected ConfusionCell m_conf_aa = new ConfusionCell();
    protected ConfusionCell m_conf_ab = new ConfusionCell();
    protected ConfusionCell m_conf_ba = new ConfusionCell();
    protected ConfusionCell m_conf_bb = new ConfusionCell();
    
    // Cost matrix stuff
    protected JLabel m_cost_predictedA = new JLabel("Predicted (a)", SwingConstants.RIGHT);
    protected JLabel m_cost_predictedB = new JLabel("Predicted (b)", SwingConstants.RIGHT);
    protected JLabel m_cost_actualA = new JLabel(" Actual (a)");
    protected JLabel m_cost_actualB = new JLabel(" Actual (b)");
    protected JTextField m_cost_aa = new JTextField("0.0", 5);
    protected JTextField m_cost_ab = new JTextField("1.0", 5);
    protected JTextField m_cost_ba = new JTextField("1.0", 5);
    protected JTextField m_cost_bb = new JTextField("0.0" ,5);
    protected JButton m_maximizeCB = new JButton("Maximize Cost/Benefit");
    protected JButton m_minimizeCB = new JButton("Minimize Cost/Benefit");
    protected JRadioButton m_costR = new JRadioButton("Cost");
    protected JRadioButton m_benefitR = new JRadioButton("Benefit");
    protected JLabel m_costBenefitL = new JLabel("Cost: ", SwingConstants.RIGHT);
    protected JLabel m_costBenefitV = new JLabel("0");
    protected JLabel m_randomV = new JLabel("0");
    protected JLabel m_gainV = new JLabel("0");
    
    protected int m_originalPopSize;
    
    /** Population text field */
    protected JTextField m_totalPopField = new JTextField(6);
    protected int m_totalPopPrevious;
    
    /** Classification accuracy */
    protected JLabel m_classificationAccV = new JLabel("-");
    
    // Only update curve & stats if values in cost matrix have changed
    protected double m_tpPrevious;
    protected double m_fpPrevious;
    protected double m_tnPrevious;
    protected double m_fnPrevious;
    
    /**
     * Inner class for handling a single cell in the confusion matrix.
     * Displays the value, value as a percentage of total population and
     * graphical depiction of percentage.
     * 
     * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
     */
    protected static class ConfusionCell extends JPanel {

      /** For serialization */
      private static final long serialVersionUID = 6148640235434494767L;
      
      private JLabel m_conf_cell = new JLabel("-", SwingConstants.RIGHT);
      JLabel m_conf_perc = new JLabel("-", SwingConstants.RIGHT);
      
      private JPanel m_percentageP;
      
      protected double m_percentage = 0;
      
      public ConfusionCell() {
        setLayout(new BorderLayout());
        setBorder(BorderFactory.createEtchedBorder());
        
        add(m_conf_cell, BorderLayout.NORTH);
        
        m_percentageP = new JPanel() {
          public void paintComponent(Graphics gx) {
            super.paintComponent(gx);
            
            if (m_percentage > 0) {
              gx.setColor(Color.BLUE);
              int height = this.getHeight();
              double width = this.getWidth();
              int barWidth = (int)(m_percentage * width); 
              gx.fillRect(0, 0, barWidth, height);
            }
          }
        };
        
        Dimension d = new Dimension(30,5);
        m_percentageP.setMinimumSize(d);
        m_percentageP.setPreferredSize(d);
        JPanel percHolder = new JPanel();
        percHolder.setLayout(new BorderLayout());
        percHolder.add(m_percentageP, BorderLayout.CENTER);
        percHolder.add(m_conf_perc, BorderLayout.EAST);
        
        add(percHolder, BorderLayout.SOUTH);
      }
      
      /**
       * Set the value of a cell.
       * 
       * @param cellValue the value of the cell
       * @param max the max (for setting value as a percentage)
       * @param scaleFactor scale the value by this amount
       * @param precision precision for the percentage value
       */
      public void setCellValue(double cellValue, double max, double scaleFactor, int precision) {
        if (!Utils.isMissingValue(cellValue)) {
          m_percentage = cellValue / max;
        } else {
          m_percentage = 0;
        }
        
        m_conf_cell.setText(Utils.doubleToString((cellValue * scaleFactor), 0));
        m_conf_perc.setText(Utils.doubleToString(m_percentage * 100.0, precision) + "%");
        
        // refresh the percentage bar
        m_percentageP.repaint();
      }
    }
    
    public AnalysisPanel() {
      setLayout(new BorderLayout());
      m_performancePanel.setShowAttBars(false);
      m_performancePanel.setShowClassPanel(false);
      m_costBenefitPanel.setShowAttBars(false);
      m_costBenefitPanel.setShowClassPanel(false);
      
      Dimension size = new Dimension(500, 400);
      m_performancePanel.setPreferredSize(size);
      m_performancePanel.setMinimumSize(size);
      
      size = new Dimension(500, 400);
      m_costBenefitPanel.setMinimumSize(size);
      m_costBenefitPanel.setPreferredSize(size);
      
      m_thresholdSlider.addChangeListener(new ChangeListener() {
        public void stateChanged(ChangeEvent e) {
          updateInfoForSliderValue((double)m_thresholdSlider.getValue() / 100.0);
        }
      });
      
      JPanel plotHolder = new JPanel();
      plotHolder.setLayout(new GridLayout(1,2));      
      plotHolder.add(m_performancePanel);
      plotHolder.add(m_costBenefitPanel);
      add(plotHolder, BorderLayout.CENTER);
      
      JPanel lowerPanel = new JPanel();
      lowerPanel.setLayout(new BorderLayout());
      
      ButtonGroup bGroup = new ButtonGroup();
      bGroup.add(m_percPop);
      bGroup.add(m_percOfTarget);
      bGroup.add(m_threshold);
      
      ButtonGroup bGroup2 = new ButtonGroup();
      bGroup2.add(m_costR);
      bGroup2.add(m_benefitR);
      ActionListener rl = new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          if (m_costR.isSelected()) {
            m_costBenefitL.setText("Cost: ");
          } else {
            m_costBenefitL.setText("Benefit: ");
          }

          double gain = Double.parseDouble(m_gainV.getText());
          gain = -gain;
          m_gainV.setText(Utils.doubleToString(gain, 2));
        }
      };
      m_costR.addActionListener(rl);
      m_benefitR.addActionListener(rl);
      m_costR.setSelected(true);
      
      m_percPop.setSelected(true);
      JPanel threshPanel = new JPanel();
      threshPanel.setLayout(new BorderLayout());
      JPanel radioHolder = new JPanel();
      radioHolder.setLayout(new FlowLayout());
      radioHolder.add(m_percPop);
      radioHolder.add(m_percOfTarget);
      radioHolder.add(m_threshold);
      threshPanel.add(radioHolder, BorderLayout.NORTH);
      threshPanel.add(m_thresholdSlider, BorderLayout.SOUTH);
      
      JPanel threshInfoPanel = new JPanel();
      threshInfoPanel.setLayout(new GridLayout(3,2));
      threshInfoPanel.add(new JLabel("% of Population: ", SwingConstants.RIGHT));
      threshInfoPanel.add(m_percPopLab);
      threshInfoPanel.add(new JLabel("% of Target: ", SwingConstants.RIGHT));
      threshInfoPanel.add(m_percOfTargetLab);
      threshInfoPanel.add(new JLabel("Score Threshold: ", SwingConstants.RIGHT));
      threshInfoPanel.add(m_thresholdLab);
      
      JPanel threshHolder = new JPanel();
      threshHolder.setBorder(BorderFactory.createTitledBorder("Threshold"));
      threshHolder.setLayout(new BorderLayout());
      threshHolder.add(threshPanel, BorderLayout.CENTER);
      threshHolder.add(threshInfoPanel, BorderLayout.EAST);
      
      lowerPanel.add(threshHolder, BorderLayout.NORTH);
      
      // holder for the two matrixes
      JPanel matrixHolder = new JPanel();
      matrixHolder.setLayout(new GridLayout(1,2));
      
      // confusion matrix
      JPanel confusionPanel = new JPanel();
      confusionPanel.setLayout(new GridLayout(3,3));
      confusionPanel.add(m_conf_predictedA);
      confusionPanel.add(m_conf_predictedB);
      confusionPanel.add(new JLabel()); // dummy
      confusionPanel.add(m_conf_aa);
      confusionPanel.add(m_conf_ab);
      confusionPanel.add(m_conf_actualA);
      confusionPanel.add(m_conf_ba);
      confusionPanel.add(m_conf_bb);
      confusionPanel.add(m_conf_actualB);
      JPanel tempHolderCA = new JPanel();
      tempHolderCA.setLayout(new BorderLayout());
      tempHolderCA.setBorder(BorderFactory.createTitledBorder("Confusion Matrix"));
      tempHolderCA.add(confusionPanel, BorderLayout.CENTER);
      
      JPanel accHolder = new JPanel();
      accHolder.setLayout(new FlowLayout(FlowLayout.LEFT));
      accHolder.add(new JLabel("Classification Accuracy: "));
      accHolder.add(m_classificationAccV);
      tempHolderCA.add(accHolder, BorderLayout.SOUTH);
      
      matrixHolder.add(tempHolderCA);
      
      // cost matrix
      JPanel costPanel = new JPanel();
      costPanel.setBorder(BorderFactory.createTitledBorder("Cost Matrix"));
      costPanel.setLayout(new BorderLayout());
      
      JPanel cmHolder = new JPanel();
      cmHolder.setLayout(new GridLayout(3, 3));
      cmHolder.add(m_cost_predictedA);      
      cmHolder.add(m_cost_predictedB);
      cmHolder.add(new JLabel()); // dummy
      cmHolder.add(m_cost_aa);
      cmHolder.add(m_cost_ab);
      cmHolder.add(m_cost_actualA);
      cmHolder.add(m_cost_ba);
      cmHolder.add(m_cost_bb);
      cmHolder.add(m_cost_actualB);
      costPanel.add(cmHolder, BorderLayout.CENTER);
      
      FocusListener fl = new FocusListener() {
        public void focusGained(FocusEvent e) {
          
        }
        
        public void focusLost(FocusEvent e) {
          if (constructCostBenefitData()) {
            try {
              m_costBenefitPanel.setMasterPlot(m_costBenefit);
              m_costBenefitPanel.validate(); m_costBenefitPanel.repaint();
            } catch (Exception ex) {
              ex.printStackTrace();
            }
            updateCostBenefit();
          }
        }
      };
      
      ActionListener al = new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          if (constructCostBenefitData()) {
            try {
              m_costBenefitPanel.setMasterPlot(m_costBenefit);
              m_costBenefitPanel.validate(); m_costBenefitPanel.repaint();
            } catch (Exception ex) {
              ex.printStackTrace();
            }
            updateCostBenefit();
          }
        }
      };
            
      m_cost_aa.addFocusListener(fl);
      m_cost_aa.addActionListener(al);
      m_cost_ab.addFocusListener(fl);
      m_cost_ab.addActionListener(al);
      m_cost_ba.addFocusListener(fl);
      m_cost_ba.addActionListener(al);
      m_cost_bb.addFocusListener(fl);
      m_cost_bb.addActionListener(al);
      
      m_totalPopField.addFocusListener(fl);
      m_totalPopField.addActionListener(al);
      
      JPanel cbHolder = new JPanel();
      cbHolder.setLayout(new BorderLayout());
      JPanel tempP = new JPanel();
      tempP.setLayout(new GridLayout(3, 2));
      tempP.add(m_costBenefitL);
      tempP.add(m_costBenefitV);
      tempP.add(new JLabel("Random: ", SwingConstants.RIGHT));
      tempP.add(m_randomV);
      tempP.add(new JLabel("Gain: ", SwingConstants.RIGHT));
      tempP.add(m_gainV);
      cbHolder.add(tempP, BorderLayout.NORTH);
      JPanel butHolder = new JPanel();
      butHolder.setLayout(new GridLayout(2, 1));
      butHolder.add(m_maximizeCB);
      butHolder.add(m_minimizeCB);
      m_maximizeCB.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          findMaxMinCB(true);
        }
      });
      
      m_minimizeCB.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          findMaxMinCB(false);
        }
      });
      
      cbHolder.add(butHolder, BorderLayout.SOUTH);
      costPanel.add(cbHolder, BorderLayout.EAST);
      
      JPanel popCBR = new JPanel();
      popCBR.setLayout(new GridLayout(1, 2));
      JPanel popHolder = new JPanel();
      popHolder.setLayout(new FlowLayout(FlowLayout.LEFT));
      popHolder.add(new JLabel("Total Population: "));
      popHolder.add(m_totalPopField);
      
      JPanel radioHolder2 = new JPanel();
      radioHolder2.setLayout(new FlowLayout(FlowLayout.RIGHT));
      radioHolder2.add(m_costR);
      radioHolder2.add(m_benefitR);
      popCBR.add(popHolder);
      popCBR.add(radioHolder2);
      
      costPanel.add(popCBR, BorderLayout.SOUTH);
      
      matrixHolder.add(costPanel);
      
      
      lowerPanel.add(matrixHolder, BorderLayout.SOUTH);
      


//      popAccHolder.add(popHolder);
      
      //popAccHolder.add(accHolder);
      
      /*JPanel lowerPanel2 = new JPanel();
      lowerPanel2.setLayout(new BorderLayout());
      lowerPanel2.add(lowerPanel, BorderLayout.NORTH);
      lowerPanel2.add(popAccHolder, BorderLayout.SOUTH); */
      
      add(lowerPanel, BorderLayout.SOUTH);
      
    }
    
    private void findMaxMinCB(boolean max) {
      double maxMin = (max) 
      ? Double.NEGATIVE_INFINITY 
          : Double.POSITIVE_INFINITY;
      
      Instances cBCurve = m_costBenefit.getPlotInstances();
      int maxMinIndex = 0;
      
      for (int i = 0; i < cBCurve.numInstances(); i++) {
        Instance current = cBCurve.instance(i);
        if (max) {
          if (current.value(1) > maxMin) {
            maxMin = current.value(1);
            maxMinIndex = i;
          }
        } else {
          if (current.value(1) < maxMin) {
            maxMin = current.value(1);
            maxMinIndex = i;
          }
        }
      }
      
      
      // set the slider to the correct position
      int indexOfSampleSize = 
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.SAMPLE_SIZE_NAME).index();
      int indexOfPercOfTarget = 
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.RECALL_NAME).index();
      int indexOfThreshold =
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.THRESHOLD_NAME).index();
      int indexOfMetric;
      
      if (m_percPop.isSelected()) {
        indexOfMetric = indexOfSampleSize;           
      } else if (m_percOfTarget.isSelected()) {
        indexOfMetric = indexOfPercOfTarget;
      } else {
        indexOfMetric = indexOfThreshold;
      }
      
      double valueOfMetric = m_masterPlot.getPlotInstances().instance(maxMinIndex).value(indexOfMetric);
      valueOfMetric *= 100.0;
      
      // set the approximate location of the slider
      m_thresholdSlider.setValue((int)valueOfMetric);
      
      // make sure the actual values relate to the true min/max rather
      // than being off due to slider location error.
      updateInfoGivenIndex(maxMinIndex);
    }
    
    private void updateCostBenefit() {
      double value = (double)m_thresholdSlider.getValue() / 100.0;
      Instances plotInstances = m_masterPlot.getPlotInstances();
      int indexOfSampleSize = 
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.SAMPLE_SIZE_NAME).index();
      int indexOfPercOfTarget = 
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.RECALL_NAME).index();
      int indexOfThreshold =
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.THRESHOLD_NAME).index();
      int indexOfMetric;
      
      if (m_percPop.isSelected()) {
        indexOfMetric = indexOfSampleSize;           
      } else if (m_percOfTarget.isSelected()) {
        indexOfMetric = indexOfPercOfTarget;
      } else {
        indexOfMetric = indexOfThreshold;
      }
      
      int index = findIndexForValue(value, plotInstances, indexOfMetric);
      updateCBRandomGainInfo(index);
    }
    
    private void updateCBRandomGainInfo(int index) {
      double requestedPopSize = m_originalPopSize;
      try {
        requestedPopSize = Double.parseDouble(m_totalPopField.getText());
      } catch (NumberFormatException e) {}
      double scaleFactor = requestedPopSize / m_originalPopSize;
      
      double CB = m_costBenefit.
        getPlotInstances().instance(index).value(1);
      m_costBenefitV.setText(Utils.doubleToString(CB,2));
      
      double totalRandomCB = 0.0;
      Instance first = m_masterPlot.getPlotInstances().instance(0);
      double totalPos = first.value(m_masterPlot.getPlotInstances().
          attribute(ThresholdCurve.TRUE_POS_NAME).index()) * scaleFactor;
      double totalNeg = first.value(m_masterPlot.getPlotInstances().
          attribute(ThresholdCurve.FALSE_POS_NAME)) * scaleFactor;

      double posInSample = (totalPos * (Double.parseDouble(m_percPopLab.getText()) / 100.0));
      double negInSample = (totalNeg * (Double.parseDouble(m_percPopLab.getText()) / 100.0));
      double posOutSample = totalPos - posInSample;
      double negOutSample = totalNeg - negInSample;
      
      double tpCost = 0.0;
      try {
        tpCost = Double.parseDouble(m_cost_aa.getText());
      } catch (NumberFormatException n) {}
      double fpCost = 0.0;
      try {
        fpCost = Double.parseDouble(m_cost_ba.getText());
      } catch (NumberFormatException n) {}
      double tnCost = 0.0;
      try {
        tnCost = Double.parseDouble(m_cost_bb.getText());
      } catch (NumberFormatException n) {}
      double fnCost = 0.0;
      try {
        fnCost = Double.parseDouble(m_cost_ab.getText());
      } catch (NumberFormatException n) {}
            
      totalRandomCB += posInSample * tpCost;
      totalRandomCB += negInSample * fpCost;
      totalRandomCB += posOutSample * fnCost;
      totalRandomCB += negOutSample * tnCost;
      
      m_randomV.setText(Utils.doubleToString(totalRandomCB, 2));
      double gain = (m_costR.isSelected()) 
      ? totalRandomCB - CB
          : CB - totalRandomCB;
      m_gainV.setText(Utils.doubleToString(gain, 2));
      
      // update classification rate
      Instance currentInst = m_masterPlot.getPlotInstances().instance(index);
      double tp = currentInst.value(m_masterPlot.getPlotInstances().
          attribute(ThresholdCurve.TRUE_POS_NAME).index());
      double tn = currentInst.value(m_masterPlot.getPlotInstances().
          attribute(ThresholdCurve.TRUE_NEG_NAME).index());
      m_classificationAccV.
        setText(Utils.doubleToString((tp + tn) / (totalPos + totalNeg) * 100.0, 4) + "%");      
    }
    
    private void updateInfoGivenIndex(int index) {
      Instances plotInstances = m_masterPlot.getPlotInstances();
      int indexOfSampleSize = 
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.SAMPLE_SIZE_NAME).index();
      int indexOfPercOfTarget = 
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.RECALL_NAME).index();
      int indexOfThreshold =
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.THRESHOLD_NAME).index();
      
      // update labels
      m_percPopLab.setText(Utils.
          doubleToString(100.0 * plotInstances.instance(index).value(indexOfSampleSize), 4));
      m_percOfTargetLab.setText(Utils.doubleToString(
          100.0 * plotInstances.instance(index).value(indexOfPercOfTarget), 4));
      m_thresholdLab.setText(Utils.doubleToString(plotInstances.instance(index).value(indexOfThreshold), 4));
      /*if (m_percPop.isSelected()) {
        m_percPopLab.setText(Utils.doubleToString(100.0 * value, 4));
      } else if (m_percOfTarget.isSelected()) {
        m_percOfTargetLab.setText(Utils.doubleToString(100.0 * value, 4));
      } else {
        m_thresholdLab.setText(Utils.doubleToString(value, 4));
      }*/
      
      // Update the highlighted point on the graphs */
      if (m_previousShapeIndex >= 0) {
        m_shapeSizes[m_previousShapeIndex] = 1;
      }
     
      m_shapeSizes[index] = 10;
      m_previousShapeIndex = index;
      
      // Update the confusion matrix
//      double totalInstances = 
      int tp = plotInstances.attribute(ThresholdCurve.TRUE_POS_NAME).index();
      int fp = plotInstances.attribute(ThresholdCurve.FALSE_POS_NAME).index();
      int tn = plotInstances.attribute(ThresholdCurve.TRUE_NEG_NAME).index();
      int fn = plotInstances.attribute(ThresholdCurve.FALSE_NEG_NAME).index();
      Instance temp = plotInstances.instance(index);
      double totalInstances = temp.value(tp) + temp.value(fp) + temp.value(tn) + temp.value(fn);
      // get the value out of the total pop field (if possible)
      double requestedPopSize = totalInstances;
      try {
        requestedPopSize = Double.parseDouble(m_totalPopField.getText());
      } catch (NumberFormatException e) {}
      
      m_conf_aa.setCellValue(temp.value(tp), totalInstances, 
          requestedPopSize / totalInstances, 2);
      m_conf_ab.setCellValue(temp.value(fn), totalInstances, 
          requestedPopSize / totalInstances, 2);
      m_conf_ba.setCellValue(temp.value(fp), totalInstances, 
          requestedPopSize / totalInstances, 2);
      m_conf_bb.setCellValue(temp.value(tn), totalInstances, 
            requestedPopSize / totalInstances, 2);
      
      updateCBRandomGainInfo(index);
      
      repaint();
    }
    
    private void updateInfoForSliderValue(double value) {
      int indexOfSampleSize = 
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.SAMPLE_SIZE_NAME).index();
      int indexOfPercOfTarget = 
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.RECALL_NAME).index();
      int indexOfThreshold =
        m_masterPlot.getPlotInstances().attribute(ThresholdCurve.THRESHOLD_NAME).index();
      int indexOfMetric;
      
      if (m_percPop.isSelected()) {
        indexOfMetric = indexOfSampleSize;           
      } else if (m_percOfTarget.isSelected()) {
        indexOfMetric = indexOfPercOfTarget;
      } else {
        indexOfMetric = indexOfThreshold;
      }
      
      Instances plotInstances = m_masterPlot.getPlotInstances();
      int index = findIndexForValue(value, plotInstances, indexOfMetric);
      updateInfoGivenIndex(index);
    }
    
    private int findIndexForValue(double value, Instances plotInstances, int indexOfMetric) {
      // binary search
      // threshold curve is sorted ascending in the threshold (thus
      // descending for recall and pop size)
      int index = -1;
      int lower = 0;
      int upper = plotInstances.numInstances() - 1;
      int mid = (upper - lower) / 2;
      boolean done = false;
      while (!done) {
        if (upper - lower <= 1) {
          
          // choose the one closest to the value
          double comp1 = plotInstances.instance(upper).value(indexOfMetric);
          double comp2 = plotInstances.instance(lower).value(indexOfMetric);
          if (Math.abs(comp1 - value) < Math.abs(comp2 - value)) {
            index = upper;
          } else {
            index = lower;
          }
          
          break;
        }
        double comparisonVal = plotInstances.instance(mid).value(indexOfMetric);
        if (value > comparisonVal) {
          if (m_threshold.isSelected()) {
            lower = mid;
            mid += (upper - lower) / 2;
          } else {
            upper = mid;
            mid -= (upper - lower) / 2;
          }
        } else if (value < comparisonVal) {
          if (m_threshold.isSelected()) {
            upper = mid;
            mid -= (upper - lower) / 2;
          } else {
            lower = mid;
            mid += (upper - lower) / 2;
          }
        } else {
          index = mid;
          done = true;
        }
      }
      
      // now check for ties in the appropriate direction
      if (!m_threshold.isSelected()) {
        while (index + 1 < plotInstances.numInstances()) {
          if (plotInstances.instance(index + 1).value(indexOfMetric) == 
            plotInstances.instance(index).value(indexOfMetric)) {
            index++;
          } else {
            break;
          }
        }
      } else {
        while (index - 1 >= 0) {
          if (plotInstances.instance(index - 1).value(indexOfMetric) == 
            plotInstances.instance(index).value(indexOfMetric)) {
            index--;
          } else {
            break;
          } 
        }
      }
      return index;
    }
    
    /**
     * Set the threshold data for the panel to use.
     * 
     * @param data PlotData2D object encapsulating the threshold data.
     * @param classAtt the class attribute from the original data used to generate
     * the threshold data.
     * @throws Exception if something goes wrong.
     */
    public synchronized void setDataSet(PlotData2D data, Attribute classAtt) throws Exception {      
      // make a copy of the PlotData2D object
      m_masterPlot = new PlotData2D(data.getPlotInstances());
      boolean[] connectPoints = new boolean[m_masterPlot.getPlotInstances().numInstances()];
      for (int i = 1; i < connectPoints.length; i++) {
        connectPoints[i] = true;
      }
      m_masterPlot.setConnectPoints(connectPoints);

      m_masterPlot.m_alwaysDisplayPointsOfThisSize = 10;
      setClassForConfusionMatrix(classAtt);
      m_performancePanel.setMasterPlot(m_masterPlot);
      m_performancePanel.validate(); m_performancePanel.repaint();

      m_shapeSizes = new int[m_masterPlot.getPlotInstances().numInstances()];
      for (int i = 0; i < m_shapeSizes.length; i++) {
        m_shapeSizes[i] = 1;
      }
      m_masterPlot.setShapeSize(m_shapeSizes);
      constructCostBenefitData();
      m_costBenefitPanel.setMasterPlot(m_costBenefit);
      m_costBenefitPanel.validate(); m_costBenefitPanel.repaint();

      m_totalPopPrevious = 0;
      m_fpPrevious = 0;
      m_tpPrevious = 0;
      m_tnPrevious = 0;
      m_fnPrevious = 0;
      m_previousShapeIndex = -1;

      // set the total population size
      Instance first = m_masterPlot.getPlotInstances().instance(0);
      double totalPos = first.value(m_masterPlot.getPlotInstances().
          attribute(ThresholdCurve.TRUE_POS_NAME).index());
      double totalNeg = first.value(m_masterPlot.getPlotInstances().
          attribute(ThresholdCurve.FALSE_POS_NAME));
      m_originalPopSize = (int)(totalPos + totalNeg);
      m_totalPopField.setText("" + m_originalPopSize);

      m_performancePanel.setYIndex(5);
      m_performancePanel.setXIndex(10);
      m_costBenefitPanel.setXIndex(0);
      m_costBenefitPanel.setYIndex(1);
      //      System.err.println(m_masterPlot.getPlotInstances());
      updateInfoForSliderValue((double)m_thresholdSlider.getValue() / 100.0);
    }
    
    private void setClassForConfusionMatrix(Attribute classAtt) {
      m_classAttribute = classAtt;
      m_conf_actualA.setText(" Actual (a): " + classAtt.value(0));
      m_conf_actualA.setToolTipText(classAtt.value(0));
      String negClasses = "";
      for (int i = 1; i < classAtt.numValues(); i++) {
        negClasses += classAtt.value(i);
        if (i < classAtt.numValues() - 1) {
          negClasses += ",";
        }
      }
      m_conf_actualB.setText(" Actual (b): " + negClasses);
      m_conf_actualB.setToolTipText(negClasses);
    }
    
    private boolean constructCostBenefitData() {
      double tpCost = 0.0;
      try {
        tpCost = Double.parseDouble(m_cost_aa.getText());
      } catch (NumberFormatException n) {}
      double fpCost = 0.0;
      try {
        fpCost = Double.parseDouble(m_cost_ba.getText());
      } catch (NumberFormatException n) {}
      double tnCost = 0.0;
      try {
        tnCost = Double.parseDouble(m_cost_bb.getText());
      } catch (NumberFormatException n) {}
      double fnCost = 0.0;
      try {
        fnCost = Double.parseDouble(m_cost_ab.getText());
      } catch (NumberFormatException n) {}
      
      double requestedPopSize = m_originalPopSize;
      try {
        requestedPopSize = Double.parseDouble(m_totalPopField.getText());
      } catch (NumberFormatException e) {}
      
      double scaleFactor = 1.0;
      if (m_originalPopSize != 0) {
        scaleFactor = requestedPopSize / m_originalPopSize;
      }
      
      if (tpCost == m_tpPrevious && fpCost == m_fpPrevious &&
          tnCost == m_tnPrevious && fnCost == m_fnPrevious &&
          requestedPopSize == m_totalPopPrevious) {
        return false;
      }
      
      // First construct some Instances for the curve
      FastVector fv = new FastVector();
      fv.addElement(new Attribute("Sample Size"));
      fv.addElement(new Attribute("Cost/Benefit"));
      Instances costBenefitI = new Instances("Cost/Benefit Curve", fv, 100);
      
      // process the performance data to make this curve
      Instances performanceI = m_masterPlot.getPlotInstances();
      
      for (int i = 0; i < performanceI.numInstances(); i++) {
        Instance current = performanceI.instance(i);
        
        double[] vals = new double[2];
        vals[0] = current.value(10); // sample size
        vals[1] = (current.value(0) * tpCost
            + current.value(1) * fnCost
            + current.value(2) * fpCost
            + current.value(3) * tnCost) * scaleFactor;
        Instance newInst = new DenseInstance(1.0, vals);
        costBenefitI.add(newInst);
      }
      
      costBenefitI.compactify();
      
      // now set up the plot data
      m_costBenefit = new PlotData2D(costBenefitI);
      m_costBenefit.m_alwaysDisplayPointsOfThisSize = 10;
      m_costBenefit.setPlotName("Cost/benefit curve");
      boolean[] connectPoints = new boolean[costBenefitI.numInstances()];
      
      for (int i = 0; i < connectPoints.length; i++) {
        connectPoints[i] = true;
      }
      try {
        m_costBenefit.setConnectPoints(connectPoints);
        m_costBenefit.setShapeSize(m_shapeSizes);
      } catch (Exception ex) {
        // ignore
      }
      
      m_tpPrevious = tpCost;
      m_fpPrevious = fpCost;
      m_tnPrevious = tnCost;
      m_fnPrevious = fnCost;
      
      return true;
    }
  }
  
  /**
   * Constructor.
   */
  public CostBenefitAnalysis() {
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

  /**
   * Accept a threshold data event and set up the visualization.
   * @param e a threshold data event
   */
  public void acceptDataSet(ThresholdDataEvent e) {
    try {
      setCurveData(e.getDataSet(), e.getClassAttribute());
    } catch (Exception ex) {
      System.err.println("[CostBenefitAnalysis] Problem setting up visualization.");
      ex.printStackTrace();
    }
  }
  
  /**
   * Set the threshold curve data to use.
   * 
   * @param curveData a PlotData2D object set up with the curve data.
   * @param origClassAtt the class attribute from the original data used to
   * generate the curve.
   * @throws Exception if somthing goes wrong during the setup process.
   */
  public void setCurveData(PlotData2D curveData, Attribute origClassAtt) 
    throws Exception {
    if (m_analysisPanel == null) {
      m_analysisPanel = new AnalysisPanel();
    }
    m_analysisPanel.setDataSet(curveData, origClassAtt);
  }

  public BeanVisual getVisual() {
    return m_visual;
  }

  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"DefaultDataVisualizer.gif",
        BeanVisual.ICON_PATH+"DefaultDataVisualizer_animated.gif");
  }

  public Enumeration enumerateRequests() {
    Vector newVector = new Vector(0);
    if (m_analysisPanel != null) {
      if (m_analysisPanel.m_masterPlot != null) {
        newVector.addElement("Show analysis");
      }
    }
    return newVector.elements();
  }

  public void performRequest(String request) {
    if (request.compareTo("Show analysis") == 0) {
      try {
        // popup visualize panel
        if (!m_framePoppedUp) {
          m_framePoppedUp = true;

          final javax.swing.JFrame jf = 
            new javax.swing.JFrame("Cost/Benefit Analysis");
          jf.setSize(1000,600);
          jf.getContentPane().setLayout(new BorderLayout());
          jf.getContentPane().add(m_analysisPanel, BorderLayout.CENTER);
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
    } else {
      throw new IllegalArgumentException(request
          + " not supported (Cost/Benefit Analysis");
    }
  }

  public void addVetoableChangeListener(String name, VetoableChangeListener vcl) {
    m_bcSupport.addVetoableChangeListener(name, vcl);
  }

  public BeanContext getBeanContext() {
    return m_beanContext;
  }

  public void removeVetoableChangeListener(String name,
      VetoableChangeListener vcl) {
    m_bcSupport.removeVetoableChangeListener(name, vcl);
  }
  
  protected void appearanceFinal() {
    removeAll();
    setLayout(new BorderLayout());
    setUpFinal();
  }
  
  protected void setUpFinal() {
    if (m_analysisPanel == null) {
      m_analysisPanel = new AnalysisPanel();
    }
    add(m_analysisPanel, BorderLayout.CENTER);
  }
  
  protected void appearanceDesign() {
    removeAll();
    m_visual = new BeanVisual("CostBenefitAnalysis", 
                              BeanVisual.ICON_PATH+"ModelPerformanceChart.gif",
                              BeanVisual.ICON_PATH
                              +"ModelPerformanceChart_animated.gif");
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
  }

  public void setBeanContext(BeanContext bc) throws PropertyVetoException {
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
   * Returns true if, at this time, 
   * the object will accept a connection via the named event
   *
   * @param eventName the name of the event in question
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {
    return (m_listenee == null);
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
   * Notify this object that it has been deregistered as a listener with
   * a source for named event. This object is responsible
   * for recording this fact.
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public void disconnectionNotification(String eventName, Object source) {
    if (m_listenee == source) {
      m_listenee = null;
    }
    
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
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return false;
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
   * Set a logger
   *
   * @param logger a <code>weka.gui.Logger</code> value
   */
  public void setLog(Logger logger) {
    // we don't need to do any logging    
  }

  /**
   * Stop any processing that the bean might be doing.
   */
  public void stop() {
    // nothing to do here
  }
    
  public static void main(String[] args) {
    try {
      Instances train = new Instances(new java.io.BufferedReader(new java.io.FileReader(args[0])));
      train.setClassIndex(train.numAttributes() - 1);
      weka.classifiers.evaluation.ThresholdCurve tc = 
        new weka.classifiers.evaluation.ThresholdCurve();
      weka.classifiers.evaluation.EvaluationUtils eu = 
        new weka.classifiers.evaluation.EvaluationUtils();
      //weka.classifiers.Classifier classifier = new weka.classifiers.functions.Logistic();
      weka.classifiers.Classifier classifier = new weka.classifiers.bayes.NaiveBayes();
      FastVector predictions = new FastVector();
      eu.setSeed(1);
      predictions.appendElements(eu.getCVPredictions(classifier, train, 10));
      Instances result = tc.getCurve(predictions, 0);
      PlotData2D pd = new PlotData2D(result);
      pd.m_alwaysDisplayPointsOfThisSize = 10;

      boolean[] connectPoints = new boolean[result.numInstances()];
      for (int i = 1; i < connectPoints.length; i++) {
        connectPoints[i] = true;
      }
      pd.setConnectPoints(connectPoints);
      final javax.swing.JFrame jf = 
        new javax.swing.JFrame("CostBenefitTest");
      jf.setSize(1000,600);
      //jf.pack();
      jf.getContentPane().setLayout(new BorderLayout());
      final CostBenefitAnalysis.AnalysisPanel analysisPanel = 
        new CostBenefitAnalysis.AnalysisPanel();
      
      jf.getContentPane().add(analysisPanel, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
        public void windowClosing(java.awt.event.WindowEvent e) {
          jf.dispose();
          System.exit(0);
        }
      });
      
      jf.setVisible(true);
      
      analysisPanel.setDataSet(pd, train.classAttribute());
      
    } catch (Exception ex) {
      ex.printStackTrace();
    }
 
  }
}
