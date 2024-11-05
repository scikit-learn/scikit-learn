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
 *    AttributeSummaryPanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui;

import weka.core.Attribute;
import weka.core.AttributeStats;
import weka.core.Instances;
import weka.core.Utils;
import weka.core.labelranking.PreferenceAttribute;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.core.labelranking.RankUtilities;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;

/** 
 * This panel displays summary statistics about an attribute: name, type
 * number/% of missing/unique values, number of distinct values. For
 * numeric attributes gives some other stats (mean/std dev), for nominal
 * attributes gives counts for each attribute value.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5297 $
 */
public class AttributeSummaryPanel 
  extends JPanel {
  
  /** for serialization */
  static final long serialVersionUID = -5434987925737735880L;

  /** Message shown when no instances have been loaded and no attribute set */
  protected static final String NO_SOURCE = "None";

  /** Displays the name of the relation */
  protected JLabel m_AttributeNameLab = new JLabel(NO_SOURCE);
  
  /** Displays the type of attribute */
  protected JLabel m_AttributeTypeLab = new JLabel(NO_SOURCE);
  
  /** Displays the number of missing values */
  protected JLabel m_MissingLab = new JLabel(NO_SOURCE);
    
  /** Displays the number of unique values */
  protected JLabel m_UniqueLab = new JLabel(NO_SOURCE);
    
  /** Displays the number of distinct values */
  protected JLabel m_DistinctLab = new JLabel(NO_SOURCE);

  /** Displays other stats in a table */
  protected JTable m_StatsTable = new JTable() {
    /** for serialization */
    private static final long serialVersionUID = 7165142874670048578L;

    /**
     * returns always false, since it's just information for the user
     * 
     * @param row	the row
     * @param column	the column
     * @return		always false, i.e., the whole table is not editable
     */
    public boolean isCellEditable(int row, int column) {
      return false;
    }
  };
  
  /** The instances we're playing with */
  protected Instances m_Instances;

  /** Cached stats on the attributes we've summarized so far */
  protected AttributeStats [] m_AttributeStats;
  
  /** Do all instances have the same weight */
  protected boolean m_allEqualWeights = true;
  
  /**
   * Creates the instances panel with no initial instances.
   */
  public AttributeSummaryPanel() {

    JPanel simple = new JPanel();
    GridBagLayout gbL = new GridBagLayout();
    simple.setLayout(gbL);
    JLabel lab = new JLabel("Name:", SwingConstants.RIGHT);
    lab.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
    GridBagConstraints gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.EAST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 0;     gbC.gridx = 0;
    gbL.setConstraints(lab, gbC);
    simple.add(lab);
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.WEST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 0;     gbC.gridx = 1;
    gbC.weightx = 100; gbC.gridwidth = 3;
    gbL.setConstraints(m_AttributeNameLab, gbC);
    simple.add(m_AttributeNameLab);
    m_AttributeNameLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 10));
    
    lab = new JLabel("Type:", SwingConstants.RIGHT);
    lab.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.EAST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 0;     gbC.gridx = 4;
    gbL.setConstraints(lab, gbC);
    simple.add(lab);
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.WEST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 0;     gbC.gridx = 5;
    gbC.weightx = 100;
    gbL.setConstraints(m_AttributeTypeLab, gbC);
    simple.add(m_AttributeTypeLab);
    m_AttributeTypeLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 10));

    // Put into a separate panel?
    lab = new JLabel("Missing:", SwingConstants.RIGHT);
    lab.setBorder(BorderFactory.createEmptyBorder(0, 10, 5, 0));
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.EAST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 1;     gbC.gridx = 0;
    gbL.setConstraints(lab, gbC);
    simple.add(lab);
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.WEST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 1;     gbC.gridx = 1;
    gbC.weightx = 100;
    gbL.setConstraints(m_MissingLab, gbC);
    simple.add(m_MissingLab);
    m_MissingLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 5, 10));

    lab = new JLabel("Distinct:", SwingConstants.RIGHT);
    lab.setBorder(BorderFactory.createEmptyBorder(0, 10, 5, 0));
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.EAST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 1;     gbC.gridx = 2;
    gbL.setConstraints(lab, gbC);
    simple.add(lab);
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.WEST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 1;     gbC.gridx = 3;
    gbC.weightx = 100;
    gbL.setConstraints(m_DistinctLab, gbC);
    simple.add(m_DistinctLab);
    m_DistinctLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 5, 10));

    lab = new JLabel("Unique:", SwingConstants.RIGHT);
    lab.setBorder(BorderFactory.createEmptyBorder(0, 10, 5, 0));
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.EAST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 1;     gbC.gridx = 4;
    gbL.setConstraints(lab, gbC);
    simple.add(lab);
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.WEST;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 1;     gbC.gridx = 5;
    gbC.weightx = 100;
    gbL.setConstraints(m_UniqueLab, gbC);
    simple.add(m_UniqueLab);
    m_UniqueLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 5, 10));
    
    setLayout(new BorderLayout());
    add(simple, BorderLayout.NORTH);
    add(new JScrollPane(m_StatsTable), BorderLayout.CENTER);
    m_StatsTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
  }

  /**
   * Tells the panel to use a new set of instances.
   *
   * @param inst a set of Instances
   */
  public void setInstances(Instances inst) {
    
    m_Instances = inst;
    m_AttributeStats = new AttributeStats [inst.numAttributes()];
    m_AttributeNameLab.setText(NO_SOURCE);
    m_AttributeTypeLab.setText(NO_SOURCE);
    m_MissingLab.setText(NO_SOURCE);
    m_UniqueLab.setText(NO_SOURCE);
    m_DistinctLab.setText(NO_SOURCE);
    m_StatsTable.setModel(new DefaultTableModel());
    
    m_allEqualWeights = true;
    double w = m_Instances.instance(0).weight();
    for (int i = 1; i < m_Instances.numInstances(); i++) {
      if (m_Instances.instance(i).weight() != w) {
        m_allEqualWeights = false;
        break;
      }
    }
  }

  /**
   * Sets the attribute that statistics will be displayed for.
   *
   * @param index the index of the attribute to display
   */
  public void setAttribute(final int index) {

    setHeader(index);
    if (m_AttributeStats[index] == null) {
      Thread t = new Thread() {
	public void run() {
	  m_AttributeStats[index] = m_Instances
	  .attributeStats(index);
	  SwingUtilities.invokeLater(new Runnable() {
	    public void run() {
	      setDerived(index);
	      m_StatsTable.sizeColumnsToFit(-1);
	      m_StatsTable.revalidate();
	      m_StatsTable.repaint();
	    }
	  });
	}
      };
      t.setPriority(Thread.MIN_PRIORITY);
      t.start();
    } else {
      setDerived(index);
    }
  }
  
  /**
   * Sets the gui elements for fields that are stored in the AttributeStats
   * structure.
   * 
   * @param index	the index of the attribute
   */
  protected void setDerived(int index) {
    
    AttributeStats as = m_AttributeStats[index];
    long percent = Math.round(100.0 * as.missingCount / as.totalCount);
    m_MissingLab.setText("" + as.missingCount + " (" + percent + "%)");
    percent = Math.round(100.0 * as.uniqueCount / as.totalCount);
    m_UniqueLab.setText("" + as.uniqueCount + " (" + percent + "%)");
    m_DistinctLab.setText("" + as.distinctCount);
    setTable(as, index);
  }

  /**
   * Creates a tablemodel for the attribute being displayed
   * 
   * @param as		the attribute statistics
   * @param index	the index of the attribute
   */
  protected void setTable(AttributeStats as, int index) {

  if(m_Instances.attribute(index).isRanking()){
    if (as.nominalCounts != null) {
    	//RANKING BEGIN
  	  if(m_Instances.labels.size()==0)
		  	m_Instances.setLabels(RankUtilities.labels);
      Attribute att = m_Instances.attribute(index);
      Object[] colNames;
      if(!(att instanceof PreferenceAttribute || m_Instances.instance(0) instanceof PreferenceDenseInstance)){
    	colNames = new Object[] {"No.", "Label", "Count", "Weight"};
      }
      else{
    	  colNames = new Object[m_Instances.getLabels().size()+1];
    	  colNames[0]=">";
    	  for(int i=1; i<m_Instances.getLabels().size()+1; i++){
    		  colNames[i] = m_Instances.getLabels().get(i-1)[0];
    	  }
//    	colNames = new Object[]{"No.", "Rank", "Count", "Weight"};	  
      }

      Object [][] data;
      
      if(!(att instanceof PreferenceAttribute || m_Instances.instance(0) instanceof PreferenceDenseInstance)){
    	  data = new Object [as.nominalCounts.length][4];
    	  for (int i = 0; i < as.nominalCounts.length; i++) {
    		  data[i][0] = new Integer(i + 1);
    		  data[i][1] = att.value(i);
    		  data[i][2] = new Integer(as.nominalCounts[i]);
    		  data[i][3] = new Double(Utils.doubleToString(as.nominalWeights[i], 3));
    	  }
      }
      else{
    	  //Here, we count the amount of pairwise rankings in our data set,
    	  //displaying them in the PreprocessPanel.

    	  int[][] rankMatrix = new int[m_Instances.labels.size()][m_Instances.labels.size()];
    	  
    	  for(int b=0; b<m_Instances.numInstances(); b++){
    		  PreferenceDenseInstance current = (PreferenceDenseInstance) m_Instances.get(b);
    		  for(int c=0; ;c++){
    			  if(current.getHashMap().get(c)!=null){
    				  int[][] pairPrefs = current.getHashMap().get(c);
    				  for(int d=0; d<rankMatrix.length; d++){
    					  for(int e=0; e<rankMatrix.length; e++){
    						  if(pairPrefs[d][e]==1){
    							  rankMatrix[d][e]++;
    						  }
    					  }
    				  }
    				  break;
    			  }
    		  }
    	  }
    	  //putting our values from rankMatrix and the weights of all pairs
    	  //into the summary table.
    	  int border = rankMatrix.length;

    	  data = new Object [border][border+1];
    	  for (int i = 0; i < border; i++){
    		  for(int j=0; j < border+1; j++){ 
    			  if(j==0)
    				  data[i][j] = m_Instances.getLabels().get(i)[0].trim();
    			  else{
    				  data[i][j] = rankMatrix[i][j-1];
    			  }
    		  }
    	  }
    	  m_Instances.rankMatrix = rankMatrix;
      }
      //RANKING END
      m_StatsTable.setModel(new DefaultTableModel(data, colNames));
      m_StatsTable.getColumnModel().getColumn(0).setMaxWidth(60);
      DefaultTableCellRenderer tempR = new DefaultTableCellRenderer();
      tempR.setHorizontalAlignment(JLabel.RIGHT);
      m_StatsTable.getColumnModel().getColumn(0).setCellRenderer(tempR);
    } else if (as.numericStats != null) {
      Object [] colNames = {"Statistic", "Value"};
      Object [][] data = new Object [4][2];
      data[0][0] = "Minimum"; data[0][1] = Utils.doubleToString(as.numericStats.min, 3);
      data[1][0] = "Maximum"; data[1][1] = Utils.doubleToString(as.numericStats.max, 3);
      data[2][0] = "Mean" + ((!m_allEqualWeights) ? " (weighted)" : "");    
      data[2][1] = Utils.doubleToString(as.numericStats.mean, 3);
      data[3][0] = "StdDev" + ((!m_allEqualWeights) ? " (weighted)" : "");  
      data[3][1] = Utils.doubleToString(as.numericStats.stdDev, 3);
      m_StatsTable.setModel(new DefaultTableModel(data, colNames));
    } else {
      m_StatsTable.setModel(new DefaultTableModel());
    }
    m_StatsTable.getColumnModel().setColumnMargin(4);
  }
  else{
	  if (as.nominalCounts != null) {
	      Attribute att = m_Instances.attribute(index);
	      Object [] colNames = {"No.", "Label", "Count", "Weight"};
	      Object [][] data = new Object [as.nominalCounts.length][4];
	      for (int i = 0; i < as.nominalCounts.length; i++) {
		data[i][0] = new Integer(i + 1);
		data[i][1] = att.value(i);
		data[i][2] = new Integer(as.nominalCounts[i]);
		data[i][3] = new Double(Utils.doubleToString(as.nominalWeights[i], 3));
	      }
	      m_StatsTable.setModel(new DefaultTableModel(data, colNames));
	      m_StatsTable.getColumnModel().getColumn(0).setMaxWidth(60);
	      DefaultTableCellRenderer tempR = new DefaultTableCellRenderer();
	      tempR.setHorizontalAlignment(JLabel.RIGHT);
	      m_StatsTable.getColumnModel().getColumn(0).setCellRenderer(tempR);
	    } else if (as.numericStats != null) {
	      Object [] colNames = {"Statistic", "Value"};
	      Object [][] data = new Object [4][2];
	      data[0][0] = "Minimum"; data[0][1] = Utils.doubleToString(as.numericStats.min, 3);
	      data[1][0] = "Maximum"; data[1][1] = Utils.doubleToString(as.numericStats.max, 3);
	      data[2][0] = "Mean" + ((!m_allEqualWeights) ? " (weighted)" : "");    
	      data[2][1] = Utils.doubleToString(as.numericStats.mean, 3);
	      data[3][0] = "StdDev" + ((!m_allEqualWeights) ? " (weighted)" : "");  
	      data[3][1] = Utils.doubleToString(as.numericStats.stdDev, 3);
	      m_StatsTable.setModel(new DefaultTableModel(data, colNames));
	    } else {
	      m_StatsTable.setModel(new DefaultTableModel());
	    }
	    m_StatsTable.getColumnModel().setColumnMargin(4);
  }
  
  }
  
  /**
   * Sets the labels for fields we can determine just from the instance
   * header.
   * 
   * @param index	the index of the attribute
   */
  protected void setHeader(int index) {
    
    Attribute att = m_Instances.attribute(index);
    m_AttributeNameLab.setText(att.name());
    switch (att.type()) {
    case Attribute.NOMINAL:
    case PreferenceAttribute.RANKING:
    	if(att instanceof PreferenceAttribute || m_Instances.instance(0) instanceof PreferenceDenseInstance)m_AttributeTypeLab.setText("Ranking");
    	else
    		m_AttributeTypeLab.setText("Nominal");
      break;
    case Attribute.NUMERIC:
      m_AttributeTypeLab.setText("Numeric");
      break;
    case Attribute.STRING:
      m_AttributeTypeLab.setText("String");
      break;
    case Attribute.DATE:
      m_AttributeTypeLab.setText("Date");
      break;
    case Attribute.RELATIONAL:
      m_AttributeTypeLab.setText("Relational");
      break;
    default:
      m_AttributeTypeLab.setText("Unknown");
      break;
    }
    m_MissingLab.setText("...");
    m_UniqueLab.setText("...");
    m_DistinctLab.setText("...");
  }

  /**
   * Tests out the attribute summary panel from the command line.
   *
   * @param args optional name of dataset to load
   */
  public static void main(String [] args) {

    try {
      final javax.swing.JFrame jf = new javax.swing.JFrame("Attribute Panel");
      jf.getContentPane().setLayout(new BorderLayout());
      final AttributeSummaryPanel p = new AttributeSummaryPanel();
      p.setBorder(BorderFactory.createTitledBorder("Attribute"));
      jf.getContentPane().add(p, BorderLayout.CENTER);
      final javax.swing.JComboBox j = new javax.swing.JComboBox();
      j.setEnabled(false);
      j.addActionListener(new java.awt.event.ActionListener() {
	public void actionPerformed(java.awt.event.ActionEvent e) {
	  p.setAttribute(j.getSelectedIndex());
	}
      });
      jf.getContentPane().add(j, BorderLayout.NORTH);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	public void windowClosing(java.awt.event.WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
      if (args.length == 1) {
	java.io.Reader r = new java.io.BufferedReader(
			   new java.io.FileReader(args[0]));
	Instances inst = new Instances(r);
	p.setInstances(inst);
	p.setAttribute(0);
	String [] names = new String [inst.numAttributes()];
	for (int i = 0; i < names.length; i++) {
	  names[i] = inst.attribute(i).name();
	}
	j.setModel(new javax.swing.DefaultComboBoxModel(names));
	j.setEnabled(true);
      }
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
