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
 *    AttributeListPanel.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import weka.core.Instances;

import java.awt.BorderLayout;
import java.awt.Dimension;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumnModel;

/**
 * Creates a panel that displays the attributes contained in a set of
 * instances, letting the user select a single attribute for inspection.
 *
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.3 $
 */
public class AttributeListPanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -2030706987910400362L;

  /**
   * A table model that looks at the names of attributes.
   */
  class AttributeTableModel
    extends AbstractTableModel {

    /** for serialization */
    private static final long serialVersionUID = -7345701953670327707L;

    /** The instances who's attribute structure we are reporting */
    protected Instances m_Instances;
    
    /**
     * Creates the tablemodel with the given set of instances.
     *
     * @param instances the initial set of Instances
     */
    public AttributeTableModel(Instances instances) {

      setInstances(instances);
    }

    /**
     * Sets the tablemodel to look at a new set of instances.
     *
     * @param instances the new set of Instances.
     */
    public void setInstances(Instances instances) {

      m_Instances = instances;
    }
    
    /**
     * Gets the number of attributes.
     *
     * @return the number of attributes.
     */
    public int getRowCount() {
      
      return m_Instances.numAttributes();
    }
    
    /**
     * Gets the number of columns: 2
     *
     * @return 2
     */
    public int getColumnCount() {
      
      return 2;
    }
    
    /**
     * Gets a table cell
     *
     * @param row the row index
     * @param column the column index
     * @return the value at row, column
     */
    public Object getValueAt(int row, int column) {
      
      switch (column) {
      case 0:
	return new Integer(row + 1);
      case 1:
	return m_Instances.attribute(row).name();
      default:
	return null;
      }
    }
    
    /**
     * Gets the name for a column.
     *
     * @param column the column index.
     * @return the name of the column.
     */
    public String getColumnName(int column) {
      
      switch (column) {
      case 0:
	return new String("No.");
      case 1:
	return new String("Name");
      default:
	return null;
      }
    }
    
    /**
     * Gets the class of elements in a column.
     *
     * @param col the column index.
     * @return the class of elements in the column.
     */
    public Class getColumnClass(int col) {
      return getValueAt(0, col).getClass();
    }

    /**
     * Returns false
     *
     * @param row ignored
     * @param col ignored
     * @return false
     */
    public boolean isCellEditable(int row, int col) {

      return false;
    }
  }

  /** The table displaying attribute names */
  protected JTable m_Table = new JTable();

  /** The table model containing attribute names */
  protected AttributeTableModel m_Model;
  
  /**
   * Creates the attribute selection panel with no initial instances.
   */
  public AttributeListPanel() {

    m_Table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    m_Table.setColumnSelectionAllowed(false); 
    m_Table.setPreferredScrollableViewportSize(new Dimension(250, 150));

    setLayout(new BorderLayout());
    add(new JScrollPane(m_Table), BorderLayout.CENTER);
  }

  /**
   * Sets the instances who's attribute names will be displayed.
   *
   * @param newInstances the new set of instances
   */
  public void setInstances(Instances newInstances) {

    if (m_Model == null) {
      m_Model = new AttributeTableModel(newInstances);
      m_Table.setModel(m_Model);
      TableColumnModel tcm = m_Table.getColumnModel();
      tcm.getColumn(0).setMaxWidth(60);
      tcm.getColumn(1).setMinWidth(100);
    } else {
      m_Model.setInstances(newInstances);
    }
    m_Table.sizeColumnsToFit(-1);
    m_Table.revalidate();
    m_Table.repaint();
  }
  
  /**
   * Gets the selection model used by the table.
   *
   * @return a value of type 'ListSelectionModel'
   */
  public ListSelectionModel getSelectionModel() {

    return m_Table.getSelectionModel();
  }
  
  /**
   * Tests the attribute list panel from the command line.
   *
   * @param args must contain the name of an arff file to load.
   */
  public static void main(String[] args) {

    try {
      if (args.length == 0) {
	throw new Exception("supply the name of an arff file");
      }
      Instances i = new Instances(new java.io.BufferedReader(
				  new java.io.FileReader(args[0])));
      AttributeListPanel asp = new AttributeListPanel();
      final javax.swing.JFrame jf =
	new javax.swing.JFrame("Attribute List Panel");
      jf.getContentPane().setLayout(new BorderLayout());
      jf.getContentPane().add(asp, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	public void windowClosing(java.awt.event.WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
      asp.setInstances(i);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
  
} // AttributeListPanel
