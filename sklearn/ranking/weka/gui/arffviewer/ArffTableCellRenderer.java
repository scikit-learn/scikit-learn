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
 * ArffTableCellRenderer.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.arffviewer;

import weka.core.Attribute;
import java.awt.Color;
import java.awt.Component;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.table.DefaultTableCellRenderer;

/**
 * Handles the background colors for missing values differently than the
 * DefaultTableCellRenderer.
 *
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.4 $ 
 */

public class ArffTableCellRenderer 
  extends DefaultTableCellRenderer {
  
  /** for serialization */
  static final long serialVersionUID = 9195794493301191171L;
  
  /** the color for missing values */
  private Color           missingColor;
  /** the color for selected missing values */
  private Color           missingColorSelected;
  /** the color for highlighted values */
  private Color           highlightColor;
  /** the color for selected highlighted values */
  private Color           highlightColorSelected;
  
  /**
   * initializes the Renderer with a standard color
   */
  public ArffTableCellRenderer() {
    this( new Color(223, 223, 223), 
        new Color(192, 192, 192) );
  }
  
  /**
   * initializes the Renderer with the given colors
   * 
   * @param missingColor		the color for missing values
   * @param missingColorSelected	the color selected missing values
   */
  public ArffTableCellRenderer( Color missingColor, 
      Color missingColorSelected ) {
    this( missingColor,
        missingColorSelected,
        Color.RED,
        Color.RED.darker() );
  }
  
  /**
   * initializes the Renderer with the given colors
   * 
   * @param missingColor		the color for missing values
   * @param missingColorSelected	the color selected missing values
   * @param highlightColor		the color for highlighted values
   * @param highlightColorSelected	the color selected highlighted values
   */
  public ArffTableCellRenderer( Color missingColor, 
      Color missingColorSelected,
      Color highlightColor,
      Color highlightColorSelected ) {
    super();
    
    this.missingColor           = missingColor;
    this.missingColorSelected   = missingColorSelected;
    this.highlightColor         = highlightColor;
    this.highlightColorSelected = highlightColorSelected;
  }
  
  /**
   * Returns the default table cell renderer.
   * stuff for the header is taken from <a href="http://www.chka.de/swing/table/faq.html">here</a>
   * 
   * @param table		the table this object belongs to
   * @param value		the actual cell value
   * @param isSelected		whether the cell is selected
   * @param hasFocus		whether the cell has the focus
   * @param row			the row in the table
   * @param column		the column in the table
   * @return			the rendering component
   */
  public Component getTableCellRendererComponent(
      JTable table, Object value, boolean isSelected, 
      boolean hasFocus, int row, int column ) {
    ArffSortedTableModel            model;
    Component                  result;
    String                     searchString;
    boolean                    found;
    
    result = super.getTableCellRendererComponent(
        table, value, isSelected, hasFocus, row, column);
    
    // search
    if (table instanceof ArffTable)
      searchString = ((ArffTable) table).getSearchString();
    else
      searchString = null;
    if ( (searchString != null) && (!searchString.equals("")) )
      found = (searchString.equals(value.toString()));
    else
      found = false;
    
    if (table.getModel() instanceof ArffSortedTableModel) {
      model = (ArffSortedTableModel) table.getModel();
      // normal cell
      if (row >= 0) {
        if (model.isMissingAt(row, column)) {
          setToolTipText("missing");
          if (found) {
            if (isSelected)
              result.setBackground(highlightColorSelected);
            else
              result.setBackground(highlightColor);
          }
          else {
            if (isSelected)
              result.setBackground(missingColorSelected);
            else
              result.setBackground(missingColor);
          }
        }
        else {
          setToolTipText(null);
          if (found) {
            if (isSelected)
              result.setBackground(highlightColorSelected);
            else
              result.setBackground(highlightColor);
          }
          else {
            if (isSelected)
              result.setBackground(table.getSelectionBackground());
            else
              result.setBackground(Color.WHITE);
          }
        }
        
        // alignment
        if (model.getType(row, column) == Attribute.NUMERIC)
          setHorizontalAlignment(SwingConstants.RIGHT);
        else
          setHorizontalAlignment(SwingConstants.LEFT);
      }
      // header
      else {
        setBorder(UIManager.getBorder("TableHeader.cellBorder"));
        setHorizontalAlignment(SwingConstants.CENTER);
        if (table.getColumnModel().getSelectionModel().isSelectedIndex(column))
          result.setBackground(UIManager.getColor("TableHeader.background").darker());
        else
          result.setBackground(UIManager.getColor("TableHeader.background"));
      }
    }
    
    return result;
  }
}

