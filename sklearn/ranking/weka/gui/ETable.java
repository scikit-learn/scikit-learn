/* 
 * This software is provided under the terms 
 * of the GNU Lesser General Public License, Version 2.1. You may not use 
 * this file except in compliance with the license. If you need a copy of the license, 
 * please go to http://www.gnu.org/licenses/lgpl-2.1.txt. The Original Code is Pentaho 
 * Data Integration.  The Initial Developer is Pentaho Corporation.
 *
 * Software distributed under the GNU Lesser Public License is distributed on an "AS IS" 
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or  implied. Please refer to 
 * the license for the specific language governing your rights and limitations.*/

/*
 * VersionPackageConstraint
 * Copyright (c) 2006 Elliot Hughes 
 */

package weka.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;

import javax.swing.BorderFactory;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JViewport;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;

/**
 * A better-looking table than JTable. In particular, on Mac OS this looks
 * more like a Cocoa table than the default Aqua LAF manages.
 *
 * @author Elliott Hughes
 */
public class ETable extends JTable {
  /**
   * For serialization
   */
  private static final long serialVersionUID = -3028630226368293049L;
  
  private final Color MAC_FOCUSED_SELECTED_CELL_HORIZONTAL_LINE_COLOR = new Color(0x7daaea);
  private final Color MAC_UNFOCUSED_SELECTED_CELL_HORIZONTAL_LINE_COLOR = new Color(0xe0e0e0);

  private final Color MAC_UNFOCUSED_SELECTED_CELL_BACKGROUND_COLOR = new Color(0xc0c0c0);

  private final Color MAC_FOCUSED_UNSELECTED_VERTICAL_LINE_COLOR = new Color(0xd9d9d9);
  private final Color MAC_FOCUSED_SELECTED_VERTICAL_LINE_COLOR = new Color(0x346dbe);
  private final Color MAC_UNFOCUSED_UNSELECTED_VERTICAL_LINE_COLOR = new Color(0xd9d9d9);
  private final Color MAC_UNFOCUSED_SELECTED_VERTICAL_LINE_COLOR = new Color(0xacacac);

  private final Color MAC_OS_ALTERNATE_ROW_COLOR = new Color(0.92f, 0.95f, 0.99f);

  public ETable() {
    // Although it's the JTable default, most systems' tables don't draw a grid by default.
    // Worse, it's not easy (or possible?) for us to take over grid painting ourselves for those LAFs (Metal, for example) that do paint grids.
    // The Aqua and GTK LAFs ignore the grid settings anyway, so this causes no change there.
    setShowGrid(false);

    // Tighten the cells up, and enable the manual painting of the vertical grid lines.
    setIntercellSpacing(new Dimension());

    // Table column re-ordering is too badly implemented to enable.
    getTableHeader().setReorderingAllowed(false);

    if (System.getProperty("os.name").contains("Mac")) {
      // Work-around for Apple 4352937.
      JLabel.class.cast(getTableHeader().getDefaultRenderer()).setHorizontalAlignment(SwingConstants.LEADING);

      // Use an iTunes-style vertical-only "grid".
      setShowHorizontalLines(false);
      setShowVerticalLines(true);
    }
  }

  /**
   * Paints empty rows too, after letting the UI delegate do
   * its painting.
   */
  public void paint(Graphics g) {
    super.paint(g);
    paintEmptyRows(g);
  }

  /**
   * Paints the backgrounds of the implied empty rows when the
   * table model is insufficient to fill all the visible area
   * available to us. We don't involve cell renderers, because
   * we have no data.
   */
  protected void paintEmptyRows(Graphics g) {
    final int rowCount = getRowCount();
    final Rectangle clip = g.getClipBounds();
    final int height = clip.y + clip.height;
    if (rowCount * rowHeight < height) {
      for (int i = rowCount; i <= height/rowHeight; ++i) {
        g.setColor(colorForRow(i));
        g.fillRect(clip.x, i * rowHeight, clip.width, rowHeight);
      }

      // Mac OS' Aqua LAF never draws vertical grid lines, so we have to draw them ourselves.
      if (System.getProperty("os.name").contains("Mac") && getShowVerticalLines()) {
        g.setColor(MAC_UNFOCUSED_UNSELECTED_VERTICAL_LINE_COLOR);
        TableColumnModel columnModel = getColumnModel();
        int x = 0;
        for (int i = 0; i < columnModel.getColumnCount(); ++i) {
          TableColumn column = columnModel.getColumn(i);
          x += column.getWidth();
          g.drawLine(x - 1, rowCount * rowHeight, x - 1, height);
        }
      }
    }
  }

  /**
   * Changes the behavior of a table in a JScrollPane to be more like
   * the behavior of JList, which expands to fill the available space.
   * JTable normally restricts its size to just what's needed by its
   * model.
   */
  public boolean getScrollableTracksViewportHeight() {
    if (getParent() instanceof JViewport) {
      JViewport parent = (JViewport) getParent();
      return (parent.getHeight() > getPreferredSize().height);
    }
    return false;
  }

  /**
   * Returns the appropriate background color for the given row.
   */
  protected Color colorForRow(int row) {
    return (row % 2 == 0) ? alternateRowColor() : getBackground();
  }

  private Color alternateRowColor() {
    return UIManager.getLookAndFeel().getClass().getName().contains("GTK") ? Color.WHITE : MAC_OS_ALTERNATE_ROW_COLOR;
  }

  /**
   * Shades alternate rows in different colors.
   */
  public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
    Component c = super.prepareRenderer(renderer, row, column);
    boolean focused = hasFocus();
    boolean selected = isCellSelected(row, column);
    if (selected) {
      if (System.getProperty("os.name").contains("Mac") && focused == false) {
        // Native Mac OS renders the selection differently if the table doesn't have the focus.
        // The Mac OS LAF doesn't imitate this for us.
        c. setBackground(MAC_UNFOCUSED_SELECTED_CELL_BACKGROUND_COLOR);
        c.setForeground(UIManager.getColor("Table.foreground"));
      } else {
        c.setBackground(UIManager.getColor("Table.selectionBackground"));
        c.setForeground(UIManager.getColor("Table.selectionForeground"));
      }
    } else {
      // Outside of selected rows, we want to alternate the background color.
      c.setBackground(colorForRow(row));
      c.setForeground(UIManager.getColor("Table.foreground"));
    }

    if (c instanceof JComponent) {
      JComponent jc = (JComponent) c;

      // The Java 6 GTK LAF JCheckBox doesn't paint its background by default.
      // Sun 5043225 says this is the intended behavior, though presumably not when it's being used as a table cell renderer.
      if (UIManager.getLookAndFeel().getClass().getName().contains("GTK") && c instanceof JCheckBox) {
        jc.setOpaque(true);
      }

      if (getCellSelectionEnabled() == false && isEditing() == false) {
        if (System.getProperty("os.name").contains("Mac")) {
          // Native Mac OS doesn't draw a border on the selected cell.
          // It does however draw a horizontal line under the whole row, and a vertical line separating each column.
          fixMacOsCellRendererBorder(jc, selected, focused);
        } else {
          // FIXME: doesn't Windows have row-wide selection focus?
          // Hide the cell focus.
          jc.setBorder(null);
        }
      }

      initToolTip(jc, row, column);
    }

    return c;
  }

  private void fixMacOsCellRendererBorder(JComponent renderer, boolean selected, boolean focused) {
    Border border;
    if (selected) {
      border = BorderFactory.createMatteBorder(0, 0, 1, 0, focused ? MAC_FOCUSED_SELECTED_CELL_HORIZONTAL_LINE_COLOR : MAC_UNFOCUSED_SELECTED_CELL_HORIZONTAL_LINE_COLOR);
    } else {
      border = BorderFactory.createEmptyBorder(0, 0, 1, 0);
    }

    // Mac OS' Aqua LAF never draws vertical grid lines, so we have to draw them ourselves.
    if (getShowVerticalLines()) {
      Color verticalLineColor;
      if (focused) {
        verticalLineColor = selected ? MAC_FOCUSED_SELECTED_VERTICAL_LINE_COLOR : MAC_FOCUSED_UNSELECTED_VERTICAL_LINE_COLOR;
      } else {
        verticalLineColor = selected ? MAC_UNFOCUSED_SELECTED_VERTICAL_LINE_COLOR : MAC_UNFOCUSED_UNSELECTED_VERTICAL_LINE_COLOR;
      }
      Border verticalBorder = BorderFactory.createMatteBorder(0, 0, 0, 1, verticalLineColor);
      border = BorderFactory.createCompoundBorder(border, verticalBorder);
    }

    renderer.setBorder(border);
  }

  /**
   * Sets the component's tool tip if the component is being rendered smaller than its preferred size.
   * This means that all users automatically get tool tips on truncated text fields that show them the full value.
   */
  private void initToolTip(JComponent c, int row, int column) {
    String toolTipText = null;
    if (c.getPreferredSize().width > getCellRect(row, column, false).width) {
      toolTipText = getValueAt(row, column).toString();
    }
    c.setToolTipText(toolTipText);
  }

  /**
   * Places tool tips over the cell they correspond to. MS Outlook does this, and it works rather well.
   * Swing will automatically override our suggested location if it would cause the tool tip to go off the display.
   */
  public Point getToolTipLocation(MouseEvent e) {
    // After a tool tip has been displayed for a cell that has a tool tip, cells without tool tips will show an empty tool tip until the tool tip mode times out (or the table has a global default tool tip).
    // (ToolTipManager.checkForTipChange considers a non-null result from getToolTipText *or* a non-null result from getToolTipLocation as implying that the tool tip should be displayed. This seems like a bug, but that's the way it is.)
    if (getToolTipText(e) == null) {
      return null;
    }
    final int row = rowAtPoint(e.getPoint());
    final int column = columnAtPoint(e.getPoint());
    if (row == -1 || column == -1) {
      return null;
    }
    return getCellRect(row, column, false).getLocation();
  }

  /**
   * Improve the appearance of of a table in a JScrollPane on Mac OS, where there's otherwise an unsightly hole.
   */
  @Override
  protected void configureEnclosingScrollPane() {
    super.configureEnclosingScrollPane();

    if (System.getProperty("os.name").contains("Mac") == false) {
      return;
    }

    Container p = getParent();
    if (p instanceof JViewport) {
      Container gp = p.getParent();
      if (gp instanceof JScrollPane) {
        JScrollPane scrollPane = (JScrollPane)gp;
        // Make certain we are the viewPort's view and not, for
        // example, the rowHeaderView of the scrollPane -
        // an implementor of fixed columns might do this.
        JViewport viewport = scrollPane.getViewport();
        if (viewport == null || viewport.getView() != this) {
          return;
        }

        // JTable copy & paste above this point; our code below.

        // Put a dummy header in the upper-right corner.
        final Component renderer = new JTableHeader().getDefaultRenderer().getTableCellRendererComponent(null, "", false, false, -1, 0);
        JPanel panel = new JPanel(new BorderLayout());
        panel.add(renderer, BorderLayout.CENTER);
        scrollPane.setCorner(JScrollPane.UPPER_RIGHT_CORNER, panel);
      }
    }
  }
}

