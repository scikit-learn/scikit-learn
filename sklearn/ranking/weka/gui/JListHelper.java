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
 * JListHelper.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui;

import javax.swing.DefaultListModel;
import javax.swing.JList;

/**
 * A helper class for JList GUI elements with DefaultListModel or 
 * derived models.
 *
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.1 $
 * @see     JList
 * @see     DefaultListModel
 */
public class JListHelper {
  
  /** moves items up */
  public final static int MOVE_UP = 0;

  /** moves items down */
  public final static int MOVE_DOWN = 1;
  
  /**
   * moves the selected items by a certain amount of items in a given direction
   *
   * @param list        the JList to work on
   * @param moveby      the number of items to move by
   * @param direction   the direction to move in
   * @see               #MOVE_UP
   * @see               #MOVE_DOWN
   */
  protected static void moveItems(JList list, int moveby, int direction) {
    int[]               indices;
    int                 i;
    Object              o;
    DefaultListModel    model;

    model = (DefaultListModel) list.getModel();

    switch (direction) {
      case MOVE_UP:
        indices = list.getSelectedIndices();
        for (i = 0; i < indices.length; i++) {
          if (indices[i] == 0)
            continue;
          o = model.remove(indices[i]);
          indices[i] -= moveby;
          model.insertElementAt(o, indices[i]);
        }
        list.setSelectedIndices(indices);
        break;

      case MOVE_DOWN:
        indices = list.getSelectedIndices();
        for (i = indices.length - 1; i >= 0; i--) {
          if (indices[i] == model.getSize() - 1)
            continue;
          o = model.remove(indices[i]);
          indices[i] += moveby;
          model.insertElementAt(o, indices[i]);
        }
        list.setSelectedIndices(indices);
        break;

      default:
        System.err.println(
            JListHelper.class.getName() + ": direction '" 
            + direction + "' is unknown!");
    }
  }

  /**
   * moves the selected items up by 1
   *
   * @param list        the JList to work on
   */
  public static void moveUp(JList list) {
    if (canMoveUp(list))
      moveItems(list, 1, MOVE_UP);
  }

  /**
   * moves the selected item down by 1
   *
   * @param list        the JList to work on
   */
  public static void moveDown(JList list) {
    if (canMoveDown(list))
      moveItems(list, 1, MOVE_DOWN);
  }

  /**
   * moves the selected items to the top
   *
   * @param list        the JList to work on
   */
  public static void moveTop(JList list) {
    int[]     indices;
    int       diff;

    if (canMoveUp(list)) {
      indices = list.getSelectedIndices();
      diff    = indices[0];
      moveItems(list, diff, MOVE_UP);
    }
  }

  /**
   * moves the selected items to the end
   *
   * @param list        the JList to work on
   */
  public static void moveBottom(JList list) {
    int[]     indices;
    int       diff;

    if (canMoveDown(list)) {
      indices = list.getSelectedIndices();
      diff    = list.getModel().getSize() - 1 - indices[indices.length - 1];
      moveItems(list, diff, MOVE_DOWN);
    }
  }

  /**
   * checks whether the selected items can be moved up
   *
   * @param list        the JList to work on
   */
  public static boolean canMoveUp(JList list) {
    boolean   result;
    int[]     indices;

    result = false;
    
    indices = list.getSelectedIndices();
    if (indices.length > 0) {
      if (indices[0] > 0)
        result = true;
    }

    return result;
  }

  /**
   * checks whether the selected items can be moved down
   *
   * @param list        the JList to work on
   */
  public static boolean canMoveDown(JList list) {
    boolean   result;
    int[]     indices;

    result = false;
    
    indices = list.getSelectedIndices();
    if (indices.length > 0) {
      if (indices[indices.length - 1] < list.getModel().getSize() - 1)
        result = true;
    }

    return result;
  }
}
