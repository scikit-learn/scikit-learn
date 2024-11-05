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
 *    TreeDisplayEvent.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.treevisualizer;


/**
 * An event containing the user selection from the tree display
 *
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public class TreeDisplayEvent {
  public static final int NO_COMMAND = 0;
  public static final int ADD_CHILDREN = 1;
  public static final int REMOVE_CHILDREN = 2;

  /** States that the user has accepted the tree. */
  public static final int ACCEPT = 3;

  /** Asks for another learning scheme to classify this node. */
  public static final int CLASSIFY_CHILD = 4;

  /** Command to remove instances from this node and send them to the 
   * VisualizePanel. */
  public static final int SEND_INSTANCES = 5;

  /** The int representing the action. */
  private int m_command;

  /** The id string for the node to alter. */
  private String m_nodeId;

  /**
   * Constructs an event with the specified command
   * and what the command is applied to.
   * @param ar The event type.
   * @param id The id string for the node to perform the action on.
   */
  public TreeDisplayEvent(int ar, String id) {
    m_command = 0;
    if (ar == 1 || ar == 2 || ar == 3 || ar == 4 || ar == 5) {
      //then command is good
      m_command = ar;
    }
    m_nodeId = id;
  }

  /**
   * @return The command.
   */
  public int getCommand() {
    return m_command;
  }
  
  /**
   * @return The id of the node.
   */
  public String getID() {
    return m_nodeId;
  }
}





