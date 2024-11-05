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
 *    NodePlace.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.treevisualizer;

/**
 * This is an interface for classes that wish to take a node structure and 
 * arrange them
 *
 * @author Malcolm F Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public interface NodePlace {
 
  /**
   * The function to call to postion the tree that starts at Node r
   *
   * @param r The top of the tree.
   */
   void place(Node r);
  
} 
