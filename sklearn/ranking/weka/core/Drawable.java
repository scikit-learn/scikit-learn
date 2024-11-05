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
 *    Drawable.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

/** 
 * Interface to something that can be drawn as a graph.
 *
 * @author Ashraf M. Kibriya(amk14@cs.waikato.ac.nz), Eibe Frank(eibe@cs.waikato.ac.nz)
 * @version $Revision: 5961 $
 */
public interface Drawable {

  int NOT_DRAWABLE = 0, TREE = 1, BayesNet = 2, Newick = 3;

  /**
   * Returns the type of graph representing
   * the object.
   *
   * @return the type of graph representing the object
   */
  int graphType();

  /**
   * Returns a string that describes a graph representing
   * the object. The string should be in XMLBIF ver.
   * 0.3 format if the graph is a BayesNet, otherwise
   * it should be in dotty format.
   *
   * @return the graph described by a string
   * @exception Exception if the graph can't be computed
   */
  String graph() throws Exception;
}








