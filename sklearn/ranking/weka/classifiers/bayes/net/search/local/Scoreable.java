/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * Scoreable.java
 * Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 * 
 */
package weka.classifiers.bayes.net.search.local;

/**
 * Interface for allowing to score a classifier
 * 
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 1.4 $
 */
public interface Scoreable {

  /**
   * score types
   */
  int BAYES = 0;
  int BDeu = 1;
  int MDL = 2;
  int ENTROPY = 3;
  int AIC = 4;

  /**
   * Returns log-score
   * 
   * @param nType the score type
   * @return the log-score
   */
  double logScore(int nType, int nCardinality);
}    // interface Scoreable




