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
 * Stemmer.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.stemmers;

import weka.core.RevisionHandler;

import java.io.Serializable;

/**
 * Interface for all stemming algorithms.
 *
 * @author    FracPete (fracpete at waikato dot ac dot nz)
 * @version   $Revision: 5953 $
 */
public interface Stemmer 
  extends Serializable, RevisionHandler {

  /**
   * Stems the given word and returns the stemmed version
   *
   * @param word      the unstemmed word
   * @return          the stemmed word
   */
  public String stem(String word);
}
