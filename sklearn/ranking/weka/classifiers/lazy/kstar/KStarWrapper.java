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

/**
 *    KStarWrapper.java
 *    Copyright (C) 1995 Univeristy of Waikato
 *    Java port to Weka by Abdelaziz Mahoui (am14@cs.waikato.ac.nz).
 *
 */

package weka.classifiers.lazy.kstar;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

/*
 * @author Len Trigg (len@reeltwo.com)
 * @author Abdelaziz Mahoui (am14@cs.waikato.ac.nz)
 * @version $Revision 1.0 $
 */
public class KStarWrapper
  implements RevisionHandler {

  /** used/reused to hold the sphere size */
  public double sphere = 0.0;

  /** used/reused to hold the actual entropy */
  public double actEntropy = 0.0;

  /** used/reused to hold the random entropy */
  public double randEntropy = 0.0;

  /** used/reused to hold the average transformation probability */
  public double avgProb = 0.0;

  /** used/reused to hold the smallest transformation probability */
  public double minProb = 0.0;
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.7 $");
  }

}
