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
 *    BatchClustererListener.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.util.EventListener;

/**
 * Interface to something that can process a BatchClustererEvent
 *
 * @author <a href="mailto:mutter@cs.waikato.ac.nz">MStefan Mutter</a>
 * @version $Revision: 1.2 $
 * @see EventListener
 */
public interface BatchClustererListener extends EventListener {
  
  /**
   * Accept a BatchClustererEvent
   *
   * @param e a <code>BatchClustererEvent</code> value
   */
  void acceptClusterer(BatchClustererEvent e);
}
