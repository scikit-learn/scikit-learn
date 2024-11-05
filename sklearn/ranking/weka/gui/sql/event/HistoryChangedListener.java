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
 * HistoryChangedListener.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui.sql.event;

import java.util.EventListener;
import java.util.EventObject;

/**
 * A listener for changes in a history.
 *
 * @author      FracPete (fracpete at waikato dot ac dot nz)
 * @version     $Revision: 1.1 $
 */

public interface HistoryChangedListener extends EventListener {
  /**
   * This method gets called when a history is modified.
   */
  public void historyChanged(HistoryChangedEvent evt);
}
