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
 * ScriptExecutionEvent.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.scripting.event;

import weka.gui.scripting.Script;

import java.util.EventObject;

/**
 * Event that gets sent when a script is executed.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5142 $
 */
public class ScriptExecutionEvent
  extends EventObject {

  /** for serialization. */
  private static final long serialVersionUID = -8357216611114356632L;

  /**
   * Defines the type of event.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5142 $
   */
  public enum Type {
    /** started execution. */
    STARTED,
    /** finished normal. */
    FINISHED,
    /** finished with error. */
    ERROR,
    /** got stopped by user. */
    STOPPED
  }
  
  /** the type of event. */
  protected Type m_Type;
  
  /** optional additional information. */
  protected Object m_Additional;
  
  /**
   * Initializes the event.
   * 
   * @param source	the script that triggered the event
   * @param type	the type of finish
   */
  public ScriptExecutionEvent(Script source, Type type) {
    this(source, type, null);
  }
  
  /**
   * Initializes the event.
   * 
   * @param source	the script that triggered the event
   * @param type	the type of finish
   * @param additional	additional information, can be null
   */
  public ScriptExecutionEvent(Script source, Type type, Object additional) {
    super(source);
    
    m_Type       = type;
    m_Additional = additional;
  }
  
  /**
   * Returns the script that triggered the event.
   * 
   * @return		the script
   */
  public Script getScript() {
    return (Script) getSource();
  }
  
  /**
   * Returns the type of event.
   * 
   * @return		the type
   */
  public Type getType() {
    return m_Type;
  }
  
  /**
   * Returns whether additional information is available.
   * 
   * @return		true if additional information is available
   * @see		#getAdditional()
   */
  public boolean hasAdditional() {
    return (m_Additional != null);
  }
  
  /**
   * Returns the additional information.
   * 
   * @return		the additional information, can be null
   */
  public Object getAdditional() {
    return m_Additional;
  }
}
