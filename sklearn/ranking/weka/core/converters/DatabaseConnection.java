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
 *    DatabaseConnection.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.converters;

import weka.core.RevisionUtils;
import weka.experiment.DatabaseUtils;

import java.sql.DatabaseMetaData;
import java.sql.SQLException;

/**
 * Connects to a database.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Stefan Mutter (mutter@cs.waikato.ac.nz)
 * @version $Revision: 5953 $
 */
public class DatabaseConnection 
  extends DatabaseUtils {

  /** for serialization */
  static final long serialVersionUID = 1673169848863178695L;
  
  /**
   * Sets up the database drivers
   *
   * @throws Exception if an error occurs
   */
  public DatabaseConnection() throws Exception {
    super();
  }

  /** 
   * Check if the property checkUpperCaseNames in the DatabaseUtils file is 
   * set to true or false.
   *
   * @return  	true if the property checkUpperCaseNames in the DatabaseUtils 
   * 		file is set to true, false otherwise.
   */
  public boolean getUpperCase(){
    return m_checkForUpperCaseNames;
  }
  
  /**
   * Gets meta data for the database connection object.
   *
   * @return the meta data.
   * @throws Exception if an error occurs
   */
  public DatabaseMetaData getMetaData() throws Exception{
    if (!isConnected())
      throw new IllegalStateException("Not connected, please connect first!");
    
    return m_Connection.getMetaData();
  }
  
  /**
   * Dewtermines if the current query retrieves a result set or updates a table
   *
   * @return the update count (-1 if the query retrieves a result set).
   * @throws SQLException if an error occurs
   */
  public int getUpdateCount() throws SQLException {
    if (!isConnected())
      throw new IllegalStateException("Not connected, please connect first!");
    
    return m_PreparedStatement.getUpdateCount();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }
}
