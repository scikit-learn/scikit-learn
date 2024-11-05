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
 *    BeanVisual.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Point;
import java.awt.Toolkit;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * BeanVisual encapsulates icons and label for a given bean. Has methods
 * to load icons, set label text and toggle between static and animated
 * versions of a bean's icon.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.10 $
 * @since 1.0
 * @see JPanel
 * @see Serializable
 */
public class BeanVisual
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -6677473561687129614L;

  public static final String ICON_PATH="weka/gui/beans/icons/";

  public static final int NORTH_CONNECTOR = 0;
  public static final int SOUTH_CONNECTOR = 1;
  public static final int EAST_CONNECTOR = 2;
  public static final int WEST_CONNECTOR = 3;

  /**
   * Holds name (including path) of the static icon
   */
  protected String m_iconPath;

  /**
   * Holds name (including path) of the animated icon
   */
  protected String m_animatedIconPath;

  /**
   * ImageIcons for the icons. Is transient because for some reason
   * animated gifs cease to be animated after restoring from serialization.
   * Icons are re-loaded from source after deserialization
   */
  protected transient ImageIcon m_icon;
  protected transient ImageIcon m_animatedIcon;

  /**
   * Name for the bean
   */
  protected String m_visualName;

  protected JLabel m_visualLabel;
  
  /**
   * Container for the icon
   */
  //  protected IconHolder m_visualHolder;

  //  protected JLabel m_textLabel;
  private boolean m_stationary = true;

  private PropertyChangeSupport m_pcs = new PropertyChangeSupport(this);
  
  private boolean m_displayConnectors = false;
  private Color m_connectorColor = Color.blue;

  /**
   * Constructor
   *
   * @param visualName name for the bean
   * @param iconPath path to the icon file
   * @param animatedIconPath path to the animated icon file
   */
  public BeanVisual(String visualName, String iconPath, 
		    String animatedIconPath) {

    loadIcons(iconPath, animatedIconPath);
    m_visualName = visualName;
    //    m_textLabel = new JLabel(m_visualName, JLabel.CENTER);
    m_visualLabel = new JLabel(m_icon);

    setLayout(new BorderLayout());
   
    //    m_visualHolder = new IconHolder(m_visualLabel);
    
    add(m_visualLabel, BorderLayout.CENTER);
    Dimension d = m_visualLabel.getPreferredSize();
    //      this.setSize((int)d.getWidth()+50, (int)d.getHeight()+50);
    Dimension d2 = new Dimension((int)d.getWidth() + 10, 
				 (int)d.getHeight() + 10);
    setMinimumSize(d2);
    setPreferredSize(d2);
    setMaximumSize(d2);
  }

  /**
   * Reduce this BeanVisual's icon size by the given factor
   *
   * @param factor the factor by which to reduce the icon size by
   */
  public void scale(int factor) {
    if (m_icon != null) {
      removeAll();
      Image pic = m_icon.getImage();
      int width = m_icon.getIconWidth();
      int height = m_icon.getIconHeight();
      int reduction = width / factor;
      width -= reduction;
      height -= reduction;
      pic = pic.getScaledInstance(width, height, Image.SCALE_SMOOTH);
      m_icon = new ImageIcon(pic);
      m_visualLabel = new JLabel(m_icon);
      add(m_visualLabel, BorderLayout.CENTER);
      Dimension d = m_visualLabel.getPreferredSize();
      //      this.setSize((int)d.getWidth()+50, (int)d.getHeight()+50);
      Dimension d2 = new Dimension((int)d.getWidth() + 10, 
				   (int)d.getHeight() + 10);
      setMinimumSize(d2);
      setPreferredSize(d2);
      setMaximumSize(d2);   
    }
  }

  /**
   * Loads static and animated versions of a beans icons. These are
   * assumed to be defined in the system resource location (i.e. in the
   * CLASSPATH). If the named icons do not exist, no changes to the
   * visual appearance is made. Since default icons for generic
   * types of beans (eg. DataSource, Classifier etc)
   * are assumed to exist, it allows developers to add custom icons for
   * for specific instantiations of these beans 
   * (eg. J48, DiscretizeFilter etc) at their leisure.
   *
   * @param iconPath path to
   * @param animatedIconPath a <code>String</code> value
   */
  public boolean loadIcons(String iconPath, String animatedIconPath) {
    boolean success = true;
    //    java.net.URL imageURL = ClassLoader.getSystemResource(iconPath);
    java.net.URL imageURL = this.getClass().getClassLoader().getResource(iconPath);
    if (imageURL == null) {
      //      System.err.println("Warning: unable to load "+iconPath);
    } else {
      Image pic = Toolkit.getDefaultToolkit().
	getImage(imageURL);

      m_icon = new ImageIcon(pic);
      if (m_visualLabel != null) {
	m_visualLabel.setIcon(m_icon);
      }
    }
    
    //    imageURL = ClassLoader.getSystemResource(animatedIconPath);
    imageURL = this.getClass().getClassLoader().getResource(animatedIconPath);
    if (imageURL == null) {
      //      System.err.println("Warning: unable to load "+animatedIconPath);
      success = false;
    } else {
      Image pic2 = Toolkit.getDefaultToolkit().
	getImage(imageURL);
      m_animatedIcon = new ImageIcon(pic2);
    }
    m_iconPath = iconPath;
    m_animatedIconPath = animatedIconPath;
    return success;
  }

  /**
   * Set the label for the visual. Informs any property change listeners
   *
   * @param text the label
   */
  public void setText(String text) {
    m_visualName = text;
    //    m_textLabel.setText(m_visualName);
    m_pcs.firePropertyChange("label",null,null);
  }

  /**
   * Get the visual's label
   *
   * @return a <code>String</code> value
   */
  public String getText() {
    return m_visualName;
  }

  /**
   * Set the static version of the icon
   *
   */
  public void setStatic() {
    m_visualLabel.setIcon(m_icon);
  }

  /**
   * Set the animated version of the icon
   *
   */
  public void setAnimated() {
    m_visualLabel.setIcon(m_animatedIcon);
  }

  /**
   * Returns the coordinates of the closest "connector" point to the
   * supplied point. Coordinates are in the parent containers coordinate
   * space.
   *
   * @param pt the reference point
   * @return the closest connector point
   */
  public Point getClosestConnectorPoint(Point pt) {
    int sourceX = getParent().getX();
    int sourceY = getParent().getY();
    int sourceWidth = getWidth();
    int sourceHeight = getHeight();
    int sourceMidX = sourceX + (sourceWidth / 2);
    int sourceMidY = sourceY + (sourceHeight / 2);
    int x = (int)pt.getX();
    int y = (int)pt.getY();
    
    Point closest = new Point();
    int cx = (Math.abs(x - sourceMidX) < Math.abs(y - sourceMidY)) ? 
      sourceMidX :
      ((x < sourceMidX) ? sourceX : sourceX + sourceWidth);
    int cy = (Math.abs(y - sourceMidY) < Math.abs(x - sourceMidX)) ? 
      sourceMidY :
      ((y < sourceMidY) ? sourceY : sourceY + sourceHeight) ;
    closest.setLocation(cx, cy);
    return closest;
  }

  /**
   * Returns the coordinates of the connector point given a compass point
   *
   * @param compassPoint a compass point
   * @return a <code>Point</code> value
   */
  public Point getConnectorPoint(int compassPoint) {
    int sourceX = getParent().getX();
    int sourceY = getParent().getY();
    int sourceWidth = getWidth();
    int sourceHeight = getHeight();
    int sourceMidX = sourceX + (sourceWidth / 2);
    int sourceMidY = sourceY + (sourceHeight / 2);

    switch (compassPoint) {
    case NORTH_CONNECTOR : return new Point(sourceMidX, sourceY);
    case SOUTH_CONNECTOR : return new Point(sourceMidX, sourceY+sourceHeight);
    case WEST_CONNECTOR :  return new Point(sourceX, sourceMidY);
    case EAST_CONNECTOR :  return new Point(sourceX+sourceWidth, sourceMidY);
    default : System.err.println("Unrecognised connectorPoint (BeanVisual)");
    }
    return new Point(sourceX, sourceY);
  }

  /**
   * Returns the static icon
   *
   * @return an <code>ImageIcon</code> value
   */
  public ImageIcon getStaticIcon() {
    return m_icon;
  }

  /**
   * Returns the animated icon
   *
   * @return an <code>ImageIcon</code> value
   */
  public ImageIcon getAnimatedIcon() {
    return m_animatedIcon;
  }

  /**
   * returns the path for the icon
   * 
   * @return the path for the icon
   */
  public String getIconPath() {
    return m_iconPath;
  }

  /**
   * returns the path for the animated icon
   * 
   * @return the path for the animated icon
   */
  public String getAnimatedIconPath() {
    return m_animatedIconPath;
  }

  /**
   * Turn on/off the connector points
   *
   * @param dc a <code>boolean</code> value
   */
  public void setDisplayConnectors(boolean dc) {
    //    m_visualHolder.setDisplayConnectors(dc);
    m_displayConnectors = dc;
    m_connectorColor = Color.blue;
    repaint();
  }

  /**
   * Turn on/off the connector points
   *
   * @param dc a <code>boolean</code> value
   * @param c the Color to use
   */
  public void setDisplayConnectors(boolean dc,
                                   Color c) {
    setDisplayConnectors(dc);
    m_connectorColor = c;
  }

  /**
   * Add a listener for property change events
   *
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void addPropertyChangeListener(PropertyChangeListener pcl) {
    m_pcs.addPropertyChangeListener(pcl);
  }

  /**
   * Remove a property change listener
   *
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void removePropertyChangeListener(PropertyChangeListener pcl) {
    m_pcs.removePropertyChangeListener(pcl);
  }

  public void paintComponent(Graphics gx) {
    super.paintComponent(gx);
    if (m_displayConnectors) {
      gx.setColor(m_connectorColor);
      
      int midx = (int)(this.getWidth() / 2.0);
      int midy = (int)(this.getHeight() / 2.0);
      gx.fillOval(midx-2, 0, 5, 5);
      gx.fillOval(midx-2, this.getHeight()-5, 5, 5);
      gx.fillOval(0, midy-2, 5, 5);
      gx.fillOval(this.getWidth()-5, midy-2, 5, 5);
    }
  }

  /**
   * Overides default read object in order to reload icons.
   * This is necessary because for some strange reason animated
   * gifs stop being animated after being serialized/deserialized.
   *
   * @param ois an <code>ObjectInputStream</code> value
   * @exception IOException if an error occurs
   * @exception ClassNotFoundException if an error occurs
   */
  private void readObject(ObjectInputStream ois) 
    throws IOException, ClassNotFoundException {
    try {
      ois.defaultReadObject();
      remove(m_visualLabel);
      m_visualLabel = new JLabel(m_icon);
      loadIcons(m_iconPath, m_animatedIconPath);
      add(m_visualLabel, BorderLayout.CENTER);
      Dimension d = m_visualLabel.getPreferredSize();
      Dimension d2 = new Dimension((int)d.getWidth() + 10, 
				   (int)d.getHeight() + 10);
      setMinimumSize(d2);
      setPreferredSize(d2);
      setMaximumSize(d2);   
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
