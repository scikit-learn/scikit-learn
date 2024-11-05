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
 *    BeanInstance.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.beans.Beans;
import java.io.Serializable;
import java.util.Vector;

import javax.swing.JComponent;

/**
 * Class that manages a set of beans.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version  $Revision: 4829 $
 * @since 1.0
 */
public class BeanInstance
  implements Serializable {

  /** for serialization */
  private static final long serialVersionUID = -7575653109025406342L;

  /**
   * class variable holding all the beans
   */
  private static Vector COMPONENTS = new Vector();

  public static final int IDLE = 0;
  public static final int BEAN_EXECUTING = 1;
  
  /**
   * Holds the bean encapsulated in this instance
   */
  private Object m_bean;
  private int m_x;
  private int m_y;


  /**
   * Reset the list of beans
   */
  public static void reset(JComponent container) {
    // remove beans from container if necessary
    removeAllBeansFromContainer(container);
    COMPONENTS = new Vector();
  }

  /**
   * Removes all beans from containing component
   *
   * @param container a <code>JComponent</code> value
   */
  public static void removeAllBeansFromContainer(JComponent container) {
    if (container != null) {
      if (COMPONENTS != null) {
	for (int i = 0; i < COMPONENTS.size(); i++) {
	  BeanInstance tempInstance = (BeanInstance)COMPONENTS.elementAt(i);
	  Object tempBean = tempInstance.getBean();
	  if (Beans.isInstanceOf(tempBean, JComponent.class)) {
	    container.remove((JComponent)tempBean);
	  }
	}
      }
      container.revalidate();
    }
  }

  /**
   * Adds all beans to the supplied component
   *
   * @param container a <code>JComponent</code> value
   */
  public static void addAllBeansToContainer(JComponent container) {
    if (container != null) {
      if (COMPONENTS != null) {
	for (int i = 0; i < COMPONENTS.size(); i++) {
	  BeanInstance tempInstance = (BeanInstance)COMPONENTS.elementAt(i);
	  Object tempBean = tempInstance.getBean();
	  if (Beans.isInstanceOf(tempBean, JComponent.class)) {
	    container.add((JComponent)tempBean);
	  }
	}
      }
      container.revalidate();
    }
  }

  /**
   * Return the list of displayed beans
   *
   * @return a vector of beans
   */
  public static Vector getBeanInstances() {
    return COMPONENTS;
  }

  /**
   * Describe <code>setBeanInstances</code> method here.
   *
   * @param beanInstances a <code>Vector</code> value
   * @param container a <code>JComponent</code> value
   */
  public static void setBeanInstances(Vector beanInstances, 
				      JComponent container) {
    reset(container);
    
    if (container != null) {
      for (int i = 0; i < beanInstances.size(); i++) {
	Object bean = ((BeanInstance)beanInstances.elementAt(i)).getBean();
	if (Beans.isInstanceOf(bean, JComponent.class)) {
	  container.add((JComponent)bean);
	}
      }
      container.revalidate();
      container.repaint();
    }
    COMPONENTS = beanInstances;
  }

  /**
   * Renders the textual labels for the beans.
   *
   * @param gx a <code>Graphics</code> object on which to render
   * the labels
   */
  public static void paintLabels(Graphics gx) {
    gx.setFont(new Font(null, Font.PLAIN, 9));
    FontMetrics fm = gx.getFontMetrics();
    int hf = fm.getAscent();
    for (int i = 0; i < COMPONENTS.size(); i++) {
      BeanInstance bi = (BeanInstance)COMPONENTS.elementAt(i);
      if (!(bi.getBean() instanceof Visible)) {
	continue;
      }
      int cx = bi.getX(); int cy = bi.getY();
      int width = ((JComponent)bi.getBean()).getWidth();
      int height = ((JComponent)bi.getBean()).getHeight();
      String label = ((Visible)bi.getBean()).getVisual().getText();
      int labelwidth = fm.stringWidth(label);
      if (labelwidth < width) {
	gx.drawString(label, (cx+(width/2)) - (labelwidth / 2), cy+height+hf+2);
      } else {
	// split label

	// find mid point
	int mid = label.length() / 2;
	// look for split point closest to the mid
	int closest = label.length();
	int closestI = -1;
	for (int z = 0; z < label.length(); z++) {
	  if (label.charAt(z) < 'a') {
	    if (Math.abs(mid - z) < closest) {
	      closest = Math.abs(mid - z);
	      closestI = z;
	    }
	  }
	}
	if (closestI != -1) {
	  String left = label.substring(0, closestI);
	  String right = label.substring(closestI, label.length());
	  if (left.length() > 1 && right.length() > 1) {
	    gx.drawString(left, (cx+(width/2)) - (fm.stringWidth(left) / 2), 
			  cy+height+(hf * 1)+2);
	    gx.drawString(right, (cx+(width/2)) - (fm.stringWidth(right) / 2), 
			  cy+height+(hf * 2)+2);
	  } else {
	    gx.drawString(label, (cx+(width/2)) - (fm.stringWidth(label) / 2), 
			  cy+height+(hf * 1)+2);
	  }
	} else {
	  gx.drawString(label, (cx+(width/2)) - (fm.stringWidth(label) / 2), 
			cy+height+(hf * 1)+2);
	}
      }
    }
  }

  /**
   * Looks for a bean (if any) whose bounds contain the supplied point
   *
   * @param p a point
   * @return a bean that contains the supplied point or null if no bean
   * contains the point
   */
  public static BeanInstance findInstance(Point p) {
    Rectangle tempBounds = new Rectangle();
    for (int i=0; i < COMPONENTS.size(); i++) {
      
      BeanInstance t = (BeanInstance)COMPONENTS.elementAt(i);
      JComponent temp = (JComponent)t.getBean();
				      
      tempBounds = temp.getBounds(tempBounds);
      if (tempBounds.contains(p)) {
	return t;
      }
    }
    return null;
  }

  /**
   * Looks for all beans (if any) located within
   * the supplied bounding box. Also adjusts the
   * bounding box to be a tight fit around all
   * contained beans
   *
   * @param boundingBox the bounding rectangle
   * @return a Vector of BeanInstances
   */
  public static Vector findInstances(Rectangle boundingBox) {
    Graphics gx = null;
    FontMetrics fm = null;
    
    int centerX, centerY;
    int startX, startY, endX, endY;
    startX = (int)boundingBox.getX();
    startY = (int)boundingBox.getY();
    endX = (int)boundingBox.getMaxX();
    endY = (int)boundingBox.getMaxY();
    int minX = Integer.MAX_VALUE;
    int minY = Integer.MAX_VALUE;
    int maxX = Integer.MIN_VALUE;
    int maxY = Integer.MIN_VALUE;
    Vector result = new Vector();
    for (int i = 0; i < COMPONENTS.size(); i++) {
      BeanInstance t = (BeanInstance)COMPONENTS.elementAt(i);
      centerX = t.getX() + (t.getWidth()/2);
      centerY = t.getY() + (t.getHeight()/2);
      if (boundingBox.contains(centerX, centerY)) {
	result.addElement(t);


	// adjust bounding box stuff
//	int hf = 0;
	if (gx == null) {
	  gx = ((JComponent)t.getBean()).getGraphics();
	  gx.setFont(new Font(null, Font.PLAIN, 9));
	  fm = gx.getFontMetrics();
//	  hf = fm.getAscent();
	}
	String label = "";
	if (t.getBean() instanceof Visible) {
	  label = ((Visible)t.getBean()).getVisual().getText();
	}
	int labelwidth = fm.stringWidth(label);
	int heightMultiplier = (labelwidth > t.getWidth())
	? 2
	: 1;
	/*if (label.length() == 0) {
	  heightMultiplier = 0;
	}*/
	int brx = 0;
	int blx = 0;
	if (centerX - (labelwidth / 2) - 2 < t.getX()) {
	  blx = (centerX - (labelwidth / 2) - 2);
	  brx = centerX + (labelwidth / 2) + 2;
	} else {
	  blx = t.getX() - 2;
	  brx = t.getX() + t.getWidth() + 2;
	}

	if (blx < minX) {
	  minX = blx;
	}
	if (brx > maxX) {
	  maxX = brx;
	}
	if (t.getY() - 2 < minY) {
	  minY = t.getY() - 2;
	}
	if (t.getY() + t.getHeight() + 2 > maxY) {
	  maxY = t.getY() + t.getHeight()  + 2;
	}
      }
    }
    boundingBox.setBounds(minX, minY, maxX - minX, maxY - minY);
    
    return result;
  }

  /**
   * Creates a new <code>BeanInstance</code> instance.
   *
   * @param container a <code>JComponent</code> to add the bean to
   * @param bean the bean to add
   * @param x the x coordinate of the bean
   * @param y the y coordinate of the bean
   */
  public BeanInstance(JComponent container, Object bean, int x, int y) {
    m_bean = bean;
    m_x = x;
    m_y = y;
    addBean(container);
  }

  /**
   * Creates a new <code>BeanInstance</code> instance given the fully
   * qualified name of the bean
   *
   * @param container a <code>JComponent</code> to add the bean to
   * @param beanName the fully qualified name of the bean
   * @param x the x coordinate of the bean
   * @param y th y coordinate of the bean
   */
  public BeanInstance(JComponent container, String beanName, int x, int y) {
    m_x = x;
    m_y = y;
    
    // try and instantiate the named component
    try {
      m_bean = Beans.instantiate(null, beanName);
    } catch (Exception ex) {
      ex.printStackTrace();
      return;
    }

    addBean(container);
  }

  /**
   * Remove this bean from the list of beans and from the containing component
   *
   * @param container the <code>JComponent</code> that holds the bean
   */
  public void removeBean(JComponent container) {
    for (int i = 0; i < COMPONENTS.size(); i++) {
      if ((BeanInstance)COMPONENTS.elementAt(i) == this) {
	System.err.println("Removing bean");
	COMPONENTS.removeElementAt(i);
      }
    }
    if (container != null) {
      container.remove((JComponent)m_bean);
      container.revalidate();
      container.repaint();
    }
  }

  /**
   * Adds this bean to the global list of beans and
   * to the supplied container. The constructor
   * calls this method, so a client should not need
   * to unless they have called removeBean and then
   * wish to have it added again.
   *
   * @param container the Component on which this
   * BeanInstance will be displayed
   */
  public void addBean(JComponent container) {

    // do nothing if we are already in the list
    if (COMPONENTS.contains(this)) {
      return;
    }

    // Ignore invisible components
    if (!Beans.isInstanceOf(m_bean, JComponent.class)) {
      System.err.println("Component is invisible!");
      return;
    }
    
    COMPONENTS.addElement(this);
    
    // Position and layout the component
    JComponent c = (JComponent)m_bean;
    Dimension d = c.getPreferredSize();
    int dx = (int)(d.getWidth() / 2);
    int dy = (int)(d.getHeight() / 2);
    m_x -= dx;
    m_y -= dy;
    c.setLocation(m_x, m_y);
    //    c.doLayout();
    c.validate();
    //    bp.addBean(c);
    //    c.repaint();
    if (container != null) {
      container.add(c);
      container.revalidate();
    }
  }

  /**
   * Gets the bean encapsulated in this instance
   *
   * @return an <code>Object</code> value
   */
  public Object getBean() {
    return m_bean;
  }

  /**
   * Gets the x coordinate of this bean
   *
   * @return an <code>int</code> value
   */
  public int getX() {
    return m_x;
  }

  /**
   * Gets the y coordinate of this bean
   *
   * @return an <code>int</code> value
   */
  public int getY() {
    return m_y;
  }

  /**
   * Gets the width of this bean
   *
   * @return an <code>int</code> value
   */
  public int getWidth() {
    return ((JComponent)m_bean).getWidth();
  }

  /**
   * Gets the height of this bean
   *
   * @return an <code>int</code> value
   */
  public int getHeight() {
    return ((JComponent)m_bean).getHeight();
  }
 
  /**
   * Set the x and y coordinates of this bean
   *
   * @param newX the x coordinate
   * @param newY the y coordinate
   */
  public void setXY(int newX, int newY) {
    setX(newX);
    setY(newY);
    if (getBean() instanceof MetaBean) {
      ((MetaBean)getBean()).shiftBeans(this, false);
    }
  }

  /**
   * Sets the x coordinate of this bean
   *
   * @param newX an <code>int</code> value
   */
  public void setX(int newX) {
    m_x = newX;
    ((JComponent)m_bean).setLocation(m_x, m_y);
    ((JComponent)m_bean).validate();
  }

  /**
   * Sets the y coordinate of this bean
   *
   * @param newY an <code>int</code> value
   */
  public void setY(int newY) {
    m_y = newY;
    ((JComponent)m_bean).setLocation(m_x, m_y);
    ((JComponent)m_bean).validate();
  }
}
