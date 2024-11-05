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
 * CheckBoxList.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.gui;

import java.awt.Component;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Vector;

import javax.swing.DefaultListModel;
import javax.swing.JCheckBox;
import javax.swing.JList;
import javax.swing.ListCellRenderer;
import javax.swing.ListModel;

/**
 * An extended JList that contains CheckBoxes. If necessary a CheckBoxListItem
 * wrapper is added around the displayed object in any of the Model methods, 
 * e.g., addElement. For methods returning objects the opposite takes place, 
 * the wrapper is removed and only the payload object is returned.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.3 $
 */
public class CheckBoxList
  extends JList {

  /** for serialization */
  private static final long serialVersionUID = -4359573373359270258L;
  
  /**
   * represents an item in the CheckBoxListModel
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 1.3 $
   */
  protected class CheckBoxListItem {
    
    /** whether item is checked or not */
    private boolean m_Checked = false;
    
    /** the actual object */
    private Object m_Content = null;
    
    /**
     * initializes the item with the given object and initially unchecked
     * 
     * @param o		the content object
     */
    public CheckBoxListItem(Object o) {
      this(o, false);
    }
    
    /**
     * initializes the item with the given object and whether it's checked
     * initially
     * 
     * @param o		the content object
     * @param checked	whether the item should be checked initially
     */
    public CheckBoxListItem(Object o, boolean checked) {
      m_Checked = checked;
      m_Content = o;
    }
    
    /**
     * returns the content object
     */
    public Object getContent() {
      return m_Content;
    }
    
    /**
     * sets the checked state of the item
     */
    public void setChecked(boolean value) {
      m_Checked = value;
    }
    
    /**
     * returns the checked state of the item
     */
    public boolean getChecked() {
      return m_Checked;
    }
    
    /**
     * returns the string representation of the content object
     */
    public String toString() {
      return m_Content.toString();
    }
    
    /**
     * returns true if the "payload" objects of the current and the given
     * CheckBoxListItem are the same.
     * 
     * @param o		the CheckBoxListItem to check
     * @throws IllegalArgumentException if the provided object is not a CheckBoxListItem
     */
    public boolean equals(Object o) {
      if (!(o instanceof CheckBoxListItem))
	throw new IllegalArgumentException("Must be a CheckBoxListItem!");
      
      return getContent().equals(((CheckBoxListItem) o).getContent());
    }
  }
  
  /**
   * A specialized model.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 1.3 $
   */
  public class CheckBoxListModel
    extends DefaultListModel {

    /** for serialization */
    private static final long serialVersionUID = 7772455499540273507L;
    
    /**
     * initializes the model with no data.
     */
    public CheckBoxListModel() {
      super();
    }
    
    /**
     * Constructs a CheckBoxListModel from an array of objects and then applies 
     * setModel to it.
     * 
     * @param listData	the data to use
     */
    public CheckBoxListModel(Object[] listData) {
      for (int i = 0; i < listData.length; i++)
        addElement(listData[i]);
    }
    
    /**
     * Constructs a CheckBoxListModel from a Vector and then applies setModel 
     * to it.
     */
    public CheckBoxListModel(Vector listData) {
      for (int i = 0; i < listData.size(); i++)
        addElement(listData.get(i));
    }
    
    /**
     * Inserts the specified element at the specified position in this list.
     * 
     * @param index	index at which the specified element is to be inserted
     * @param element	element to be inserted
     */
    public void add(int index, Object element) {
      if (!(element instanceof CheckBoxListItem))
	super.add(index, new CheckBoxListItem(element));
      else
	super.add(index, element);
    }
    
    /**
     * Adds the specified component to the end of this list.
     * 
     * @param obj 	the component to be added
     */
    public void addElement(Object obj) {
      if (!(obj instanceof CheckBoxListItem))
	super.addElement(new CheckBoxListItem(obj));
      else
	super.addElement(obj);
    }
    
    /**
     * Tests whether the specified object is a component in this list.
     * 
     * @param elem	the element to check
     * @return		true if the element is in the list
     */
    public boolean contains(Object elem) {
      if (!(elem instanceof CheckBoxListItem))
	return super.contains(new CheckBoxListItem(elem));
      else
	return super.contains(elem);
    }
    
    /**
     * Copies the components of this list into the specified array.
     * 
     * @param anArray	the array into which the components get copied
     * @throws IndexOutOfBoundsException if the array is not big enough
     */
    public void copyInto(Object[] anArray) {
      if (anArray.length < getSize())
	throw new IndexOutOfBoundsException("Array not big enough!");
      
      for (int i = 0; i < getSize(); i++)
	anArray[i] = ((CheckBoxListItem) getElementAt(i)).getContent();
    }
    
    /**
     * Returns the component at the specified index. Throws an 
     * ArrayIndexOutOfBoundsException if the index is negative or not less 
     * than the size of the list.
     * 
     * @param index	an index into this list
     * @return 		the component at the specified index
     * @throws ArrayIndexOutOfBoundsException
     */
    public Object elementAt(int index) {
      return ((CheckBoxListItem) super.elementAt(index)).getContent();
    }
    
    /**
     * Returns the first component of this list. Throws a 
     * NoSuchElementException if this vector has no components.
     * 
     * @return		the first component of this list
     * @throws NoSuchElementException
     */
    public Object firstElement() {
      return ((CheckBoxListItem) super.firstElement()).getContent();
    }
    
    /**
     * Returns the element at the specified position in this list.
     * 
     * @param index of element to return
     * @throws ArrayIndexOutOfBoundsException
     */
    public Object get(int index) {
      return ((CheckBoxListItem) super.get(index)).getContent();
    }
    
    /**
     * Returns the component at the specified index.
     * 
     * @param index 	an index into this list
     * @return 		the component at the specified index 
     * @throws ArrayIndexOutOfBoundsException
     */
    public Object getElementAt(int index) {
      return ((CheckBoxListItem) super.getElementAt(index)).getContent();
    }
    
    /**
     * Searches for the first occurrence of elem.
     * 
     * @param elem	an object
     * @return 		the index of the first occurrence of the argument in this list; 
     * 			returns -1 if the object is not found
     */
    public int indexOf(Object elem) {
      if (!(elem instanceof CheckBoxListItem))
	return super.indexOf(new CheckBoxListItem(elem));
      else
	return super.indexOf(elem);
    }
    
    /**
     * Searches for the first occurrence of elem, beginning the search at index.
     * 
     * @param elem 	the desired component
     * @param index	the index from which to begin searching
     * @return		the index where the first occurrence of elem  is found after index; 
     * 			returns -1  if the elem is not found in the list
     */
    public int indexOf(Object elem, int index) {
      if (!(elem instanceof CheckBoxListItem))
	return super.indexOf(new CheckBoxListItem(elem), index);
      else
	return super.indexOf(elem, index);
    }
    
    /**
     * Inserts the specified object as a component in this list at the 
     * specified index.
     * 
     * @param obj	the component to insert
     * @param index	where to insert the new component
     * @throws ArrayIndexOutOfBoundsException 
     */
    public void insertElementAt(Object obj, int index) {
      if (!(obj instanceof CheckBoxListItem))
	super.insertElementAt(new CheckBoxListItem(obj), index);
      else
	super.insertElementAt(obj, index);
    }
    
    /**
     * Returns the last component of the list. Throws a NoSuchElementException 
     * if this vector has no components.
     * 
     * @return 		the last component of the list
     * @throws NoSuchElementException
     */
    public Object lastElement() {
      return ((CheckBoxListItem) super.lastElement()).getContent();
    }
    
    /**
     * Returns the index of the last occurrence of elem.
     * 
     * @param elem	the desired component
     * @return		the index of the last occurrence of elem  in the list; 
     * 			returns -1 if the object is not found
     */
    public int lastIndexOf(Object elem) {
      if (!(elem instanceof CheckBoxListItem))
	return super.lastIndexOf(new CheckBoxListItem(elem));
      else
	return super.lastIndexOf(elem);
    }
    
    /**
     * Searches backwards for elem, starting from the specified index, 
     * and returns an index to it.
     * 
     * @param elem	the desired component
     * @param index	the index to start searching from
     * @return		the index of the last occurrence of the elem in this 
     * 			list at position less than index; returns -1 if the 
     * 			object is not found
     */
    public int lastIndexOf(Object elem, int index) {
      if (!(elem instanceof CheckBoxListItem))
	return super.lastIndexOf(new CheckBoxListItem(elem), index);
      else
	return super.lastIndexOf(elem, index);
    }
    
    /**
     * Removes the element at the specified position in this list. Returns the 
     * element that was removed from the list.
     * 
     * @param index	the index of the element to removed
     * @throws ArrayIndexOutOfBoundsException
     */
    public Object remove(int index) {
      return ((CheckBoxListItem) super.remove(index)).getContent();
    }
    
    /**
     * Removes the first (lowest-indexed) occurrence of the argument from this 
     * list.
     * 
     * @param obj	the component to be removed
     * @return		true if the argument was a component of this list; 
     * 			false otherwise
     */
    public boolean removeElement(Object obj) {
      if (!(obj instanceof CheckBoxListItem))
	return super.removeElement(new CheckBoxListItem(obj));
      else
	return super.removeElement(obj);
    }
    
    /**
     * Replaces the element at the specified position in this list with the 
     * specified element.
     * 
     * @param index	index of element to replace
     * @param element	element to be stored at the specified position
     * @throws ArrayIndexOutOfBoundsException 
     */
    public Object set(int index, Object element) {
      if (!(element instanceof CheckBoxListItem))
	return ((CheckBoxListItem) super.set(index, new CheckBoxListItem(element))).getContent();
      else
	return ((CheckBoxListItem) super.set(index, element)).getContent();
    }
    
    /**
     * Sets the component at the specified index of this list to be the 
     * specified object. The previous component at that position is discarded.
     * 
     * @param obj	what the component is to be set to
     * @param index	the specified index
     * @throws ArrayIndexOutOfBoundsException
     */
    public void setElementAt(Object obj, int index) {
      if (!(obj instanceof CheckBoxListItem))
	super.setElementAt(new CheckBoxListItem(obj), index);
      else
	super.setElementAt(obj, index);
    }
    
    /**
     * Returns an array containing all of the elements in this list in the 
     * correct order.
     * 
     * @return 		an array containing the elements of the list
     */
    public Object[] toArray() {
      Object[]		result;
      Object[]		internal;
      int		i;
      
      internal = super.toArray();
      result   = new Object[internal.length];
      
      for (i = 0; i < internal.length; i++)
	result[i] = ((CheckBoxListItem) internal[i]).getContent();
      
      return result;
    }
    
    /**
     * returns the checked state of the element at the given index
     * 
     * @param index	the index of the element to return the checked state for
     * @return		the checked state of the specifed element
     */
    public boolean getChecked(int index) {
      return ((CheckBoxListItem) super.getElementAt(index)).getChecked();
    }
    
    /**
     * sets the checked state of the element at the given index
     * 
     * @param index	the index of the element to set the checked state for
     * @param checked	the new checked state
     */
    public void setChecked(int index, boolean checked) {
      ((CheckBoxListItem) super.getElementAt(index)).setChecked(checked);
    }
  }
  
  /**
   * A specialized CellRenderer for the CheckBoxList
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 1.3 $
   * @see CheckBoxList
   */
  public class CheckBoxListRenderer 
    extends JCheckBox 
    implements ListCellRenderer {

    /** for serialization */
    private static final long serialVersionUID = 1059591605858524586L;
  
    /**
     * Return a component that has been configured to display the specified 
     * value.
     * 
     * @param list	The JList we're painting.
     * @param value	The value returned by list.getModel().getElementAt(index).
     * @param index	The cells index.
     * @param isSelected	True if the specified cell was selected.
     * @param cellHasFocus	True if the specified cell has the focus.
     * @return 		A component whose paint() method will render the 
     * 			specified value.
     */
    public Component getListCellRendererComponent(
	JList list,
        Object value,
        int index,
        boolean isSelected,
        boolean cellHasFocus) {
      
      setText(value.toString());
      setSelected(((CheckBoxList) list).getChecked(index));
      setBackground(isSelected ? list.getSelectionBackground() : list.getBackground());
      setForeground(isSelected ? list.getSelectionForeground() : list.getForeground());
      setFocusPainted(false);
      
      return this;
    }
  }
  
  /**
   * initializes the list with an empty CheckBoxListModel
   */
  public CheckBoxList() {
    this(null);
  }
  
  /**
   * initializes the list with the given CheckBoxListModel
   * 
   * @param model	the model to initialize with
   */
  public CheckBoxList(CheckBoxListModel model) {
    super();
    
    if (model == null)
      model = this.new CheckBoxListModel();
    
    setModel(model);
    setCellRenderer(new CheckBoxListRenderer());
    
    addMouseListener(new MouseAdapter() {
      public void mousePressed(MouseEvent e) {
	int index = locationToIndex(e.getPoint());
	
	if (index != -1) {
	  setChecked(index, !getChecked(index));
	  repaint();
	}
      }
    });
    
    addKeyListener(new KeyAdapter() {
      public void keyTyped(KeyEvent e) {
        if ( (e.getKeyChar() == ' ') && (e.getModifiers() == 0) ) {
          int index = getSelectedIndex();
          setChecked(index, !getChecked(index));
          e.consume();
	  repaint();
        }
      }
    });
  }
  
  /**
   * sets the model - must be an instance of CheckBoxListModel
   * 
   * @param model			the model to use
   * @throws IllegalArgumentException 	if the model is not an instance of
   * 					CheckBoxListModel
   * @see CheckBoxListModel
   */
  public void setModel(ListModel model) {
    if (!(model instanceof CheckBoxListModel))
      throw new IllegalArgumentException("Model must be an instance of CheckBoxListModel!");
    
    super.setModel(model);
  }
  
  /**
   * Constructs a CheckBoxListModel from an array of objects and then applies 
   * setModel to it.
   * 
   * @param listData	the data to use
   */
  public void setListData(Object[] listData) {
    setModel(new CheckBoxListModel(listData));
  }
  
  /**
   * Constructs a CheckBoxListModel from a Vector and then applies setModel 
   * to it.
   */
  public void setListData(Vector listData) {
    setModel(new CheckBoxListModel(listData));
  }
  
  /**
   * returns the checked state of the element at the given index
   * 
   * @param index	the index of the element to return the checked state for
   * @return		the checked state of the specifed element
   */
  public boolean getChecked(int index) {
    return ((CheckBoxListModel) getModel()).getChecked(index);
  }
  
  /**
   * sets the checked state of the element at the given index
   * 
   * @param index	the index of the element to set the checked state for
   * @param checked	the new checked state
   */
  public void setChecked(int index, boolean checked) {
    ((CheckBoxListModel) getModel()).setChecked(index, checked);
  }
  
  /**
   * returns an array with the indices of all checked items
   * 
   * @return		the indices of all items that are currently checked
   */
  public int[] getCheckedIndices() {
    Vector	list;
    int[]	result;
    int		i;
    
    // traverse over model
    list = new Vector();
    for (i = 0; i < getModel().getSize(); i++) {
      if (getChecked(i))
	list.add(new Integer(i));
    }
    
    // generate result array
    result = new int[list.size()];
    for (i = 0; i < list.size(); i++) {
      result[i] = ((Integer) list.get(i)).intValue();
    }
    
    return result;
  }
}
